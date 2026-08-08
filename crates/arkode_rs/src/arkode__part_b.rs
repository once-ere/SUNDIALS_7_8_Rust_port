/* =================================================================
   FRAGMENT: src/arkode/arkode.c, PART B (every function whose C
   definition begins at line 2000 or later).  This file contains ONLY
   function definitions; `arkode.rs` supplies the module doc comment,
   the `use` statements and any module-scope constants.  Every constant
   used below (ZERO, ONE, TWO, FOUR, TENTH, ONEPSM, FUZZ_FACTOR,
   H0_LBFACTOR, H0_UBFACTOR, H0_BIAS, H0_ITERS, ARK_* return codes,
   PREDICT_AGAIN/CONV_FAIL/TRY_AGAIN/..., MSG_ARK_*) comes from the
   frozen contract `crate::arkode_impl`.

   Reference build: SUNDIALS_LOGGING_LEVEL = 2, so every SUNLogInfo /
   SUNLogInfoIf / SUNLogDebug / SUNLogExtraDebug call site in the C is
   omitted at translation time; ARK_WARNING messages (none in this
   part) would still go through arkProcessError.
   ================================================================= */

/*---------------------------------------------------------------
  arkInitialSetup

  This routine performs all necessary items to prepare ARKODE for
  the first internal step after initialization, reinitialization,
  a reset() call, or a resize() call, including:
  - input consistency checks
  - (re)initializes the stepper
  - computes error and residual weights
  - (re)initialize the interpolation structure
  - checks for valid initial step input or estimates first step
  - checks for approach to tstop
  - checks for root near t0
  ---------------------------------------------------------------*/
pub fn arkInitialSetup(ark_mem: &ARKodeMem, tout: sunrealtype) -> i32 {
    /* Is tout too close to tn? */
    let (tcur, uround) = {
        let m = ark_mem.borrow();
        (m.tcur, m.uround)
    };
    let tdist = SUNRabs(tout - tcur);
    let tround = uround * SUNMAX(SUNRabs(tcur), SUNRabs(tout));

    if tdist == ZERO || tdist < TWO * tround {
        arkProcessError(
            Some(ark_mem),
            ARK_TOO_CLOSE,
            line!() as i32,
            "arkInitialSetup",
            file!(),
            MSG_ARK_TOO_CLOSE,
        );
        return ARK_TOO_CLOSE;
    }

    /* Check that user has supplied an initial step size if fixedstep mode is on */
    let (fixedstep, hin) = {
        let m = ark_mem.borrow();
        (m.fixedstep, m.hin)
    };
    if fixedstep && (hin == ZERO) {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!() as i32,
            "arkInitialSetup",
            file!(),
            "Fixed step mode enabled, but no step size set",
        );
        return ARK_ILL_INPUT;
    }

    /* Perform additional N_Vector checks here, now that ARKODE has been
    fully configured by the user */
    if !arkCheckNvectorOptional(ark_mem) {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!() as i32,
            "arkInitialSetup",
            file!(),
            MSG_ARK_BAD_NVECTOR,
        );
        return ARK_ILL_INPUT;
    }

    /* Test input tstop for legality (correct direction of integration) */
    if ark_mem.borrow().tstopset {
        let (h, tstop, tcur) = {
            let m = ark_mem.borrow();
            (m.h, m.tstop, m.tcur)
        };
        let htmp = if h == ZERO { tout - tcur } else { h };
        if (tstop - tcur) * htmp <= ZERO {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!() as i32,
                "arkInitialSetup",
                file!(),
                &MSG_ARK_BAD_TSTOP(tstop, tcur),
            );
            return ARK_ILL_INPUT;
        }
    }

    /* Check to see if y0 satisfies constraints */
    let constraints = ark_mem.borrow().constraints.clone();
    if let Some(constraints) = constraints {
        let (yn, tempv1) = {
            let m = ark_mem.borrow();
            (
                m.yn.clone().expect("yn allocated"),
                m.tempv1.clone().expect("tempv1 allocated"),
            )
        };
        let conOK = N_VConstrMask(&constraints, &yn, &tempv1);
        if !conOK {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!() as i32,
                "arkInitialSetup",
                file!(),
                MSG_ARK_Y0_FAIL_CONSTR,
            );
            return ARK_ILL_INPUT;
        }
    }

    /* Load initial error weights.

    C: `ark_mem->efun(ark_mem->yn, ark_mem->ewt, ark_mem->e_data)`, where
    `e_data` aliases `ark_mem` for the built-in weight functions and
    `ark_mem->user_data` when the user supplied `efun`.  A `Box` token
    cannot alias (deviation class 6), so the user-efun case passes the
    CURRENT `user_data` box and the built-in case passes `e_data` (which
    holds a boxed `ARKodeMem` handle clone).  The box is taken out of the
    mem around the call and restored on every path. */
    let (yn, ewt) = {
        let m = ark_mem.borrow();
        (
            m.yn.clone().expect("yn allocated"),
            m.ewt.clone().expect("ewt allocated"),
        )
    };
    let (efun, user_efun) = {
        let m = ark_mem.borrow();
        (m.efun, m.user_efun)
    };
    let efun = efun.expect("efun set");
    let retval = if user_efun {
        let mut data = ark_mem.borrow_mut().user_data.take();
        let r = efun(&yn, &ewt, &mut data);
        ark_mem.borrow_mut().user_data = data;
        r
    } else {
        let mut data = ark_mem.borrow_mut().e_data.take();
        let r = efun(&yn, &ewt, &mut data);
        ark_mem.borrow_mut().e_data = data;
        r
    };
    if retval != 0 {
        if ark_mem.borrow().itol == ARK_WF {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!() as i32,
                "arkInitialSetup",
                file!(),
                MSG_ARK_EWT_FAIL,
            );
        } else {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!() as i32,
                "arkInitialSetup",
                file!(),
                MSG_ARK_BAD_EWT,
            );
        }
        return ARK_ILL_INPUT;
    }

    /* Set up the time stepper module if not done so already */
    if !ark_mem.borrow().preallocated {
        let step_init = ark_mem.borrow().step_init;
        let step_init = match step_init {
            None => {
                arkProcessError(
                    Some(ark_mem),
                    ARK_ILL_INPUT,
                    line!() as i32,
                    "arkInitialSetup",
                    file!(),
                    "Time stepper module is missing",
                );
                return ARK_ILL_INPUT;
            }
            Some(f) => f,
        };
        let init_type = ark_mem.borrow().init_type;
        let retval = step_init(ark_mem, init_type);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                retval,
                line!() as i32,
                "arkInitialSetup",
                file!(),
                "Error in initialization of time stepper module",
            );
            return retval;
        }
    }

    /* Load initial residual weights */
    if ark_mem.borrow().rwt_is_ewt {
        /* update pointer to ewt */
        let ewt = ark_mem.borrow().ewt.clone();
        ark_mem.borrow_mut().rwt = ewt;
    } else {
        let (yn, rwt) = {
            let m = ark_mem.borrow();
            (
                m.yn.clone().expect("yn allocated"),
                m.rwt.clone().expect("rwt allocated"),
            )
        };
        let (rfun, user_rfun) = {
            let m = ark_mem.borrow();
            (m.rfun, m.user_rfun)
        };
        let rfun = rfun.expect("rfun set");
        let retval = if user_rfun {
            let mut data = ark_mem.borrow_mut().user_data.take();
            let r = rfun(&yn, &rwt, &mut data);
            ark_mem.borrow_mut().user_data = data;
            r
        } else {
            let mut data = ark_mem.borrow_mut().r_data.take();
            let r = rfun(&yn, &rwt, &mut data);
            ark_mem.borrow_mut().r_data = data;
            r
        };
        if retval != 0 {
            if ark_mem.borrow().itol == ARK_WF {
                arkProcessError(
                    Some(ark_mem),
                    ARK_ILL_INPUT,
                    line!() as i32,
                    "arkInitialSetup",
                    file!(),
                    MSG_ARK_RWT_FAIL,
                );
            } else {
                arkProcessError(
                    Some(ark_mem),
                    ARK_ILL_INPUT,
                    line!() as i32,
                    "arkInitialSetup",
                    file!(),
                    MSG_ARK_BAD_RWT,
                );
            }
            return ARK_ILL_INPUT;
        }
    }

    /* Create default interpolation module (if needed) */
    let (interp_type, interp_present) = {
        let m = ark_mem.borrow();
        (m.interp_type, m.interp.is_some())
    };
    if interp_type != ARK_INTERP_NONE && !interp_present {
        let interp_degree = ark_mem.borrow().interp_degree;
        let interp = if interp_type == ARK_INTERP_LAGRANGE {
            crate::arkode_interp::arkInterpCreate_Lagrange(ark_mem, interp_degree)
        } else {
            crate::arkode_interp::arkInterpCreate_Hermite(ark_mem, interp_degree)
        };
        let is_none = interp.is_none();
        ark_mem.borrow_mut().interp = interp;
        if is_none {
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_FAIL,
                line!() as i32,
                "arkInitialSetup",
                file!(),
                "Unable to allocate interpolation module",
            );
            return ARK_MEM_FAIL;
        }
    }

    /* Fill initial interpolation data (if needed) */
    let interp = ark_mem.borrow().interp.clone();
    if let Some(interp) = interp.as_ref() {
        /* Stepper init may have limited the interpolation degree */
        let interp_degree = ark_mem.borrow().interp_degree;
        if crate::arkode_interp::arkInterpSetDegree(ark_mem, interp, interp_degree) != 0 {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!() as i32,
                "arkInitialSetup",
                file!(),
                "Unable to update interpolation polynomial degree",
            );
            return ARK_ILL_INPUT;
        }

        let tcur = ark_mem.borrow().tcur;
        if crate::arkode_interp::arkInterpInit(ark_mem, interp, tcur) != 0 {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!() as i32,
                "arkInitialSetup",
                file!(),
                "Unable to initialize interpolation module",
            );
            return ARK_ILL_INPUT;
        }
    }

    /* Check if the configuration requires interpolation */
    let (root_present, interp_present, tstopinterp) = {
        let m = ark_mem.borrow();
        (m.root_mem.is_some(), m.interp.is_some(), m.tstopinterp)
    };
    if root_present && !interp_present {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!() as i32,
            "arkInitialSetup",
            file!(),
            "Rootfinding requires an interpolation module",
        );
        return ARK_ILL_INPUT;
    }

    if tstopinterp && !interp_present {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!() as i32,
            "arkInitialSetup",
            file!(),
            "Stop time interpolation requires an interpolation module",
        );
        return ARK_ILL_INPUT;
    }

    /* Call stepper-provided initial step size estimation routine to fill
    ark_mem->hin, if applicable. */
    let (h0u, hin, fixedstep, step_H0) = {
        let m = ark_mem.borrow();
        (m.h0u, m.hin, m.fixedstep, m.step_H0)
    };
    if h0u == ZERO && hin == ZERO && !fixedstep && step_H0.is_some() {
        let step_H0 = step_H0.expect("step_H0 set");
        /* C passes `&(ark_mem->hin)` straight into the stepper; the port
        copies the field out, calls, and writes it back on every path
        (binding invariant B). */
        let mut hin_out = ark_mem.borrow().hin;
        let retval = step_H0(ark_mem, tout, &mut hin_out);
        ark_mem.borrow_mut().hin = hin_out;
        if retval != 0 {
            arkProcessError(
                Some(ark_mem),
                ARK_STEP_H0_FAIL,
                line!() as i32,
                "arkInitialSetup",
                file!(),
                "Failure in timestepping module h0 calculation",
            );
            return ARK_STEP_H0_FAIL;
        }
    }

    /* If fullrhs will be called (to estimate initial step, explicit steppers, Hermite
    interpolation module, and possibly (but not always) arkRootCheck1), then
    ensure that it is provided, and space is allocated for fn.  Otherwise,
    we should free ark_mem->fn if it is allocated. */
    let (call_fullrhs, h0u, hin, root_present) = {
        let m = ark_mem.borrow();
        (m.call_fullrhs, m.h0u, m.hin, m.root_mem.is_some())
    };
    if call_fullrhs || (h0u == ZERO && hin == ZERO) || root_present {
        if ark_mem.borrow().step_fullrhs.is_none() {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!() as i32,
                "arkInitialSetup",
                file!(),
                MSG_ARK_MISSING_FULLRHS,
            );
            return ARK_ILL_INPUT;
        }

        let yn = ark_mem.borrow().yn.clone().expect("yn allocated");
        let mut fn_ = ark_mem.borrow_mut().fn_.take();
        let allocOK = arkAllocVec(ark_mem, &yn, &mut fn_);
        ark_mem.borrow_mut().fn_ = fn_;
        if !allocOK {
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_FAIL,
                line!() as i32,
                "arkInitialSetup",
                file!(),
                MSG_ARK_MEM_FAIL,
            );
            return ARK_MEM_FAIL;
        }
    } else if ark_mem.borrow().fn_.is_some() {
        let mut fn_ = ark_mem.borrow_mut().fn_.take();
        arkFreeVec(ark_mem, &mut fn_);
        ark_mem.borrow_mut().fn_ = fn_;
    }

    /* initialization complete */
    ark_mem.borrow_mut().initialized = SUNTRUE;

    /* Set initial step size */
    if ark_mem.borrow().h0u == ZERO {
        /* Check input h for validity */
        {
            let mut guard = ark_mem.borrow_mut();
            let m = &mut *guard;
            m.h = m.hin;
        }
        let (h, tcur) = {
            let m = ark_mem.borrow();
            (m.h, m.tcur)
        };
        if (h != ZERO) && ((tout - tcur) * h < ZERO) {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!() as i32,
                "arkInitialSetup",
                file!(),
                MSG_ARK_BAD_H0,
            );
            return ARK_ILL_INPUT;
        }

        /* Estimate initial h if not set */
        if h == ZERO {
            /* If necessary, temporarily set h as it is used to compute the tolerance
            in a potential mass matrix solve when computing the full rhs */
            {
                let mut guard = ark_mem.borrow_mut();
                let m = &mut *guard;
                m.h = SUNRabs(tout - m.tcur);
                if m.h == ZERO {
                    m.h = ONE;
                }
            }

            /* Estimate the first step size */
            let mut tout_hin = tout;
            let (tstopset, tstop, tcur) = {
                let m = ark_mem.borrow();
                (m.tstopset, m.tstop, m.tcur)
            };
            if tstopset && (tout - tcur) * (tout - tstop) > ZERO {
                tout_hin = tstop;
            }
            let hflag = arkHin(ark_mem, tout_hin);
            if hflag != ARK_SUCCESS {
                let istate = arkHandleFailure(ark_mem, hflag);
                return istate;
            }

            /* Use first step growth factor for estimated h */
            {
                let mut guard = ark_mem.borrow_mut();
                let m = &mut *guard;
                let hadapt_mem = m.hadapt_mem.as_mut().expect("hadapt_mem allocated");
                let etamx1 = hadapt_mem.etamx1;
                hadapt_mem.etamax = etamx1;
            }
        } else if ark_mem.borrow().nst == 0 {
            /* Use first step growth factor for user defined h */
            let mut guard = ark_mem.borrow_mut();
            let m = &mut *guard;
            let hadapt_mem = m.hadapt_mem.as_mut().expect("hadapt_mem allocated");
            let etamx1 = hadapt_mem.etamx1;
            hadapt_mem.etamax = etamx1;
        } else {
            /* Use standard growth factor (e.g., for reset) */
            let mut guard = ark_mem.borrow_mut();
            let m = &mut *guard;
            let hadapt_mem = m.hadapt_mem.as_mut().expect("hadapt_mem allocated");
            let growth = hadapt_mem.growth;
            hadapt_mem.etamax = growth;
        }

        /* Enforce step size bounds */
        {
            let mut guard = ark_mem.borrow_mut();
            let m = &mut *guard;
            let rh = SUNRabs(m.h) * m.hmax_inv;
            if rh > ONE {
                m.h /= rh;
            }
            let habs = SUNRabs(m.h);
            if habs < m.hmin {
                let scale = m.hmin / habs;
                m.h *= scale;
            }
        }

        /* Check for approach to tstop */
        {
            let mut guard = ark_mem.borrow_mut();
            let m = &mut *guard;
            if m.tstopset && ((m.tcur + m.h - m.tstop) * m.h > ZERO) {
                m.h = (m.tstop - m.tcur) * (ONE - FOUR * m.uround);
            }
        }

        /* Set initial time step factors */
        {
            let mut guard = ark_mem.borrow_mut();
            let m = &mut *guard;
            m.h0u = m.h;
            m.eta = ONE;
            m.hprime = m.h;
        }
    } else {
        /* If next step would overtake tstop, adjust stepsize */
        let mut guard = ark_mem.borrow_mut();
        let m = &mut *guard;
        if m.tstopset && ((m.tcur + m.hprime - m.tstop) * m.h > ZERO) {
            m.hprime = (m.tstop - m.tcur) * (ONE - FOUR * m.uround);
            m.eta = m.hprime / m.h;
        }
    }

    /* Check for zeros of root function g at and near t0. */
    let nrtfn = ark_mem
        .borrow()
        .root_mem
        .as_ref()
        .map(|r| r.nrtfn)
        .unwrap_or(0);
    if ark_mem.borrow().root_mem.is_some() && nrtfn > 0 {
        let retval = crate::arkode_root::arkRootCheck1(ark_mem);
        if retval != ARK_SUCCESS {
            return retval;
        }
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStopTests

  This routine performs relevant stopping tests:
  - check for root in last step
  - check if we passed tstop
  - check if we passed tout (NORMAL mode)
  - check if current tn was returned (ONE_STEP mode)
  - check if we are close to tstop
  (adjust step size if needed)
  ---------------------------------------------------------------*/
pub fn arkStopTests(
    ark_mem: &ARKodeMem,
    tout: sunrealtype,
    yout: &N_Vector,
    tret: &mut sunrealtype,
    itask: i32,
    ier: &mut i32,
) -> i32 {
    /* Estimate an infinitesimal time interval to be used as
    a roundoff for time quantities (based on current time
    and step size) */
    let troundoff = {
        let m = ark_mem.borrow();
        FUZZ_FACTOR * m.uround * (SUNRabs(m.tcur) + SUNRabs(m.h))
    };

    /* First, check for a root in the last step taken, other than the
    last root found, if any.  If itask = ARK_ONE_STEP and y(tn) was not
    returned because of an intervening root, return y(tn) now.     */
    if ark_mem.borrow().root_mem.is_some() {
        let nrtfn = ark_mem
            .borrow()
            .root_mem
            .as_ref()
            .expect("root_mem allocated")
            .nrtfn;
        if nrtfn > 0 {
            /* Shortcut to roots found in previous step */
            let irfndp = ark_mem
                .borrow()
                .root_mem
                .as_ref()
                .expect("root_mem allocated")
                .irfnd;

            /* If the full RHS was not computed in the last call to arkCompleteStep
            and roots were found in the previous step, then compute the full rhs
            for possible use in arkRootCheck2 (not always necessary) */
            let fn_is_current = ark_mem.borrow().fn_is_current;
            if !fn_is_current && irfndp != 0 {
                let (step_fullrhs, tn, yn, fn_) = {
                    let m = ark_mem.borrow();
                    (
                        m.step_fullrhs,
                        m.tn,
                        m.yn.clone().expect("yn allocated"),
                        m.fn_.clone().expect("fn allocated"),
                    )
                };
                let step_fullrhs = step_fullrhs.expect("step_fullrhs set");
                let retval = step_fullrhs(ark_mem, tn, &yn, &fn_, ARK_FULLRHS_END);
                if retval != 0 {
                    /* NOTE: upstream C passes MSG_ARK_RHSFUNC_FAILED (which carries a
                    SUN_FORMAT_G conversion) with NO argument here -- undefined
                    behavior in C.  The port supplies ark_mem->tcur, the value every
                    other MSG_ARK_RHSFUNC_FAILED call site uses. */
                    let tcur = ark_mem.borrow().tcur;
                    arkProcessError(
                        Some(ark_mem),
                        ARK_RHSFUNC_FAIL,
                        line!() as i32,
                        "arkStopTests",
                        file!(),
                        &MSG_ARK_RHSFUNC_FAILED(tcur),
                    );
                    *ier = ARK_RHSFUNC_FAIL;
                    return 1;
                }
                ark_mem.borrow_mut().fn_is_current = SUNTRUE;
            }

            let retval = crate::arkode_root::arkRootCheck2(ark_mem);

            if retval == CLOSERT {
                let tlo = ark_mem
                    .borrow()
                    .root_mem
                    .as_ref()
                    .expect("root_mem allocated")
                    .tlo;
                arkProcessError(
                    Some(ark_mem),
                    ARK_ILL_INPUT,
                    line!() as i32,
                    "arkStopTests",
                    file!(),
                    &MSG_ARK_CLOSE_ROOTS(tlo),
                );
                *ier = ARK_ILL_INPUT;
                return 1;
            } else if retval == ARK_RTFUNC_FAIL {
                let tlo = ark_mem
                    .borrow()
                    .root_mem
                    .as_ref()
                    .expect("root_mem allocated")
                    .tlo;
                arkProcessError(
                    Some(ark_mem),
                    ARK_RTFUNC_FAIL,
                    line!() as i32,
                    "arkStopTests",
                    file!(),
                    &MSG_ARK_RTFUNC_FAILED(tlo),
                );
                *ier = ARK_RTFUNC_FAIL;
                return 1;
            } else if retval == RTFOUND {
                let tlo = ark_mem
                    .borrow()
                    .root_mem
                    .as_ref()
                    .expect("root_mem allocated")
                    .tlo;
                ark_mem.borrow_mut().tretlast = tlo;
                *tret = tlo;
                *ier = ARK_ROOT_RETURN;
                return 1;
            }

            /* If tn is distinct from tretlast (within roundoff),
            check remaining interval for roots */
            let (tcur, tretlast) = {
                let m = ark_mem.borrow();
                (m.tcur, m.tretlast)
            };
            if SUNRabs(tcur - tretlast) > troundoff {
                let retval = crate::arkode_root::arkRootCheck3(ark_mem, tout, itask);

                if retval == ARK_SUCCESS {
                    /* no root found */
                    ark_mem
                        .borrow_mut()
                        .root_mem
                        .as_mut()
                        .expect("root_mem allocated")
                        .irfnd = 0;
                    if (irfndp == 1) && (itask == ARK_ONE_STEP) {
                        let (tcur, yn) = {
                            let m = ark_mem.borrow();
                            (m.tcur, m.yn.clone().expect("yn allocated"))
                        };
                        ark_mem.borrow_mut().tretlast = tcur;
                        *tret = tcur;
                        N_VScale(ONE, &yn, yout);
                        *ier = ARK_SUCCESS;
                        return 1;
                    }
                } else if retval == RTFOUND {
                    /* a new root was found */
                    let tlo = {
                        let mut m = ark_mem.borrow_mut();
                        let root_mem = m.root_mem.as_mut().expect("root_mem allocated");
                        root_mem.irfnd = 1;
                        root_mem.tlo
                    };
                    ark_mem.borrow_mut().tretlast = tlo;
                    *tret = tlo;
                    *ier = ARK_ROOT_RETURN;
                    return 1;
                } else if retval == ARK_RTFUNC_FAIL {
                    /* g failed */
                    let tlo = ark_mem
                        .borrow()
                        .root_mem
                        .as_ref()
                        .expect("root_mem allocated")
                        .tlo;
                    arkProcessError(
                        Some(ark_mem),
                        ARK_RTFUNC_FAIL,
                        line!() as i32,
                        "arkStopTests",
                        file!(),
                        &MSG_ARK_RTFUNC_FAILED(tlo),
                    );
                    *ier = ARK_RTFUNC_FAIL;
                    return 1;
                }
            }
        } /* end of root stop check */
    }

    /* Test for tn at tstop or near tstop */
    if ark_mem.borrow().tstopset {
        let (tcur, tstop, h) = {
            let m = ark_mem.borrow();
            (m.tcur, m.tstop, m.h)
        };
        if SUNRabs(tcur - tstop) <= troundoff {
            /* Ensure tout >= tstop, otherwise check for tout return below */
            if (tout - tstop) * h >= ZERO || SUNRabs(tout - tstop) <= troundoff {
                let (tstopinterp, interp_present) = {
                    let m = ark_mem.borrow();
                    (m.tstopinterp, m.interp.is_some())
                };
                if tstopinterp && interp_present {
                    *ier = ARKodeGetDky(ark_mem, tstop, 0, yout);
                    if *ier != ARK_SUCCESS {
                        arkProcessError(
                            Some(ark_mem),
                            ARK_ILL_INPUT,
                            line!() as i32,
                            "arkStopTests",
                            file!(),
                            &MSG_ARK_BAD_TSTOP(tstop, tcur),
                        );
                        *ier = ARK_ILL_INPUT;
                        return 1;
                    }
                } else {
                    let yn = ark_mem.borrow().yn.clone().expect("yn allocated");
                    N_VScale(ONE, &yn, yout);
                }
                {
                    let mut m = ark_mem.borrow_mut();
                    m.tretlast = tstop;
                    m.tstopset = SUNFALSE;
                }
                *tret = tstop;
                *ier = ARK_TSTOP_RETURN;
                return 1;
            }
        }
        /* If next step would overtake tstop, adjust stepsize */
        else {
            let mut guard = ark_mem.borrow_mut();
            let m = &mut *guard;
            if (m.tcur + m.hprime - m.tstop) * m.h > ZERO {
                m.hprime = (m.tstop - m.tcur) * (ONE - FOUR * m.uround);
                m.eta = m.hprime / m.h;
            }
        }
    }

    /* In ARK_NORMAL mode, test if tout was reached */
    let (tcur, h) = {
        let m = ark_mem.borrow();
        (m.tcur, m.h)
    };
    if (itask == ARK_NORMAL) && ((tcur - tout) * h >= ZERO) {
        if ark_mem.borrow().interp.is_some() {
            *ier = ARKodeGetDky(ark_mem, tout, 0, yout);
            if *ier != ARK_SUCCESS {
                arkProcessError(
                    Some(ark_mem),
                    ARK_ILL_INPUT,
                    line!() as i32,
                    "arkStopTests",
                    file!(),
                    &MSG_ARK_BAD_TOUT(tout),
                );
                *ier = ARK_ILL_INPUT;
                return 1;
            }
            ark_mem.borrow_mut().tretlast = tout;
            *tret = tout;
        } else {
            let yn = ark_mem.borrow().yn.clone().expect("yn allocated");
            N_VScale(ONE, &yn, yout);
            ark_mem.borrow_mut().tretlast = tcur;
            *tret = tcur;
        }
        *ier = ARK_SUCCESS;
        return 1;
    }

    /* In ARK_ONE_STEP mode, test if tn was returned */
    let tretlast = ark_mem.borrow().tretlast;
    if itask == ARK_ONE_STEP && SUNRabs(tcur - tretlast) > troundoff {
        let yn = ark_mem.borrow().yn.clone().expect("yn allocated");
        ark_mem.borrow_mut().tretlast = tcur;
        *tret = tcur;
        N_VScale(ONE, &yn, yout);
        *ier = ARK_SUCCESS;
        return 1;
    }

    0
}

/*---------------------------------------------------------------
  arkHin

  This routine computes a tentative initial step size h0.
  Note that here tout is either the value passed to ARKodeEvolve
  at the first call or the value of tstop (if tstop is enabled and
  it is closer to t0=tn than tout). If the RHS function fails
  unrecoverably, arkHin returns ARK_RHSFUNC_FAIL. If the RHS
  function fails recoverably too many times and recovery is not
  possible, arkHin returns ARK_REPTD_RHSFUNC_ERR. Otherwise, arkHin
  sets h to the chosen value h0 and returns ARK_SUCCESS.

  The algorithm used seeks to find h0 as a solution of
  (WRMS norm of (h0^2 ydd / 2)) = 1,
  where ydd = estimated second derivative of y.

  We start with an initial estimate equal to the geometric mean
  of the lower and upper bounds on the step size.

  Loop up to H0_ITERS times to find h0.
  Stop if new and previous values differ by a factor < 2.
  Stop if hnew/hg > 2 after one iteration, as this probably
  means that the ydd value is bad because of cancellation error.

  For each new proposed hg, we allow H0_ITERS attempts to
  resolve a possible recoverable failure from f() by reducing
  the proposed stepsize by a factor of 0.2. If a legal stepsize
  still cannot be found, fall back on a previous value if
  possible, or else return ARK_REPTD_RHSFUNC_ERR.

  Finally, we apply a bias (0.5) and verify that h0 is within
  bounds.
  ---------------------------------------------------------------*/
pub fn arkHin(ark_mem: &ARKodeMem, tout: sunrealtype) -> i32 {
    /* arkInitialSetup checks for tdiff = 0 or < 2 * troundoff */
    let (tcur, uround) = {
        let m = ark_mem.borrow();
        (m.tcur, m.uround)
    };
    let tdiff = tout - tcur;
    let sign: i32 = if tdiff > ZERO { 1 } else { -1 };
    let tdist = SUNRabs(tdiff);
    let tround = uround * SUNMAX(SUNRabs(tcur), SUNRabs(tout));

    /* call full RHS if needed */
    if !ark_mem.borrow().fn_is_current {
        /* NOTE: The step size (h) is used in setting the tolerance in a potential
        mass matrix solve when computing the full RHS. Before calling arkHin, h
        is set to |tout - tcur| or 1 and so we do not need to guard against
        h == 0 here before calling the full RHS. */
        let (step_fullrhs, tn, yn, fn_) = {
            let m = ark_mem.borrow();
            (
                m.step_fullrhs,
                m.tn,
                m.yn.clone().expect("yn allocated"),
                m.fn_.clone().expect("fn allocated"),
            )
        };
        let step_fullrhs = step_fullrhs.expect("step_fullrhs set");
        let retval = step_fullrhs(ark_mem, tn, &yn, &fn_, ARK_FULLRHS_START);
        if retval != 0 {
            return ARK_RHSFUNC_FAIL;
        }
        ark_mem.borrow_mut().fn_is_current = SUNTRUE;
    }

    /* Set lower and upper bounds on h0, and take geometric mean
    as first trial value.
    Exit with this value if the bounds cross each other. */
    let hlb = H0_LBFACTOR * tround;
    let hub = arkUpperBoundH0(ark_mem, tdist);

    let mut hg = SUNRsqrt(hlb * hub);

    if hub < hlb {
        if sign == -1 {
            ark_mem.borrow_mut().h = -hg;
        } else {
            ark_mem.borrow_mut().h = hg;
        }
        return ARK_SUCCESS;
    }

    /* Outer loop */
    let mut hs = hg; /* safeguard against 'uninitialized variable' warning */
    let mut hnew = ZERO;
    let mut yddnrm = ZERO;
    for count1 in 1..=H0_ITERS {
        /* Attempts to estimate ydd */
        let mut hgOK = SUNFALSE;

        for _count2 in 1..=H0_ITERS {
            let hgs = hg * sign as sunrealtype;
            let retval = arkYddNorm(ark_mem, hgs, &mut yddnrm);
            /* If f() failed unrecoverably, give up */
            if retval < 0 {
                return ARK_RHSFUNC_FAIL;
            }
            /* If successful, we can use ydd */
            if retval == ARK_SUCCESS {
                hgOK = SUNTRUE;
                break;
            }
            /* f() failed recoverably; cut step size and test it again */
            hg *= 0.2;
        }

        /* If f() failed recoverably H0_ITERS times */
        if !hgOK {
            /* Exit if this is the first or second pass. No recovery possible */
            if count1 <= 2 {
                return ARK_REPTD_RHSFUNC_ERR;
            }
            /* We have a fall-back option. The value hs is a previous hnew which
            passed through f(). Use it and break */
            hnew = hs;
            break;
        }

        /* The proposed step size is feasible. Save it. */
        hs = hg;

        /* Propose new step size */
        hnew = if yddnrm * hub * hub > TWO {
            SUNRsqrt(TWO / yddnrm)
        } else {
            SUNRsqrt(hg * hub)
        };

        /* If last pass, stop now with hnew */
        if count1 == H0_ITERS {
            break;
        }

        let hrat = hnew / hg;

        /* Accept hnew if it does not differ from hg by more than a factor of 2 */
        if (hrat > HALF) && (hrat < TWO) {
            break;
        }

        /* After one pass, if ydd seems to be bad, use fall-back value. */
        if (count1 > 1) && (hrat > TWO) {
            hnew = hg;
            break;
        }

        /* Send this value back through f() */
        hg = hnew;
    }

    /* Apply bounds, bias factor, and attach sign */
    let mut h0 = H0_BIAS * hnew;
    if h0 < hlb {
        h0 = hlb;
    }
    if h0 > hub {
        h0 = hub;
    }
    if sign == -1 {
        h0 = -h0;
    }
    ark_mem.borrow_mut().h = h0;

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkUpperBoundH0

  This routine sets an upper bound on abs(h0) based on
  tdist = tn - t0 and the values of y[i]/y'[i].
  ---------------------------------------------------------------*/
pub fn arkUpperBoundH0(ark_mem: &ARKodeMem, tdist: sunrealtype) -> sunrealtype {
    /* Bound based on |y0|/|y0'| -- allow at most an increase of
     * H0_UBFACTOR in y0 (based on a forward Euler step). The weight
     * factor is used as a safeguard against zero components in y0. */
    let (temp1, temp2, yn, fn_) = {
        let m = ark_mem.borrow();
        (
            m.tempv1.clone().expect("tempv1 allocated"),
            m.tempv2.clone().expect("tempv2 allocated"),
            m.yn.clone().expect("yn allocated"),
            m.fn_.clone().expect("fn allocated"),
        )
    };

    N_VAbs(&yn, &temp2);

    /* C: ark_mem->efun(ark_mem->yn, temp1, ark_mem->e_data); return ignored */
    {
        let (efun, user_efun) = {
            let m = ark_mem.borrow();
            (m.efun, m.user_efun)
        };
        let efun = efun.expect("efun set");
        if user_efun {
            let mut data = ark_mem.borrow_mut().user_data.take();
            let _ = efun(&yn, &temp1, &mut data);
            ark_mem.borrow_mut().user_data = data;
        } else {
            let mut data = ark_mem.borrow_mut().e_data.take();
            let _ = efun(&yn, &temp1, &mut data);
            ark_mem.borrow_mut().e_data = data;
        }
    }

    N_VInv(&temp1, &temp1);
    N_VLinearSum(H0_UBFACTOR, &temp2, ONE, &temp1, &temp1);

    N_VAbs(&fn_, &temp2);

    N_VDiv(&temp2, &temp1, &temp1);
    let hub_inv = N_VMaxNorm(&temp1);

    /* bound based on tdist -- allow at most a step of magnitude
     * H0_UBFACTOR * tdist */
    let mut hub = H0_UBFACTOR * tdist;

    /* Use the smaller of the two */
    if hub * hub_inv > ONE {
        hub = ONE / hub_inv;
    }

    hub
}

/*---------------------------------------------------------------
  arkYddNorm

  This routine computes an estimate of the second derivative of y
  using a difference quotient, and returns its WRMS norm.
  ---------------------------------------------------------------*/
pub fn arkYddNorm(ark_mem: &ARKodeMem, hg: sunrealtype, yddnrm: &mut sunrealtype) -> i32 {
    let (fn_, yn, ycur, tempv1, ewt, tcur) = {
        let m = ark_mem.borrow();
        (
            m.fn_.clone().expect("fn allocated"),
            m.yn.clone().expect("yn allocated"),
            m.ycur.clone().expect("ycur set"),
            m.tempv1.clone().expect("tempv1 allocated"),
            m.ewt.clone().expect("ewt allocated"),
            m.tcur,
        )
    };

    /* increment y with a multiple of f */
    N_VLinearSum(hg, &fn_, ONE, &yn, &ycur);

    /* compute y', via the ODE RHS routine */
    let step_fullrhs = ark_mem.borrow().step_fullrhs.expect("step_fullrhs set");
    let retval = step_fullrhs(ark_mem, tcur + hg, &ycur, &tempv1, ARK_FULLRHS_OTHER);
    if retval != 0 {
        return ARK_RHSFUNC_FAIL;
    }

    /* difference new f and original f to estimate y'' */
    N_VLinearSum(ONE / hg, &tempv1, -ONE / hg, &fn_, &tempv1);

    /* reset ycur to equal yn (unnecessary?) */
    N_VScale(ONE, &yn, &ycur);

    /* compute norm of y'' */
    *yddnrm = N_VWrmsNorm(&tempv1, &ewt);

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkCompleteStep

  This routine performs various update operations when the step
  solution is complete.  It is assumed that the timestepper
  module has stored the time-evolved solution in ark_mem->ycur,
  and the step that gave rise to this solution in ark_mem->h.
  We update the current time (tn), the current solution (yn),
  increment the overall step counter nst, record the values hold
  and tnew, allow for user-provided postprocessing, and update
  the interpolation structure.
  ---------------------------------------------------------------*/
pub fn arkCompleteStep(ark_mem: &ARKodeMem, dsm: sunrealtype) -> i32 {
    /* Set current time to the end of the step (in case the last
    stage time does not coincide with the step solution time).
    If tstop is enabled, it is possible for tn + h to be past
    tstop by roundoff, and in that case, we reset tn (after
    incrementing by h) to tstop. */

    /* During long-time integration, roundoff can creep into tcur.
    Compensated summation fixes this but with increased cost, so it is optional. */
    {
        let mut guard = ark_mem.borrow_mut();
        let m = &mut *guard;
        if m.use_compensated_sums {
            sundials_core::sundials_utils::sunCompensatedSum(
                m.tn,
                m.h,
                &mut m.tcur,
                &mut m.terr,
            );
        } else {
            m.tcur = m.tn + m.h;
        }

        if m.tstopset {
            let troundoff = FUZZ_FACTOR * m.uround * (SUNRabs(m.tcur) + SUNRabs(m.h));
            if SUNRabs(m.tcur - m.tstop) <= troundoff {
                m.tcur = m.tstop;
            }
        }

        /* store this step's contribution to accumulated temporal error */
        if m.AccumErrorType != ARK_ACCUMERROR_NONE {
            if m.AccumErrorType == ARK_ACCUMERROR_MAX {
                m.AccumError = SUNMAX(dsm, m.AccumError);
            } else if m.AccumErrorType == ARK_ACCUMERROR_SUM {
                m.AccumError += dsm;
            } else
            /* ARK_ACCUMERROR_AVG */
            {
                m.AccumError += dsm * m.h;
            }
        }
    }

    /* call the user-supplied post-step function (if supplied) */
    let PostStepFn = ark_mem.borrow().PostStepFn;
    if let Some(PostStepFn) = PostStepFn {
        let (tcur, ycur, nst) = {
            let m = ark_mem.borrow();
            (m.tcur, m.ycur.clone().expect("ycur set"), m.nst)
        };
        let mut user_data = ark_mem.borrow_mut().user_data.take();
        let retval = PostStepFn(tcur, &ycur, nst, &mut user_data);
        ark_mem.borrow_mut().user_data = user_data;
        if retval != 0 {
            return ARK_POSTSTEPFN_FAIL;
        }
    }

    /* update interpolation structure

    NOTE: This must be called before updating yn with ycur as the interpolation
    module may need to save tn, yn from the start of this step. */
    let interp = ark_mem.borrow().interp.clone();
    if let Some(interp) = interp.as_ref() {
        let tcur = ark_mem.borrow().tcur;
        let retval = crate::arkode_interp::arkInterpUpdate(ark_mem, interp, tcur);
        if retval != ARK_SUCCESS {
            return retval;
        }
    }

    /* update yn to current solution */
    let (ycur, yn) = {
        let m = ark_mem.borrow();
        (
            m.ycur.clone().expect("ycur set"),
            m.yn.clone().expect("yn allocated"),
        )
    };
    N_VScale(ONE, &ycur, &yn);
    ark_mem.borrow_mut().fn_is_current = SUNFALSE;

    /* Notify time step controller object of successful step */
    let hcontroller = ark_mem
        .borrow()
        .hadapt_mem
        .as_ref()
        .expect("hadapt_mem allocated")
        .hcontroller
        .clone();
    if let Some(hcontroller) = hcontroller.as_ref() {
        let h = ark_mem.borrow().h;
        let retval =
            sundials_core::sundials_adaptcontroller::SUNAdaptController_UpdateH(hcontroller, h, dsm);
        if retval != sundials_core::sundials_errors::SUN_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_CONTROLLER_ERR,
                line!() as i32,
                "arkCompleteStep",
                file!(),
                "Failure updating controller object",
            );
            return ARK_CONTROLLER_ERR;
        }
    }

    /* update scalar quantities */
    {
        let mut guard = ark_mem.borrow_mut();
        let m = &mut *guard;
        m.nst += 1;
        m.checkpoint_step_idx += 1;
        m.hold = m.h;
        m.tn = m.tcur;
        m.hprime = m.h * m.eta;

        /* Reset growth factor for subsequent time step */
        let hadapt_mem = m.hadapt_mem.as_mut().expect("hadapt_mem allocated");
        let growth = hadapt_mem.growth;
        hadapt_mem.etamax = growth;
    }

    /* Turn off flag indicating initial step and first stage */
    {
        let mut m = ark_mem.borrow_mut();
        m.initsetup = SUNFALSE;
        m.firststage = SUNFALSE;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkHandleFailure

  This routine prints error messages for all cases of failure by
  arkHin and ark_step. It returns to ARKODE the value that ARKODE
  is to return to the user.
  ---------------------------------------------------------------*/
pub fn arkHandleFailure(ark_mem: &ARKodeMem, flag: i32) -> i32 {
    let (tcur, h) = {
        let m = ark_mem.borrow();
        (m.tcur, m.h)
    };

    /* Depending on flag, print error message and return error flag */
    match flag {
        ARK_ERR_FAILURE => {
            arkProcessError(
                Some(ark_mem),
                ARK_ERR_FAILURE,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                &MSG_ARK_ERR_FAILS(tcur, h),
            );
        }
        ARK_CONV_FAILURE => {
            arkProcessError(
                Some(ark_mem),
                ARK_CONV_FAILURE,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                &MSG_ARK_CONV_FAILS(tcur, h),
            );
        }
        ARK_LSETUP_FAIL => {
            arkProcessError(
                Some(ark_mem),
                ARK_LSETUP_FAIL,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                &MSG_ARK_SETUP_FAILED(tcur),
            );
        }
        ARK_LSOLVE_FAIL => {
            arkProcessError(
                Some(ark_mem),
                ARK_LSOLVE_FAIL,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                &MSG_ARK_SOLVE_FAILED(tcur),
            );
        }
        ARK_RHSFUNC_FAIL => {
            arkProcessError(
                Some(ark_mem),
                ARK_RHSFUNC_FAIL,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                &MSG_ARK_RHSFUNC_FAILED(tcur),
            );
        }
        ARK_UNREC_RHSFUNC_ERR => {
            arkProcessError(
                Some(ark_mem),
                ARK_UNREC_RHSFUNC_ERR,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                &MSG_ARK_RHSFUNC_UNREC(tcur),
            );
        }
        ARK_REPTD_RHSFUNC_ERR => {
            arkProcessError(
                Some(ark_mem),
                ARK_REPTD_RHSFUNC_ERR,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                &MSG_ARK_RHSFUNC_REPTD(tcur),
            );
        }
        ARK_RTFUNC_FAIL => {
            arkProcessError(
                Some(ark_mem),
                ARK_RTFUNC_FAIL,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                &MSG_ARK_RTFUNC_FAILED(tcur),
            );
        }
        ARK_TOO_CLOSE => {
            arkProcessError(
                Some(ark_mem),
                ARK_TOO_CLOSE,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                MSG_ARK_TOO_CLOSE,
            );
        }
        ARK_CONSTR_FAIL => {
            arkProcessError(
                Some(ark_mem),
                ARK_CONSTR_FAIL,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                &MSG_ARK_FAILED_CONSTR(tcur),
            );
        }
        ARK_MASSSOLVE_FAIL => {
            arkProcessError(
                Some(ark_mem),
                ARK_MASSSOLVE_FAIL,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                MSG_ARK_MASSSOLVE_FAIL,
            );
        }
        ARK_NLS_SETUP_FAIL => {
            arkProcessError(
                Some(ark_mem),
                ARK_NLS_SETUP_FAIL,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                &format!(
                    "At t = {} the nonlinear solver setup failed unrecoverably",
                    sundials_core::sundials_utils::sun_format_g(tcur)
                ),
            );
        }
        ARK_VECTOROP_ERR => {
            arkProcessError(
                Some(ark_mem),
                ARK_VECTOROP_ERR,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                &MSG_ARK_VECTOROP_ERR(tcur),
            );
        }
        ARK_INNERSTEP_FAIL => {
            arkProcessError(
                Some(ark_mem),
                ARK_INNERSTEP_FAIL,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                &MSG_ARK_INNERSTEP_FAILED(tcur),
            );
        }
        ARK_NLS_OP_ERR => {
            arkProcessError(
                Some(ark_mem),
                ARK_NLS_OP_ERR,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                &MSG_ARK_NLS_FAIL(tcur),
            );
        }
        ARK_USER_PREDICT_FAIL => {
            arkProcessError(
                Some(ark_mem),
                ARK_USER_PREDICT_FAIL,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                &MSG_ARK_USER_PREDICT_FAIL(tcur),
            );
        }
        ARK_POSTPROCESS_STEP_FAIL => {
            arkProcessError(
                Some(ark_mem),
                ARK_POSTPROCESS_STEP_FAIL,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                &MSG_ARK_POSTPROCESS_STEP_FAIL(tcur),
            );
        }
        ARK_POSTPROCESS_STAGE_FAIL => {
            arkProcessError(
                Some(ark_mem),
                ARK_POSTPROCESS_STAGE_FAIL,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                &MSG_ARK_POSTPROCESS_STAGE_FAIL(tcur),
            );
        }
        ARK_PRESTEPFN_FAIL => {
            arkProcessError(
                Some(ark_mem),
                ARK_PRESTEPFN_FAIL,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                &MSG_ARK_PRESTEPFN_FAIL(tcur),
            );
        }
        ARK_POSTSTEPFN_FAIL => {
            arkProcessError(
                Some(ark_mem),
                ARK_POSTSTEPFN_FAIL,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                &MSG_ARK_POSTSTEPFN_FAIL(tcur),
            );
        }
        ARK_PRERHSFN_FAIL => {
            arkProcessError(
                Some(ark_mem),
                ARK_PRERHSFN_FAIL,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                &MSG_ARK_PRERHSFN_FAIL(tcur),
            );
        }
        ARK_INTERP_FAIL => {
            arkProcessError(
                Some(ark_mem),
                ARK_INTERP_FAIL,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                &format!(
                    "At t = {} the interpolation module failed unrecoverably",
                    sundials_core::sundials_utils::sun_format_g(tcur)
                ),
            );
        }
        ARK_INVALID_TABLE => {
            arkProcessError(
                Some(ark_mem),
                ARK_INVALID_TABLE,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                "ARKODE was provided an invalid method table",
            );
        }
        ARK_RELAX_FAIL => {
            arkProcessError(
                Some(ark_mem),
                ARK_RELAX_FAIL,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                &format!(
                    "At t = {} the relaxation module failed",
                    sundials_core::sundials_utils::sun_format_g(tcur)
                ),
            );
        }
        ARK_RELAX_MEM_NULL => {
            arkProcessError(
                Some(ark_mem),
                ARK_RELAX_MEM_NULL,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                "The ARKODE relaxation module memory is NULL",
            );
        }
        ARK_RELAX_FUNC_FAIL => {
            arkProcessError(
                Some(ark_mem),
                ARK_RELAX_FUNC_FAIL,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                "The relaxation function failed unrecoverably",
            );
        }
        ARK_RELAX_JAC_FAIL => {
            arkProcessError(
                Some(ark_mem),
                ARK_RELAX_JAC_FAIL,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                "The relaxation Jacobian failed unrecoverably",
            );
        }
        ARK_ADJ_RECOMPUTE_FAIL => {
            arkProcessError(
                Some(ark_mem),
                ARK_ADJ_RECOMPUTE_FAIL,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                "The forward recomputation of step failed unrecoverably",
            );
        }
        ARK_ADJ_CHECKPOINT_FAIL => {
            arkProcessError(
                Some(ark_mem),
                ARK_ADJ_CHECKPOINT_FAIL,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                "A checkpoint operation failed unrecoverably",
            );
        }
        ARK_SUNADJSTEPPER_ERR => {
            arkProcessError(
                Some(ark_mem),
                ARK_SUNADJSTEPPER_ERR,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                "A SUNAdjStepper operation failed unrecoverably",
            );
        }
        ARK_DOMEIG_FAIL => {
            arkProcessError(
                Some(ark_mem),
                ARK_DOMEIG_FAIL,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                "The dominant eigenvalue function failed unrecoverably",
            );
        }
        ARK_MAX_STAGE_LIMIT_FAIL => {
            arkProcessError(
                Some(ark_mem),
                ARK_MAX_STAGE_LIMIT_FAIL,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                "The max stage limit failed unrecoverably",
            );
        }
        ARK_SUNSTEPPER_ERR => {
            arkProcessError(
                Some(ark_mem),
                ARK_SUNSTEPPER_ERR,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                "An inner SUNStepper error occurred",
            );
        }
        _ => {
            /* This return should never happen */
            arkProcessError(
                Some(ark_mem),
                ARK_UNRECOGNIZED_ERROR,
                line!() as i32,
                "arkHandleFailure",
                file!(),
                "ARKODE encountered an unrecognized error. Please report this to the Sundials developers at sundials-users@llnl.gov",
            );
            return ARK_UNRECOGNIZED_ERROR;
        }
    }

    flag
}

/*---------------------------------------------------------------
  arkEwtSetSS

  This routine is responsible for setting the error weight vector
  ewt as follows:

  ewt[i] = 1 / (reltol * SUNRabs(ycur[i]) + abstol), i=0,...,neq-1

  When the absolute tolerance is zero, it tests for non-positive
  components before inverting. arkEwtSetSS returns 0 if ewt is
  successfully set to a positive vector and -1 otherwise. In the
  latter case, ewt is considered undefined.
  ---------------------------------------------------------------*/
pub fn arkEwtSetSS(
    ycur: &N_Vector,
    weight: &N_Vector,
    arkode_mem: &mut Option<Box<dyn std::any::Any>>,
) -> i32 {
    /* arkode_mem points to ark_mem here (a boxed ARKodeMem handle clone;
    C's cast of a NULL/foreign pointer is UB -> deterministic panic) */
    let ark_mem = arkode_mem
        .as_ref()
        .and_then(|b| b.downcast_ref::<ARKodeMem>())
        .cloned()
        .expect("arkEwtSetSS data holds ARKodeMem");

    let (tempv1, reltol, Sabstol, atolmin0) = {
        let m = ark_mem.borrow();
        (
            m.tempv1.clone().expect("tempv1 allocated"),
            m.reltol,
            m.Sabstol,
            m.atolmin0,
        )
    };
    N_VAbs(ycur, &tempv1);
    N_VScale(reltol, &tempv1, &tempv1);
    N_VAddConst(&tempv1, Sabstol, &tempv1);
    if atolmin0 && N_VMin(&tempv1) <= ZERO {
        return -1;
    }
    N_VInv(&tempv1, weight);
    0
}

/*---------------------------------------------------------------
  arkEwtSetSV

  This routine is responsible for setting the error weight vector
  ewt as follows:

  ewt[i] = 1 / (reltol * SUNRabs(ycur[i]) + abstol[i]), i=0,...,neq-1

  When any absolute tolerance is zero, it tests for non-positive
  components before inverting. arkEwtSetSV returns 0 if ewt is
  successfully set to a positive vector and -1 otherwise. In the
  latter case, ewt is considered undefined.
  ---------------------------------------------------------------*/
pub fn arkEwtSetSV(
    ycur: &N_Vector,
    weight: &N_Vector,
    arkode_mem: &mut Option<Box<dyn std::any::Any>>,
) -> i32 {
    let ark_mem = arkode_mem
        .as_ref()
        .and_then(|b| b.downcast_ref::<ARKodeMem>())
        .cloned()
        .expect("arkEwtSetSV data holds ARKodeMem");

    let (tempv1, reltol, Vabstol, atolmin0) = {
        let m = ark_mem.borrow();
        (
            m.tempv1.clone().expect("tempv1 allocated"),
            m.reltol,
            m.Vabstol.clone().expect("Vabstol allocated"),
            m.atolmin0,
        )
    };
    N_VAbs(ycur, &tempv1);
    N_VLinearSum(reltol, &tempv1, ONE, &Vabstol, &tempv1);
    if atolmin0 && N_VMin(&tempv1) <= ZERO {
        return -1;
    }
    N_VInv(&tempv1, weight);
    0
}

/*---------------------------------------------------------------
  arkEwtSetSmallReal

  This routine is responsible for setting the error weight vector
  ewt as follows:

  ewt[i] = SUN_SMALL_REAL

  This is routine is only used with explicit time stepping with
  a fixed step size to avoid a potential too much error return
  to the user.
  ---------------------------------------------------------------*/
pub fn arkEwtSetSmallReal(
    _ycur: &N_Vector, /* SUNDIALS_MAYBE_UNUSED in C */
    weight: &N_Vector,
    _arkode_mem: &mut Option<Box<dyn std::any::Any>>, /* SUNDIALS_MAYBE_UNUSED in C */
) -> i32 {
    N_VConst(SUN_SMALL_REAL, weight);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkRwtSetSS

  This routine sets rwt as described above in the case tol_type = ARK_SS.
  When the absolute tolerance is zero, it tests for non-positive
  components before inverting. arkRwtSetSS returns 0 if rwt is
  successfully set to a positive vector and -1 otherwise. In the
  latter case, rwt is considered undefined.
  ---------------------------------------------------------------*/
pub fn arkRwtSetSS(ark_mem: &ARKodeMem, My: &N_Vector, weight: &N_Vector) -> i32 {
    let (tempv1, reltol, SRabstol, Ratolmin0) = {
        let m = ark_mem.borrow();
        (
            m.tempv1.clone().expect("tempv1 allocated"),
            m.reltol,
            m.SRabstol,
            m.Ratolmin0,
        )
    };
    N_VAbs(My, &tempv1);
    N_VScale(reltol, &tempv1, &tempv1);
    N_VAddConst(&tempv1, SRabstol, &tempv1);
    if Ratolmin0 && N_VMin(&tempv1) <= ZERO {
        return -1;
    }
    N_VInv(&tempv1, weight);
    0
}

/*---------------------------------------------------------------
  arkRwtSetSV

  This routine sets rwt as described above in the case tol_type = ARK_SV.
  When any absolute tolerance is zero, it tests for non-positive
  components before inverting. arkRwtSetSV returns 0 if rwt is
  successfully set to a positive vector and -1 otherwise. In the
  latter case, rwt is considered undefined.
  ---------------------------------------------------------------*/
pub fn arkRwtSetSV(ark_mem: &ARKodeMem, My: &N_Vector, weight: &N_Vector) -> i32 {
    let (tempv1, reltol, VRabstol, Ratolmin0) = {
        let m = ark_mem.borrow();
        (
            m.tempv1.clone().expect("tempv1 allocated"),
            m.reltol,
            m.VRabstol.clone().expect("VRabstol allocated"),
            m.Ratolmin0,
        )
    };
    N_VAbs(My, &tempv1);
    N_VLinearSum(reltol, &tempv1, ONE, &VRabstol, &tempv1);
    if Ratolmin0 && N_VMin(&tempv1) <= ZERO {
        return -1;
    }
    N_VInv(&tempv1, weight);
    0
}

/*---------------------------------------------------------------
  arkPredict_MaximumOrder

  This routine predicts the nonlinear implicit stage solution
  using the ARKode interpolation module.  This uses the
  highest-degree interpolant supported by the module (stored
  in the interpolation module).
  ---------------------------------------------------------------*/
pub fn arkPredict_MaximumOrder(
    ark_mem: &ARKodeMem,
    tau: sunrealtype,
    yguess: &N_Vector,
) -> i32 {
    /* verify that ark_mem and interpolation structure are provided.
    The C `ark_mem == NULL` guard cannot fire: the Rust handle is a
    `&ARKodeMem`. */
    let interp = ark_mem.borrow().interp.clone();
    let interp = match interp {
        None => {
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_NULL,
                line!() as i32,
                "arkPredict_MaximumOrder",
                file!(),
                "ARKodeInterpMem structure is NULL",
            );
            return ARK_MEM_NULL;
        }
        Some(i) => i,
    };

    /* call the interpolation module to do the work */
    crate::arkode_interp::arkInterpEvaluate(
        ark_mem,
        &interp,
        tau,
        0,
        ARK_INTERP_MAX_DEGREE,
        yguess,
    )
}

/*---------------------------------------------------------------
  arkPredict_VariableOrder

  This routine predicts the nonlinear implicit stage solution
  using the ARKODE interpolation module.  The degree of the
  interpolant is based on the level of extrapolation outside the
  preceding time step.
  ---------------------------------------------------------------*/
pub fn arkPredict_VariableOrder(
    ark_mem: &ARKodeMem,
    tau: sunrealtype,
    yguess: &N_Vector,
) -> i32 {
    let ord: i32;
    let tau_tol: sunrealtype = HALF;
    let tau_tol2: sunrealtype = 0.75;

    /* verify that ark_mem and interpolation structure are provided */
    let interp = ark_mem.borrow().interp.clone();
    let interp = match interp {
        None => {
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_NULL,
                line!() as i32,
                "arkPredict_VariableOrder",
                file!(),
                "ARKodeInterpMem structure is NULL",
            );
            return ARK_MEM_NULL;
        }
        Some(i) => i,
    };

    /* set the polynomial order based on tau input */
    if tau <= tau_tol {
        ord = 3;
    } else if tau <= tau_tol2 {
        ord = 2;
    } else {
        ord = 1;
    }

    /* call the interpolation module to do the work */
    crate::arkode_interp::arkInterpEvaluate(ark_mem, &interp, tau, 0, ord, yguess)
}

/*---------------------------------------------------------------
  arkPredict_CutoffOrder

  This routine predicts the nonlinear implicit stage solution
  using the ARKODE interpolation module.  If the level of
  extrapolation is small enough, it uses the maximum degree
  polynomial available (stored in the interpolation module
  structure); otherwise it uses a linear polynomial.
  ---------------------------------------------------------------*/
pub fn arkPredict_CutoffOrder(ark_mem: &ARKodeMem, tau: sunrealtype, yguess: &N_Vector) -> i32 {
    let ord: i32;
    let tau_tol: sunrealtype = HALF;

    /* verify that ark_mem and interpolation structure are provided */
    let interp = ark_mem.borrow().interp.clone();
    let interp = match interp {
        None => {
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_NULL,
                line!() as i32,
                "arkPredict_CutoffOrder",
                file!(),
                "ARKodeInterpMem structure is NULL",
            );
            return ARK_MEM_NULL;
        }
        Some(i) => i,
    };

    /* set the polynomial order based on tau input */
    if tau <= tau_tol {
        ord = ARK_INTERP_MAX_DEGREE;
    } else {
        ord = 1;
    }

    /* call the interpolation module to do the work */
    crate::arkode_interp::arkInterpEvaluate(ark_mem, &interp, tau, 0, ord, yguess)
}

/*---------------------------------------------------------------
  arkPredict_Bootstrap

  This routine predicts the nonlinear implicit stage solution
  using a quadratic Hermite interpolating polynomial, based on
  the data {y_n, f(t_n,y_n), f(t_n+hj,z_j)}.

  Note: we assume that ftemp = f(t_n+hj,z_j) can be computed via
     N_VLinearCombination(nvec, cvals, Xvecs, ftemp),
  i.e. the inputs cvals[0:nvec-1] and Xvecs[0:nvec-1] may be
  combined to form f(t_n+hj,z_j).

  PORT NOTE (call-site contract): C requires the caller's `cvals` and
  `Xvecs` scratch arrays to hold at least `nvec + 2` slots (steppers size
  them with `nfusedopvecs`).  `cvals` is therefore `&mut [sunrealtype]`
  and MUST already be that long.  `Xvecs` is the locked "handle scratch
  rebuilt on demand" `Vec<N_Vector>` (an `N_Vector` array cannot be left
  uninitialized in safe Rust), so it is taken as `&mut Vec<N_Vector>` and
  grown here to `nvec + 2` if the caller pushed only `nvec` handles; every
  slot 0..nvec+2 is written before use, so the filler is never observable.
  The in-place forward shift is transcribed literally, including its
  self-overwriting behavior for `i >= 2` (unreachable: `nvec <= 2` at both
  upstream call sites).
  ---------------------------------------------------------------*/
pub fn arkPredict_Bootstrap(
    ark_mem: &ARKodeMem,
    hj: sunrealtype,
    tau: sunrealtype,
    nvec: i32,
    cvals: &mut [sunrealtype],
    Xvecs: &mut Vec<N_Vector>,
    yguess: &N_Vector,
) -> i32 {
    /* verify that ark_mem and interpolation structure are provided */
    if ark_mem.borrow().interp.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_NULL,
            line!() as i32,
            "arkPredict_Bootstrap",
            file!(),
            "ARKodeInterpMem structure is NULL",
        );
        return ARK_MEM_NULL;
    }

    let (yn, fn_) = {
        let m = ark_mem.borrow();
        (
            m.yn.clone().expect("yn allocated"),
            m.fn_.clone().expect("fn allocated"),
        )
    };

    /* set coefficients for Hermite interpolant */
    let a0 = ONE;
    let a2 = tau * tau / TWO / hj;
    let a1 = tau - a2;

    /* set arrays for fused vector operation; shift inputs for
    f(t_n+hj,z_j) to end of queue */
    let n = nvec as usize;
    if Xvecs.len() < n + 2 {
        Xvecs.resize(n + 2, yn.clone());
    }
    for i in 0..n {
        cvals[2 + i] = a2 * cvals[i];
        Xvecs[2 + i] = Xvecs[i].clone();
    }
    cvals[0] = a0;
    Xvecs[0] = yn;
    cvals[1] = a1;
    Xvecs[1] = fn_;

    /* call fused vector operation to compute prediction */
    let retval = N_VLinearCombination(nvec + 2, cvals, &Xvecs[..], yguess);
    if retval != 0 {
        return ARK_VECTOROP_ERR;
    }
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkCheckConvergence

  This routine checks the return flag from the time-stepper's
  "step" routine for algebraic solver convergence issues.

  Returns ARK_SUCCESS (0) if successful, PREDICT_AGAIN (>0)
  on a recoverable convergence failure, or a relevant
  nonrecoverable failure flag (<0).
  --------------------------------------------------------------*/
pub fn arkCheckConvergence(ark_mem: &ARKodeMem, nflagPtr: &mut i32, ncfPtr: &mut i32) -> i32 {
    /* If nonlinear solver succeeded, return with ARK_SUCCESS */
    if *nflagPtr == ARK_SUCCESS {
        return ARK_SUCCESS;
    }
    /* Returns with an ARK_RETRY_STEP flag occur at a stage well before
    any algebraic solvers are involved. On the other hand,
    the arkCheckConvergence function handles the results from algebraic
    solvers, which never take place with an ARK_RETRY_STEP flag.
    Therefore, we immediately return from arkCheckConvergence,
    as it is irrelevant in the case of an ARK_RETRY_STEP */
    if *nflagPtr == ARK_RETRY_STEP {
        return ARK_RETRY_STEP;
    }

    /* The nonlinear soln. failed; increment ncfn */
    ark_mem.borrow_mut().ncfn += 1;

    /* If fixed time stepping, then return with convergence failure */
    if ark_mem.borrow().fixedstep {
        return ARK_CONV_FAILURE;
    }

    /* Otherwise, access adaptivity structure */
    if ark_mem.borrow().hadapt_mem.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_NULL,
            line!() as i32,
            "arkCheckConvergence",
            file!(),
            MSG_ARKADAPT_NO_MEM,
        );
        return ARK_MEM_NULL;
    }

    /* Return if lsetup, lsolve, or rhs failed unrecoverably */
    if *nflagPtr < 0 {
        if *nflagPtr == ARK_LSETUP_FAIL {
            return ARK_LSETUP_FAIL;
        } else if *nflagPtr == ARK_LSOLVE_FAIL {
            return ARK_LSOLVE_FAIL;
        } else if *nflagPtr == ARK_RHSFUNC_FAIL {
            return ARK_RHSFUNC_FAIL;
        } else {
            return ARK_NLS_OP_ERR;
        }
    }

    /* At this point, nflag = CONV_FAIL or RHSFUNC_RECVR; increment ncf */
    *ncfPtr += 1;
    {
        let mut guard = ark_mem.borrow_mut();
        let m = &mut *guard;
        let hadapt_mem = m.hadapt_mem.as_mut().expect("hadapt_mem allocated");
        hadapt_mem.etamax = ONE;
    }

    /* If we had maxncf failures, or if |h| = hmin,
    return ARK_CONV_FAILURE or ARK_REPTD_RHSFUNC_ERR. */
    let (maxncf, h, hmin) = {
        let m = ark_mem.borrow();
        (m.maxncf, m.h, m.hmin)
    };
    if (*ncfPtr == maxncf) || (SUNRabs(h) <= hmin * ONEPSM) {
        if *nflagPtr == CONV_FAIL {
            return ARK_CONV_FAILURE;
        }
        if *nflagPtr == RHSFUNC_RECVR {
            return ARK_REPTD_RHSFUNC_ERR;
        }
    }

    /* Reduce step size due to convergence failure */
    {
        let mut guard = ark_mem.borrow_mut();
        let m = &mut *guard;
        let etacf = m.hadapt_mem.as_ref().expect("hadapt_mem allocated").etacf;
        m.eta = etacf;
    }

    /* Signal for Jacobian/preconditioner setup */
    *nflagPtr = PREV_CONV_FAIL;

    /* Return to reattempt the step */
    PREDICT_AGAIN
}

/*---------------------------------------------------------------
  arkCheckConstraints

  This routine determines if the constraints of the problem
  are satisfied by the proposed step

  Returns ARK_SUCCESS if successful, otherwise CONSTR_RECVR
  --------------------------------------------------------------*/
pub fn arkCheckConstraints(ark_mem: &ARKodeMem, constrfails: &mut i32, nflag: &mut i32) -> i32 {
    let (mm, tmp, constraints, ycur, yn) = {
        let m = ark_mem.borrow();
        (
            m.tempv4.clone().expect("tempv4 allocated"),
            m.tempv3.clone().expect("tempv3 allocated"),
            m.constraints.clone().expect("constraints set"),
            m.ycur.clone().expect("ycur set"),
            m.yn.clone().expect("yn allocated"),
        )
    };

    /* Check constraints and get mask vector mm for where constraints failed */
    let constraintsPassed = N_VConstrMask(&constraints, &ycur, &mm);
    if constraintsPassed {
        return ARK_SUCCESS;
    }

    /* Constraints not met */

    /* Update total fails and fails in current step */
    ark_mem.borrow_mut().nconstrfails += 1;
    *constrfails += 1;

    /* Return with error if reached max fails in a step */
    if *constrfails == ark_mem.borrow().maxconstrfails {
        return ARK_CONSTR_FAIL;
    }

    /* Return with error if using fixed step sizes */
    if ark_mem.borrow().fixedstep {
        return ARK_CONSTR_FAIL;
    }

    /* Return with error if |h| == hmin */
    let (h, hmin) = {
        let m = ark_mem.borrow();
        (m.h, m.hmin)
    };
    if SUNRabs(h) <= hmin * ONEPSM {
        return ARK_CONSTR_FAIL;
    }

    /* Reduce h by computing eta = h'/h */
    N_VLinearSum(ONE, &yn, -ONE, &ycur, &tmp);
    N_VProd(&mm, &tmp, &tmp);
    let eta = 0.9 * N_VMinQuotient(&yn, &tmp);
    let eta = SUNMAX(eta, TENTH);
    ark_mem.borrow_mut().eta = eta;

    /* Signal for Jacobian/preconditioner setup */
    *nflag = PREV_CONV_FAIL;

    /* Return to reattempt the step */
    CONSTR_RECVR
}

/*---------------------------------------------------------------
  arkCheckTemporalError

  This routine performs the local error test for the method.
  The weighted local error norm dsm is passed in.  This value is
  used to predict the next step to attempt based on dsm.
  The test dsm <= 1 is made, and if this fails then additional
  checks are performed based on the number of successive error
  test failures.

  Returns ARK_SUCCESS if the test passes.

  If the test fails:
    - if maxnef error test failures have occurred or if
      SUNRabs(h) = hmin, we return ARK_ERR_FAILURE.
    - otherwise: set *nflagPtr to PREV_ERR_FAIL, and
      return TRY_AGAIN.
  --------------------------------------------------------------*/
pub fn arkCheckTemporalError(
    ark_mem: &ARKodeMem,
    nflagPtr: &mut i32,
    nefPtr: &mut i32,
    dsm: sunrealtype,
) -> i32 {
    /* Access hadapt_mem structure */
    if ark_mem.borrow().hadapt_mem.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_NULL,
            line!() as i32,
            "arkCheckTemporalError",
            file!(),
            MSG_ARKADAPT_NO_MEM,
        );
        return ARK_MEM_NULL;
    }

    /* consider change of step size for next step attempt (may be
    larger/smaller than current step, depending on dsm) */
    let (tn, h, ycur) = {
        let m = ark_mem.borrow();
        (m.tn, m.h, m.ycur.clone().expect("ycur set"))
    };
    let ttmp = if dsm <= ONE { tn + h } else { tn };
    let retval = crate::arkode_adapt::arkAdapt(ark_mem, &ycur, ttmp, h, dsm);
    if retval != ARK_SUCCESS {
        return ARK_ERR_FAILURE;
    }

    /* if we've made it here then no nonrecoverable failures occurred; someone above
    has recommended an 'eta' value for the next step -- enforce bounds on that value
    and set upcoming step size */
    {
        let mut guard = ark_mem.borrow_mut();
        let m = &mut *guard;
        let etamax = m.hadapt_mem.as_ref().expect("hadapt_mem allocated").etamax;
        m.eta = SUNMIN(m.eta, etamax);
        m.eta = SUNMAX(m.eta, m.hmin / SUNRabs(m.h));
        let denom = SUNMAX(ONE, SUNRabs(m.h) * m.hmax_inv * m.eta);
        m.eta /= denom;
    }

    /* If est. local error norm dsm passes test, return ARK_SUCCESS */
    if dsm <= ONE {
        return ARK_SUCCESS;
    }

    /* Test failed; increment counters, set nflag */
    *nefPtr += 1;
    ark_mem.borrow_mut().netf += 1;
    *nflagPtr = PREV_ERR_FAIL;

    /* At maxnef failures, return ARK_ERR_FAILURE */
    if *nefPtr == ark_mem.borrow().maxnef {
        return ARK_ERR_FAILURE;
    }

    /* Set etamax=1 to prevent step size increase at end of this step */
    {
        let mut guard = ark_mem.borrow_mut();
        let m = &mut *guard;
        let hadapt_mem = m.hadapt_mem.as_mut().expect("hadapt_mem allocated");
        hadapt_mem.etamax = ONE;
    }

    /* Enforce failure bounds on eta */
    {
        let mut guard = ark_mem.borrow_mut();
        let m = &mut *guard;
        let (small_nef, etamxf) = {
            let hadapt_mem = m.hadapt_mem.as_ref().expect("hadapt_mem allocated");
            (hadapt_mem.small_nef, hadapt_mem.etamxf)
        };
        if *nefPtr >= small_nef {
            m.eta = SUNMIN(m.eta, etamxf);
        }

        /* Enforce min/max step bounds once again due to adjustments above */
        let etamax = m.hadapt_mem.as_ref().expect("hadapt_mem allocated").etamax;
        m.eta = SUNMIN(m.eta, etamax);
        m.eta = SUNMAX(m.eta, m.hmin / SUNRabs(m.h));
        let denom = SUNMAX(ONE, SUNRabs(m.h) * m.hmax_inv * m.eta);
        m.eta /= denom;
    }

    TRY_AGAIN
}

/*---------------------------------------------------------------
  arkAllocVec and arkAllocVecArray:

  These routines allocate (respectively) single vector or a vector
  array based on a template vector.  If the target vector or vector
  array already exists it is left alone; otherwise it is allocated
  by cloning the input vector.

  This routine also updates the optional outputs lrw and liw, which
  are (respectively) the lengths of the overall ARKODE real and
  integer work spaces.

  SUNTRUE is returned if the allocation is successful (or if the
  target vector or vector array already exists) otherwise SUNFALSE
  is returned.

  PORT NOTE: C passes `&ark_mem->ewt` etc. -- an interior pointer into
  the mem that these routines also mutate.  Rust call sites must
  `Option::take` the field out of the mem, call, and store the result
  back (the failure path in C leaves `*v == NULL`, so restoring `None`
  is equivalent).
  ---------------------------------------------------------------*/
pub fn arkAllocVec(
    ark_mem: &ARKodeMem,
    tmpl: &N_Vector,
    v: &mut Option<N_Vector>,
) -> sunbooleantype {
    /* return failure if N_VClone or N_VDestroy is not implemented */
    {
        let ops = tmpl.ops.borrow();
        if ops.nvclone.is_none() || ops.nvdestroy.is_none() {
            return SUNFALSE;
        }
    }

    /* allocate the new vector if necessary */
    if v.is_none() {
        *v = N_VClone(tmpl);
        if v.is_none() {
            arkFreeVectors(ark_mem);
            return SUNFALSE;
        } else {
            let mut guard = ark_mem.borrow_mut();
            let m = &mut *guard;
            m.lrw += m.lrw1;
            m.liw += m.liw1;
        }
    }
    SUNTRUE
}

pub fn arkAllocVecArray(
    count: i32,
    tmpl: &N_Vector,
    v: &mut Option<Vec<N_Vector>>,
    lrw1: sunindextype,
    lrw: &mut i64,
    liw1: sunindextype,
    liw: &mut i64,
) -> sunbooleantype {
    /* allocate the new vector array if necessary */
    if v.is_none() {
        *v = N_VCloneVectorArray(count, tmpl);
        if v.is_none() {
            return SUNFALSE;
        }
        *lrw += count as i64 * lrw1;
        *liw += count as i64 * liw1;
    }
    SUNTRUE
}

/*---------------------------------------------------------------
  arkFreeVec and arkFreeVecArray:

  These routines (respectively) free a single vector or a vector
  array. If the target vector or vector array is already NULL it
  is left alone; otherwise it is freed and the optional outputs
  lrw and liw are updated accordingly.
  ---------------------------------------------------------------*/
pub fn arkFreeVec(ark_mem: &ARKodeMem, v: &mut Option<N_Vector>) {
    if v.is_some() {
        N_VDestroy(v.take().expect("vector present"));
        let mut guard = ark_mem.borrow_mut();
        let m = &mut *guard;
        m.lrw -= m.lrw1;
        m.liw -= m.liw1;
    }
}

pub fn arkFreeVecArray(
    count: i32,
    v: &mut Option<Vec<N_Vector>>,
    lrw1: sunindextype,
    lrw: &mut i64,
    liw1: sunindextype,
    liw: &mut i64,
) {
    if v.is_some() {
        N_VDestroyVectorArray(v.take().expect("vector array present"), count);
        *lrw -= count as i64 * lrw1;
        *liw -= count as i64 * liw1;
    }
}

/*---------------------------------------------------------------
  arkResizeVec and arkResizeVecArray:

  This routines (respectively) resize a single vector or a vector
  array based on a template vector. If the ARKVecResizeFn function
  is non-NULL, then it calls that routine to perform the resize;
  otherwise it deallocates and reallocates the target vector or
  vector array based on the template vector. These routines also
  updates the optional outputs lrw and liw, which are
  (respectively) the lengths of the overall ARKODE real and
  integer work spaces.

  SUNTRUE is returned if the resize is successful otherwise
  SUNFALSE is returned.
  ---------------------------------------------------------------*/
pub fn arkResizeVec(
    ark_mem: &ARKodeMem,
    resize: Option<ARKVecResizeFn>,
    resize_data: &mut Option<Box<dyn std::any::Any>>,
    lrw_diff: sunindextype,
    liw_diff: sunindextype,
    tmpl: &N_Vector,
    v: &mut Option<N_Vector>,
) -> sunbooleantype {
    if v.is_some() {
        match resize {
            None => {
                N_VDestroy(v.take().expect("vector present"));
                *v = None;
                *v = N_VClone(tmpl);
                if v.is_none() {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_MEM_FAIL,
                        line!() as i32,
                        "arkResizeVec",
                        file!(),
                        "Unable to clone vector",
                    );
                    return SUNFALSE;
                }
            }
            Some(resize) => {
                let vv = v.as_ref().expect("vector present").clone();
                if resize(&vv, tmpl, resize_data) != 0 {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_MEM_FAIL,
                        line!() as i32,
                        "arkResizeVec",
                        file!(),
                        MSG_ARK_RESIZE_FAIL,
                    );
                    return SUNFALSE;
                }
            }
        }
        let mut guard = ark_mem.borrow_mut();
        let m = &mut *guard;
        m.lrw += lrw_diff;
        m.liw += liw_diff;
    }
    SUNTRUE
}

pub fn arkResizeVecArray(
    resize: Option<ARKVecResizeFn>,
    resize_data: &mut Option<Box<dyn std::any::Any>>,
    count: i32,
    tmpl: &N_Vector,
    v: &mut Option<Vec<N_Vector>>,
    lrw_diff: sunindextype,
    lrw: &mut i64,
    liw_diff: sunindextype,
    liw: &mut i64,
) -> sunbooleantype {
    if v.is_some() {
        match resize {
            None => {
                N_VDestroyVectorArray(v.take().expect("vector array present"), count);
                *v = None;
                *v = N_VCloneVectorArray(count, tmpl);
                if v.is_none() {
                    return SUNFALSE;
                }
            }
            Some(resize) => {
                for i in 0..count as usize {
                    let vi = v.as_ref().expect("vector array present")[i].clone();
                    if resize(&vi, tmpl, resize_data) != 0 {
                        return SUNFALSE;
                    }
                }
            }
        }
        *lrw += count as i64 * lrw_diff;
        *liw += count as i64 * liw_diff;
    }
    SUNTRUE
}

/*---------------------------------------------------------------
  arkAllocVectors:

  This routine allocates the ARKODE vectors ewt, yn, tempv* and
  ftemp. If any of these vectors already exist, they are left
  alone. Otherwise, it will allocate each vector by cloning the
  input vector. This routine also updates the optional outputs
  lrw and liw, which are (respectively) the lengths of the real
  and integer work spaces.

  If all memory allocations are successful, arkAllocVectors
  returns SUNTRUE, otherwise it returns SUNFALSE.
  ---------------------------------------------------------------*/
pub fn arkAllocVectors(ark_mem: &ARKodeMem, tmpl: &N_Vector) -> sunbooleantype {
    /* Allocate ewt if needed */
    let mut v = ark_mem.borrow_mut().ewt.take();
    let ok = arkAllocVec(ark_mem, tmpl, &mut v);
    ark_mem.borrow_mut().ewt = v;
    if !ok {
        return SUNFALSE;
    }

    /* Set rwt to point at ewt */
    if ark_mem.borrow().rwt_is_ewt {
        let ewt = ark_mem.borrow().ewt.clone();
        ark_mem.borrow_mut().rwt = ewt;
    }

    /* Allocate yn if needed */
    let mut v = ark_mem.borrow_mut().yn.take();
    let ok = arkAllocVec(ark_mem, tmpl, &mut v);
    ark_mem.borrow_mut().yn = v;
    if !ok {
        return SUNFALSE;
    }

    /* Allocate tempv1 if needed */
    let mut v = ark_mem.borrow_mut().tempv1.take();
    let ok = arkAllocVec(ark_mem, tmpl, &mut v);
    ark_mem.borrow_mut().tempv1 = v;
    if !ok {
        return SUNFALSE;
    }

    /* Allocate tempv2 if needed */
    let mut v = ark_mem.borrow_mut().tempv2.take();
    let ok = arkAllocVec(ark_mem, tmpl, &mut v);
    ark_mem.borrow_mut().tempv2 = v;
    if !ok {
        return SUNFALSE;
    }

    /* Allocate tempv3 if needed */
    let mut v = ark_mem.borrow_mut().tempv3.take();
    let ok = arkAllocVec(ark_mem, tmpl, &mut v);
    ark_mem.borrow_mut().tempv3 = v;
    if !ok {
        return SUNFALSE;
    }

    /* Allocate tempv4 if needed */
    let mut v = ark_mem.borrow_mut().tempv4.take();
    let ok = arkAllocVec(ark_mem, tmpl, &mut v);
    ark_mem.borrow_mut().tempv4 = v;
    if !ok {
        return SUNFALSE;
    }

    SUNTRUE
}

/*---------------------------------------------------------------
  arkResizeVectors:

  This routine resizes all ARKODE vectors if they exist,
  otherwise they are left alone. If a resize function is provided
  it is called to resize the vectors otherwise the vector is
  freed and a new vector is created by cloning in input vector.
  This routine also updates the optional outputs lrw and liw,
  which are (respectively) the lengths of the real and integer
  work spaces.

  If all memory allocations are successful, arkResizeVectors
  returns SUNTRUE, otherwise it returns SUNFALSE.
  ---------------------------------------------------------------*/
pub fn arkResizeVectors(
    ark_mem: &ARKodeMem,
    resize: Option<ARKVecResizeFn>,
    resize_data: &mut Option<Box<dyn std::any::Any>>,
    lrw_diff: sunindextype,
    liw_diff: sunindextype,
    tmpl: &N_Vector,
) -> sunbooleantype {
    /* Vabstol */
    let mut v = ark_mem.borrow_mut().Vabstol.take();
    let ok = arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, tmpl, &mut v);
    ark_mem.borrow_mut().Vabstol = v;
    if !ok {
        return SUNFALSE;
    }

    /* VRabstol */
    let mut v = ark_mem.borrow_mut().VRabstol.take();
    let ok = arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, tmpl, &mut v);
    ark_mem.borrow_mut().VRabstol = v;
    if !ok {
        return SUNFALSE;
    }

    /* ewt */
    let mut v = ark_mem.borrow_mut().ewt.take();
    let ok = arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, tmpl, &mut v);
    ark_mem.borrow_mut().ewt = v;
    if !ok {
        return SUNFALSE;
    }

    /* rwt  */
    if ark_mem.borrow().rwt_is_ewt {
        /* update pointer to ewt */
        let ewt = ark_mem.borrow().ewt.clone();
        ark_mem.borrow_mut().rwt = ewt;
    } else {
        /* resize if distinct from ewt */
        let mut v = ark_mem.borrow_mut().rwt.take();
        let ok = arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, tmpl, &mut v);
        ark_mem.borrow_mut().rwt = v;
        if !ok {
            return SUNFALSE;
        }
    }

    /* yn */
    let mut v = ark_mem.borrow_mut().yn.take();
    let ok = arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, tmpl, &mut v);
    ark_mem.borrow_mut().yn = v;
    if !ok {
        return SUNFALSE;
    }

    /* fn */
    let mut v = ark_mem.borrow_mut().fn_.take();
    let ok = arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, tmpl, &mut v);
    ark_mem.borrow_mut().fn_ = v;
    if !ok {
        return SUNFALSE;
    }

    /* tempv* */
    let mut v = ark_mem.borrow_mut().tempv1.take();
    let ok = arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, tmpl, &mut v);
    ark_mem.borrow_mut().tempv1 = v;
    if !ok {
        return SUNFALSE;
    }

    let mut v = ark_mem.borrow_mut().tempv2.take();
    let ok = arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, tmpl, &mut v);
    ark_mem.borrow_mut().tempv2 = v;
    if !ok {
        return SUNFALSE;
    }

    let mut v = ark_mem.borrow_mut().tempv3.take();
    let ok = arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, tmpl, &mut v);
    ark_mem.borrow_mut().tempv3 = v;
    if !ok {
        return SUNFALSE;
    }

    let mut v = ark_mem.borrow_mut().tempv4.take();
    let ok = arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, tmpl, &mut v);
    ark_mem.borrow_mut().tempv4 = v;
    if !ok {
        return SUNFALSE;
    }

    let mut v = ark_mem.borrow_mut().tempv5.take();
    let ok = arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, tmpl, &mut v);
    ark_mem.borrow_mut().tempv5 = v;
    if !ok {
        return SUNFALSE;
    }

    SUNTRUE
}

/*---------------------------------------------------------------
  arkFreeVectors

  This routine frees the ARKODE vectors allocated in both
  arkAllocVectors and arkAllocRKVectors.

  PORT NOTE: exactly as in C, `rwt` is NOT cleared when it aliases
  `ewt` (C leaves a dangling alias; the port leaves the extra handle
  clone, which keeps the buffer alive until `rwt` is overwritten -- not
  observable, and `lrw`/`liw` accounting is unchanged).
  ---------------------------------------------------------------*/
pub fn arkFreeVectors(ark_mem: &ARKodeMem) {
    let mut v = ark_mem.borrow_mut().ewt.take();
    arkFreeVec(ark_mem, &mut v);
    ark_mem.borrow_mut().ewt = v;

    if !ark_mem.borrow().rwt_is_ewt {
        let mut v = ark_mem.borrow_mut().rwt.take();
        arkFreeVec(ark_mem, &mut v);
        ark_mem.borrow_mut().rwt = v;
    }

    let mut v = ark_mem.borrow_mut().tempv1.take();
    arkFreeVec(ark_mem, &mut v);
    ark_mem.borrow_mut().tempv1 = v;

    let mut v = ark_mem.borrow_mut().tempv2.take();
    arkFreeVec(ark_mem, &mut v);
    ark_mem.borrow_mut().tempv2 = v;

    let mut v = ark_mem.borrow_mut().tempv3.take();
    arkFreeVec(ark_mem, &mut v);
    ark_mem.borrow_mut().tempv3 = v;

    let mut v = ark_mem.borrow_mut().tempv4.take();
    arkFreeVec(ark_mem, &mut v);
    ark_mem.borrow_mut().tempv4 = v;

    let mut v = ark_mem.borrow_mut().tempv5.take();
    arkFreeVec(ark_mem, &mut v);
    ark_mem.borrow_mut().tempv5 = v;

    let mut v = ark_mem.borrow_mut().yn.take();
    arkFreeVec(ark_mem, &mut v);
    ark_mem.borrow_mut().yn = v;

    let mut v = ark_mem.borrow_mut().fn_.take();
    arkFreeVec(ark_mem, &mut v);
    ark_mem.borrow_mut().fn_ = v;

    let mut v = ark_mem.borrow_mut().Vabstol.take();
    arkFreeVec(ark_mem, &mut v);
    ark_mem.borrow_mut().Vabstol = v;

    let mut v = ark_mem.borrow_mut().constraints.take();
    arkFreeVec(ark_mem, &mut v);
    ark_mem.borrow_mut().constraints = v;
}

/*---------------------------------------------------------------
  arkAccessHAdaptMem:

  Shortcut routine to unpack ark_mem and hadapt_mem structures from
  void* pointer.  If either is missing it returns ARK_MEM_NULL.

  PORT NOTE: `ARKodeHAdaptMem` is a `Box` owned by `ark_mem` and cannot
  be handed out, so -- exactly like `step_getlinmem` in the frozen
  contract -- this becomes a PRESENCE CHECK.  On `ARK_SUCCESS` the
  caller reaches the record through
  `ark_mem.borrow[_mut]().hadapt_mem.as_[mut_]ref().expect(...)`.
  The C `arkode_mem == NULL` branch (the only user of `fname`) cannot
  fire because the Rust handle is a `&ARKodeMem`.
  ---------------------------------------------------------------*/
pub fn arkAccessHAdaptMem(ark_mem: &ARKodeMem, fname: &str) -> i32 {
    let _ = fname; /* used only by C's unreachable NULL-handle branch */
    if ark_mem.borrow().hadapt_mem.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_NULL,
            line!() as i32,
            "arkAccessHAdaptMem",
            file!(),
            MSG_ARKADAPT_NO_MEM,
        );
        return ARK_MEM_NULL;
    }
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  Utility routines for ARKODE to serve as an MRIStepInnerStepper
  ---------------------------------------------------------------*/

/*------------------------------------------------------------------------------
  ark_MRIStepInnerEvolve

  Implementation of MRIStepInnerStepperEvolveFn to advance the inner (fast)
  ODE IVP.  Since the raw return value from an MRIStepInnerStepper is
  meaningless, aside from whether it is 0 (success), >0 (recoverable failure),
  and <0 (unrecoverable failure), we map various ARKODE return values
  accordingly.
  ----------------------------------------------------------------------------*/
pub fn ark_MRIStepInnerEvolve(
    stepper: &crate::arkode_mristep::MRIStepInnerStepper,
    _t0: sunrealtype, /* SUNDIALS_MAYBE_UNUSED in C */
    tout: sunrealtype,
    y: &N_Vector,
) -> i32 {
    /* extract the ARKODE memory struct */
    let mut arkode_mem: Option<ARKodeMem> = None;
    let retval =
        crate::arkode_mristep::MRIStepInnerStepper_GetContentAs::<ARKodeMem>(
            stepper,
            &mut arkode_mem,
        );
    if retval != ARK_SUCCESS {
        return -1;
    }
    let ark_mem = match arkode_mem {
        None => {
            arkProcessError(
                None,
                ARK_MEM_NULL,
                line!() as i32,
                "ark_MRIStepInnerEvolve",
                file!(),
                MSG_ARK_NO_MEM,
            );
            return -1;
        }
        Some(m) => m,
    };

    /* get the forcing data */
    let mut tshift: sunrealtype = ZERO;
    let mut tscale: sunrealtype = ZERO;
    let mut forcing: Vec<N_Vector> = Vec::new();
    let mut nforcing: i32 = 0;
    let retval = crate::arkode_mristep::MRIStepInnerStepper_GetForcingData(
        stepper,
        &mut tshift,
        &mut tscale,
        &mut forcing,
        &mut nforcing,
    );
    if retval != ARK_SUCCESS {
        return -1;
    }

    /* set the inner forcing data */
    let step_setforcing = ark_mem.borrow().step_setforcing.expect("step_setforcing set");
    let retval = step_setforcing(&ark_mem, tshift, tscale, &forcing, nforcing);
    if retval != ARK_SUCCESS {
        return -1;
    }

    /* set the stop time */
    let retval = crate::arkode_io::ARKodeSetStopTime(&ark_mem, tout);
    if retval != ARK_SUCCESS {
        return -1;
    }

    /* evolve inner ODE, consider all positive return values as 'success' */
    let mut tret: sunrealtype = ZERO;
    let mut retval = ARKodeEvolve(&ark_mem, tout, y, &mut tret, ARK_NORMAL);
    if retval > 0 {
        retval = 0;
    }

    /* set a recoverable failure for a few ARKODE failure modes;
    on other ARKODE errors return with an unrecoverable failure */
    if retval < 0 {
        if (retval == ARK_TOO_MUCH_WORK)
            || (retval == ARK_CONV_FAILURE)
            || (retval == ARK_ERR_FAILURE)
        {
            retval = 1;
        } else {
            return -1;
        }
    }

    /* disable inner forcing */
    let step_setforcing = ark_mem.borrow().step_setforcing.expect("step_setforcing set");
    if step_setforcing(&ark_mem, ZERO, ONE, &[], 0) != ARK_SUCCESS {
        return -1;
    }

    retval
}

/*------------------------------------------------------------------------------
  ark_MRIStepInnerFullRhs

  Implementation of MRIStepInnerStepperFullRhsFn to compute the full inner
  (fast) ODE IVP RHS.
  ----------------------------------------------------------------------------*/
pub fn ark_MRIStepInnerFullRhs(
    stepper: &crate::arkode_mristep::MRIStepInnerStepper,
    t: sunrealtype,
    y: &N_Vector,
    f: &N_Vector,
    mode: i32,
) -> i32 {
    let mut arkode_mem: Option<ARKodeMem> = None;
    let retval =
        crate::arkode_mristep::MRIStepInnerStepper_GetContentAs::<ARKodeMem>(
            stepper,
            &mut arkode_mem,
        );
    if retval != ARK_SUCCESS {
        return -1;
    }
    let ark_mem = match arkode_mem {
        None => {
            arkProcessError(
                None,
                ARK_MEM_NULL,
                line!() as i32,
                "ark_MRIStepInnerFullRhs",
                file!(),
                MSG_ARK_NO_MEM,
            );
            return -1;
        }
        Some(m) => m,
    };
    let step_fullrhs = ark_mem.borrow().step_fullrhs.expect("step_fullrhs set");
    let retval = step_fullrhs(&ark_mem, t, y, f, mode);
    if retval == ARK_SUCCESS {
        return 0;
    }
    -1
}

/*------------------------------------------------------------------------------
  ark_MRIStepInnerReset

  Implementation of MRIStepInnerStepperResetFn to reset the inner (fast) stepper
  state.
  ----------------------------------------------------------------------------*/
pub fn ark_MRIStepInnerReset(
    stepper: &crate::arkode_mristep::MRIStepInnerStepper,
    tR: sunrealtype,
    yR: &N_Vector,
) -> i32 {
    let mut arkode_mem: Option<ARKodeMem> = None;
    let retval =
        crate::arkode_mristep::MRIStepInnerStepper_GetContentAs::<ARKodeMem>(
            stepper,
            &mut arkode_mem,
        );
    if retval != ARK_SUCCESS {
        return -1;
    }
    /* C hands the (possibly NULL) void* straight to ARKodeReset, which then
    returns ARK_MEM_NULL -> -1; the Rust handle model reaches the same
    result without the intermediate call. */
    let ark_mem = match arkode_mem {
        None => return -1,
        Some(m) => m,
    };
    let retval = ARKodeReset(&ark_mem, tR, yR);
    if retval == ARK_SUCCESS {
        return 0;
    }
    -1
}

/*------------------------------------------------------------------------------
  ark_MRIStepInnerGetAccumulatedError

  Implementation of MRIStepInnerGetAccumulatedError to retrieve the accumulated
  temporal error estimate from the inner (fast) stepper.
  ----------------------------------------------------------------------------*/
pub fn ark_MRIStepInnerGetAccumulatedError(
    stepper: &crate::arkode_mristep::MRIStepInnerStepper,
    accum_error: &mut sunrealtype,
) -> i32 {
    let mut arkode_mem: Option<ARKodeMem> = None;
    let retval =
        crate::arkode_mristep::MRIStepInnerStepper_GetContentAs::<ARKodeMem>(
            stepper,
            &mut arkode_mem,
        );
    if retval != ARK_SUCCESS {
        return -1;
    }
    let ark_mem = match arkode_mem {
        None => return -1,
        Some(m) => m,
    };
    let retval = crate::arkode_io::ARKodeGetAccumulatedError(&ark_mem, accum_error);
    if retval == ARK_SUCCESS {
        return 0;
    }
    if retval > 0 {
        return 1;
    }
    -1
}

/*------------------------------------------------------------------------------
  ark_MRIStepInnerResetAccumulatedError

  Implementation of MRIStepInnerResetAccumulatedError to reset the accumulated
  temporal error estimator in the inner (fast) stepper.
  ----------------------------------------------------------------------------*/
pub fn ark_MRIStepInnerResetAccumulatedError(
    stepper: &crate::arkode_mristep::MRIStepInnerStepper,
) -> i32 {
    let mut arkode_mem: Option<ARKodeMem> = None;
    let retval =
        crate::arkode_mristep::MRIStepInnerStepper_GetContentAs::<ARKodeMem>(
            stepper,
            &mut arkode_mem,
        );
    if retval != ARK_SUCCESS {
        return -1;
    }
    let ark_mem = match arkode_mem {
        None => return -1,
        Some(m) => m,
    };
    let retval = crate::arkode_io::ARKodeResetAccumulatedError(&ark_mem);
    if retval == ARK_SUCCESS {
        return 0;
    }
    -1
}

/*------------------------------------------------------------------------------
  ark_MRIStepInnerSetRTol

  Implementation of MRIStepInnerSetRTol to set a relative tolerance for the
  upcoming evolution using the inner (fast) stepper.
  ----------------------------------------------------------------------------*/
pub fn ark_MRIStepInnerSetRTol(
    stepper: &crate::arkode_mristep::MRIStepInnerStepper,
    rtol: sunrealtype,
) -> i32 {
    let mut arkode_mem: Option<ARKodeMem> = None;
    let retval =
        crate::arkode_mristep::MRIStepInnerStepper_GetContentAs::<ARKodeMem>(
            stepper,
            &mut arkode_mem,
        );
    if retval != ARK_SUCCESS {
        return -1;
    }
    let ark_mem = match arkode_mem {
        None => {
            arkProcessError(
                None,
                ARK_MEM_NULL,
                line!() as i32,
                "ark_MRIStepInnerSetRTol",
                file!(),
                MSG_ARK_NO_MEM,
            );
            return -1;
        }
        Some(m) => m,
    };
    if rtol > ZERO {
        ark_mem.borrow_mut().reltol = rtol;
        0
    } else {
        -1
    }
}
