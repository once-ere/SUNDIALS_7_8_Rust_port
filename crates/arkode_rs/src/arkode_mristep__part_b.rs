/*---------------------------------------------------------------
  mriStep_TakeStepMERK:

  This routine performs a single MERK step.

  Both the vectors ark_mem->yn and ark_mem->ycur hold the previous
  time-step solution on input, and the vector ark_mem->ycur should
  hold the result of this step on output.

  If timestep adaptivity is enabled, this routine also computes
  the error estimate y-ytilde, where ytilde is the
  embedded solution, and the norm weights come from ark_ewt.
  This estimate is stored in ark_mem->tempv1, in case the calling
  routine wishes to examine the error locations.

  The output variable dsmPtr should contain a scalar-valued
  estimate of the temporal error from this step, ||y-ytilde||_WRMS
  if timestep adaptivity is enabled; otherwise it should be 0.

  The input/output variable nflagPtr is used to gauge convergence
  of any algebraic solvers within the step.  At the start of a new
  time step, this will initially have the value FIRST_CALL.  On
  return from this function, nflagPtr should have a value:
            0 => algebraic solve completed successfully
           >0 => solve did not converge at this step size
                 (but may with a smaller stepsize)
           <0 => solve encountered an unrecoverable failure
  Since the fast-scale evolution could be considered a different
  type of "algebraic solver", we similarly report any fast-scale
  evolution error as a recoverable nflagPtr value.

  The return value from this routine is:
            0 => step completed successfully
           >0 => step encountered recoverable failure;
                 reduce step and retry (if possible)
           <0 => step encountered unrecoverable failure
  ---------------------------------------------------------------*/
pub fn mriStep_TakeStepMERK(
    ark_mem: &ARKodeMem,
    dsmPtr: &mut sunrealtype,
    nflagPtr: &mut i32,
) -> i32 {
    /* access the MRIStep mem structure */
    let retval = mriStep_AccessStepMem(ark_mem, "mriStep_TakeStepMERK");
    if retval != ARK_SUCCESS {
        return retval;
    }

    /* initialize algebraic solver convergence flag to success;
       error estimate to zero */
    *nflagPtr = ARK_SUCCESS;
    *dsmPtr = ZERO;

    /* set N_Vector shortcuts */
    let (ytilde, ytemp) = {
        let m = ark_mem.borrow();
        (
            m.tempv4.clone().expect("tempv4"),
            m.tempv2.clone().expect("tempv2"),
        )
    };

    /* initial time for step */
    let mut t0 = ark_mem.borrow().tn;

    /* initialize the current stage index */
    {
        let mut step_mem = mriStep_mem_mut(ark_mem);
        step_mem.cur_stage = 0;
        step_mem.istage = step_mem.cur_stage;
    }

    /* handles that do not change during the step */
    let stepper = mriStep_mem_mut(ark_mem).stepper.clone().expect("stepper");
    let mric = mriStep_mem_mut(ark_mem).MRIC.clone().expect("MRIC");
    let stages = mriStep_mem_mut(ark_mem).stages;

    /* if MRI adaptivity is enabled: reset fast accumulated error,
       and send appropriate control parameter to the fast integrator */
    let hcontroller = ark_mem
        .borrow()
        .hadapt_mem
        .as_ref()
        .expect("hadapt_mem")
        .hcontroller
        .clone();
    let adapt_type = match hcontroller.as_ref() {
        Some(C) => SUNAdaptController_GetType(C),
        None => SUN_ADAPTCONTROLLER_NONE,
    };
    let mut need_inner_dsm = SUNFALSE;
    if adapt_type == SUN_ADAPTCONTROLLER_MRI_H_TOL {
        need_inner_dsm = SUNTRUE;
        mriStep_mem_mut(ark_mem).inner_dsm = ZERO;
        let retval = mriStepInnerStepper_ResetAccumulatedError(&stepper);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_INNERSTEP_FAIL,
                line!() as i32,
                "mriStep_TakeStepMERK",
                file!(),
                "Unable to reset the inner stepper error estimate",
            );
            return ARK_INNERSTEP_FAIL;
        }
        let inner_rtol_factor = mriStep_mem_mut(ark_mem).inner_rtol_factor;
        let reltol = ark_mem.borrow().reltol;
        let retval = mriStepInnerStepper_SetRTol(&stepper, inner_rtol_factor * reltol);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_INNERSTEP_FAIL,
                line!() as i32,
                "mriStep_TakeStepMERK",
                file!(),
                "Unable to set the inner stepper tolerance",
            );
            return ARK_INNERSTEP_FAIL;
        }
    }

    /* for adaptive computations, reset the inner integrator to the beginning of this step */
    let fixedstep = ark_mem.borrow().fixedstep;
    if !fixedstep {
        let (tcur, ycur) = {
            let m = ark_mem.borrow();
            (m.tcur, m.ycur.clone().expect("ycur"))
        };
        let retval = mriStepInnerStepper_Reset(&stepper, tcur, &ycur);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_INNERSTEP_FAIL,
                line!() as i32,
                "mriStep_TakeStepMERK",
                file!(),
                "Unable to reset the inner stepper",
            );
            return ARK_INNERSTEP_FAIL;
        }
    }

    /* Evaluate the slow RHS function if needed. NOTE: we decide between calling the
       full RHS function (if ark_mem->fn is non-NULL and MRIStep is not an inner
       integrator) versus just updating the stored value of Fse[0]. In either case,
       we use ARK_FULLRHS_START mode because MERK methods do not evaluate Fse at the
       end of the time step (so nothing can be leveraged). */
    let nested_mri = {
        let step_mem = mriStep_mem_mut(ark_mem);
        step_mem.expforcing || step_mem.impforcing
    };
    let fn_is_null = ark_mem.borrow().fn_.is_none();
    let fn_is_current = ark_mem.borrow().fn_is_current;
    if fn_is_null || nested_mri {
        let (tcur, ycur) = {
            let m = ark_mem.borrow();
            (m.tcur, m.ycur.clone().expect("ycur"))
        };
        let retval = mriStep_UpdateF0(ark_mem, tcur, &ycur, ARK_FULLRHS_START);
        if retval != 0 {
            return ARK_RHSFUNC_FAIL;
        }
    } else if !fn_is_null && !fn_is_current {
        let (tcur, ycur, fn_) = {
            let m = ark_mem.borrow();
            (
                m.tcur,
                m.ycur.clone().expect("ycur"),
                m.fn_.clone().expect("fn"),
            )
        };
        let retval = mriStep_FullRHS(ark_mem, tcur, &ycur, &fn_, ARK_FULLRHS_START);
        if retval != 0 {
            return ARK_RHSFUNC_FAIL;
        }
    }
    ark_mem.borrow_mut().fn_is_current = SUNTRUE;

    /* The first stage is the previous time-step solution, so its RHS
       is the [already-computed] slow RHS from the start of the step */

    /* Loop over stage groups */
    let ngroup = mric.borrow().ngroup;
    for ig in 0..ngroup {
        /* Find the lowest stage number in this group. The stages in a group are not
           necessarily in increasing order e.g., in MERK43 stage 3 is before stage 2
           in time. Since all the stages in a group share the same forcing vectors
           and the tables must be lower triangular, only stages up to one less than
           the lowest stage index in the group can be used in the forcing. Using the
           lowest stage number in the group prevents unintentionally including stage
           RHS values that have not been computed yet. */
        let mut lowest_stage: i32;
        {
            let C = mric.borrow();
            lowest_stage = C.group[ig as usize][0];
            for il in 1..C.stages {
                if C.group[ig as usize][il as usize] < 0 {
                    break;
                }
                lowest_stage = SUNMIN(lowest_stage, C.group[ig as usize][il as usize]);
            }
        }

        /* Set up fast RHS for this stage group */
        let (tn, h) = {
            let m = ark_mem.borrow();
            (m.tn, m.h)
        };
        let retval = mriStep_ComputeInnerForcing(ark_mem, lowest_stage, tn, tn + h);
        if retval != ARK_SUCCESS {
            return retval;
        }

        /* Set initial condition for this stage group (all but first group) */
        if ig > 0 {
            let (yn, ycur) = {
                let m = ark_mem.borrow();
                (m.yn.clone().expect("yn"), m.ycur.clone().expect("ycur"))
            };
            N_VScale(ONE, &yn, &ycur);
        }
        t0 = ark_mem.borrow().tn;

        /* Evolve fast IVP over each subinterval in stage group */
        for is in 0..stages {
            /* Get stage index from group; skip to the next group if
               we've reached the end of this one */
            let stage = mric.borrow().group[ig as usize][is as usize];
            {
                let mut step_mem = mriStep_mem_mut(ark_mem);
                step_mem.cur_stage = stage;
                step_mem.istage = step_mem.cur_stage;
            }
            if stage < 0 {
                break;
            }
            let mut nextstage = -1;
            if stage < stages {
                nextstage = mric.borrow().group[ig as usize][(is + 1) as usize];
            }

            /* Determine if this is an "embedding" or "solution" stage */
            let mut embedding = SUNFALSE;
            let mut solution = SUNFALSE;
            let ngroup = mric.borrow().ngroup;
            if ig == ngroup - 2 {
                if (stage >= 0) && (nextstage < 0) {
                    embedding = SUNTRUE;
                }
            }
            if ig == ngroup - 1 {
                if (stage >= 0) && (nextstage < 0) {
                    solution = SUNTRUE;
                }
            }

            /* Skip the embedding if we're using fixed time-stepping and
               temporal error estimation is disabled */
            let (fixedstep, accum_type) = {
                let m = ark_mem.borrow();
                (m.fixedstep, m.AccumErrorType)
            };
            if fixedstep && embedding && (accum_type == ARK_ACCUMERROR_NONE) {
                break;
            }

            /* Set current stage abscissa */
            let cstage = if stage >= stages {
                ONE
            } else {
                mric.borrow().c[stage as usize]
            };

            /* Set desired output time for subinterval */
            let (tn, h) = {
                let m = ark_mem.borrow();
                (m.tn, m.h)
            };
            let tf = tn + cstage * h;

            /* Reset the inner stepper on the first stage within all but the
               first stage group due to "stage-restart" structure */
            if (stage > 1) && (is == 0) {
                let ycur = ark_mem.borrow().ycur.clone().expect("ycur");
                let retval = mriStepInnerStepper_Reset(&stepper, t0, &ycur);
                if retval != ARK_SUCCESS {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_INNERSTEP_FAIL,
                        line!() as i32,
                        "mriStep_TakeStepMERK",
                        file!(),
                        "Unable to reset the inner stepper",
                    );
                    return ARK_INNERSTEP_FAIL;
                }
            }

            /* Evolve fast IVP for this stage, potentially get inner dsm on all
               non-embedding stages */
            let ycur = ark_mem.borrow().ycur.clone().expect("ycur");
            let retval = mriStep_StageERKFast(
                ark_mem,
                t0,
                tf,
                &ycur,
                &ytemp,
                need_inner_dsm && !embedding,
            );
            if retval != ARK_SUCCESS {
                *nflagPtr = CONV_FAIL;
                return retval;
            }

            /* Update "initial time" for next stage in group */
            t0 = tf;

            /* set current stage time for postprocessing and RHS calls */
            ark_mem.borrow_mut().tcur = tf;

            /* apply user-supplied stage postprocessing function (if supplied),
               and reset the inner integrator with the modified stage solution */
            let (PostProcessStageFn, PostProcessStepFn) = {
                let m = ark_mem.borrow();
                (m.PostProcessStageFn, m.PostProcessStepFn)
            };
            if !solution && !embedding && PostProcessStageFn.is_some() {
                let (tcur, ycur) = {
                    let m = ark_mem.borrow();
                    (m.tcur, m.ycur.clone().expect("ycur"))
                };
                let mut user_data = ark_mem.borrow_mut().user_data.take();
                let retval =
                    PostProcessStageFn.expect("PostProcessStageFn")(tcur, &ycur, &mut user_data);
                ark_mem.borrow_mut().user_data = user_data;
                if retval != 0 {
                    return ARK_POSTPROCESS_STAGE_FAIL;
                }

                let retval = mriStepInnerStepper_Reset(&stepper, tcur, &ycur);
                if retval != ARK_SUCCESS {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_INNERSTEP_FAIL,
                        line!() as i32,
                        "mriStep_TakeStepMERK",
                        file!(),
                        "Unable to reset the inner stepper",
                    );
                    return ARK_INNERSTEP_FAIL;
                }
            } else if solution && PostProcessStepFn.is_some() {
                let (tcur, ycur) = {
                    let m = ark_mem.borrow();
                    (m.tcur, m.ycur.clone().expect("ycur"))
                };
                let mut user_data = ark_mem.borrow_mut().user_data.take();
                let retval =
                    PostProcessStepFn.expect("PostProcessStepFn")(tcur, &ycur, &mut user_data);
                ark_mem.borrow_mut().user_data = user_data;
                if retval != 0 {
                    return ARK_POSTPROCESS_STEP_FAIL;
                }

                let retval = mriStepInnerStepper_Reset(&stepper, tcur, &ycur);
                if retval != ARK_SUCCESS {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_INNERSTEP_FAIL,
                        line!() as i32,
                        "mriStep_TakeStepMERK",
                        file!(),
                        "Unable to reset the inner stepper",
                    );
                    return ARK_INNERSTEP_FAIL;
                }
            }

            /* Compute updated slow RHS (except for final solution or embedding) */
            if !solution && !embedding {
                /* call the user-supplied pre-RHS function (if supplied) */
                let PreRhsFn = ark_mem.borrow().PreRhsFn;
                if let Some(PreRhsFn) = PreRhsFn {
                    let (tcur, ycur) = {
                        let m = ark_mem.borrow();
                        (m.tcur, m.ycur.clone().expect("ycur"))
                    };
                    let mut user_data = ark_mem.borrow_mut().user_data.take();
                    let retval = PreRhsFn(tcur, &ycur, &mut user_data);
                    ark_mem.borrow_mut().user_data = user_data;
                    if retval != 0 {
                        return ARK_PRERHSFN_FAIL;
                    }
                }

                /* store explicit slow rhs */
                let (tcur, ycur) = {
                    let m = ark_mem.borrow();
                    (m.tcur, m.ycur.clone().expect("ycur"))
                };
                let (fse, Fse_stage) = {
                    let step_mem = mriStep_mem_mut(ark_mem);
                    (
                        step_mem.fse.expect("fse"),
                        step_mem.Fse[stage as usize].clone(),
                    )
                };
                let mut user_data = ark_mem.borrow_mut().user_data.take();
                let retval = fse(tcur, &ycur, &Fse_stage, &mut user_data);
                ark_mem.borrow_mut().user_data = user_data;
                mriStep_mem_mut(ark_mem).nfse += 1;

                if retval < 0 {
                    return ARK_RHSFUNC_FAIL;
                }
                if retval > 0 {
                    return ARK_UNREC_RHSFUNC_ERR;
                }

                /* Add external forcing to Fse[stage], if applicable */
                let expforcing = mriStep_mem_mut(ark_mem).expforcing;
                if expforcing {
                    let mut cvals: Vec<sunrealtype> = Vec::new();
                    let mut Xvecs: Vec<N_Vector> = Vec::new();
                    cvals.push(ONE);
                    Xvecs.push(Fse_stage.clone());
                    let mut nvec = 1;
                    {
                        let step_mem = mriStep_mem_mut(ark_mem);
                        mriStep_ApplyForcing(
                            &step_mem,
                            tcur,
                            ONE,
                            &mut nvec,
                            &mut cvals,
                            &mut Xvecs,
                        );
                    }
                    N_VLinearCombination(nvec, &cvals, &Xvecs, &Fse_stage);
                }
            }

            /* If this is the embedding stage, archive solution for error estimation */
            if embedding {
                let ycur = ark_mem.borrow().ycur.clone().expect("ycur");
                N_VScale(ONE, &ycur, &ytilde);
            }
        } /* loop over stages */
    } /* loop over stage groups */

    /* if temporal error estimation is enabled: compute estimate via difference between
       step solution and embedding, store in ark_mem->tempv1, and store norm in dsmPtr */
    let (fixedstep, accum_type) = {
        let m = ark_mem.borrow();
        (m.fixedstep, m.AccumErrorType)
    };
    if !fixedstep || (accum_type != ARK_ACCUMERROR_NONE) {
        let (ycur, tempv1, ewt) = {
            let m = ark_mem.borrow();
            (
                m.ycur.clone().expect("ycur"),
                m.tempv1.clone().expect("tempv1"),
                m.ewt.clone().expect("ewt"),
            )
        };
        N_VLinearSum(ONE, &ytilde, -ONE, &ycur, &tempv1);
        *dsmPtr = N_VWrmsNorm(&tempv1, &ewt);
    }

    ARK_SUCCESS
}

/*===============================================================
  Internal utility routines
  ===============================================================*/

/*---------------------------------------------------------------
  mriStep_AccessARKODEStepMem:

  Shortcut routine to unpack ark_mem and step_mem structures from
  void* pointer.  If either is missing it returns ARK_MEM_NULL.

  Port note (frozen seam spec, section 3): handles are never NULL in
  Rust, so the C out-params `ARKodeMem* ark_mem` /
  `ARKodeMRIStepMem* step_mem` disappear and this collapses to the
  step-memory PRESENCE CHECK.  Use `mriStep_mem_mut(ark_mem)` at each
  use site to reach the record itself.
  ---------------------------------------------------------------*/
pub fn mriStep_AccessARKODEStepMem(arkode_mem: &ARKodeMem, fname: &str) -> i32 {
    /* access ARKodeMem structure: `&ARKodeMem` is never NULL */

    /* access ARKodeMRIStepMem structure */
    if arkode_mem.borrow().step_mem.is_none() {
        arkProcessError(
            Some(arkode_mem),
            ARK_MEM_NULL,
            line!() as i32,
            fname,
            file!(),
            MSG_MRISTEP_NO_MEM,
        );
        return ARK_MEM_NULL;
    }
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_AccessStepMem:

  Shortcut routine to unpack step_mem structure from ark_mem.
  If missing it returns ARK_MEM_NULL.

  Port note: presence check only (see mriStep_AccessARKODEStepMem).
  ---------------------------------------------------------------*/
pub fn mriStep_AccessStepMem(ark_mem: &ARKodeMem, fname: &str) -> i32 {
    /* access ARKodeMRIStepMem structure */
    if ark_mem.borrow().step_mem.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_NULL,
            line!() as i32,
            fname,
            file!(),
            MSG_MRISTEP_NO_MEM,
        );
        return ARK_MEM_NULL;
    }
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_SetCoupling

  This routine determines the MRI method to use, based on the
  desired accuracy and fixed/adaptive time stepping choice.
  ---------------------------------------------------------------*/
pub fn mriStep_SetCoupling(ark_mem: &ARKodeMem) -> i32 {
    let mut Cliw: sunindextype = 0;
    let mut Clrw: sunindextype = 0;
    let mut table_id: i32 = ARKODE_MRI_NONE;

    /* access ARKodeMRIStepMem structure */
    if ark_mem.borrow().step_mem.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_NULL,
            line!() as i32,
            "mriStep_SetCoupling",
            file!(),
            MSG_MRISTEP_NO_MEM,
        );
        return ARK_MEM_NULL;
    }

    /* if coupling has already been specified, just return */
    let have_coupling = mriStep_mem_mut(ark_mem).MRIC.is_some();
    if have_coupling {
        return ARK_SUCCESS;
    }

    let (implicit_rhs, explicit_rhs, q) = {
        let step_mem = mriStep_mem_mut(ark_mem);
        (step_mem.implicit_rhs, step_mem.explicit_rhs, step_mem.q)
    };
    let fixedstep = ark_mem.borrow().fixedstep;

    /* select method based on order and type */
    if fixedstep
    /**** fixed-step methods ****/
    {
        if implicit_rhs && explicit_rhs
        /**** ImEx methods ****/
        {
            match q {
                1 => table_id = MRISTEP_DEFAULT_IMEX_SD_1,
                2 => table_id = MRISTEP_DEFAULT_IMEX_SD_2,
                3 => table_id = MRISTEP_DEFAULT_IMEX_SD_3,
                4 => table_id = MRISTEP_DEFAULT_IMEX_SD_4,
                _ => {}
            }
        } else if implicit_rhs
        /**** implicit methods ****/
        {
            match q {
                1 => table_id = MRISTEP_DEFAULT_IMPL_SD_1,
                2 => table_id = MRISTEP_DEFAULT_IMPL_SD_2,
                3 => table_id = MRISTEP_DEFAULT_IMPL_SD_3,
                4 => table_id = MRISTEP_DEFAULT_IMPL_SD_4,
                _ => {}
            }
        } else
        /**** explicit methods ****/
        {
            match q {
                1 => table_id = MRISTEP_DEFAULT_EXPL_1,
                2 => table_id = MRISTEP_DEFAULT_EXPL_2,
                3 => table_id = MRISTEP_DEFAULT_EXPL_3,
                4 => table_id = MRISTEP_DEFAULT_EXPL_4,
                5 => table_id = MRISTEP_DEFAULT_EXPL_5_AD,
                _ => {}
            }
        }
    } else
    /**** adaptive methods ****/
    {
        if implicit_rhs && explicit_rhs
        /**** ImEx methods ****/
        {
            match q {
                2 => table_id = MRISTEP_DEFAULT_IMEX_SD_2_AD,
                3 => table_id = MRISTEP_DEFAULT_IMEX_SD_3_AD,
                4 => table_id = MRISTEP_DEFAULT_IMEX_SD_4_AD,
                _ => {}
            }
        } else if implicit_rhs
        /**** implicit methods ****/
        {
            match q {
                2 => table_id = MRISTEP_DEFAULT_IMPL_SD_2,
                3 => table_id = MRISTEP_DEFAULT_IMPL_SD_3,
                4 => table_id = MRISTEP_DEFAULT_IMPL_SD_4,
                _ => {}
            }
        } else
        /**** explicit methods ****/
        {
            match q {
                2 => table_id = MRISTEP_DEFAULT_EXPL_2_AD,
                3 => table_id = MRISTEP_DEFAULT_EXPL_3_AD,
                4 => table_id = MRISTEP_DEFAULT_EXPL_4_AD,
                5 => table_id = MRISTEP_DEFAULT_EXPL_5_AD,
                _ => {}
            }
        }
    }
    if table_id == ARKODE_MRI_NONE {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!() as i32,
            "mriStep_SetCoupling",
            file!(),
            "No MRI method is available for the requested configuration.",
        );
        return ARK_ILL_INPUT;
    }

    mriStep_mem_mut(ark_mem).MRIC = MRIStepCoupling_LoadTable(table_id);
    let MRIC = mriStep_mem_mut(ark_mem).MRIC.clone();
    if MRIC.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_INVALID_TABLE,
            line!() as i32,
            "mriStep_SetCoupling",
            file!(),
            "An error occurred in constructing coupling table.",
        );
        return ARK_INVALID_TABLE;
    }
    let MRIC = MRIC.expect("MRIC");

    /* note coupling structure space requirements */
    MRIStepCoupling_Space(&MRIC, &mut Cliw, &mut Clrw);
    {
        let mut m = ark_mem.borrow_mut();
        m.liw += Cliw;
        m.lrw += Clrw;
    }

    /* set [redundant] stored values for stage numbers and
       method/embedding orders */
    {
        let C = MRIC.borrow();
        let mut step_mem = mriStep_mem_mut(ark_mem);
        step_mem.stages = C.stages;
        step_mem.q = C.q;
        step_mem.p = C.p;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_CheckCoupling

  This routine runs through the MRI coupling structure to ensure
  that it meets all necessary requirements, including:
    sorted abscissae, with c[0] = 0 and c[end] = 1
    lower-triangular (i.e., ERK or DIRK)
    all DIRK stages are solve-decoupled [temporarily]
    method order q > 0 (all)
    stages > 0 (all)

  Returns ARK_SUCCESS if it passes, ARK_INVALID_TABLE otherwise.
  ---------------------------------------------------------------*/
pub fn mriStep_CheckCoupling(ark_mem: &ARKodeMem) -> i32 {
    let mut okay: sunbooleantype;
    let mut Gabs: sunrealtype;
    let mut Wabs: sunrealtype;
    let tol: sunrealtype = 100.0 * SUN_UNIT_ROUNDOFF;

    /* access ARKodeMRIStepMem structure */
    if ark_mem.borrow().step_mem.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_NULL,
            line!() as i32,
            "mriStep_CheckCoupling",
            file!(),
            MSG_MRISTEP_NO_MEM,
        );
        return ARK_MEM_NULL;
    }

    let (mric, implicit_rhs, explicit_rhs) = {
        let step_mem = mriStep_mem_mut(ark_mem);
        (
            step_mem.MRIC.clone().expect("MRIC"),
            step_mem.implicit_rhs,
            step_mem.explicit_rhs,
        )
    };
    let fixedstep = ark_mem.borrow().fixedstep;

    let C = mric.borrow();

    /* check that stages > 0 */
    if C.stages < 1 {
        arkProcessError(
            Some(ark_mem),
            ARK_INVALID_TABLE,
            line!() as i32,
            "mriStep_CheckCoupling",
            file!(),
            "stages < 1!",
        );
        return ARK_INVALID_TABLE;
    }

    /* check that method order q > 0 */
    if C.q < 1 {
        arkProcessError(
            Some(ark_mem),
            ARK_INVALID_TABLE,
            line!() as i32,
            "mriStep_CheckCoupling",
            file!(),
            "method order < 1",
        );
        return ARK_INVALID_TABLE;
    }

    /* check that embedding order p > 0 (if adaptive) */
    if (C.p < 1) && (!fixedstep) {
        arkProcessError(
            Some(ark_mem),
            ARK_INVALID_TABLE,
            line!() as i32,
            "mriStep_CheckCoupling",
            file!(),
            "embedding order < 1, but ARKodeSetFixedStep was not called",
        );
        return ARK_INVALID_TABLE;
    }

    /* Check that coupling table has compatible type */
    if implicit_rhs && explicit_rhs && (C.type_ != MRISTEP_IMEX) && (C.type_ != MRISTEP_SR) {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!() as i32,
            "mriStep_CheckCoupling",
            file!(),
            "Invalid coupling table for an IMEX problem!",
        );
        return ARK_ILL_INPUT;
    }
    if explicit_rhs
        && (C.type_ != MRISTEP_EXPLICIT)
        && (C.type_ != MRISTEP_IMEX)
        && (C.type_ != MRISTEP_MERK)
        && (C.type_ != MRISTEP_SR)
    {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!() as i32,
            "mriStep_CheckCoupling",
            file!(),
            "Invalid coupling table for an explicit problem!",
        );
        return ARK_ILL_INPUT;
    }
    if implicit_rhs
        && (C.type_ != MRISTEP_IMPLICIT)
        && (C.type_ != MRISTEP_IMEX)
        && (C.type_ != MRISTEP_SR)
    {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!() as i32,
            "mriStep_CheckCoupling",
            file!(),
            "Invalid coupling table for an implicit problem!",
        );
        return ARK_ILL_INPUT;
    }

    /* Check that the matrices are defined appropriately */
    if (C.type_ == MRISTEP_IMEX) || (C.type_ == MRISTEP_SR) {
        /* ImEx */
        if C.W.is_empty() || C.G.is_empty() {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!() as i32,
                "mriStep_CheckCoupling",
                file!(),
                "Invalid coupling table for an IMEX problem!",
            );
            return ARK_ILL_INPUT;
        }
    } else if (C.type_ == MRISTEP_EXPLICIT) || (C.type_ == MRISTEP_MERK) {
        /* Explicit */
        if C.W.is_empty() || !C.G.is_empty() {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!() as i32,
                "mriStep_CheckCoupling",
                file!(),
                "Invalid coupling table for an explicit problem!",
            );
            return ARK_ILL_INPUT;
        }
    } else if C.type_ == MRISTEP_IMPLICIT {
        /* Implicit */
        if !C.W.is_empty() || C.G.is_empty() {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!() as i32,
                "mriStep_CheckCoupling",
                file!(),
                "Invalid coupling table for an implicit problem!",
            );
            return ARK_ILL_INPUT;
        }
    }

    /* Check that W tables are strictly lower triangular */
    if !C.W.is_empty() {
        Wabs = 0.0;
        for k in 0..C.nmat {
            for i in 0..C.stages {
                for j in i..C.stages {
                    Wabs += SUNRabs(C.W[k as usize][i as usize][j as usize]);
                }
            }
        }
        if Wabs > tol {
            arkProcessError(
                Some(ark_mem),
                ARK_INVALID_TABLE,
                line!() as i32,
                "mriStep_CheckCoupling",
                file!(),
                "Coupling can be up to ERK (at most)!",
            );
            return ARK_INVALID_TABLE;
        }
    }

    /* Check that G tables are lower triangular */
    if !C.G.is_empty() {
        Gabs = 0.0;
        for k in 0..C.nmat {
            for i in 0..C.stages {
                for j in (i + 1)..C.stages {
                    Gabs += SUNRabs(C.G[k as usize][i as usize][j as usize]);
                }
            }
        }
        if Gabs > tol {
            arkProcessError(
                Some(ark_mem),
                ARK_INVALID_TABLE,
                line!() as i32,
                "mriStep_CheckCoupling",
                file!(),
                "Coupling can be up to DIRK (at most)!",
            );
            return ARK_INVALID_TABLE;
        }
    }

    /* Check that MERK "groups" are structured appropriately */
    if C.type_ == MRISTEP_MERK {
        let mut group_counter: Vec<i32> = vec![0; (C.stages + 1) as usize];
        for i in 0..C.ngroup {
            for j in 0..C.stages {
                let k = C.group[i as usize][j as usize];
                if k == -1 {
                    break;
                }
                if (k < 0) || (k > C.stages) {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_INVALID_TABLE,
                        line!() as i32,
                        "mriStep_CheckCoupling",
                        file!(),
                        "Invalid MERK group index!",
                    );
                    return ARK_INVALID_TABLE;
                }
                group_counter[k as usize] += 1;
            }
        }
        for i in 1..=C.stages {
            if (group_counter[i as usize] == 0) || (group_counter[i as usize] > 1) {
                arkProcessError(
                    Some(ark_mem),
                    ARK_INVALID_TABLE,
                    line!() as i32,
                    "mriStep_CheckCoupling",
                    file!(),
                    "Duplicated/missing stages from MERK groups!",
                );
                return ARK_INVALID_TABLE;
            }
        }
    }

    /* Check that no stage has MRISTAGE_DIRK_FAST type (for now) */
    let stages = C.stages;
    drop(C);
    okay = SUNTRUE;
    for i in 0..stages {
        if mriStepCoupling_GetStageType(&mric, i) == MRISTAGE_DIRK_FAST {
            okay = SUNFALSE;
        }
    }
    if !okay {
        arkProcessError(
            Some(ark_mem),
            ARK_INVALID_TABLE,
            line!() as i32,
            "mriStep_CheckCoupling",
            file!(),
            "solve-coupled DIRK stages not currently supported",
        );
        return ARK_INVALID_TABLE;
    }
    let C = mric.borrow();

    /* check that MRI-GARK stage times are sorted */
    if (C.type_ == MRISTEP_IMPLICIT) || (C.type_ == MRISTEP_EXPLICIT) || (C.type_ == MRISTEP_IMEX)
    {
        okay = SUNTRUE;
        for i in 1..C.stages {
            if (C.c[i as usize] - C.c[(i - 1) as usize]) < -tol {
                okay = SUNFALSE;
            }
        }
        if !okay {
            arkProcessError(
                Some(ark_mem),
                ARK_INVALID_TABLE,
                line!() as i32,
                "mriStep_CheckCoupling",
                file!(),
                "Stage times must be sorted.",
            );
            return ARK_INVALID_TABLE;
        }
    }

    /* check that the first stage is just the old step solution */
    Gabs = SUNRabs(C.c[0]);
    for k in 0..C.nmat {
        for j in 0..C.stages {
            if !C.W.is_empty() {
                Gabs += SUNRabs(C.W[k as usize][0][j as usize]);
            }
            if !C.G.is_empty() {
                Gabs += SUNRabs(C.G[k as usize][0][j as usize]);
            }
        }
    }
    if Gabs > tol {
        arkProcessError(
            Some(ark_mem),
            ARK_INVALID_TABLE,
            line!() as i32,
            "mriStep_CheckCoupling",
            file!(),
            "First stage must equal old solution.",
        );
        return ARK_INVALID_TABLE;
    }

    /* check that the last stage is at the final time */
    if SUNRabs(ONE - C.c[(C.stages - 1) as usize]) > tol {
        arkProcessError(
            Some(ark_mem),
            ARK_INVALID_TABLE,
            line!() as i32,
            "mriStep_CheckCoupling",
            file!(),
            "Final stage time must be equal 1.",
        );
        return ARK_INVALID_TABLE;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_StageERKFast

  This routine performs a single MRI stage, is, with explicit
  slow time scale and fast time scale that requires evolution.

  On input, ycur is the initial condition for the fast IVP at t0.
  On output, ycur is the solution of the fast IVP at tf.
  The vector ytemp is only used if temporal adaptivity is enabled,
  and the fast error is not provided by the fast integrator.

  get_inner_dsm indicates whether this stage is one that should
  accumulate an inner temporal error estimate.

  Port note: the C `ARKodeMRIStepMem step_mem` parameter is dropped;
  the record is re-acquired granularly through `mriStep_mem_mut`,
  because this routine re-enters `ark_mem` (inner-stepper evolve,
  user callbacks, arkProcessError) and no borrow may be live then.
  ---------------------------------------------------------------*/
pub fn mriStep_StageERKFast(
    ark_mem: &ARKodeMem,
    t0: sunrealtype,
    tf: sunrealtype,
    ycur: &N_Vector,
    ytemp: &N_Vector,
    get_inner_dsm: sunbooleantype,
) -> i32 {
    let _ = ytemp; /* SUNDIALS_MAYBE_UNUSED */

    let stepper = mriStep_mem_mut(ark_mem).stepper.clone().expect("stepper");

    /* pre inner evolve function (if supplied) */
    let pre_inner_evolve = mriStep_mem_mut(ark_mem).pre_inner_evolve;
    if let Some(pre_inner_evolve) = pre_inner_evolve {
        let forcing = stepper.forcing.borrow().clone();
        let nforcing = *stepper.nforcing.borrow();
        let mut user_data = ark_mem.borrow_mut().user_data.take();
        let retval = pre_inner_evolve(t0, &forcing, nforcing, &mut user_data);
        ark_mem.borrow_mut().user_data = user_data;
        if retval != 0 {
            return ARK_OUTERTOINNER_FAIL;
        }
    }

    /* Get the adaptivity type (if applicable) */
    let adapt_type = if get_inner_dsm {
        let hcontroller = ark_mem
            .borrow()
            .hadapt_mem
            .as_ref()
            .expect("hadapt_mem")
            .hcontroller
            .clone();
        match hcontroller.as_ref() {
            Some(C) => SUNAdaptController_GetType(C),
            None => SUN_ADAPTCONTROLLER_NONE,
        }
    } else {
        SUN_ADAPTCONTROLLER_NONE
    };

    /* advance inner method in time */
    let retval = mriStepInnerStepper_Evolve(&stepper, t0, tf, ycur);

    if retval < 0 {
        arkProcessError(
            Some(ark_mem),
            ARK_INNERSTEP_FAIL,
            line!() as i32,
            "mriStep_StageERKFast",
            file!(),
            "Failure when evolving the inner stepper",
        );
        return ARK_INNERSTEP_FAIL;
    }
    if retval > 0 {
        /* increment stepper-specific counter, and decrement ARKODE-level nonlinear
           solver counter (since that will be incremented automatically by ARKODE).
           Return with "TRY_AGAIN" which should cause ARKODE to cut the step size
           and retry the step. */
        mriStep_mem_mut(ark_mem).inner_fails += 1;
        ark_mem.borrow_mut().ncfn -= 1;
        return TRY_AGAIN;
    }

    /* for normal stages (i.e., not the embedding) with MRI adaptivity enabled, get an
       estimate for the fast time scale error */
    if get_inner_dsm {
        /* if the fast integrator uses adaptive steps, retrieve the error estimate */
        if adapt_type == SUN_ADAPTCONTROLLER_MRI_H_TOL {
            /* C passes &step_mem->inner_dsm; mirror the write back into the field */
            let mut inner_dsm: sunrealtype = mriStep_mem_mut(ark_mem).inner_dsm;
            let retval = mriStepInnerStepper_GetAccumulatedError(&stepper, &mut inner_dsm);
            mriStep_mem_mut(ark_mem).inner_dsm = inner_dsm;
            if retval != ARK_SUCCESS {
                arkProcessError(
                    Some(ark_mem),
                    ARK_INNERSTEP_FAIL,
                    line!() as i32,
                    "mriStep_StageERKFast",
                    file!(),
                    "Unable to get accumulated error from the inner stepper",
                );
                return ARK_INNERSTEP_FAIL;
            }

            /* scale the error estimate by 1/rtol to account for different inner/outer tolerances */
            let reltol = ark_mem.borrow().reltol;
            mriStep_mem_mut(ark_mem).inner_dsm /= reltol;
        }
    }

    /* post inner evolve function (if supplied) */
    let post_inner_evolve = mriStep_mem_mut(ark_mem).post_inner_evolve;
    if let Some(post_inner_evolve) = post_inner_evolve {
        let mut user_data = ark_mem.borrow_mut().user_data.take();
        let retval = post_inner_evolve(tf, ycur, &mut user_data);
        ark_mem.borrow_mut().user_data = user_data;
        if retval != 0 {
            return ARK_INNERTOOUTER_FAIL;
        }
    }

    /* return with success */
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_StageERKNoFast

  This routine performs a single MRI stage with explicit slow
  time scale only (no fast time scale evolution).

  Port note: the C `ARKodeMRIStepMem step_mem` parameter is dropped
  (see mriStep_StageERKFast).
  ---------------------------------------------------------------*/
pub fn mriStep_StageERKNoFast(ark_mem: &ARKodeMem, is: i32) -> i32 {
    let mric = mriStep_mem_mut(ark_mem).MRIC.clone().expect("MRIC");

    /* determine effective ERK coefficients (store in Ae_row and Ai_row) */
    let retval = {
        let mut guard = mriStep_mem_mut(ark_mem);
        let step_mem = &mut *guard;
        mriStep_RKCoeffs(
            &mric,
            is,
            &step_mem.stage_map,
            &mut step_mem.Ae_row,
            &mut step_mem.Ai_row,
        )
    };
    if retval != ARK_SUCCESS {
        return retval;
    }

    /* call fused vector operation to perform ERK update -- bound on
       j needs "SUNMIN" to handle the case of an "embedding" stage */
    let (h, ycur) = {
        let m = ark_mem.borrow();
        (m.h, m.ycur.clone().expect("ycur"))
    };
    let mut cvals: Vec<sunrealtype> = Vec::new();
    let mut Xvecs: Vec<N_Vector> = Vec::new();
    cvals.push(ONE);
    Xvecs.push(ycur.clone());
    let mut nvec = 1;
    {
        let step_mem = mriStep_mem_mut(ark_mem);
        for j in 0..SUNMIN(is, step_mem.stages) {
            if step_mem.explicit_rhs && step_mem.stage_map[j as usize] > -1 {
                cvals.push(h * step_mem.Ae_row[step_mem.stage_map[j as usize] as usize]);
                Xvecs.push(step_mem.Fse[step_mem.stage_map[j as usize] as usize].clone());
                nvec += 1;
            }
            if step_mem.implicit_rhs && step_mem.stage_map[j as usize] > -1 {
                cvals.push(h * step_mem.Ai_row[step_mem.stage_map[j as usize] as usize]);
                Xvecs.push(step_mem.Fsi[step_mem.stage_map[j as usize] as usize].clone());
                nvec += 1;
            }
        }
    }
    /* Is there a case where we have an explicit update with Fsi? */

    let retval = N_VLinearCombination(nvec, &cvals, &Xvecs, &ycur);
    if retval != 0 {
        return ARK_VECTOROP_ERR;
    }
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_StageDIRKFast

  This routine performs a single stage of a "solve coupled"
  MRI method, i.e. a stage that is DIRK on the slow time scale
  and involves evolution of the fast time scale, in a
  fully-coupled fashion.

  Port note: the C `ARKodeMRIStepMem step_mem` parameter is dropped
  (see mriStep_StageERKFast).
  ---------------------------------------------------------------*/
pub fn mriStep_StageDIRKFast(ark_mem: &ARKodeMem, is: i32, nflagPtr: &mut i32) -> i32 {
    let _ = is; /* SUNDIALS_MAYBE_UNUSED */
    let _ = nflagPtr; /* SUNDIALS_MAYBE_UNUSED */

    /* this is not currently implemented */
    arkProcessError(
        Some(ark_mem),
        ARK_INVALID_TABLE,
        line!() as i32,
        "mriStep_StageDIRKFast",
        file!(),
        "This routine is not yet implemented.",
    );
    ARK_INVALID_TABLE
}

/*---------------------------------------------------------------
  mriStep_StageDIRKNoFast

  This routine performs a single MRI stage with implicit slow
  time scale only (no fast time scale evolution).

  Port note: the C `ARKodeMRIStepMem step_mem` parameter is dropped
  (see mriStep_StageERKFast).
  ---------------------------------------------------------------*/
pub fn mriStep_StageDIRKNoFast(ark_mem: &ARKodeMem, is: i32, nflagPtr: &mut i32) -> i32 {
    /* store current stage index (for an "embedded" stage, subtract 1) */
    let istage = {
        let mut step_mem = mriStep_mem_mut(ark_mem);
        step_mem.istage = if is == step_mem.stages { is - 1 } else { is };
        step_mem.istage
    };

    /* Call predictor for current stage solution (result placed in zpred) */
    let zpred = mriStep_mem_mut(ark_mem).zpred.clone().expect("zpred");
    let retval = mriStep_Predict(ark_mem, istage, &zpred);
    if retval != ARK_SUCCESS {
        return retval;
    }

    /* If a user-supplied predictor routine is provided, call that here
       Note that mriStep_Predict is *still* called, so this user-supplied
       routine can just "clean up" the built-in prediction, if desired. */
    let stage_predict = mriStep_mem_mut(ark_mem).stage_predict;
    if let Some(stage_predict) = stage_predict {
        let tcur = ark_mem.borrow().tcur;
        let mut user_data = ark_mem.borrow_mut().user_data.take();
        let retval = stage_predict(tcur, &zpred, &mut user_data);
        ark_mem.borrow_mut().user_data = user_data;
        if retval < 0 {
            return ARK_USER_PREDICT_FAIL;
        }
        if retval > 0 {
            return TRY_AGAIN;
        }
    }

    /* determine effective DIRK coefficients (store in cvals) */
    let mric = mriStep_mem_mut(ark_mem).MRIC.clone().expect("MRIC");
    let retval = {
        let mut guard = mriStep_mem_mut(ark_mem);
        let step_mem = &mut *guard;
        mriStep_RKCoeffs(
            &mric,
            is,
            &step_mem.stage_map,
            &mut step_mem.Ae_row,
            &mut step_mem.Ai_row,
        )
    };
    if retval != ARK_SUCCESS {
        return retval;
    }

    /* Set up data for evaluation of DIRK stage residual (data stored in sdata) */
    let retval = mriStep_StageSetup(ark_mem);
    if retval != ARK_SUCCESS {
        return retval;
    }

    /* perform implicit solve (result is stored in ark_mem->ycur); return
       with positive value on anything but success */
    *nflagPtr = mriStep_Nls(ark_mem, *nflagPtr);
    if *nflagPtr != ARK_SUCCESS {
        return TRY_AGAIN;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_ComputeInnerForcing

  Constructs the 'coefficient' vectors for the forcing polynomial
  for a 'fast' outer MRI-GARK stage i:

  p_i(theta) = sum_{k=0}^{n-1} forcing[k] * theta^k

  where theta = (t - t0) / (tf-t0) is the mapped 'time' for
  each 'fast' MRIStep evolution, with:
  * t0 -- the start of this outer MRIStep stage
  * tf-t0, the temporal width of this MRIStep stage
  * n -- shorthand for MRIC->nmat

  Defining cdiff = (tf-t0)/h, explicit and solve-decoupled
  implicit or IMEX MRI-based methods define this forcing polynomial
  for each outer stage i > 0:

  p_i(theta) = w_i,0(theta) * fse_0 + ... + w_i,{i-1}(theta) * fse_{i-1}
             + g_i,0(theta) * fsi_0 + ... + g_i,{i-1}(theta) * fsi_{i-1}

  where

  w_i,j(theta) = w_0,i,j + w_1,i,j * theta + ... + w_n,i,j * theta^{n-1},
  w_k,i,j = 1/cdiff * MRIC->W[k][i][j]

  and

  g_i,j(theta) = g_0,i,j + g_1,i,j * theta + ... + g_n,i,j * theta^{n-1},
  g_k,i,j = 1/cdiff * MRIC->G[k][i][j]

  Converting to the appropriate form, we have

  p_i(theta) = ( w_0,i,0 * fse_0 + ... + w_0,i,{i-1} * fse_{i-1} +
                 g_0,i,0 * fsi_0 + ... + g_0,i,{i-1} * fsi_{i-1} ) * theta^0
             + ( w_1,i,0 * fse_0 + ... + w_1,i,{i-1} * fse_{i-1} +
                 g_1,i,0 * fsi_0 + ... + g_1,i,{i-1} * fsi_{i-1} ) * theta^1
                                    .
                                    .
                                    .
             + ( w_n,i,0 * fse_0 + ... + w_n,i,{i-1} * fse_{i-1} +
                 g_n,i,0 * fsi_0 + ... + g_n,i,{i-1} * fsi_{i-1} ) * theta^{n-1}

  Thus we define the forcing vectors for k = 0,...,nmat - 1

  forcing[k] = w_k,i,0 * fse_0 + ... + w_k,i,{i-1} * fse_{i-1}
             + g_k,i,0 * fsi_0 + ... + g_k,i,{i-1} * fsi_{i-1}

             = 1 / cdiff *
               ( W[k][i][0] * fse_0 + ... + W[k][i][i-1] * fse_{i-1} +
               ( G[k][i][0] * fsi_0 + ... + G[k][i][i-1] * fsi_{i-1} )

  We may use an identical formula for MERK methods, so long as we set t0=tn,
  tf=tn+h, stage_map[j]=j (identity map), and implicit_rhs=SUNFALSE.
  With this configuration: tf-t0=h, theta = (t-tn)/h, and cdiff=1.  MERK methods
  define the forcing polynomial for each outer stage i > 0 as:

  p_i(theta) = w_i,0(theta) * fse_0 + ... + w_i,{i-1}(theta) * fse_{i-1}

  where

  w_i,j(theta) = w_0,i,j + w_1,i,j * theta + ... + w_n,i,j * theta^{n-1},
  w_k,i,j = MRIC->W[k][i][j]

  which is equivalent to the formula above.

  We may use a similar formula for MRISR methods, so long as we set t0=tn,
  tf=tn+h*ci, stage_map[j]=j (identity map), and implicit_rhs=SUNFALSE.
  With this configuration: tf-t0=ci*h, theta = (t-tn)/(ci*h), and cdiff=1/ci.
  MRISR methods define the forcing polynomial for each outer stage i > 0 as:

  p_i(theta) = w_i,0(theta) * fs_0 + ... + w_i,{i-1}(theta) * fs_{i-1}

  where fs_j = fse_j + fsi_j and

  w_i,j(theta) = w_0,i,j + w_1,i,j * theta + ... + w_n,i,j * theta^{n-1},
  w_k,i,j = 1/ci * MRIC->W[k][i][j]

  which is equivalent to the formula above, so long as the stage RHS vectors
  Fse[j] are repurposed to instead store (fse_j + fsi_j).

  This routine additionally returns a success/failure flag:
     ARK_SUCCESS -- successful evaluation

  Port note: the C `ARKodeMRIStepMem step_mem` parameter is dropped
  (see mriStep_StageERKFast); `cvals`/`Xvecs` are function-local
  rebuilds of the step_mem scratch arrays (locked house pattern).
  ---------------------------------------------------------------*/
pub fn mriStep_ComputeInnerForcing(
    ark_mem: &ARKodeMem,
    stage: i32,
    t0: sunrealtype,
    tf: sunrealtype,
) -> i32 {
    let (stepper, mric, mut implicit_rhs, mut explicit_rhs, stages) = {
        let step_mem = mriStep_mem_mut(ark_mem);
        (
            step_mem.stepper.clone().expect("stepper"),
            step_mem.MRIC.clone().expect("MRIC"),
            step_mem.implicit_rhs,
            step_mem.explicit_rhs,
            step_mem.stages,
        )
    };

    /* Set inner forcing time normalization constants */
    *stepper.tshift.borrow_mut() = t0;
    *stepper.tscale.borrow_mut() = tf - t0;

    /* Adjust implicit/explicit RHS flags for MRISR methods, since these
       ignore the G coefficients in the forcing function */
    let is_mrisr = mric.borrow().type_ == MRISTEP_SR;
    if is_mrisr {
        implicit_rhs = SUNFALSE;
        explicit_rhs = SUNTRUE;
    }

    /* compute inner forcing vectors (assumes cdiff != 0) */
    let mut Xvecs: Vec<N_Vector> = Vec::new();
    {
        let step_mem = mriStep_mem_mut(ark_mem);
        for j in 0..SUNMIN(stage, stages) {
            if explicit_rhs && step_mem.stage_map[j as usize] > -1 {
                Xvecs.push(step_mem.Fse[step_mem.stage_map[j as usize] as usize].clone());
            }
            if implicit_rhs && step_mem.stage_map[j as usize] > -1 {
                Xvecs.push(step_mem.Fsi[step_mem.stage_map[j as usize] as usize].clone());
            }
        }
    }

    let nmat = mric.borrow().nmat;
    let rcdiff = ark_mem.borrow().h / (tf - t0);

    for k in 0..nmat {
        let mut cvals: Vec<sunrealtype> = Vec::new();
        let mut nstore = 0;
        {
            let C = mric.borrow();
            let step_mem = mriStep_mem_mut(ark_mem);
            for j in 0..SUNMIN(stage, stages) {
                if step_mem.stage_map[j as usize] > -1 {
                    if explicit_rhs && implicit_rhs {
                        /* ImEx */
                        cvals.push(rcdiff * C.W[k as usize][stage as usize][j as usize]);
                        nstore += 1;
                        cvals.push(rcdiff * C.G[k as usize][stage as usize][j as usize]);
                        nstore += 1;
                    } else if explicit_rhs {
                        /* explicit only */
                        cvals.push(rcdiff * C.W[k as usize][stage as usize][j as usize]);
                        nstore += 1;
                    } else {
                        /* implicit only */
                        cvals.push(rcdiff * C.G[k as usize][stage as usize][j as usize]);
                        nstore += 1;
                    }
                }
            }
        }

        let forcing_k = stepper.forcing.borrow()[k as usize].clone();
        let retval = N_VLinearCombination(nstore, &cvals, &Xvecs, &forcing_k);
        if retval != 0 {
            return ARK_VECTOROP_ERR;
        }
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  Compute/return the effective RK coefficients for a "nofast"
  stage.  We may assume that "A" has already been allocated.
  ---------------------------------------------------------------*/
pub fn mriStep_RKCoeffs(
    MRIC: &MRIStepCoupling,
    is: i32,
    stage_map: &[i32],
    Ae_row: &mut [sunrealtype],
    Ai_row: &mut [sunrealtype],
) -> i32 {
    let C = MRIC.borrow();

    if is < 1 || is > C.stages || stage_map.is_empty() || Ae_row.is_empty() || Ai_row.is_empty() {
        return ARK_INVALID_TABLE;
    }

    /* initialize RK coefficient array */
    for j in 0..C.stages {
        Ae_row[j as usize] = ZERO;
        Ai_row[j as usize] = ZERO;
    }

    /* compute RK coefficients -- note that bounds on j need
       "SUNMIN" to handle the case of an "embedding" stage */
    for k in 0..C.nmat {
        let kconst = ONE / (k as sunrealtype + ONE);
        if !C.W.is_empty() {
            for j in 0..SUNMIN(is, C.stages - 1) {
                if stage_map[j as usize] > -1 {
                    Ae_row[stage_map[j as usize] as usize] +=
                        C.W[k as usize][is as usize][j as usize] * kconst;
                }
            }
        }
        if !C.G.is_empty() {
            for j in 0..=SUNMIN(is, C.stages - 1) {
                if stage_map[j as usize] > -1 {
                    Ai_row[stage_map[j as usize] as usize] +=
                        C.G[k as usize][is as usize][j as usize] * kconst;
                }
            }
        }
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_Predict

  This routine computes the prediction for a specific internal
  stage solution, storing the result in yguess.  The
  prediction is done using the interpolation structure in
  extrapolation mode, hence stages "far" from the previous time
  interval are predicted using lower order polynomials than the
  "nearby" stages.
  ---------------------------------------------------------------*/
pub fn mriStep_Predict(ark_mem: &ARKodeMem, istage: i32, yguess: &N_Vector) -> i32 {
    /* access ARKodeMRIStepMem structure */
    if ark_mem.borrow().step_mem.is_none() {
        arkProcessError(
            None,
            ARK_MEM_NULL,
            line!() as i32,
            "mriStep_Predict",
            file!(),
            MSG_MRISTEP_NO_MEM,
        );
        return ARK_MEM_NULL;
    }

    let predictor = mriStep_mem_mut(ark_mem).predictor;

    /* verify that interpolation structure is provided */
    let no_interp = ark_mem.borrow().interp.is_none();
    if no_interp && (predictor > 0) {
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_NULL,
            line!() as i32,
            "mriStep_Predict",
            file!(),
            "Interpolation structure is NULL",
        );
        return ARK_MEM_NULL;
    }

    /* local shortcuts for use with fused vector operations */
    let mut cvals: Vec<sunrealtype> = Vec::new();
    let mut Xvecs: Vec<N_Vector> = Vec::new();

    /* if the first step (or if resized), use initial condition as guess */
    let initsetup = ark_mem.borrow().initsetup;
    if initsetup {
        let yn = ark_mem.borrow().yn.clone().expect("yn");
        N_VScale(ONE, &yn, yguess);
        return ARK_SUCCESS;
    }

    let mric = mriStep_mem_mut(ark_mem).MRIC.clone().expect("MRIC");

    /* set evaluation time tau as relative shift from previous successful time */
    let mut tau = {
        let m = ark_mem.borrow();
        mric.borrow().c[istage as usize] * m.h / m.hold
    };

    /* use requested predictor formula */
    match predictor {
        1 => {
            /***** Interpolatory Predictor 1 -- all to max order *****/
            let retval = arkPredict_MaximumOrder(ark_mem, tau, yguess);
            if retval != ARK_ILL_INPUT {
                return retval;
            }
        }

        2 => {
            /***** Interpolatory Predictor 2 -- decrease order w/ increasing level of extrapolation *****/
            let retval = arkPredict_VariableOrder(ark_mem, tau, yguess);
            if retval != ARK_ILL_INPUT {
                return retval;
            }
        }

        3 => {
            /***** Cutoff predictor: max order interpolatory output for stages "close"
                   to previous step, first-order predictor for subsequent stages *****/
            let retval = arkPredict_CutoffOrder(ark_mem, tau, yguess);
            if retval != ARK_ILL_INPUT {
                return retval;
            }
        }

        4 => {
            /***** Bootstrap predictor: if any previous stage in step has nonzero c_i,
                   construct a quadratic Hermite interpolant for prediction; otherwise
                   use the trivial predictor.  The actual calculations are performed in
                   arkPredict_Bootstrap, but here we need to determine the appropriate
                   stage, c_j, to use. *****/

            /* determine if any previous stages in step meet criteria */
            let mut jstage: i32 = -1;
            {
                let C = mric.borrow();
                for i in 0..istage {
                    jstage = if C.c[i as usize] != ZERO { i } else { jstage };
                }
            }

            /* if using the trivial predictor, break */
            if jstage != -1 {
                /* find the "optimal" previous stage to use */
                {
                    let C = mric.borrow();
                    let step_mem = mriStep_mem_mut(ark_mem);
                    for i in 0..istage {
                        if (C.c[i as usize] > C.c[jstage as usize])
                            && (C.c[i as usize] != ZERO)
                            && step_mem.stage_map[i as usize] > -1
                        {
                            jstage = i;
                        }
                    }
                }

                /* set stage time, stage RHS and interpolation values */
                let ark_h = ark_mem.borrow().h;
                let h = ark_h * mric.borrow().c[jstage as usize];
                tau = ark_h * mric.borrow().c[istage as usize];
                let mut nvec = 0;
                {
                    let step_mem = mriStep_mem_mut(ark_mem);
                    if step_mem.implicit_rhs {
                        /* Implicit piece */
                        cvals.push(ONE);
                        Xvecs.push(
                            step_mem.Fsi[step_mem.stage_map[jstage as usize] as usize].clone(),
                        );
                        nvec += 1;
                    }
                    if step_mem.explicit_rhs {
                        /* Explicit piece */
                        cvals.push(ONE);
                        Xvecs.push(
                            step_mem.Fse[step_mem.stage_map[jstage as usize] as usize].clone(),
                        );
                        nvec += 1;
                    }
                }

                /* call predictor routine */
                let retval = arkPredict_Bootstrap(ark_mem, h, tau, nvec, &cvals, &Xvecs, yguess);
                if retval != ARK_ILL_INPUT {
                    return retval;
                }
            }
        }

        _ => {}
    }

    /* if we made it here, use the trivial predictor (previous step solution) */
    let yn = ark_mem.borrow().yn.clone().expect("yn");
    N_VScale(ONE, &yn, yguess);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_StageSetup

  This routine sets up the stage data for computing the
  solve-decoupled MRI stage residual, along with the step- and
  method-related factors gamma, gammap and gamrat.

  At the ith stage, we compute the residual vector for
  z=z_i=zp+zc:
    r = z - z_{i-1} - h*sum_{j=0}^{i} A(i,j)*F(z_j)
    r = (zp + zc) - z_{i-1} - h*sum_{j=0}^{i} A(i,j)*F(z_j)
    r = (zc - gamma*F(z)) - data,
  where data = (z_{i-1} - zp + h*sum_{j=0}^{i-1} A(i,j)*F(z_j))
  corresponds to existing information.  This routine computes
  this 'data' vector and stores in step_mem->sdata.

  Note: on input, this row A(i,:) is already stored in rkcoeffs.
  ---------------------------------------------------------------*/
pub fn mriStep_StageSetup(ark_mem: &ARKodeMem) -> i32 {
    /* access ARKodeMRIStepMem structure */
    if ark_mem.borrow().step_mem.is_none() {
        arkProcessError(
            None,
            ARK_MEM_NULL,
            line!() as i32,
            "mriStep_StageSetup",
            file!(),
            MSG_MRISTEP_NO_MEM,
        );
        return ARK_MEM_NULL;
    }

    let (h, firststage, ycur) = {
        let m = ark_mem.borrow();
        (m.h, m.firststage, m.ycur.clone().expect("ycur"))
    };

    /* Set shortcut to current stage index */
    let i = mriStep_mem_mut(ark_mem).istage;

    /* local shortcuts for fused vector operations */
    let mut cvals: Vec<sunrealtype> = Vec::new();
    let mut Xvecs: Vec<N_Vector> = Vec::new();

    let sdata;
    let mut nvec;
    {
        let mut step_mem = mriStep_mem_mut(ark_mem);

        /* Update gamma (if the method contains an implicit component) */
        step_mem.gamma = h * step_mem.Ai_row[step_mem.stage_map[i as usize] as usize];

        if firststage {
            step_mem.gammap = step_mem.gamma;
        }
        step_mem.gamrat = if firststage {
            ONE
        } else {
            step_mem.gamma / step_mem.gammap
        };

        /* set cvals and Xvecs for setting stage data */
        cvals.push(ONE);
        Xvecs.push(ycur.clone());
        cvals.push(-ONE);
        Xvecs.push(step_mem.zpred.clone().expect("zpred"));
        nvec = 2;

        for j in 0..i {
            if step_mem.explicit_rhs && step_mem.stage_map[j as usize] > -1 {
                cvals.push(h * step_mem.Ae_row[step_mem.stage_map[j as usize] as usize]);
                Xvecs.push(step_mem.Fse[step_mem.stage_map[j as usize] as usize].clone());
                nvec += 1;
            }
            if step_mem.implicit_rhs && step_mem.stage_map[j as usize] > -1 {
                cvals.push(h * step_mem.Ai_row[step_mem.stage_map[j as usize] as usize]);
                Xvecs.push(step_mem.Fsi[step_mem.stage_map[j as usize] as usize].clone());
                nvec += 1;
            }
        }

        sdata = step_mem.sdata.clone().expect("sdata");
    }

    /* call fused vector operation to do the work */
    let retval = N_VLinearCombination(nvec, &cvals, &Xvecs, &sdata);
    if retval != 0 {
        return ARK_VECTOROP_ERR;
    }

    /* return with success */
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_SlowRHS:

  Wrapper routine to call the user-supplied slow RHS functions,
  f(t,y) = fse(t,y) + fsi(t,y), with API matching
  ARKTimestepFullRHSFn.  This is only used to determine an
  initial slow time-step size to use when one is not specified
  by the user (i.e., mode should correspond with
  ARK_FULLRHS_START.
  ---------------------------------------------------------------*/
pub fn mriStep_SlowRHS(
    ark_mem: &ARKodeMem,
    t: sunrealtype,
    y: &N_Vector,
    f: &N_Vector,
    mode: i32,
) -> i32 {
    let _ = mode; /* SUNDIALS_MAYBE_UNUSED */

    /* access ARKodeMRIStepMem structure */
    let retval = mriStep_AccessStepMem(ark_mem, "mriStep_SlowRHS");
    if retval != ARK_SUCCESS {
        return retval;
    }

    /* call the user-supplied pre-RHS function (if supplied) */
    let PreRhsFn = ark_mem.borrow().PreRhsFn;
    if let Some(PreRhsFn) = PreRhsFn {
        let mut user_data = ark_mem.borrow_mut().user_data.take();
        let retval = PreRhsFn(t, y, &mut user_data);
        ark_mem.borrow_mut().user_data = user_data;
        if retval != 0 {
            return ARK_PRERHSFN_FAIL;
        }
    }

    let (implicit_rhs, explicit_rhs) = {
        let step_mem = mriStep_mem_mut(ark_mem);
        (step_mem.implicit_rhs, step_mem.explicit_rhs)
    };

    /* call fsi if the problem has an implicit component */
    if implicit_rhs {
        let (fsi, Fsi0) = {
            let step_mem = mriStep_mem_mut(ark_mem);
            (step_mem.fsi.expect("fsi"), step_mem.Fsi[0].clone())
        };
        let mut user_data = ark_mem.borrow_mut().user_data.take();
        let retval = fsi(t, y, &Fsi0, &mut user_data);
        ark_mem.borrow_mut().user_data = user_data;
        {
            let mut step_mem = mriStep_mem_mut(ark_mem);
            step_mem.nfsi += 1;
            step_mem.fsi_is_current = SUNTRUE;
        }
        if retval != 0 {
            arkProcessError(
                Some(ark_mem),
                ARK_RHSFUNC_FAIL,
                line!() as i32,
                "mriStep_SlowRHS",
                file!(),
                &MSG_ARK_RHSFUNC_FAILED(t),
            );
            return ARK_RHSFUNC_FAIL;
        }

        /* Add external forcing, if applicable */
        let impforcing = mriStep_mem_mut(ark_mem).impforcing;
        if impforcing {
            let mut cvals: Vec<sunrealtype> = Vec::new();
            let mut Xvecs: Vec<N_Vector> = Vec::new();
            cvals.push(ONE);
            Xvecs.push(Fsi0.clone());
            let mut nvec = 1;
            {
                let step_mem = mriStep_mem_mut(ark_mem);
                mriStep_ApplyForcing(&step_mem, t, ONE, &mut nvec, &mut cvals, &mut Xvecs);
            }
            N_VLinearCombination(nvec, &cvals, &Xvecs, &Fsi0);
        }
    }

    /* call fse if the problem has an explicit component */
    if explicit_rhs {
        let (fse, Fse0) = {
            let step_mem = mriStep_mem_mut(ark_mem);
            (step_mem.fse.expect("fse"), step_mem.Fse[0].clone())
        };
        let mut user_data = ark_mem.borrow_mut().user_data.take();
        let retval = fse(t, y, &Fse0, &mut user_data);
        ark_mem.borrow_mut().user_data = user_data;
        {
            let mut step_mem = mriStep_mem_mut(ark_mem);
            step_mem.nfse += 1;
            step_mem.fse_is_current = SUNTRUE;
        }
        if retval != 0 {
            arkProcessError(
                Some(ark_mem),
                ARK_RHSFUNC_FAIL,
                line!() as i32,
                "mriStep_SlowRHS",
                file!(),
                &MSG_ARK_RHSFUNC_FAILED(t),
            );
            return ARK_RHSFUNC_FAIL;
        }

        /* Add external forcing, if applicable */
        let expforcing = mriStep_mem_mut(ark_mem).expforcing;
        if expforcing {
            let mut cvals: Vec<sunrealtype> = Vec::new();
            let mut Xvecs: Vec<N_Vector> = Vec::new();
            cvals.push(ONE);
            Xvecs.push(Fse0.clone());
            let mut nvec = 1;
            {
                let step_mem = mriStep_mem_mut(ark_mem);
                mriStep_ApplyForcing(&step_mem, t, ONE, &mut nvec, &mut cvals, &mut Xvecs);
            }
            N_VLinearCombination(nvec, &cvals, &Xvecs, &Fse0);
        }
    }

    /* combine RHS vectors into output */
    if explicit_rhs && implicit_rhs
    /* ImEx */
    {
        let (Fse0, Fsi0) = {
            let step_mem = mriStep_mem_mut(ark_mem);
            (step_mem.Fse[0].clone(), step_mem.Fsi[0].clone())
        };
        N_VLinearSum(ONE, &Fse0, ONE, &Fsi0, f);
    } else if implicit_rhs {
        let Fsi0 = mriStep_mem_mut(ark_mem).Fsi[0].clone();
        N_VScale(ONE, &Fsi0, f);
    } else {
        let Fse0 = mriStep_mem_mut(ark_mem).Fse[0].clone();
        N_VScale(ONE, &Fse0, f);
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_Hin

  This routine computes a tentative initial step size h0.  This
  employs the same safeguards as ARKODE's arkHin utility routine,
  but employs a simpler algorithm that estimates the first step
  such that an explicit Euler step (for only the slow RHS
  routine(s)) would be within user-specified tolerances of the
  initial condition.
  ---------------------------------------------------------------*/
pub fn mriStep_Hin(
    ark_mem: &ARKodeMem,
    tcur: sunrealtype,
    tout: sunrealtype,
    fcur: &N_Vector,
    h: &mut sunrealtype,
) -> i32 {
    /* If tout is too close to tn, give up */
    let tdiff = tout - tcur;
    if tdiff == ZERO {
        return ARK_TOO_CLOSE;
    }
    let sign: i32 = if tdiff > ZERO { 1 } else { -1 };
    let tdist = SUNRabs(tdiff);
    let tround = ark_mem.borrow().uround * SUNMAX(SUNRabs(tcur), SUNRabs(tout));
    if tdist < TWO * tround {
        return ARK_TOO_CLOSE;
    }

    /* h0 should bound the change due to a forward Euler step, and
       include safeguard against "too-small" ||f(t0,y0)||: */
    let ewt = ark_mem.borrow().ewt.clone().expect("ewt");
    let fnorm = N_VWrmsNorm(fcur, &ewt) / H0_BIAS;
    let h0_inv = SUNMAX(ONE / H0_UBFACTOR / tdist, fnorm);
    *h = (sign as sunrealtype) / h0_inv;
    ARK_SUCCESS
}

/*===============================================================
  User-callable functions for a custom inner integrator
  ===============================================================*/

/// C `MRIStepInnerStepper_Create(SUNContext sunctx, MRIStepInnerStepper* stepper)`.
///
/// The C `!sunctx` guard is unreachable through `&SUNContext`.
pub fn MRIStepInnerStepper_Create(
    sunctx: &SUNContext,
    stepper: &mut Option<MRIStepInnerStepper>,
) -> i32 {
    *stepper = None;

    /* malloc + memset(0) of the record and of its ops table */
    *stepper = Some(Rc::new(_MRIStepInnerStepper {
        content: RefCell::new(None),
        python: RefCell::new(None),
        ops: RefCell::new(MRIStepInnerStepper_Ops::default()),
        sunctx: sunctx.clone(),
        forcing: RefCell::new(Vec::new()),
        nforcing: RefCell::new(0),
        nforcing_allocated: RefCell::new(0),
        /* initialize stepper data */
        last_flag: RefCell::new(ARK_SUCCESS),
        tshift: RefCell::new(ZERO),
        tscale: RefCell::new(ZERO),
        vals: RefCell::new(Vec::new()),
        vecs: RefCell::new(Vec::new()),
        lrw1: RefCell::new(0),
        liw1: RefCell::new(0),
        lrw: RefCell::new(0),
        liw: RefCell::new(0),
    }));

    ARK_SUCCESS
}

pub fn MRIStepInnerStepper_CreateFromSUNStepper(
    sunstepper: &SUNStepper,
    stepper: &mut Option<MRIStepInnerStepper>,
) -> i32 {
    let sunctx = sunstepper.sunctx.borrow().clone();
    let retval = MRIStepInnerStepper_Create(&sunctx, stepper);
    if retval != ARK_SUCCESS {
        return retval;
    }

    let this = stepper.clone().expect("stepper");

    let retval = MRIStepInnerStepper_SetContent(&this, Some(Box::new(sunstepper.clone())));
    if retval != ARK_SUCCESS {
        return retval;
    }

    let retval = MRIStepInnerStepper_SetEvolveFn(&this, Some(mriStepInnerStepper_EvolveSUNStepper));
    if retval != ARK_SUCCESS {
        return retval;
    }

    let retval =
        MRIStepInnerStepper_SetFullRhsFn(&this, Some(mriStepInnerStepper_FullRhsSUNStepper));
    if retval != ARK_SUCCESS {
        return retval;
    }

    let retval = MRIStepInnerStepper_SetResetFn(&this, Some(mriStepInnerStepper_ResetSUNStepper));
    if retval != ARK_SUCCESS {
        return retval;
    }

    ARK_SUCCESS
}

/// C `MRIStepInnerStepper_Free(MRIStepInnerStepper* stepper)`.
///
/// Dropping the `Rc` replaces C's `free(ops)` / `free(*stepper)`; storage
/// survives while other clones of the handle do (C would leave those
/// dangling).
pub fn MRIStepInnerStepper_Free(stepper: &mut Option<MRIStepInnerStepper>) -> i32 {
    if stepper.is_none() {
        return ARK_SUCCESS;
    }

    {
        let this = stepper.as_ref().expect("stepper");

        /* free the inner forcing and fused op workspace vector */
        mriStepInnerStepper_FreeVecs(this);

        /* free operations structure: released together with the handle */

        /* free python data (SUNDIALS_ENABLE_PYTHON not built) */
        *this.python.borrow_mut() = None;
    }

    /* free inner stepper mem */
    *stepper = None;

    ARK_SUCCESS
}

/// C `MRIStepInnerStepper_SetContent(stepper, void* content)`; `None` is
/// C's `NULL`.
pub fn MRIStepInnerStepper_SetContent(
    stepper: &MRIStepInnerStepper,
    content: Option<Box<dyn Any>>,
) -> i32 {
    /* C `stepper == NULL` guard is unreachable through `&MRIStepInnerStepper` */
    *stepper.content.borrow_mut() = content;

    ARK_SUCCESS
}

/// C `MRIStepInnerStepper_GetContent(stepper, void** content)`.
///
/// A safe-Rust `Box<dyn Any>` token cannot be aliased, so the stored box is
/// SWAPPED with `content` (deviation class 6, as `SUNStepper_GetContent`):
/// the caller MUST hand it back on every return path before anything else
/// touches the stepper's content. Implementation modules should instead use
/// [`MRIStepInnerStepper_GetContentAs`], which clones the handle exactly as
/// C's pointer copy does.
pub fn MRIStepInnerStepper_GetContent(
    stepper: &MRIStepInnerStepper,
    content: &mut Option<Box<dyn Any>>,
) -> i32 {
    std::mem::swap(&mut *stepper.content.borrow_mut(), content);

    ARK_SUCCESS
}

/// Port-only, borrow-safe companion to [`MRIStepInnerStepper_GetContent`]
/// for the common case where the C `void* content` is a SUNDIALS handle
/// (the ARKODE case: `ARKodeMem`; the SUNStepper case: `SUNStepper`). The
/// stepper keeps its content; nothing has to be handed back. A content type
/// mismatch is C UB (a bad cast) and panics here (deviation class 5).
pub fn MRIStepInnerStepper_GetContentAs<T: Any + Clone>(
    stepper: &MRIStepInnerStepper,
    content: &mut Option<T>,
) -> i32 {
    *content = Some(
        stepper
            .content
            .borrow()
            .as_ref()
            .expect("MRIStepInnerStepper content")
            .downcast_ref::<T>()
            .expect("MRIStepInnerStepper content")
            .clone(),
    );

    ARK_SUCCESS
}

pub fn MRIStepInnerStepper_SetEvolveFn(
    stepper: &MRIStepInnerStepper,
    fn_: Option<MRIStepInnerEvolveFn>,
) -> i32 {
    /* C `stepper == NULL` / `stepper->ops == NULL` guards are unreachable */
    stepper.ops.borrow_mut().evolve = fn_;

    ARK_SUCCESS
}

pub fn MRIStepInnerStepper_SetFullRhsFn(
    stepper: &MRIStepInnerStepper,
    fn_: Option<MRIStepInnerFullRhsFn>,
) -> i32 {
    stepper.ops.borrow_mut().fullrhs = fn_;

    ARK_SUCCESS
}

pub fn MRIStepInnerStepper_SetResetFn(
    stepper: &MRIStepInnerStepper,
    fn_: Option<MRIStepInnerResetFn>,
) -> i32 {
    stepper.ops.borrow_mut().reset = fn_;

    ARK_SUCCESS
}

pub fn MRIStepInnerStepper_SetAccumulatedErrorGetFn(
    stepper: &MRIStepInnerStepper,
    fn_: Option<MRIStepInnerGetAccumulatedError>,
) -> i32 {
    stepper.ops.borrow_mut().geterror = fn_;

    ARK_SUCCESS
}

pub fn MRIStepInnerStepper_SetAccumulatedErrorResetFn(
    stepper: &MRIStepInnerStepper,
    fn_: Option<MRIStepInnerResetAccumulatedError>,
) -> i32 {
    stepper.ops.borrow_mut().reseterror = fn_;

    ARK_SUCCESS
}

pub fn MRIStepInnerStepper_SetRTolFn(
    stepper: &MRIStepInnerStepper,
    fn_: Option<MRIStepInnerSetRTol>,
) -> i32 {
    stepper.ops.borrow_mut().setrtol = fn_;

    ARK_SUCCESS
}

pub fn MRIStepInnerStepper_AddForcing(
    stepper: &MRIStepInnerStepper,
    t: sunrealtype,
    f: &N_Vector,
) -> i32 {
    /* C `stepper == NULL` guard is unreachable through the handle type */

    /* `vals`/`vecs` are rebuilt as locals here (an N_Vector array cannot be
    left uninitialised in safe Rust); the values, the `nvec` argument and
    therefore the arithmetic are identical to C's in-place scratch. */
    let mut vals: Vec<sunrealtype> = Vec::new();
    let mut vecs: Vec<N_Vector> = Vec::new();

    /* always append the constant forcing term */
    vals.push(ONE);
    vecs.push(f.clone());

    /* compute normalized time tau and initialize tau^i */
    let tau = (t - *stepper.tshift.borrow()) / (*stepper.tscale.borrow());
    let mut taui = ONE;

    let nforcing = *stepper.nforcing.borrow();
    for i in 0..nforcing {
        vals.push(taui);
        vecs.push(stepper.forcing.borrow()[i as usize].clone());
        taui *= tau;
    }

    N_VLinearCombination(nforcing + 1, &vals, &vecs, f);

    ARK_SUCCESS
}

/// C `MRIStepInnerStepper_GetForcingData(stepper, tshift, tscale, N_Vector**
/// forcing, nforcing)`. The C out-param hands back the internal array
/// pointer; the port hands back a `Vec` of clones of the same `N_Vector`
/// handles (C pointer copies), so the vectors themselves still alias.
pub fn MRIStepInnerStepper_GetForcingData(
    stepper: &MRIStepInnerStepper,
    tshift: &mut sunrealtype,
    tscale: &mut sunrealtype,
    forcing: &mut Vec<N_Vector>,
    nforcing: &mut i32,
) -> i32 {
    *tshift = *stepper.tshift.borrow();
    *tscale = *stepper.tscale.borrow();
    *forcing = stepper.forcing.borrow().clone();
    *nforcing = *stepper.nforcing.borrow();

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  Internal inner integrator functions
  ---------------------------------------------------------------*/

/* Check for required operations */
pub fn mriStepInnerStepper_HasRequiredOps(stepper: &MRIStepInnerStepper) -> i32 {
    /* C NULL guards on `stepper` and `stepper->ops` are unreachable */

    if stepper.ops.borrow().evolve.is_some() {
        ARK_SUCCESS
    } else {
        ARK_ILL_INPUT
    }
}

/* Check whether stepper supports fast/slow tolerance adaptivity */
pub fn mriStepInnerStepper_SupportsRTolAdaptivity(stepper: &MRIStepInnerStepper) -> sunbooleantype {
    let ops = stepper.ops.borrow();
    if ops.geterror.is_some() && ops.reseterror.is_some() && ops.setrtol.is_some() {
        SUNTRUE
    } else {
        SUNFALSE
    }
}

/* Evolve the inner (fast) ODE */
pub fn mriStepInnerStepper_Evolve(
    stepper: &MRIStepInnerStepper,
    t0: sunrealtype,
    tout: sunrealtype,
    y: &N_Vector,
) -> i32 {
    let evolve = stepper.ops.borrow().evolve;
    if evolve.is_none() {
        return ARK_ILL_INPUT;
    }
    let evolve = evolve.expect("evolve");

    let last_flag = evolve(stepper, t0, tout, y);
    *stepper.last_flag.borrow_mut() = last_flag;

    last_flag
}

pub fn mriStepInnerStepper_EvolveSUNStepper(
    stepper: &MRIStepInnerStepper,
    t0: sunrealtype,
    tout: sunrealtype,
    y: &N_Vector,
) -> i32 {
    let _ = t0; /* SUNDIALS_MAYBE_UNUSED */

    let sunstepper: SUNStepper = stepper
        .content
        .borrow()
        .as_ref()
        .expect("SUNStepper content")
        .downcast_ref::<SUNStepper>()
        .expect("SUNStepper content")
        .clone();
    let mut tret: sunrealtype = ZERO;

    let (tshift, tscale, forcing, nforcing) = (
        *stepper.tshift.borrow(),
        *stepper.tscale.borrow(),
        stepper.forcing.borrow().clone(),
        *stepper.nforcing.borrow(),
    );
    let mut err = SUNStepper_SetForcing(&sunstepper, tshift, tscale, &forcing, nforcing);
    *stepper.last_flag.borrow_mut() = *sunstepper.last_flag.borrow();
    if err != SUN_SUCCESS {
        return ARK_SUNSTEPPER_ERR;
    }

    err = SUNStepper_SetStopTime(&sunstepper, tout);
    *stepper.last_flag.borrow_mut() = *sunstepper.last_flag.borrow();
    if err != SUN_SUCCESS {
        return ARK_SUNSTEPPER_ERR;
    }

    err = SUNStepper_Evolve(&sunstepper, tout, y, &mut tret);
    *stepper.last_flag.borrow_mut() = *sunstepper.last_flag.borrow();
    if err != SUN_SUCCESS {
        return ARK_SUNSTEPPER_ERR;
    }

    err = SUNStepper_SetForcing(&sunstepper, ZERO, ONE, &[], 0);
    *stepper.last_flag.borrow_mut() = *sunstepper.last_flag.borrow();
    if err != SUN_SUCCESS {
        return ARK_SUNSTEPPER_ERR;
    }

    ARK_SUCCESS
}

/* Compute the full RHS for inner (fast) time scale TODO(DJG): This function can
   be made optional when fullrhs is not called unconditionally by the ARKODE
   infrastructure e.g., in arkInitialSetup, arkYddNorm, and arkCompleteStep. */
pub fn mriStepInnerStepper_FullRhs(
    stepper: &MRIStepInnerStepper,
    t: sunrealtype,
    y: &N_Vector,
    f: &N_Vector,
    mode: i32,
) -> i32 {
    let fullrhs = stepper.ops.borrow().fullrhs;
    if fullrhs.is_none() {
        return ARK_ILL_INPUT;
    }
    let fullrhs = fullrhs.expect("fullrhs");

    let last_flag = fullrhs(stepper, t, y, f, mode);
    *stepper.last_flag.borrow_mut() = last_flag;
    last_flag
}

pub fn mriStepInnerStepper_FullRhsSUNStepper(
    stepper: &MRIStepInnerStepper,
    t: sunrealtype,
    y: &N_Vector,
    f: &N_Vector,
    ark_mode: i32,
) -> i32 {
    let sunstepper: SUNStepper = stepper
        .content
        .borrow()
        .as_ref()
        .expect("SUNStepper content")
        .downcast_ref::<SUNStepper>()
        .expect("SUNStepper content")
        .clone();

    let mode: SUNFullRhsMode = match ark_mode {
        ARK_FULLRHS_START => SUN_FULLRHS_START,
        ARK_FULLRHS_END => SUN_FULLRHS_END,
        _ => SUN_FULLRHS_OTHER,
    };

    let err = SUNStepper_FullRhs(&sunstepper, t, y, f, mode);
    *stepper.last_flag.borrow_mut() = *sunstepper.last_flag.borrow();
    if err != SUN_SUCCESS {
        return ARK_SUNSTEPPER_ERR;
    }
    ARK_SUCCESS
}

/* Reset the inner (fast) stepper state */
pub fn mriStepInnerStepper_Reset(
    stepper: &MRIStepInnerStepper,
    tR: sunrealtype,
    yR: &N_Vector,
) -> i32 {
    let reset = stepper.ops.borrow().reset;

    if let Some(reset) = reset {
        let last_flag = reset(stepper, tR, yR);
        *stepper.last_flag.borrow_mut() = last_flag;
        last_flag
    } else {
        /* assume stepper uses input state and does not need to be reset */
        ARK_SUCCESS
    }
}

/* Gets the inner (fast) stepper accumulated error */
pub fn mriStepInnerStepper_GetAccumulatedError(
    stepper: &MRIStepInnerStepper,
    accum_error: &mut sunrealtype,
) -> i32 {
    let geterror = stepper.ops.borrow().geterror;

    if let Some(geterror) = geterror {
        let last_flag = geterror(stepper, accum_error);
        *stepper.last_flag.borrow_mut() = last_flag;
        last_flag
    } else {
        ARK_INNERSTEP_FAIL
    }
}

/* Resets the inner (fast) stepper accumulated error */
pub fn mriStepInnerStepper_ResetAccumulatedError(stepper: &MRIStepInnerStepper) -> i32 {
    /* NOTE: upstream tests `ops->geterror` here but calls `ops->reseterror`;
    the quirk is preserved (a set geterror with an unset reseterror is a NULL
    call in C and a panic here -- deviation class 5). */
    let (geterror, reseterror) = {
        let ops = stepper.ops.borrow();
        (ops.geterror, ops.reseterror)
    };

    if geterror.is_some() {
        let last_flag = reseterror.expect("reseterror")(stepper);
        *stepper.last_flag.borrow_mut() = last_flag;
        last_flag
    } else {
        /* assume stepper provides exact solution and needs no reset */
        ARK_SUCCESS
    }
}

/* Sets the inner (fast) stepper relative tolerance scaling factor */
pub fn mriStepInnerStepper_SetRTol(stepper: &MRIStepInnerStepper, rtol: sunrealtype) -> i32 {
    let setrtol = stepper.ops.borrow().setrtol;

    if let Some(setrtol) = setrtol {
        let last_flag = setrtol(stepper, rtol);
        *stepper.last_flag.borrow_mut() = last_flag;
        last_flag
    } else {
        /* assume stepper provides exact solution */
        ARK_SUCCESS
    }
}

pub fn mriStepInnerStepper_ResetSUNStepper(
    stepper: &MRIStepInnerStepper,
    tR: sunrealtype,
    yR: &N_Vector,
) -> i32 {
    let sunstepper: SUNStepper = stepper
        .content
        .borrow()
        .as_ref()
        .expect("SUNStepper content")
        .downcast_ref::<SUNStepper>()
        .expect("SUNStepper content")
        .clone();
    let err = SUNStepper_Reset(&sunstepper, tR, yR);
    *stepper.last_flag.borrow_mut() = *sunstepper.last_flag.borrow();
    if err != SUN_SUCCESS {
        return ARK_SUNSTEPPER_ERR;
    }
    ARK_SUCCESS
}

/* Allocate MRI forcing and fused op workspace vectors if necessary */
pub fn mriStepInnerStepper_AllocVecs(
    stepper: &MRIStepInnerStepper,
    count: i32,
    tmpl: &N_Vector,
) -> i32 {
    let mut lrw1: sunindextype = 0;
    let mut liw1: sunindextype = 0;

    /* Set space requirements for one N_Vector */
    let has_nvspace = tmpl.ops.borrow().nvspace.is_some();
    if has_nvspace {
        N_VSpace(tmpl, &mut lrw1, &mut liw1);
    } else {
        lrw1 = 0;
        liw1 = 0;
    }
    *stepper.lrw1.borrow_mut() = lrw1;
    *stepper.liw1.borrow_mut() = liw1;

    /* Set the number of forcing vectors and allocate vectors */
    *stepper.nforcing.borrow_mut() = count;

    let nforcing_allocated = *stepper.nforcing_allocated.borrow();
    let nforcing = *stepper.nforcing.borrow();
    if nforcing_allocated < nforcing {
        let mut forcing = std::mem::take(&mut *stepper.forcing.borrow_mut());
        let mut lrw = *stepper.lrw.borrow();
        let mut liw = *stepper.liw.borrow();
        if nforcing_allocated != 0 {
            arkFreeVecArray(
                nforcing_allocated,
                &mut forcing,
                lrw1,
                &mut lrw,
                liw1,
                &mut liw,
            );
        }
        let ok = arkAllocVecArray(nforcing, tmpl, &mut forcing, lrw1, &mut lrw, liw1, &mut liw);
        *stepper.forcing.borrow_mut() = forcing;
        *stepper.lrw.borrow_mut() = lrw;
        *stepper.liw.borrow_mut() = liw;
        if !ok {
            mriStepInnerStepper_FreeVecs(stepper);
            return ARK_MEM_FAIL;
        }
        *stepper.nforcing_allocated.borrow_mut() = nforcing;
    }

    /* Allocate fused operation workspace arrays. `vecs` is N_Vector handle
    scratch that MRIStepInnerStepper_AddForcing rebuilds on demand (an
    N_Vector array cannot be left uninitialised in safe Rust), so only
    `vals` is materialised; the C NULL-return failure branches are
    unreachable because Vec allocation aborts rather than returning NULL. */
    let vals_empty = stepper.vals.borrow().is_empty();
    if vals_empty {
        *stepper.vals.borrow_mut() = vec![ZERO; (count + 1) as usize];
    }

    ARK_SUCCESS
}

/* Resize MRI forcing and fused op workspace vectors if necessary */
pub fn mriStepInnerStepper_Resize(
    stepper: &MRIStepInnerStepper,
    resize: Option<ARKVecResizeFn>,
    resize_data: &mut Option<Box<dyn Any>>,
    lrw_diff: sunindextype,
    liw_diff: sunindextype,
    tmpl: &N_Vector,
) -> i32 {
    let nforcing_allocated = *stepper.nforcing_allocated.borrow();
    let mut forcing = std::mem::take(&mut *stepper.forcing.borrow_mut());
    let mut lrw = *stepper.lrw.borrow();
    let mut liw = *stepper.liw.borrow();

    let ok = arkResizeVecArray(
        resize,
        resize_data,
        nforcing_allocated,
        tmpl,
        &mut forcing,
        lrw_diff,
        &mut lrw,
        liw_diff,
        &mut liw,
    );

    *stepper.forcing.borrow_mut() = forcing;
    *stepper.lrw.borrow_mut() = lrw;
    *stepper.liw.borrow_mut() = liw;

    if !ok {
        return ARK_MEM_FAIL;
    }

    ARK_SUCCESS
}

/* Free MRI forcing and fused op workspace vectors if necessary */
pub fn mriStepInnerStepper_FreeVecs(stepper: &MRIStepInnerStepper) -> i32 {
    let nforcing_allocated = *stepper.nforcing_allocated.borrow();
    let lrw1 = *stepper.lrw1.borrow();
    let liw1 = *stepper.liw1.borrow();
    let mut forcing = std::mem::take(&mut *stepper.forcing.borrow_mut());
    let mut lrw = *stepper.lrw.borrow();
    let mut liw = *stepper.liw.borrow();

    arkFreeVecArray(
        nforcing_allocated,
        &mut forcing,
        lrw1,
        &mut lrw,
        liw1,
        &mut liw,
    );

    *stepper.forcing.borrow_mut() = forcing;
    *stepper.lrw.borrow_mut() = lrw;
    *stepper.liw.borrow_mut() = liw;

    let vecs_alloc = !stepper.vecs.borrow().is_empty();
    if vecs_alloc {
        *stepper.vecs.borrow_mut() = Vec::new();
    }

    let vals_alloc = !stepper.vals.borrow().is_empty();
    if vals_alloc {
        *stepper.vals.borrow_mut() = Vec::new();
    }

    ARK_SUCCESS
}

/* Print forcing vectors to output file */
pub fn mriStepInnerStepper_PrintMem(stepper: &MRIStepInnerStepper, outfile: &SUNFile) {
    /* output data from the inner stepper */
    outfile.write_str("MRIStepInnerStepper Mem:\n");
    outfile.write_str(&format!(
        "MRIStepInnerStepper: inner_nforcing = {}\n",
        *stepper.nforcing.borrow()
    ));
}

/*---------------------------------------------------------------
  Utility routines for MRIStep to serve as an MRIStepInnerStepper
  ---------------------------------------------------------------*/

/*------------------------------------------------------------------------------
  mriStep_ApplyForcing

  Determines the linear combination coefficients and vectors to apply forcing
  at a given value of the independent variable (t).  This occurs through
  appending coefficients and N_Vector pointers to the underlying cvals and Xvecs
  arrays in the step_mem structure.  The dereferenced input *nvec should indicate
  the next available entry in the cvals/Xvecs arrays.  The input 's' is a
  scaling factor that should be applied to each of these coefficients.

  Port note: C appends into `step_mem->cvals` / `step_mem->Xvecs`; because an
  `N_Vector` array cannot be left uninitialised in safe Rust, every call site
  builds those two arrays as function-local `Vec`s (the locked house pattern)
  and passes them in here. `*nvec` keeps its C meaning -- the next free slot,
  i.e. the current length of both arrays -- and every C call site fills the
  arrays strictly left-to-right, so `push` is exactly C's index assignment.
  ----------------------------------------------------------------------------*/
pub fn mriStep_ApplyForcing(
    step_mem: &ARKodeMRIStepMemRec,
    t: sunrealtype,
    s: sunrealtype,
    nvec: &mut i32,
    cvals: &mut Vec<sunrealtype>,
    Xvecs: &mut Vec<N_Vector>,
) {
    /* always append the constant forcing term */
    cvals.push(s);
    Xvecs.push(step_mem.forcing[0].clone());
    *nvec += 1;

    /* compute normalized time tau and initialize tau^i */
    let tau = (t - step_mem.tshift) / (step_mem.tscale);
    let mut taui = tau;
    for i in 1..step_mem.nforcing {
        cvals.push(s * taui);
        Xvecs.push(step_mem.forcing[i as usize].clone());
        taui *= tau;
        *nvec += 1;
    }
}

/*------------------------------------------------------------------------------
  mriStep_SetInnerForcing

  Sets an array of coefficient vectors for a time-dependent external polynomial
  forcing term in the ODE RHS i.e., y' = f(t,y) + p(t). This function is
  primarily intended for using MRIStep as an inner integrator within another
  [outer] instance of MRIStep, where this instance is used to solve a
  modified ODE at a fast time scale. The polynomial is of the form

  p(t) = sum_{i = 0}^{nvecs - 1} forcing[i] * ((t - tshift) / (tscale))^i

  where tshift and tscale are used to normalize the time t (e.g., with MRIGARK
  methods).
  ----------------------------------------------------------------------------*/
pub fn mriStep_SetInnerForcing(
    ark_mem: &ARKodeMem,
    tshift: sunrealtype,
    tscale: sunrealtype,
    forcing: &[N_Vector],
    nvecs: i32,
) -> i32 {
    /* access ARKodeMRIStepMem structure */
    let retval = mriStep_AccessStepMem(ark_mem, "mriStep_SetInnerForcing");
    if retval != ARK_SUCCESS {
        return retval;
    }

    if nvecs > 0 {
        /* enable forcing, and signal that the corresponding pre-existing RHS
           vector is no longer current, since it has a stale forcing function */
        {
            let mut step_mem = mriStep_mem_mut(ark_mem);
            if step_mem.explicit_rhs {
                step_mem.expforcing = SUNTRUE;
                step_mem.impforcing = SUNFALSE;
                step_mem.fse_is_current = SUNFALSE;
            } else {
                step_mem.expforcing = SUNFALSE;
                step_mem.impforcing = SUNTRUE;
                step_mem.fsi_is_current = SUNFALSE;
            }
            step_mem.tshift = tshift;
            step_mem.tscale = tscale;
            step_mem.forcing = forcing.to_vec();
            step_mem.nforcing = nvecs;
        }

        /* Signal that any pre-existing RHS vector is no longer current, since it
           has a stale forcing function */
        ark_mem.borrow_mut().fn_is_current = SUNFALSE;

        /* If the coupling table is NULL, then mriStep_Init has not been called and
           the number of stages has not been set yet. In this case, the workspace
           arrays for fused vector operations will be re-allocated in mriStep_Init
           if necessary to account the value of nforcing. On subsequent calls we
           check if enough space has already been allocated in case nforcing has
           increased since the original allocation. */
        let mric = mriStep_mem_mut(ark_mem).MRIC.clone();
        if let Some(mric) = mric {
            let mric_stages = mric.borrow().stages;

            /* check if there are enough reusable arrays for fused operations */
            let (nfusedopvecs, allocated) = {
                let step_mem = mriStep_mem_mut(ark_mem);
                /* `Xvecs` is N_Vector handle scratch and stays empty in this port
                (call sites rebuild it as a local), so C's `Xvecs != NULL` test
                is taken from `cvals`, which C allocates and frees in lockstep
                with it -- keeping the lrw/liw accounting identical. */
                (step_mem.nfusedopvecs, !step_mem.cvals.is_empty())
            };
            if (nfusedopvecs - nvecs) < (2 * mric_stages + 2) {
                /* free current work space */
                if allocated {
                    mriStep_mem_mut(ark_mem).cvals = Vec::new();
                    ark_mem.borrow_mut().lrw -= nfusedopvecs as i64;
                }
                if allocated {
                    mriStep_mem_mut(ark_mem).Xvecs = Vec::new();
                    ark_mem.borrow_mut().liw -= nfusedopvecs as i64;
                }

                /* allocate reusable arrays for fused vector operations */
                let new_nfusedopvecs = 2 * mric_stages + 2 + nvecs;
                {
                    let mut step_mem = mriStep_mem_mut(ark_mem);
                    step_mem.nfusedopvecs = new_nfusedopvecs;
                    step_mem.cvals = vec![ZERO; new_nfusedopvecs as usize];
                }
                ark_mem.borrow_mut().lrw += new_nfusedopvecs as i64;
                {
                    let mut step_mem = mriStep_mem_mut(ark_mem);
                    step_mem.Xvecs = Vec::new();
                }
                ark_mem.borrow_mut().liw += new_nfusedopvecs as i64;
            }
        }
    } else {
        /* disable forcing */
        let mut step_mem = mriStep_mem_mut(ark_mem);
        step_mem.expforcing = SUNFALSE;
        step_mem.impforcing = SUNFALSE;
        step_mem.tshift = ZERO;
        step_mem.tscale = ONE;
        step_mem.forcing = Vec::new();
        step_mem.nforcing = 0;
    }

    ARK_SUCCESS
}
