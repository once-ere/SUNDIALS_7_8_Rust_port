//! Port of `src/arkode/arkode.c` — the main ARKODE infrastructure,
//! independent of the ARKODE time-step module, nonlinear solver, linear
//! solver and vector modules in use.
//!
//! `arkode_impl.h` (+ the constants/typedefs of `include/arkode/arkode.h`,
//! `arkode_adapt_impl.h`, `arkode_root_impl.h` and
//! `arkode_relaxation_impl.h`) is folded into `crate::arkode_impl`, which
//! also owns `arkProcessError` and every `MSG_ARK_*` message so that all
//! arkode modules share one definition.
//!
//! Reference build configuration: SUNDIALS_LOGGING_LEVEL = 2
//! (`SUNLogInfo`/`SUNLogInfoIf`/`SUNLogDebug`/`SUNLogExtraDebug*` call
//! sites omitted at translation time; `ARK_WARNING` paths are kept because
//! they print through the logger), profiling OFF
//! (`SUNDIALS_MARK_FUNCTION_BEGIN/END` omitted), error checks OFF
//! (`SUNAssert`/`SUNCheck*` are no-ops), monitoring ENABLED, serial
//! branches only. `SUNDIALS_DEBUG_PRINTVEC` is not defined, so the vector
//! dump inside `ARKodePrintMem` is dead code and is omitted.
//! `SUNDIALS_ENABLE_PYTHON` is not defined, so
//! `arkode_user_supplied_fn_table_destroy` is not called in `ARKodeFree`
//! (the `ark_mem->python = NULL` assignment is kept).
//!
//! Handle model: `ARKodeMem = Rc<RefCell<ARKodeMemRec>>`; `ark_mem->ycur`
//! is an `Rc` clone of the caller's `yout`/`y0`, so it aliases the user
//! buffer exactly as the C pointer copy does and no explicit copy-back is
//! required.
//!
//! `arkExpStab` is declared in `arkode_impl.h` but never defined anywhere
//! in the upstream C tree (dead prototype); it is therefore not ported.

use std::any::Any;
use std::cell::RefCell;
use std::rc::Rc;

use crate::arkode_adapt::{arkAdapt, arkAdaptInit, arkPrintAdaptMem};
use crate::arkode_impl::*;
use crate::arkode_interp::{
    arkInterpCreate_Hermite, arkInterpCreate_Lagrange, arkInterpEvaluate, arkInterpFree,
    arkInterpInit, arkInterpPrintMem, arkInterpResize, arkInterpSetDegree, arkInterpUpdate,
};
use crate::arkode_io::{
    ARKodeGetAccumulatedError, ARKodeResetAccumulatedError, ARKodeSetDefaults, ARKodeSetStopTime,
};
use crate::arkode_mristep::{
    MRIStepInnerStepper, MRIStepInnerStepper_Create, MRIStepInnerStepper_GetContent,
    MRIStepInnerStepper_GetForcingData, MRIStepInnerStepper_SetAccumulatedErrorGetFn,
    MRIStepInnerStepper_SetAccumulatedErrorResetFn, MRIStepInnerStepper_SetContent,
    MRIStepInnerStepper_SetEvolveFn, MRIStepInnerStepper_SetFullRhsFn,
    MRIStepInnerStepper_SetRTolFn, MRIStepInnerStepper_SetResetFn,
};
use crate::arkode_relaxation::{arkRelax, arkRelaxDestroy};
use crate::arkode_root::{
    arkPrintRootMem, arkRootCheck1, arkRootCheck2, arkRootCheck3, arkRootFree,
};

use sundials_core::sundials_adaptcontroller::*;
use sundials_core::sundials_context::SUNContext;
use sundials_core::sundials_errors::SUN_SUCCESS;
use sundials_core::sundials_math::*;
use sundials_core::sundials_nvector::*;
use sundials_core::sundials_types::*;
use sundials_core::sundials_utils::*;

/*===============================================================
  Module-local numeric constants

  `arkode.c` has no file-scope numeric `#define`s of its own: every
  named constant it uses (ZERO, TENTH, HALF, ONE, TWO, FOUR, FUZZ_FACTOR,
  H0_LBFACTOR, H0_UBFACTOR, H0_BIAS, H0_ITERS, ONEPSM, ...) comes from
  `arkode_impl.h`, which is folded into `crate::arkode_impl` and glob-
  imported above.  Per the frozen contract (section 7) those must NOT be
  redefined here.  The only bare literals in the C file
  (`SUN_RCONST(0.2)` in `arkHin`, `SUN_RCONST(0.75)` in
  `arkPredict_VariableOrder`, `SUN_RCONST(0.9)` in `arkCheckConstraints`)
  are written inline exactly where C writes them.
  ===============================================================*/

/*===============================================================
  Callback invocation helpers

  Granular borrow discipline: the `Box<dyn Any>` data token is taken out
  of the mem around every user callback call and restored afterwards on
  every path, and no mem borrow is held across the call.
  ===============================================================*/

/// Invoke the error-weight function
/// (C: `ark_mem->efun(y, ewt, ark_mem->e_data)`).
///
/// In C, `e_data` is `ark_mem` for the built-in `arkEwtSetSS`/`arkEwtSetSV`
/// and an alias of `user_data` when the user supplied `efun` through
/// `ARKodeWFtolerances` (`ARKodeSetUserData` keeps the alias in sync).
/// A `Box` cannot alias, so the port stores a boxed `ARKodeMem` handle in
/// `e_data` for the built-in case and passes the CURRENT `user_data` box
/// when `user_efun` is set (accepted deviation class 6).
fn ark_call_efun(ark_mem: &ARKodeMem, y: &N_Vector, ewt: &N_Vector) -> i32 {
    let (efun, user_efun) = {
        let m = ark_mem.borrow();
        (m.efun, m.user_efun)
    };
    let efun = efun.expect("efun set");
    if user_efun {
        let mut data = ark_mem.borrow_mut().user_data.take();
        let retval = efun(y, ewt, &mut data);
        ark_mem.borrow_mut().user_data = data;
        retval
    } else {
        let mut data = ark_mem.borrow_mut().e_data.take();
        let retval = efun(y, ewt, &mut data);
        ark_mem.borrow_mut().e_data = data;
        retval
    }
}

/// Invoke the residual-weight function
/// (C: `ark_mem->rfun(y, rwt, ark_mem->r_data)`); same `r_data` treatment
/// as `ark_call_efun`.
fn ark_call_rfun(ark_mem: &ARKodeMem, y: &N_Vector, rwt: &N_Vector) -> i32 {
    let (rfun, user_rfun) = {
        let m = ark_mem.borrow();
        (m.rfun, m.user_rfun)
    };
    let rfun = rfun.expect("rfun set");
    if user_rfun {
        let mut data = ark_mem.borrow_mut().user_data.take();
        let retval = rfun(y, rwt, &mut data);
        ark_mem.borrow_mut().user_data = data;
        retval
    } else {
        let mut data = ark_mem.borrow_mut().r_data.take();
        let retval = rfun(y, rwt, &mut data);
        ark_mem.borrow_mut().r_data = data;
        retval
    }
}

/// Invoke the user pre-step function
/// (C: `ark_mem->PreStepFn(t, y, step, attempt, ark_mem->user_data)`).
fn ark_call_prestepfn(
    ark_mem: &ARKodeMem,
    t: sunrealtype,
    y: &N_Vector,
    step: i64,
    attempt: i32,
) -> i32 {
    let f = ark_mem.borrow().PreStepFn.expect("PreStepFn set");
    let mut user_data = ark_mem.borrow_mut().user_data.take();
    let retval = f(t, y, step, attempt, &mut user_data);
    ark_mem.borrow_mut().user_data = user_data;
    retval
}

/*===============================================================
  Exported functions
  ===============================================================*/

/*---------------------------------------------------------------
  ARKodeResize:

  ARKodeResize re-initializes ARKODE's memory for a problem with a
  changing vector size.  It is assumed that the problem dynamics
  before and after the vector resize will be comparable, so that
  all time-stepping heuristics prior to calling ARKodeResize
  remain valid after the call.  If instead the dynamics should be
  re-calibrated, the ARKODE memory structure should be deleted
  with a call to ARKodeFree, and re-created with a call to
  *StepCreate.

  To aid in the vector-resize operation, the user can supply a
  vector resize function, that will take as input an N_Vector with
  the previous size, and return as output a corresponding vector
  of the new size.  If this function (of type ARKVecResizeFn) is
  not supplied (i.e. is set to NULL), then all existing N_Vectors
  will be destroyed and re-cloned from the input vector.

  In the case that the dynamical time scale should be modified
  slightly from the previous time scale, an input "hscale" is
  allowed, that will re-scale the upcoming time step by the
  specified factor.  If a value <= 0 is specified, the default of
  1.0 will be used.

  Other arguments:
  ark_mem          Existing ARKODE memory data structure.
  y0               The newly-sized solution vector, holding
                   the current dependent variable values.
  t0               The current value of the independent
                   variable.
  resize_data      User-supplied data structure that will be
                   passed to the supplied resize function.

  The return value is ARK_SUCCESS = 0 if no errors occurred, or
  a negative value otherwise.
  ---------------------------------------------------------------*/
pub fn ARKodeResize(
    arkode_mem: &ARKodeMem,
    y0: &N_Vector,
    hscale: sunrealtype,
    t0: sunrealtype,
    resize: Option<ARKVecResizeFn>,
    resize_data: &mut Option<Box<dyn Any>>,
) -> i32 {
    let mut hscale = hscale;

    /* Check ark_mem: NULL-mem check handled by the type system */
    let ark_mem = arkode_mem;

    /* Check if ark_mem was allocated */
    if !ark_mem.borrow().MallocDone {
        arkProcessError(
            Some(ark_mem),
            ARK_NO_MALLOC,
            line!() as i32,
            "ARKodeResize",
            file!(),
            MSG_ARK_NO_MALLOC,
        );
        return ARK_NO_MALLOC;
    }

    /* Check for legal input parameters (NULL y0 handled by the type
    system; the Rc clone aliases the caller's vector exactly as the C
    pointer copy does) */
    ark_mem.borrow_mut().ycur = Some(y0.clone());

    /* Copy the input parameters into ARKODE state */
    {
        let mut m = ark_mem.borrow_mut();
        m.tcur = t0;
        m.tn = t0;
    }

    /* Update time-stepping parameters */
    /*   adjust upcoming step size depending on hscale */
    if hscale <= ZERO {
        hscale = ONE;
    }
    if hscale != ONE {
        let mut m = ark_mem.borrow_mut();

        /* Encode hscale into ark_mem structure */
        m.eta = hscale;
        m.hprime *= hscale;

        /* If next step would overtake tstop, adjust stepsize */
        if m.tstopset && (m.tcur + m.hprime - m.tstop) * m.hprime > ZERO {
            m.hprime = (m.tstop - m.tcur) * (ONE - FOUR * m.uround);
            m.eta = m.hprime / m.h;
        }
    }

    /* Determine change in vector sizes */
    let mut lrw1: sunindextype = 0;
    let mut liw1: sunindextype = 0;
    if y0.ops.borrow().nvspace.is_some() {
        N_VSpace(y0, &mut lrw1, &mut liw1);
    }
    let (lrw_diff, liw_diff) = {
        let mut m = ark_mem.borrow_mut();
        let lrw_diff = lrw1 - m.lrw1;
        let liw_diff = liw1 - m.liw1;
        m.lrw1 = lrw1;
        m.liw1 = liw1;
        (lrw_diff, liw_diff)
    };

    /* Disable constraints, the user will need to set a new constraint vector for
       the updated problem size */
    {
        let mut constraints = ark_mem.borrow_mut().constraints.take();
        arkFreeVec(ark_mem, &mut constraints);
        ark_mem.borrow_mut().constraints = constraints;
    }

    /* Resize the solver vectors (using y0 as a template) */
    let resizeOK = arkResizeVectors(ark_mem, resize, resize_data, lrw_diff, liw_diff, y0);
    if !resizeOK {
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_FAIL,
            line!() as i32,
            "ARKodeResize",
            file!(),
            "Unable to resize vector",
        );
        return ARK_MEM_FAIL;
    }

    /* Resize the interpolation structure memory */
    let interp = ark_mem.borrow().interp.clone();
    if let Some(interp) = interp {
        let retval = arkInterpResize(
            ark_mem,
            &interp,
            resize,
            resize_data,
            lrw_diff,
            liw_diff,
            y0,
        );
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                retval,
                line!() as i32,
                "ARKodeResize",
                file!(),
                "Interpolation module resize failure",
            );
            return retval;
        }
    }

    /* Copy y0 into ark_yn to set the current solution */
    let yn = ark_mem.borrow().yn.clone().expect("yn allocated");
    N_VScale(ONE, y0, &yn);

    {
        let mut m = ark_mem.borrow_mut();
        m.fn_is_current = SUNFALSE;

        /* Indicate that problem needs to be initialized */
        m.initsetup = SUNTRUE;
        m.init_type = RESIZE_INIT;
        m.firststage = SUNTRUE;
    }

    /* Call the stepper-specific resize (if provided) */
    let step_resize = ark_mem.borrow().step_resize;
    if let Some(step_resize) = step_resize {
        return step_resize(ark_mem, y0, hscale, t0, resize, resize_data);
    }

    /* Problem has been successfully re-sized */
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeReset:

  This routine resets an ARKode module to solve the same
  problem from the given time with the input state (all counter
  values are retained).
  ---------------------------------------------------------------*/
pub fn ARKodeReset(arkode_mem: &ARKodeMem, tR: sunrealtype, yR: &N_Vector) -> i32 {
    /* NULL-mem check: handled by the type system */
    let ark_mem = arkode_mem;

    /* Reset main ARKODE infrastructure */
    let retval = arkInit(ark_mem, tR, yR, RESET_INIT);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            retval,
            line!() as i32,
            "ARKodeReset",
            file!(),
            "ARKode reset failure",
        );
        return retval;
    }

    /* Call stepper routine to perform remaining reset operations (if provided) */
    let step_reset = ark_mem.borrow().step_reset;
    if let Some(step_reset) = step_reset {
        return step_reset(ark_mem, tR, yR);
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSStolerances, ARKodeSVtolerances, ARKodeWFtolerances:

  These functions specify the integration tolerances. One of them
  SHOULD be called before the first call to ARKodeEvolve; otherwise
  default values of reltol=1e-4 and abstol=1e-9 will be used,
  which may be entirely incorrect for a specific problem.

  ARKodeSStolerances specifies scalar relative and absolute
  tolerances.

  ARKodeSVtolerances specifies scalar relative tolerance and a
  vector absolute tolerance (a potentially different absolute
  tolerance for each vector component).

  ARKodeWFtolerances specifies a user-provides function (of type
  ARKEwtFn) which will be called to set the error weight vector.
  ---------------------------------------------------------------*/
pub fn ARKodeSStolerances(arkode_mem: &ARKodeMem, reltol: sunrealtype, abstol: sunrealtype) -> i32 {
    /* unpack ark_mem: NULL-mem check handled by the type system */
    let ark_mem = arkode_mem;

    /* Check inputs */
    if !ark_mem.borrow().MallocDone {
        arkProcessError(
            Some(ark_mem),
            ARK_NO_MALLOC,
            line!() as i32,
            "ARKodeSStolerances",
            file!(),
            MSG_ARK_NO_MALLOC,
        );
        return ARK_NO_MALLOC;
    }
    if reltol < ZERO {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!() as i32,
            "ARKodeSStolerances",
            file!(),
            MSG_ARK_BAD_RELTOL,
        );
        return ARK_ILL_INPUT;
    }
    if abstol < ZERO {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!() as i32,
            "ARKodeSStolerances",
            file!(),
            MSG_ARK_BAD_ABSTOL,
        );
        return ARK_ILL_INPUT;
    }

    /* Ensure that vector supports N_VAddConst */
    let has_nvaddconst = {
        let tempv1 = ark_mem.borrow().tempv1.clone().expect("tempv1 allocated");
        let has = tempv1.ops.borrow().nvaddconst.is_some();
        has
    };
    if !has_nvaddconst {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!() as i32,
            "ARKodeSStolerances",
            file!(),
            "N_VAddConst unimplemented (required for scalar abstol)",
        );
        return ARK_ILL_INPUT;
    }

    let mut m = ark_mem.borrow_mut();

    /* Set flag indicating whether abstol == 0 */
    m.atolmin0 = abstol == ZERO;

    /* Copy tolerances into memory */
    m.reltol = reltol;
    m.Sabstol = abstol;
    m.itol = ARK_SS;

    /* enforce use of arkEwtSetSS */
    m.user_efun = SUNFALSE;
    m.efun = Some(arkEwtSetSS);
    /* C: e_data = ark_mem -- the built-in error-weight function reaches
    the integrator through a boxed handle clone */
    m.e_data = Some(Box::new(ark_mem.clone()));

    ARK_SUCCESS
}

pub fn ARKodeSVtolerances(arkode_mem: &ARKodeMem, reltol: sunrealtype, abstol: &N_Vector) -> i32 {
    /* unpack ark_mem: NULL-mem check handled by the type system */
    let ark_mem = arkode_mem;

    /* Check inputs */
    if !ark_mem.borrow().MallocDone {
        arkProcessError(
            Some(ark_mem),
            ARK_NO_MALLOC,
            line!() as i32,
            "ARKodeSVtolerances",
            file!(),
            MSG_ARK_NO_MALLOC,
        );
        return ARK_NO_MALLOC;
    }
    if reltol < ZERO {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!() as i32,
            "ARKodeSVtolerances",
            file!(),
            MSG_ARK_BAD_RELTOL,
        );
        return ARK_ILL_INPUT;
    }
    /* NULL abstol check: handled by the type system */
    if abstol.ops.borrow().nvmin.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!() as i32,
            "ARKodeSVtolerances",
            file!(),
            "Missing N_VMin routine from N_Vector",
        );
        return ARK_ILL_INPUT;
    }
    let abstolmin = N_VMin(abstol);
    if abstolmin < ZERO {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!() as i32,
            "ARKodeSVtolerances",
            file!(),
            MSG_ARK_BAD_ABSTOL,
        );
        return ARK_ILL_INPUT;
    }

    /* Set flag indicating whether min(abstol) == 0 */
    ark_mem.borrow_mut().atolmin0 = abstolmin == ZERO;

    /* Copy tolerances into memory */
    if !ark_mem.borrow().VabstolMallocDone {
        let ewt = ark_mem.borrow().ewt.clone().expect("ewt allocated");
        let mut Vabstol = ark_mem.borrow_mut().Vabstol.take();
        let allocOK = arkAllocVec(ark_mem, &ewt, &mut Vabstol);
        ark_mem.borrow_mut().Vabstol = Vabstol;
        if !allocOK {
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_FAIL,
                line!() as i32,
                "ARKodeSVtolerances",
                file!(),
                MSG_ARK_ARKMEM_FAIL,
            );
            return ARK_ILL_INPUT;
        }
        ark_mem.borrow_mut().VabstolMallocDone = SUNTRUE;
    }
    let Vabstol = ark_mem.borrow().Vabstol.clone().expect("Vabstol allocated");
    N_VScale(ONE, abstol, &Vabstol);

    let mut m = ark_mem.borrow_mut();
    m.reltol = reltol;
    m.itol = ARK_SV;

    /* enforce use of arkEwtSetSV */
    m.user_efun = SUNFALSE;
    m.efun = Some(arkEwtSetSV);
    /* C: e_data = ark_mem (see ARKodeSStolerances) */
    m.e_data = Some(Box::new(ark_mem.clone()));

    ARK_SUCCESS
}

pub fn ARKodeWFtolerances(arkode_mem: &ARKodeMem, efun: ARKEwtFn) -> i32 {
    /* unpack ark_mem: NULL-mem check handled by the type system */
    let ark_mem = arkode_mem;

    if !ark_mem.borrow().MallocDone {
        arkProcessError(
            Some(ark_mem),
            ARK_NO_MALLOC,
            line!() as i32,
            "ARKodeWFtolerances",
            file!(),
            MSG_ARK_NO_MALLOC,
        );
        return ARK_NO_MALLOC;
    }

    /* Copy tolerance data into memory */
    let mut m = ark_mem.borrow_mut();
    m.itol = ARK_WF;
    m.user_efun = SUNTRUE;
    m.efun = Some(efun);
    /* C: e_data = ark_mem->user_data -- a raw pointer snapshot that a
    `Box` cannot reproduce (accepted deviation class 6).  The token is
    cleared here and `ark_call_efun` passes the CURRENT `user_data` box
    whenever `user_efun` is set. */
    m.e_data = None;

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeResStolerance, ARKodeResVtolerance, ARKodeResFtolerance:

  These functions specify the absolute residual tolerance.
  Specification of the absolute residual tolerance is only
  necessary for problems with non-identity mass matrices in which
  the units of the solution vector y dramatically differ from the
  units of the ODE right-hand side f(t,y).  If this occurs, one
  of these routines SHOULD be called before the first call to
  ARKODE; otherwise the default value of rabstol=1e-9 will be
  used, which may be entirely incorrect for a specific problem.

  ARKodeResStolerances specifies a scalar residual tolerance.

  ARKodeResVtolerances specifies a vector residual tolerance
  (a potentially different absolute residual tolerance for
  each vector component).

  ARKodeResFtolerances specifies a user-provides function (of
  type ARKRwtFn) which will be called to set the residual
  weight vector.
  ---------------------------------------------------------------*/
pub fn ARKodeResStolerance(arkode_mem: &ARKodeMem, rabstol: sunrealtype) -> i32 {
    /* unpack ark_mem: NULL-mem check handled by the type system */
    let ark_mem = arkode_mem;

    /* Guard against use for time steppers that do not support mass matrices */
    if !ark_mem.borrow().step_supports_massmatrix {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!() as i32,
            "ARKodeResStolerance",
            file!(),
            "time-stepping module does not support non-identity mass matrices",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* Check inputs */
    if !ark_mem.borrow().MallocDone {
        arkProcessError(
            Some(ark_mem),
            ARK_NO_MALLOC,
            line!() as i32,
            "ARKodeResStolerance",
            file!(),
            MSG_ARK_NO_MALLOC,
        );
        return ARK_NO_MALLOC;
    }
    if rabstol < ZERO {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!() as i32,
            "ARKodeResStolerance",
            file!(),
            MSG_ARK_BAD_RABSTOL,
        );
        return ARK_ILL_INPUT;
    }

    /* Ensure that vector supports N_VAddConst */
    let has_nvaddconst = {
        let tempv1 = ark_mem.borrow().tempv1.clone().expect("tempv1 allocated");
        let has = tempv1.ops.borrow().nvaddconst.is_some();
        has
    };
    if !has_nvaddconst {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!() as i32,
            "ARKodeResStolerance",
            file!(),
            "N_VAddConst unimplemented (required for scalar rabstol)",
        );
        return ARK_ILL_INPUT;
    }

    /* Set flag indicating whether rabstol == 0 */
    ark_mem.borrow_mut().Ratolmin0 = rabstol == ZERO;

    /* Allocate space for rwt if necessary */
    if ark_mem.borrow().rwt_is_ewt {
        ark_mem.borrow_mut().rwt = None;
        let ewt = ark_mem.borrow().ewt.clone().expect("ewt allocated");
        let mut rwt = ark_mem.borrow_mut().rwt.take();
        let allocOK = arkAllocVec(ark_mem, &ewt, &mut rwt);
        ark_mem.borrow_mut().rwt = rwt;
        if !allocOK {
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_FAIL,
                line!() as i32,
                "ARKodeResStolerance",
                file!(),
                MSG_ARK_ARKMEM_FAIL,
            );
            return ARK_ILL_INPUT;
        }
        ark_mem.borrow_mut().rwt_is_ewt = SUNFALSE;
    }

    /* Copy tolerances into memory */
    let mut m = ark_mem.borrow_mut();
    m.SRabstol = rabstol;
    m.ritol = ARK_SS;

    /* enforce use of arkRwtSet
       (upstream really does clear `user_efun` and not `user_rfun` here --
       preserved verbatim) */
    m.user_efun = SUNFALSE;
    m.rfun = Some(arkRwtSet);
    /* C: r_data = ark_mem */
    m.r_data = Some(Box::new(ark_mem.clone()));

    ARK_SUCCESS
}

pub fn ARKodeResVtolerance(arkode_mem: &ARKodeMem, rabstol: &N_Vector) -> i32 {
    /* unpack ark_mem: NULL-mem check handled by the type system */
    let ark_mem = arkode_mem;

    /* Guard against use for time steppers that do not support mass matrices */
    if !ark_mem.borrow().step_supports_massmatrix {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!() as i32,
            "ARKodeResVtolerance",
            file!(),
            "time-stepping module does not support non-identity mass matrices",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* Check inputs */
    if !ark_mem.borrow().MallocDone {
        arkProcessError(
            Some(ark_mem),
            ARK_NO_MALLOC,
            line!() as i32,
            "ARKodeResVtolerance",
            file!(),
            MSG_ARK_NO_MALLOC,
        );
        return ARK_NO_MALLOC;
    }
    /* NULL rabstol check: handled by the type system */
    if rabstol.ops.borrow().nvmin.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!() as i32,
            "ARKodeResVtolerance",
            file!(),
            "Missing N_VMin routine from N_Vector",
        );
        return ARK_ILL_INPUT;
    }
    let rabstolmin = N_VMin(rabstol);
    if rabstolmin < ZERO {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!() as i32,
            "ARKodeResVtolerance",
            file!(),
            MSG_ARK_BAD_RABSTOL,
        );
        return ARK_ILL_INPUT;
    }

    /* Set flag indicating whether min(abstol) == 0 */
    ark_mem.borrow_mut().Ratolmin0 = rabstolmin == ZERO;

    /* Allocate space for rwt if necessary */
    if ark_mem.borrow().rwt_is_ewt {
        ark_mem.borrow_mut().rwt = None;
        let ewt = ark_mem.borrow().ewt.clone().expect("ewt allocated");
        let mut rwt = ark_mem.borrow_mut().rwt.take();
        let allocOK = arkAllocVec(ark_mem, &ewt, &mut rwt);
        ark_mem.borrow_mut().rwt = rwt;
        if !allocOK {
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_FAIL,
                line!() as i32,
                "ARKodeResVtolerance",
                file!(),
                MSG_ARK_ARKMEM_FAIL,
            );
            return ARK_ILL_INPUT;
        }
        ark_mem.borrow_mut().rwt_is_ewt = SUNFALSE;
    }

    /* Copy tolerances into memory */
    if !ark_mem.borrow().VRabstolMallocDone {
        let rwt = ark_mem.borrow().rwt.clone().expect("rwt allocated");
        let mut VRabstol = ark_mem.borrow_mut().VRabstol.take();
        let allocOK = arkAllocVec(ark_mem, &rwt, &mut VRabstol);
        ark_mem.borrow_mut().VRabstol = VRabstol;
        if !allocOK {
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_FAIL,
                line!() as i32,
                "ARKodeResVtolerance",
                file!(),
                MSG_ARK_ARKMEM_FAIL,
            );
            return ARK_ILL_INPUT;
        }
        ark_mem.borrow_mut().VRabstolMallocDone = SUNTRUE;
    }
    let VRabstol = ark_mem
        .borrow()
        .VRabstol
        .clone()
        .expect("VRabstol allocated");
    N_VScale(ONE, rabstol, &VRabstol);

    let mut m = ark_mem.borrow_mut();
    m.ritol = ARK_SV;

    /* enforce use of arkRwtSet (see the note in ARKodeResStolerance) */
    m.user_efun = SUNFALSE;
    m.rfun = Some(arkRwtSet);
    /* C: r_data = ark_mem */
    m.r_data = Some(Box::new(ark_mem.clone()));

    ARK_SUCCESS
}

pub fn ARKodeResFtolerance(arkode_mem: &ARKodeMem, rfun: ARKRwtFn) -> i32 {
    /* unpack ark_mem: NULL-mem check handled by the type system */
    let ark_mem = arkode_mem;

    /* Guard against use for time steppers that do not support mass matrices */
    if !ark_mem.borrow().step_supports_massmatrix {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!() as i32,
            "ARKodeResFtolerance",
            file!(),
            "time-stepping module does not support non-identity mass matrices",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    if !ark_mem.borrow().MallocDone {
        arkProcessError(
            Some(ark_mem),
            ARK_NO_MALLOC,
            line!() as i32,
            "ARKodeResFtolerance",
            file!(),
            MSG_ARK_NO_MALLOC,
        );
        return ARK_NO_MALLOC;
    }

    /* Allocate space for rwt if necessary */
    if ark_mem.borrow().rwt_is_ewt {
        ark_mem.borrow_mut().rwt = None;
        let ewt = ark_mem.borrow().ewt.clone().expect("ewt allocated");
        let mut rwt = ark_mem.borrow_mut().rwt.take();
        let allocOK = arkAllocVec(ark_mem, &ewt, &mut rwt);
        ark_mem.borrow_mut().rwt = rwt;
        if !allocOK {
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_FAIL,
                line!() as i32,
                "ARKodeResFtolerance",
                file!(),
                MSG_ARK_ARKMEM_FAIL,
            );
            return ARK_ILL_INPUT;
        }
        ark_mem.borrow_mut().rwt_is_ewt = SUNFALSE;
    }

    /* Copy tolerance data into memory */
    let mut m = ark_mem.borrow_mut();
    m.ritol = ARK_WF;
    m.user_rfun = SUNTRUE;
    m.rfun = Some(rfun);
    /* C: r_data = ark_mem->user_data -- pointer snapshot, see
    ARKodeWFtolerances (accepted deviation class 6) */
    m.r_data = None;

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeEvolve:

  This routine is the main driver of ARKODE-based integrators.

  It integrates over a time interval defined by the user, by
  calling the time step module to do internal time steps.

  The first time that ARKodeEvolve is called for a successfully
  initialized problem, it computes a tentative initial step size.

  ARKodeEvolve supports two modes as specified by itask: ARK_NORMAL and
  ARK_ONE_STEP.  In the ARK_NORMAL mode, the solver steps until
  it reaches or passes tout and then interpolates to obtain
  y(tout).  In the ARK_ONE_STEP mode, it takes one internal step
  and returns.  The behavior of both modes can be over-rided
  through user-specification of ark_tstop (through the
  *StepSetStopTime function), in which case if a solver step
  would pass tstop, the step is shortened so that it stops at
  exactly the specified stop time, and hence interpolation of
  y(tout) is not required.
  ---------------------------------------------------------------*/
pub fn ARKodeEvolve(
    arkode_mem: &ARKodeMem,
    tout: sunrealtype,
    yout: &N_Vector,
    tret: &mut sunrealtype,
    itask: i32,
) -> i32 {
    /* C leaves `istate` uninitialized; every path that leaves the internal
    step loop assigns it exactly once before its `break`, so the Rust
    declaration is likewise deferred-initialization */
    let istate: i32;

    /* Check and process inputs */

    /* Check if ark_mem exists: handled by the type system */
    let ark_mem = arkode_mem;

    /* Check if ark_mem was allocated */
    if !ark_mem.borrow().MallocDone {
        arkProcessError(
            Some(ark_mem),
            ARK_NO_MALLOC,
            line!() as i32,
            "ARKodeEvolve",
            file!(),
            MSG_ARK_NO_MALLOC,
        );
        return ARK_NO_MALLOC;
    }

    /* Check for yout != NULL (handled by the type system; ark_mem->ycur
    aliases the user's yout -- the Rc clone shares the underlying data
    exactly as the C pointer copy does) */
    ark_mem.borrow_mut().ycur = Some(yout.clone());

    /* Check for tret != NULL: handled by the type system */

    /* Check for valid itask */
    if (itask != ARK_NORMAL) && (itask != ARK_ONE_STEP) {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!() as i32,
            "ARKodeEvolve",
            file!(),
            MSG_ARK_BAD_ITASK,
        );
        return ARK_ILL_INPUT;
    }

    /* start profiler: profiling disabled in the reference build */

    /* perform first-step-specific initializations:
       - initialize tret values to initialization time
       - perform initial integrator setup  */
    if ark_mem.borrow().initsetup {
        {
            let mut m = ark_mem.borrow_mut();
            m.tretlast = m.tcur;
            *tret = m.tcur;
        }
        let retval = arkInitialSetup(ark_mem, tout);
        if retval != ARK_SUCCESS {
            return retval;
        }
    }

    /* perform stopping tests */
    if !ark_mem.borrow().initsetup {
        let mut retval: i32 = ARK_SUCCESS;
        if arkStopTests(ark_mem, tout, yout, tret, itask, &mut retval) != 0 {
            return retval;
        }
    }

    /* fill current independent variable (and optionally ycur with yn) */
    {
        let mut m = ark_mem.borrow_mut();
        m.tcur = m.tn;
    }
    {
        let (ensure_ycur, yn, ycur) = {
            let m = ark_mem.borrow();
            (m.ensure_ycur, m.yn.clone(), m.ycur.clone())
        };
        if ensure_ycur {
            N_VScale(
                ONE,
                &yn.expect("yn allocated"),
                &ycur.expect("ycur attached"),
            );
        }
    }

    /*--------------------------------------------------
      Looping point for successful internal steps

      - update the ewt/rwt vectors for upcoming step
      - check for errors (too many steps, too much
        accuracy requested, step size too small)
      - loop over attempts at a new step:
        * try to take step (via time stepper module),
          handle solver convergence or other failures
        * if the stepper requests ARK_RETRY_STEP, we
          retry the step without accumulating failures.
          A stepper should never request this multiple
          times in a row.
        * perform constraint-handling (if selected)
        * check temporal error
        * if all of the above pass, complete step by
          updating current time, solution, error &
          stepsize history arrays.
      - perform stop tests:
        * check for root in last step taken
        * check if tout was passed
        * check if close to tstop
        * check if in ONE_STEP mode (must return)
      --------------------------------------------------*/
    let mut nstloc: i64 = 0;
    loop {
        {
            let mut m = ark_mem.borrow_mut();
            m.next_h = m.h;
        }

        /* Reset and check ewt and rwt */
        if !ark_mem.borrow().initsetup {
            let (yn, ewt) = {
                let m = ark_mem.borrow();
                (
                    m.yn.clone().expect("yn allocated"),
                    m.ewt.clone().expect("ewt allocated"),
                )
            };
            let ewtsetOK = ark_call_efun(ark_mem, &yn, &ewt);
            if ewtsetOK != 0 {
                let (itol, tcur) = {
                    let m = ark_mem.borrow();
                    (m.itol, m.tcur)
                };
                if itol == ARK_WF {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_ILL_INPUT,
                        line!() as i32,
                        "ARKodeEvolve",
                        file!(),
                        &MSG_ARK_EWT_NOW_FAIL(tcur),
                    );
                } else {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_ILL_INPUT,
                        line!() as i32,
                        "ARKodeEvolve",
                        file!(),
                        &MSG_ARK_EWT_NOW_BAD(tcur),
                    );
                }

                istate = ARK_ILL_INPUT;
                ark_mem.borrow_mut().tretlast = tcur;
                *tret = tcur;
                N_VScale(ONE, &yn, yout);
                break;
            }

            if !ark_mem.borrow().rwt_is_ewt {
                let rwt = ark_mem.borrow().rwt.clone().expect("rwt allocated");
                let ewtsetOK = ark_call_rfun(ark_mem, &yn, &rwt);
                if ewtsetOK != 0 {
                    let (itol, tcur) = {
                        let m = ark_mem.borrow();
                        (m.itol, m.tcur)
                    };
                    if itol == ARK_WF {
                        arkProcessError(
                            Some(ark_mem),
                            ARK_ILL_INPUT,
                            line!() as i32,
                            "ARKodeEvolve",
                            file!(),
                            &MSG_ARK_RWT_NOW_FAIL(tcur),
                        );
                    } else {
                        arkProcessError(
                            Some(ark_mem),
                            ARK_ILL_INPUT,
                            line!() as i32,
                            "ARKodeEvolve",
                            file!(),
                            &MSG_ARK_RWT_NOW_BAD(tcur),
                        );
                    }

                    istate = ARK_ILL_INPUT;
                    ark_mem.borrow_mut().tretlast = tcur;
                    *tret = tcur;
                    N_VScale(ONE, &yn, yout);
                    break;
                }
            }
        }

        /* Check for too many steps */
        {
            let (mxstep, tcur) = {
                let m = ark_mem.borrow();
                (m.mxstep, m.tcur)
            };
            if (mxstep > 0) && (nstloc >= mxstep) {
                arkProcessError(
                    Some(ark_mem),
                    ARK_TOO_MUCH_WORK,
                    line!() as i32,
                    "ARKodeEvolve",
                    file!(),
                    &MSG_ARK_MAX_STEPS(tcur),
                );
                istate = ARK_TOO_MUCH_WORK;
                ark_mem.borrow_mut().tretlast = tcur;
                *tret = tcur;
                let yn = ark_mem.borrow().yn.clone().expect("yn allocated");
                N_VScale(ONE, &yn, yout);
                break;
            }
        }

        /* Check for too much accuracy requested */
        {
            let (yn, ewt, uround) = {
                let m = ark_mem.borrow();
                (
                    m.yn.clone().expect("yn allocated"),
                    m.ewt.clone().expect("ewt allocated"),
                    m.uround,
                )
            };
            let nrm = N_VWrmsNorm(&yn, &ewt);
            ark_mem.borrow_mut().tolsf = uround * nrm;
            let (tolsf, fixedstep, tcur) = {
                let m = ark_mem.borrow();
                (m.tolsf, m.fixedstep, m.tcur)
            };
            if tolsf > ONE && !fixedstep {
                arkProcessError(
                    Some(ark_mem),
                    ARK_TOO_MUCH_ACC,
                    line!() as i32,
                    "ARKodeEvolve",
                    file!(),
                    &MSG_ARK_TOO_MUCH_ACC(tcur),
                );
                istate = ARK_TOO_MUCH_ACC;
                ark_mem.borrow_mut().tretlast = tcur;
                *tret = tcur;
                N_VScale(ONE, &yn, yout);
                ark_mem.borrow_mut().tolsf *= TWO;
                break;
            } else {
                ark_mem.borrow_mut().tolsf = ONE;
            }
        }

        /* Check for h below roundoff level in tn */
        {
            let (tcur, h) = {
                let m = ark_mem.borrow();
                (m.tcur, m.h)
            };
            if tcur + h == tcur {
                ark_mem.borrow_mut().nhnil += 1;
                let (nhnil, mxhnil) = {
                    let m = ark_mem.borrow();
                    (m.nhnil, m.mxhnil)
                };
                if nhnil <= mxhnil {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_WARNING,
                        line!() as i32,
                        "ARKodeEvolve",
                        file!(),
                        &MSG_ARK_HNIL(tcur, h),
                    );
                }
                if nhnil == mxhnil {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_WARNING,
                        line!() as i32,
                        "ARKodeEvolve",
                        file!(),
                        MSG_ARK_HNIL_DONE,
                    );
                }
            }
        }

        /* Update parameter for upcoming step size */
        {
            let mut m = ark_mem.borrow_mut();
            if m.hprime != m.h {
                m.h = m.h * m.eta;
                m.next_h = m.h;
            }
            if m.fixedstep {
                m.h = m.hin;
                m.next_h = m.h;

                /* patch for 'fixedstep' + 'tstop' use case:
                   limit fixed step size if step would overtake tstop */
                if m.tstopset && (m.tcur + m.h - m.tstop) * m.h > ZERO {
                    m.h = (m.tstop - m.tcur) * (ONE - FOUR * m.uround);
                }
            }
        }

        /* Looping point for step attempts */
        let mut dsm: sunrealtype = ZERO;
        let mut kflag: i32 = ARK_SUCCESS;
        let mut relax_fails: i32 = 0;
        let mut nflag: i32 = FIRST_CALL;
        let mut attempts: i32 = 0;
        let mut ncf: i32 = 0;
        let mut nef: i32 = 0;
        let mut constrfails: i32 = 0;
        ark_mem.borrow_mut().last_kflag = 0;
        loop {
            /* increment attempt counters
               Note: kflag can only equal ARK_RETRY_STEP if the stepper rejected
               the current step size before performing calculations. Thus, we do
               not include those when keeping track of step "attempts". */
            if kflag != ARK_RETRY_STEP {
                attempts += 1;
                ark_mem.borrow_mut().nst_attempts += 1;
            }

            /* fill tcur with the last accepted step time */
            {
                let mut m = ark_mem.borrow_mut();
                m.tcur = m.tn;
            }

            /* call the user-supplied pre-step function (if it exists) */
            if ark_mem.borrow().PreStepFn.is_some() {
                let (ensure_ycur, tcur, nst, ycur, yn) = {
                    let m = ark_mem.borrow();
                    (m.ensure_ycur, m.tcur, m.nst, m.ycur.clone(), m.yn.clone())
                };
                let retval = if ensure_ycur {
                    ark_call_prestepfn(
                        ark_mem,
                        tcur,
                        &ycur.expect("ycur attached"),
                        nst,
                        attempts,
                    )
                } else {
                    ark_call_prestepfn(ark_mem, tcur, &yn.expect("yn allocated"), nst, attempts)
                };
                if retval != 0 {
                    return ARK_PRESTEPFN_FAIL;
                }
            }

            /* Call time stepper module to attempt a step:
                  0 => step completed successfully
                 >0 => step encountered recoverable failure; reduce step if possible
                 <0 => step encountered unrecoverable failure */
            let step = ark_mem.borrow().step.expect("step set");
            kflag = step(ark_mem, &mut dsm, &mut nflag);
            if kflag < 0 {
                break;
            }

            /* handle solver convergence failures */
            kflag = arkCheckConvergence(ark_mem, &mut nflag, &mut ncf);

            if kflag < 0 {
                break;
            }

            /* Perform relaxation:
                 - computes relaxation parameter
                 - on success, updates ycur, h, and dsm
                 - on recoverable failure, updates eta and signals to retry step
                 - on fatal error, returns negative error flag */
            if ark_mem.borrow().relax_enabled && (kflag == ARK_SUCCESS) {
                kflag = arkRelax(ark_mem, &mut relax_fails, &mut dsm);

                if kflag < 0 {
                    break;
                }
            }

            /* perform constraint-handling (if selected, and if solver check passed) */
            if ark_mem.borrow().constraints.is_some() && (kflag == ARK_SUCCESS) {
                kflag = arkCheckConstraints(ark_mem, &mut constrfails, &mut nflag);

                if kflag < 0 {
                    break;
                }
            }

            /* when fixed time-stepping is enabled, 'success' == successful stage solves
               (checked in previous block), so just enforce no step size change */
            if ark_mem.borrow().fixedstep {
                ark_mem.borrow_mut().eta = ONE;
                break;
            }

            /* check temporal error (if checks above passed) */
            if kflag == ARK_SUCCESS {
                kflag = arkCheckTemporalError(ark_mem, &mut nflag, &mut nef, dsm);

                if kflag < 0 {
                    break;
                }
            }

            /* if ignoring temporal error test result (XBraid) force step to pass */
            if ark_mem.borrow().force_pass {
                ark_mem.borrow_mut().last_kflag = kflag;
                kflag = ARK_SUCCESS;
                break;
            }

            /* break attempt loop on successful step */
            if kflag == ARK_SUCCESS {
                break;
            }

            /* unsuccessful step, if |h| = hmin, return ARK_ERR_FAILURE */
            {
                let m = ark_mem.borrow();
                if SUNRabs(m.h) <= m.hmin * ONEPSM {
                    return ARK_ERR_FAILURE;
                }
            }

            /* update h, hprime and next_h for next iteration */
            {
                let mut m = ark_mem.borrow_mut();
                m.h *= m.eta;
                m.hprime = m.h;
                m.next_h = m.hprime;

                /* reset tcur to last saved internal time before reattempting step
                   (and optionally ycur to yn ) */
                m.tcur = m.tn;
            }
            {
                let (ensure_ycur, yn, ycur) = {
                    let m = ark_mem.borrow();
                    (m.ensure_ycur, m.yn.clone(), m.ycur.clone())
                };
                if ensure_ycur {
                    N_VScale(
                        ONE,
                        &yn.expect("yn allocated"),
                        &ycur.expect("ycur attached"),
                    );
                }
            }
        } /* end looping for step attempts */

        /* If step attempt loop succeeded, complete step (update current time, solution,
           error stepsize history arrays; call user-supplied step postprocessing function)
           (added stuff from arkStep_PrepareNextStep -- revisit) */
        if kflag == ARK_SUCCESS {
            kflag = arkCompleteStep(ark_mem, dsm);
        }

        /* If step attempt loop failed, process flag and return to user */
        if kflag != ARK_SUCCESS {
            istate = arkHandleFailure(ark_mem, kflag);
            let (tcur, yn) = {
                let mut m = ark_mem.borrow_mut();
                m.tretlast = m.tcur;
                (m.tcur, m.yn.clone().expect("yn allocated"))
            };
            *tret = tcur;
            N_VScale(ONE, &yn, yout);
            break;
        }

        nstloc += 1;

        /* Check for root in last step taken. */
        if ark_mem.borrow().root_mem.is_some() {
            let nrtfn = ark_mem.borrow().root_mem.as_ref().expect("root_mem").nrtfn;
            if nrtfn > 0 {
                let retval = arkRootCheck3(ark_mem, tout, itask);
                if retval == RTFOUND {
                    /* A new root was found */
                    let tlo = {
                        let mut m = ark_mem.borrow_mut();
                        let root_mem = m.root_mem.as_mut().expect("root_mem");
                        root_mem.irfnd = 1;
                        root_mem.tlo
                    };
                    istate = ARK_ROOT_RETURN;
                    ark_mem.borrow_mut().tretlast = tlo;
                    *tret = tlo;
                    break;
                } else if retval == ARK_RTFUNC_FAIL {
                    /* g failed */
                    let tlo = ark_mem.borrow().root_mem.as_ref().expect("root_mem").tlo;
                    arkProcessError(
                        Some(ark_mem),
                        ARK_RTFUNC_FAIL,
                        line!() as i32,
                        "ARKodeEvolve",
                        file!(),
                        &MSG_ARK_RTFUNC_FAILED(tlo),
                    );
                    istate = ARK_RTFUNC_FAIL;
                    break;
                }

                /* If we are at the end of the first step and we still have
                   some event functions that are inactive, issue a warning
                   as this may indicate a user error in the implementation
                   of the root function. */
                if ark_mem.borrow().nst == 1 {
                    let (inactive_roots, mxgnull) = {
                        let m = ark_mem.borrow();
                        let root_mem = m.root_mem.as_ref().expect("root_mem");
                        let mut inactive_roots = SUNFALSE;
                        for ir in 0..root_mem.nrtfn as usize {
                            if !root_mem.gactive[ir] {
                                inactive_roots = SUNTRUE;
                                break;
                            }
                        }
                        (inactive_roots, root_mem.mxgnull)
                    };
                    if (mxgnull > 0) && inactive_roots {
                        arkProcessError(
                            Some(ark_mem),
                            ARK_WARNING,
                            line!() as i32,
                            "ARKodeEvolve",
                            file!(),
                            MSG_ARK_INACTIVE_ROOTS,
                        );
                    }
                }
            }
        }

        /* Check if tn is at tstop or near tstop */
        if ark_mem.borrow().tstopset {
            let (tcur, h, hprime, tstop, tstopinterp, has_interp, troundoff) = {
                let m = ark_mem.borrow();
                (
                    m.tcur,
                    m.h,
                    m.hprime,
                    m.tstop,
                    m.tstopinterp,
                    m.interp.is_some(),
                    FUZZ_FACTOR * m.uround * (SUNRabs(m.tcur) + SUNRabs(m.h)),
                )
            };

            if SUNRabs(tcur - tstop) <= troundoff {
                /* Ensure tout >= tstop, otherwise check for tout return below */
                if (tout - tstop) * h >= ZERO || SUNRabs(tout - tstop) <= troundoff {
                    if tstopinterp && has_interp {
                        let retval = ARKodeGetDky(ark_mem, tstop, 0, yout);
                        if retval != ARK_SUCCESS {
                            arkProcessError(
                                Some(ark_mem),
                                retval,
                                line!() as i32,
                                "ARKodeEvolve",
                                file!(),
                                &MSG_ARK_INTERPOLATION_FAIL(tstop),
                            );
                            istate = retval;
                            break;
                        }
                    } else {
                        let yn = ark_mem.borrow().yn.clone().expect("yn allocated");
                        N_VScale(ONE, &yn, yout);
                    }
                    {
                        let mut m = ark_mem.borrow_mut();
                        m.tretlast = m.tstop;
                        *tret = m.tstop;
                        m.tstopset = SUNFALSE;
                    }
                    istate = ARK_TSTOP_RETURN;
                    break;
                }
            }
            /* limit upcoming step if it will overcome tstop */
            else if (tcur + hprime - tstop) * h > ZERO {
                let mut m = ark_mem.borrow_mut();
                m.hprime = (m.tstop - m.tcur) * (ONE - FOUR * m.uround);
                m.eta = m.hprime / m.h;
            }
        }

        /* In NORMAL mode, check if tout reached */
        {
            let (tcur, h, has_interp) = {
                let m = ark_mem.borrow();
                (m.tcur, m.h, m.interp.is_some())
            };
            if (itask == ARK_NORMAL) && (tcur - tout) * h >= ZERO {
                if has_interp {
                    let retval = ARKodeGetDky(ark_mem, tout, 0, yout);
                    if retval != ARK_SUCCESS {
                        arkProcessError(
                            Some(ark_mem),
                            retval,
                            line!() as i32,
                            "ARKodeEvolve",
                            file!(),
                            &MSG_ARK_INTERPOLATION_FAIL(tout),
                        );
                        istate = retval;
                        break;
                    }
                    ark_mem.borrow_mut().tretlast = tout;
                    *tret = tout;
                } else {
                    let yn = ark_mem.borrow().yn.clone().expect("yn allocated");
                    N_VScale(ONE, &yn, yout);
                    let mut m = ark_mem.borrow_mut();
                    m.tretlast = m.tcur;
                    *tret = m.tcur;
                }
                {
                    let mut m = ark_mem.borrow_mut();
                    m.next_h = m.hprime;
                }
                istate = ARK_SUCCESS;
                break;
            }
        }

        /* In ONE_STEP mode, exit loop (arkCompleteStep already copied yn to ycur, an alias to yout) */
        if itask == ARK_ONE_STEP {
            istate = ARK_SUCCESS;
            let mut m = ark_mem.borrow_mut();
            m.tretlast = m.tcur;
            *tret = m.tcur;
            m.next_h = m.hprime;
            break;
        }
    } /* end looping for internal steps */

    /* stop profiler and return: profiling disabled in the reference build */
    istate
}

/*---------------------------------------------------------------
  ARKodeGetDky:

  This routine computes the k-th derivative of the interpolating
  polynomial at the time t and stores the result in the vector
  dky. This routine internally calls arkInterpEvaluate to perform
  the interpolation.  We have the restriction that 0 <= k <= 3.
  This routine uses an interpolating polynomial of degree
  max(deg, k), i.e. it will form a polynomial of the degree
  available by the interpolation module and/or requested by
  the user through deg, unless higher-order derivatives are
  requested.

  This function is called by ARKodeEvolve with k=0 and t=tout to
  perform interpolation of outputs, but may also be called
  indirectly by the user via time step module *StepGetDky calls.
  Note: in all cases it will be called after ark_tcur has been
  updated to correspond with the end time of the last successful
  step.
  ---------------------------------------------------------------*/
pub fn ARKodeGetDky(arkode_mem: &ARKodeMem, t: sunrealtype, k: i32, dky: &N_Vector) -> i32 {
    /* Check if ark_mem exists: handled by the type system */
    let ark_mem = arkode_mem;

    /* Check all inputs for legality (NULL dky handled by the type system) */
    let interp = ark_mem.borrow().interp.clone();
    let interp = match interp {
        None => {
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_NULL,
                line!() as i32,
                "ARKodeGetDky",
                file!(),
                "Missing interpolation structure",
            );
            return ARK_MEM_NULL;
        }
        Some(interp) => interp,
    };

    /* Allow for some slack */
    let (tcur, hold, h, uround) = {
        let m = ark_mem.borrow();
        (m.tcur, m.hold, m.h, m.uround)
    };
    let mut tfuzz = FUZZ_FACTOR * uround * (SUNRabs(tcur) + SUNRabs(hold));
    if hold < ZERO {
        tfuzz = -tfuzz;
    }
    let tp = tcur - hold - tfuzz;
    let tn1 = tcur + tfuzz;
    if (t - tp) * (t - tn1) > ZERO {
        arkProcessError(
            Some(ark_mem),
            ARK_BAD_T,
            line!() as i32,
            "ARKodeGetDky",
            file!(),
            &MSG_ARK_BAD_T(t, tcur - hold, tcur),
        );
        return ARK_BAD_T;
    }

    /* call arkInterpEvaluate to evaluate result */
    let s = (t - tcur) / h;
    let retval = arkInterpEvaluate(ark_mem, &interp, s, k, ARK_INTERP_MAX_DEGREE, dky);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            retval,
            line!() as i32,
            "ARKodeGetDky",
            file!(),
            "Error calling arkInterpEvaluate",
        );
        return retval;
    }
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeFree:

  This routine frees the ARKODE infrastructure memory.
  ---------------------------------------------------------------*/
pub fn ARKodeFree(arkode_mem: &mut Option<ARKodeMem>) {
    if arkode_mem.is_none() {
        return;
    }

    let ark_mem = arkode_mem.as_ref().expect("arkode_mem").clone();

    /* free the time-stepper module memory (if provided) */
    let step_free = ark_mem.borrow().step_free;
    if let Some(step_free) = step_free {
        step_free(&ark_mem);
    }

    /* free vector storage */
    arkFreeVectors(&ark_mem);

    /* free the time step adaptivity module */
    if ark_mem.borrow().hadapt_mem.is_some() {
        let owncontroller = ark_mem
            .borrow()
            .hadapt_mem
            .as_ref()
            .expect("hadapt_mem")
            .owncontroller;
        if owncontroller {
            let hcontroller = ark_mem
                .borrow_mut()
                .hadapt_mem
                .as_mut()
                .expect("hadapt_mem")
                .hcontroller
                .take();
            let _ = SUNAdaptController_Destroy(hcontroller);
            ark_mem
                .borrow_mut()
                .hadapt_mem
                .as_mut()
                .expect("hadapt_mem")
                .owncontroller = SUNFALSE;
        }
        ark_mem.borrow_mut().hadapt_mem = None;
    }

    /* free the interpolation module */
    let interp = ark_mem.borrow().interp.clone();
    if let Some(interp) = interp {
        arkInterpFree(&ark_mem, &interp);
        ark_mem.borrow_mut().interp = None;
    }

    /* free the root-finding module */
    if ark_mem.borrow().root_mem.is_some() {
        let _ = arkRootFree(&ark_mem);
        ark_mem.borrow_mut().root_mem = None;
    }

    /* free the relaxation module */
    if ark_mem.borrow().relax_mem.is_some() {
        let relax_mem = ark_mem.borrow_mut().relax_mem.take();
        let _ = arkRelaxDestroy(relax_mem);
    }

    /* SUNDIALS_ENABLE_PYTHON is not defined:
       arkode_user_supplied_fn_table_destroy(ark_mem->python) is not called */
    ark_mem.borrow_mut().python = None;

    /* C frees the mem struct wholesale; the Rust handle is dropped by the
    caller, so break the Rc cycles the built-in ewt/rwt data tokens create
    (e_data / r_data hold an ARKodeMem clone pointing back at this record) */
    {
        let mut m = ark_mem.borrow_mut();
        m.e_data = None;
        m.r_data = None;
    }

    *arkode_mem = None;
}

/*---------------------------------------------------------------
  ARKodePrintMem:

  This routine outputs the ark_mem structure to a specified file
  pointer.
  ---------------------------------------------------------------*/
pub fn ARKodePrintMem(arkode_mem: &ARKodeMem, outfile: &SUNFile) {
    /* Check if ark_mem exists: handled by the type system */
    let ark_mem = arkode_mem;

    /* if outfile==NULL, set it to stdout */
    let stdout_ = SUNFile::Stdout;
    let outfile = if outfile.is_null() { &stdout_ } else { outfile };

    {
        let m = ark_mem.borrow();

        /* output general values */
        outfile.write_str(&format!("itol = {}\n", m.itol));
        outfile.write_str(&format!("ritol = {}\n", m.ritol));
        outfile.write_str(&format!("mxhnil = {}\n", m.mxhnil));
        outfile.write_str(&format!("mxstep = {}\n", m.mxstep));
        outfile.write_str(&format!("lrw1 = {}\n", m.lrw1));
        outfile.write_str(&format!("liw1 = {}\n", m.liw1));
        outfile.write_str(&format!("lrw = {}\n", m.lrw));
        outfile.write_str(&format!("liw = {}\n", m.liw));
        outfile.write_str(&format!("user_efun = {}\n", m.user_efun as i32));
        outfile.write_str(&format!("tstopset = {}\n", m.tstopset as i32));
        outfile.write_str(&format!("tstopinterp = {}\n", m.tstopinterp as i32));
        outfile.write_str(&format!("tstop = {}\n", sun_format_g(m.tstop)));
        outfile.write_str(&format!(
            "VabstolMallocDone = {}\n",
            m.VabstolMallocDone as i32
        ));
        outfile.write_str(&format!("MallocDone = {}\n", m.MallocDone as i32));
        outfile.write_str(&format!("initsetup = {}\n", m.initsetup as i32));
        outfile.write_str(&format!("init_type = {}\n", m.init_type));
        outfile.write_str(&format!("firststage = {}\n", m.firststage as i32));
        outfile.write_str(&format!("uround = {}\n", sun_format_g(m.uround)));
        outfile.write_str(&format!("reltol = {}\n", sun_format_g(m.reltol)));
        outfile.write_str(&format!("Sabstol = {}\n", sun_format_g(m.Sabstol)));
        outfile.write_str(&format!("fixedstep = {}\n", m.fixedstep as i32));
        outfile.write_str(&format!("tolsf = {}\n", sun_format_g(m.tolsf)));
        outfile.write_str(&format!("call_fullrhs = {}\n", m.call_fullrhs as i32));
        outfile.write_str(&format!("do_adjoint = {}\n", m.do_adjoint as i32));
        outfile.write_str(&format!("ensure_ycur = {}\n", m.ensure_ycur as i32));

        /* output counters */
        outfile.write_str(&format!("nhnil = {}\n", m.nhnil));
        outfile.write_str(&format!("nst_attempts = {}\n", m.nst_attempts));
        outfile.write_str(&format!("nst = {}\n", m.nst));
        outfile.write_str(&format!("ncfn = {}\n", m.ncfn));
        outfile.write_str(&format!("netf = {}\n", m.netf));

        /* output time-stepping values */
        outfile.write_str(&format!("hin = {}\n", sun_format_g(m.hin)));
        outfile.write_str(&format!("h = {}\n", sun_format_g(m.h)));
        outfile.write_str(&format!("hprime = {}\n", sun_format_g(m.hprime)));
        outfile.write_str(&format!("next_h = {}\n", sun_format_g(m.next_h)));
        outfile.write_str(&format!("eta = {}\n", sun_format_g(m.eta)));
        outfile.write_str(&format!("tcur = {}\n", sun_format_g(m.tcur)));
        outfile.write_str(&format!("tretlast = {}\n", sun_format_g(m.tretlast)));
        outfile.write_str(&format!("hmin = {}\n", sun_format_g(m.hmin)));
        outfile.write_str(&format!("hmax_inv = {}\n", sun_format_g(m.hmax_inv)));
        outfile.write_str(&format!("h0u = {}\n", sun_format_g(m.h0u)));
        outfile.write_str(&format!("tn = {}\n", sun_format_g(m.tn)));
        outfile.write_str(&format!("hold = {}\n", sun_format_g(m.hold)));
        outfile.write_str(&format!("maxnef = {}\n", m.maxnef));
        outfile.write_str(&format!("maxncf = {}\n", m.maxncf));

        /* output time-stepping adaptivity structure */
        outfile.write_str("timestep adaptivity structure:\n");
        arkPrintAdaptMem(m.hadapt_mem.as_deref(), outfile);

        /* output inequality constraints quantities */
        outfile.write_str(&format!("maxconstrfails = {}\n", m.maxconstrfails));
    }

    /* output root-finding quantities */
    if ark_mem.borrow().root_mem.is_some() {
        let _ = arkPrintRootMem(ark_mem, outfile);
    }

    /* output interpolation quantities */
    let interp = ark_mem.borrow().interp.clone();
    if let Some(interp) = interp {
        arkInterpPrintMem(&interp, outfile);
    } else {
        outfile.write_str("interpolation = NULL\n");
    }

    /* SUNDIALS_DEBUG_PRINTVEC is not defined: the vector dump is omitted */

    /* Call stepper PrintMem function (if provided) */
    let step_printmem = ark_mem.borrow().step_printmem;
    if let Some(step_printmem) = step_printmem {
        step_printmem(ark_mem, outfile);
    }
}

/*------------------------------------------------------------------------------
  ARKodeCreateMRIStepInnerStepper

  Wraps an ARKODE integrator as an MRIStep inner stepper.
  ----------------------------------------------------------------------------*/

pub fn ARKodeCreateMRIStepInnerStepper(
    inner_arkode_mem: &ARKodeMem,
    stepper: &mut Option<MRIStepInnerStepper>,
) -> i32 {
    /* Check if ark_mem exists: handled by the type system */
    let ark_mem = inner_arkode_mem;

    /* return with an error if the ARKODE solver does not support forcing */
    if ark_mem.borrow().step_setforcing.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!() as i32,
            "ARKodeCreateMRIStepInnerStepper",
            file!(),
            "time-stepping module does not support forcing",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    let sunctx = ark_mem.borrow().sunctx.clone();
    let retval = MRIStepInnerStepper_Create(&sunctx, stepper);
    if retval != ARK_SUCCESS {
        return retval;
    }
    let inner = stepper.as_ref().expect("stepper created").clone();

    let retval = MRIStepInnerStepper_SetContent(&inner, Box::new(ark_mem.clone()));
    if retval != ARK_SUCCESS {
        return retval;
    }

    let retval = MRIStepInnerStepper_SetEvolveFn(&inner, Some(ark_MRIStepInnerEvolve));
    if retval != ARK_SUCCESS {
        return retval;
    }

    let retval = MRIStepInnerStepper_SetFullRhsFn(&inner, Some(ark_MRIStepInnerFullRhs));
    if retval != ARK_SUCCESS {
        return retval;
    }

    let retval = MRIStepInnerStepper_SetResetFn(&inner, Some(ark_MRIStepInnerReset));
    if retval != ARK_SUCCESS {
        return retval;
    }

    let retval = MRIStepInnerStepper_SetAccumulatedErrorGetFn(
        &inner,
        Some(ark_MRIStepInnerGetAccumulatedError),
    );
    if retval != ARK_SUCCESS {
        return retval;
    }

    let retval = MRIStepInnerStepper_SetAccumulatedErrorResetFn(
        &inner,
        Some(ark_MRIStepInnerResetAccumulatedError),
    );
    if retval != ARK_SUCCESS {
        return retval;
    }

    let retval = MRIStepInnerStepper_SetRTolFn(&inner, Some(ark_MRIStepInnerSetRTol));
    if retval != ARK_SUCCESS {
        return retval;
    }

    ARK_SUCCESS
}

/*===============================================================
  Private internal functions
  ===============================================================*/

/*---------------------------------------------------------------
  arkCreate:

  arkCreate creates an internal memory block for a problem to
  be solved by a time step module built on ARKODE.  If successful,
  arkCreate returns a pointer to the problem memory. If an
  initialization error occurs, arkCreate prints an error message
  to standard err and returns NULL.
  ---------------------------------------------------------------*/
pub fn arkCreate(sunctx: &SUNContext) -> Option<ARKodeMem> {
    /* NULL sunctx check: handled by the type system */

    /* malloc failure branch: allocation cannot fail observably in Rust.
    `ARKodeMemRec::zeroed` is C's malloc + memset(ark_mem, 0, ...) and also
    sets the context. */
    let ark_mem: ARKodeMem = Rc::new(RefCell::new(ARKodeMemRec::zeroed(sunctx.clone())));

    {
        let mut m = ark_mem.borrow_mut();

        /* Set the Python context to NULL */
        m.python = None;

        /* Set uround */
        m.uround = SUN_UNIT_ROUNDOFF;

        /* Initialize time step module to NULL */
        m.step_attachlinsol = None;
        m.step_attachmasssol = None;
        m.step_disablelsetup = None;
        m.step_disablemsetup = None;
        m.step_getlinmem = None;
        m.step_getmassmem = None;
        m.step_getimplicitrhs = None;
        m.step_mmult = None;
        m.step_getgammas = None;
        m.step_init = None;
        m.step_fullrhs = None;
        m.step = None;
        m.step_setuserdata = None;
        m.step_printallstats = None;
        m.step_writeparameters = None;
        m.step_resize = None;
        m.step_reset = None;
        m.step_free = None;
        m.step_printmem = None;
        m.step_setdefaults = None;
        m.step_computestate = None;
        m.step_setrelaxfn = None;
        m.step_setorder = None;
        m.step_setnonlinearsolver = None;
        m.step_setlinear = None;
        m.step_setnonlinear = None;
        m.step_setautonomous = None;
        m.step_setnlsrhsfn = None;
        m.step_setdeduceimplicitrhs = None;
        m.step_setnonlincrdown = None;
        m.step_setnonlinrdiv = None;
        m.step_setdeltagammamax = None;
        m.step_setlsetupfrequency = None;
        m.step_setpredictormethod = None;
        m.step_setmaxnonliniters = None;
        m.step_setnonlinconvcoef = None;
        m.step_setstagepredictfn = None;
        m.step_getnumrhsevals = None;
        m.step_setstepdirection = None;
        m.step_setoptions = None;
        m.step_getnumlinsolvsetups = None;
        m.step_H0 = None;
        m.step_setadaptcontroller = None;
        m.step_getestlocalerrors = None;
        m.step_getcurrentgamma = None;
        m.step_getnonlinearsystemdata = None;
        m.step_getnumnonlinsolviters = None;
        m.step_getnumnonlinsolvconvfails = None;
        m.step_getnonlinsolvstats = None;
        m.step_getstageindex = None;
        m.step_setforcing = None;
        m.step_mem = None;
        m.step_supports_adaptive = SUNFALSE;
        m.step_supports_implicit = SUNFALSE;
        m.step_supports_massmatrix = SUNFALSE;
        m.step_supports_relaxation = SUNFALSE;

        /* Initialize root finding variables */
        m.root_mem = None;

        /* Initialize inequality constraints variables */
        m.constraints = None;

        /* Initialize relaxation variables */
        m.relax_enabled = SUNFALSE;
        m.relax_mem = None;

        /* Initialize lrw and liw */
        m.lrw = 18;
        m.liw = 53; /* fcn/data ptr, int, long int, sunindextype, sunbooleantype */

        /* No mallocs have been done yet */
        m.VabstolMallocDone = SUNFALSE;
        m.VRabstolMallocDone = SUNFALSE;
        m.MallocDone = SUNFALSE;

        /* No user-supplied pre- or post-step functions yet */
        m.PreStepFn = None;
        m.PostStepFn = None;

        /* No user-supplied pre-RHS function yet */
        m.PreRhsFn = None;

        /* No user-supplied stage/step post-processing functions yet */
        m.PostProcessStepFn = None;
        m.PostProcessStageFn = None;

        /* No user_data pointer yet */
        m.user_data = None;
    }

    /* Allocate step adaptivity structure and note storage */
    let hadapt_mem = arkAdaptInit();
    if hadapt_mem.is_none() {
        arkProcessError(
            None,
            ARK_MEM_FAIL,
            line!() as i32,
            "arkCreate",
            file!(),
            "Allocation of step adaptivity structure failed",
        );
        let mut ark_mem = Some(ark_mem);
        ARKodeFree(&mut ark_mem);
        return None;
    }
    {
        let mut m = ark_mem.borrow_mut();
        m.hadapt_mem = hadapt_mem;
        m.lrw += ARK_ADAPT_LRW;
        m.liw += ARK_ADAPT_LIW;

        /* Initialize the interpolation structure to NULL */
        m.interp = None;
        m.interp_type = ARK_INTERP_HERMITE;
        m.interp_degree = ARK_INTERP_MAX_DEGREE;

        /* Initially, rwt should point to ewt */
        m.rwt_is_ewt = SUNTRUE;

        /* Indicate that calling the full RHS function is not required, this flag is
           updated to SUNTRUE by the interpolation module initialization function
           and/or the stepper initialization function in arkInitialSetup */
        m.call_fullrhs = SUNFALSE;

        /* Indicate that the problem needs to be initialized */
        m.initsetup = SUNTRUE;
        m.init_type = FIRST_INIT;
        m.firststage = SUNTRUE;
        m.initialized = SUNFALSE;

        /* Initial step size has not been determined yet */
        m.h = ZERO;
        m.h0u = ZERO;

        /* Accumulated error estimation strategy */
        m.AccumErrorType = ARK_ACCUMERROR_NONE;
        m.AccumError = ZERO;

        /* Default to having stepper initialize ycur during evolution */
        m.ensure_ycur = SUNFALSE;
    }

    /* Set default values for integrator and stepper optional inputs */
    let iret = ARKodeSetDefaults(&ark_mem);
    if iret != ARK_SUCCESS {
        arkProcessError(
            None,
            0,
            line!() as i32,
            "arkCreate",
            file!(),
            "Error setting default solver options",
        );
        let mut ark_mem = Some(ark_mem);
        ARKodeFree(&mut ark_mem);
        return None;
    }

    {
        let mut m = ark_mem.borrow_mut();
        m.load_checkpoint_fail = SUNFALSE;
        m.do_adjoint = SUNFALSE;
    }

    /* Return pointer to ARKODE memory block */
    Some(ark_mem)
}

/*---------------------------------------------------------------
  arkRwtSet

  This routine is responsible for setting the residual weight
  vector rwt, according to tol_type, as follows:

  (1) rwt[i] = 1 / (reltol * SUNRabs(M*ycur[i]) + rabstol), i=0,...,neq-1
      if tol_type = ARK_SS
  (2) rwt[i] = 1 / (reltol * SUNRabs(M*ycur[i]) + rabstol[i]), i=0,...,neq-1
      if tol_type = ARK_SV
  (3) unset if tol_type is any other value (occurs rwt=ewt)

  arkRwtSet returns 0 if rwt is successfully set as above to a
  positive vector and -1 otherwise. In the latter case, rwt is
  considered undefined.

  All the real work is done in the routines arkRwtSetSS, arkRwtSetSV.
  ---------------------------------------------------------------*/
pub fn arkRwtSet(y: &N_Vector, weight: &N_Vector, data: &mut Option<Box<dyn Any>>) -> i32 {
    /* data points to ark_mem here (a boxed ARKodeMem handle clone; C's cast
    of a NULL/foreign pointer is UB -> deterministic panic) */
    let ark_mem = data
        .as_mut()
        .and_then(|b| b.downcast_ref::<ARKodeMem>())
        .cloned()
        .expect("arkRwtSet data holds ARKodeMem");

    let mut flag: i32 = 0;

    /* return if rwt is just ewt */
    if ark_mem.borrow().rwt_is_ewt {
        return 0;
    }

    /* put M*y into ark_tempv1 */
    let My = ark_mem.borrow().tempv1.clone().expect("tempv1 allocated");
    let step_mmult = ark_mem.borrow().step_mmult;
    if let Some(step_mmult) = step_mmult {
        flag = step_mmult(&ark_mem, y, &My);
        if flag != ARK_SUCCESS {
            return ARK_MASSMULT_FAIL;
        }
    } else {
        /* this condition should not apply, but just in case */
        N_VScale(ONE, y, &My);
    }

    /* call appropriate routine to fill rwt */
    let ritol = ark_mem.borrow().ritol;
    match ritol {
        ARK_SS => flag = arkRwtSetSS(&ark_mem, &My, weight),
        ARK_SV => flag = arkRwtSetSV(&ark_mem, &My, weight),
        _ => {}
    }

    flag
}

/*---------------------------------------------------------------
  arkInit:

  arkInit allocates and initializes memory for a problem. All
  inputs are checked for errors. If any error occurs during
  initialization, an error flag is returned. Otherwise, it
  returns ARK_SUCCESS.

  This routine should only be called by
  (a) ARKodeReset (with the input init_type == RESET_INIT),
  (b) an ARKODE timestepper module creation routine (with
      init_type == FIRST_INIT), or
  (c) an ARKODE timestepper module re-initialization routine
      (with init_type == FIRST_INIT).
  This should never be called by the user.

  The initialization type indicates if the values of internal
  counters should be reinitialized (FIRST_INIT) or retained
  (RESET_INIT).

  This routine must be called prior to calling ARKodeEvolve
  to evolve the problem.
  ---------------------------------------------------------------*/
pub fn arkInit(ark_mem: &ARKodeMem, t0: sunrealtype, y0: &N_Vector, init_type: i32) -> i32 {
    let mut init_type = init_type;

    /* Check ark_mem: NULL-mem check handled by the type system */

    /* Check for legal input parameters (NULL y0 handled by the type
    system; the Rc clone aliases the caller's vector) */
    ark_mem.borrow_mut().ycur = Some(y0.clone());

    /* Check if reset was called before the first Evolve call */
    if init_type == RESET_INIT && !ark_mem.borrow().initialized {
        init_type = FIRST_INIT;
    }

    /* Check if allocations have been done i.e., is this first init call */
    if !ark_mem.borrow().MallocDone {
        /* Test if all required time stepper operations are implemented */
        let stepperOK = arkCheckTimestepper(ark_mem);
        if !stepperOK {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!() as i32,
                "arkInit",
                file!(),
                "Time stepper module is missing required functionality",
            );
            return ARK_ILL_INPUT;
        }

        /* Test if all required vector operations are implemented */
        let nvectorOK = arkCheckNvectorRequired(y0);
        if !nvectorOK {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!() as i32,
                "arkInit",
                file!(),
                MSG_ARK_BAD_NVECTOR,
            );
            return ARK_ILL_INPUT;
        }

        /* Set space requirements for one N_Vector */
        let mut lrw1: sunindextype = 0;
        let mut liw1: sunindextype = 0;
        if y0.ops.borrow().nvspace.is_some() {
            N_VSpace(y0, &mut lrw1, &mut liw1);
        } else {
            lrw1 = 0;
            liw1 = 0;
        }
        {
            let mut m = ark_mem.borrow_mut();
            m.lrw1 = lrw1;
            m.liw1 = liw1;
        }

        /* Allocate the solver vectors (using y0 as a template) */
        let allocOK = arkAllocVectors(ark_mem, y0);
        if !allocOK {
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_FAIL,
                line!() as i32,
                "arkInit",
                file!(),
                MSG_ARK_MEM_FAIL,
            );
            return ARK_MEM_FAIL;
        }

        /* All allocations are complete */
        ark_mem.borrow_mut().MallocDone = SUNTRUE;
    }

    /* All allocation and error checking is complete at this point */

    /* Copy the input parameters into ARKODE state */
    {
        let mut m = ark_mem.borrow_mut();
        m.tcur = t0;
        m.tn = t0;
    }

    /* Initialize yn */
    let yn = ark_mem.borrow().yn.clone().expect("yn allocated");
    N_VScale(ONE, y0, &yn);
    {
        let mut m = ark_mem.borrow_mut();
        m.fn_is_current = SUNFALSE;

        /* Clear any previous 'tstop' */
        m.tstopset = SUNFALSE;
    }

    /* Initializations on (re-)initialization call, skip on reset */
    if init_type == FIRST_INIT {
        {
            let mut m = ark_mem.borrow_mut();

            /* Counters */
            m.nst_attempts = 0;
            m.nst = 0;
            m.nhnil = 0;
            m.ncfn = 0;
            m.netf = 0;
            m.nconstrfails = 0;

            /* Initial, old, and next step sizes */
            m.h0u = ZERO;
            m.hold = ZERO;
            m.next_h = ZERO;

            /* Tolerance scale factor */
            m.tolsf = ONE;
        }

        /* Reset error controller object */
        let hcontroller = ark_mem
            .borrow()
            .hadapt_mem
            .as_ref()
            .expect("hadapt_mem")
            .hcontroller
            .clone();
        if let Some(hcontroller) = hcontroller {
            let retval = SUNAdaptController_Reset(&hcontroller);
            if retval != SUN_SUCCESS {
                arkProcessError(
                    Some(ark_mem),
                    ARK_CONTROLLER_ERR,
                    line!() as i32,
                    "arkInit",
                    file!(),
                    "Unable to reset error controller object",
                );
                return ARK_CONTROLLER_ERR;
            }
        }

        let mut m = ark_mem.borrow_mut();

        /* Adaptivity counters */
        {
            let hadapt_mem = m.hadapt_mem.as_mut().expect("hadapt_mem");
            hadapt_mem.nst_acc = 0;
            hadapt_mem.nst_exp = 0;
        }

        /* Accumulated error estimate */
        m.AccumError = ZERO;

        /* Indicate that calling the full RHS function is not required, this flag is
           updated to SUNTRUE by the interpolation module initialization function
           and/or the stepper initialization function in arkInitialSetup */
        m.call_fullrhs = SUNFALSE;

        /* Adjoint related */
        m.checkpoint_step_idx = 0;

        /* Indicate that initialization has not been done before */
        m.initialized = SUNFALSE;
    }

    /* Indicate initialization is needed */
    {
        let mut m = ark_mem.borrow_mut();
        m.initsetup = SUNTRUE;
        m.init_type = init_type;
        m.firststage = SUNTRUE;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkCheckTimestepper:

  This routine checks if all required time stepper function
  pointers have been supplied.  If any of them is missing it
  returns SUNFALSE.
  ---------------------------------------------------------------*/
pub fn arkCheckTimestepper(ark_mem: &ARKodeMem) -> sunbooleantype {
    let m = ark_mem.borrow();
    if m.step_init.is_none() || m.step.is_none() || m.step_mem.is_none() {
        return SUNFALSE;
    }
    SUNTRUE
}

/*---------------------------------------------------------------
  arkCheckNvectorRequired:

  This routine checks if all absolutely-required vector
  operations are present.  If any of them is missing it returns
  SUNFALSE.
  ---------------------------------------------------------------*/
pub fn arkCheckNvectorRequired(tmpl: &N_Vector) -> sunbooleantype {
    let ops = tmpl.ops.borrow();
    if ops.nvclone.is_none()
        || ops.nvdestroy.is_none()
        || ops.nvlinearsum.is_none()
        || ops.nvconst.is_none()
        || ops.nvdiv.is_none()
        || ops.nvscale.is_none()
        || ops.nvabs.is_none()
        || ops.nvinv.is_none()
        || ops.nvwrmsnorm.is_none()
    {
        SUNFALSE
    } else {
        SUNTRUE
    }
}

/*---------------------------------------------------------------
  arkCheckNvectorOptional:

  This routine perform conditional checks on required vector
  operations are present (i.e., if the current ARKODE
  configuration requires additional N_Vector routines).  If any
  of them is missing it returns SUNFALSE.
  ---------------------------------------------------------------*/
pub fn arkCheckNvectorOptional(ark_mem: &ARKodeMem) -> sunbooleantype {
    let (user_efun, atolmin0, user_rfun, rwt_is_ewt, Ratolmin0, h0u, hin, itol, ritol, tempv1) = {
        let m = ark_mem.borrow();
        (
            m.user_efun,
            m.atolmin0,
            m.user_rfun,
            m.rwt_is_ewt,
            m.Ratolmin0,
            m.h0u,
            m.hin,
            m.itol,
            m.ritol,
            m.tempv1.clone().expect("tempv1 allocated"),
        )
    };
    let (has_nvmin, has_nvdiv, has_nvmaxnorm, has_nvaddconst) = {
        let ops = tempv1.ops.borrow();
        (
            ops.nvmin.is_some(),
            ops.nvdiv.is_some(),
            ops.nvmaxnorm.is_some(),
            ops.nvaddconst.is_some(),
        )
    };

    /* If using a built-in routine for error/residual weights with abstol==0,
       ensure that N_VMin is available */
    if !user_efun && atolmin0 && !has_nvmin {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!() as i32,
            "arkCheckNvectorOptional",
            file!(),
            "N_VMin unimplemented (required by error-weight function)",
        );
        return SUNFALSE;
    }
    if !user_rfun && !rwt_is_ewt && Ratolmin0 && !has_nvmin {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!() as i32,
            "arkCheckNvectorOptional",
            file!(),
            "N_VMin unimplemented (required by residual-weight function)",
        );
        return SUNFALSE;
    }

    /* If the user has not specified a step size (and it will be estimated
       internally), ensure that N_VDiv and N_VMaxNorm are available */
    if (h0u == ZERO) && (hin == ZERO) && !has_nvdiv {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!() as i32,
            "arkCheckNvectorOptional",
            file!(),
            "N_VDiv unimplemented (required for initial step estimation)",
        );
        return SUNFALSE;
    }
    if (h0u == ZERO) && (hin == ZERO) && !has_nvmaxnorm {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!() as i32,
            "arkCheckNvectorOptional",
            file!(),
            "N_VMaxNorm unimplemented (required for initial step estimation)",
        );
        return SUNFALSE;
    }

    /* If using a scalar-valued absolute tolerance (for either the state or
       residual), then ensure that N_VAddConst is available */
    if (itol == ARK_SS) && !has_nvaddconst {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!() as i32,
            "arkCheckNvectorOptional",
            file!(),
            "N_VAddConst unimplemented (required for scalar abstol)",
        );
        return SUNFALSE;
    }
    if !rwt_is_ewt && (ritol == ARK_SS) && !has_nvaddconst {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!() as i32,
            "arkCheckNvectorOptional",
            file!(),
            "N_VAddConst unimplemented (required for scalar rabstol)",
        );
        return SUNFALSE;
    }

    /* If we made it here, then the vector is sufficient */
    SUNTRUE
}
