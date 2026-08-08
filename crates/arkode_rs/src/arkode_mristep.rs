//! Port of `src/arkode/arkode_mristep.c` with `src/arkode/arkode_mristep_impl.h`
//! folded in (module = C file base name; the impl header's data structures,
//! constants and error messages live here and every other MRIStep part —
//! `arkode_mristep_io.rs`, `arkode_mristep_nls.rs`,
//! `arkode_mristep_controller.rs` — reaches them through
//! `use crate::arkode_mristep::*;`).
//!
//! Reference build configuration: SUNDIALS_LOGGING_LEVEL = 2, so every
//! `SUNLogInfo` / `SUNLogInfoIf` / `SUNLogDebug` / `SUNLogExtraDebug*` call
//! site is omitted at translation time (MRIStep queues no `ARK_WARNING`
//! messages). Profiling off, error checks off (`SUNAssert` / `SUNCheck*` are
//! no-ops), monitoring on, serial branches only.
//!
//! Binding notes specific to this module:
//!  * The MRIStep content record lives BY VALUE in `ark_mem.step_mem`
//!    (`Option<Box<dyn Any>>` = C `void* step_mem`) and is reached through
//!    the single downcast helper [`mriStep_mem_mut`]. That guard IS a borrow
//!    of `ark_mem`: it is never held across `arkProcessError`, a user
//!    callback, an `N_Vector` / linear-solver / nonlinear-solver operation,
//!    an inner-stepper call, or a second borrow of the same mem.
//!  * C `step_mem->lmem` lives in `ark_mem.ark_lmem` (contract §4), so
//!    `mriStep_GetLmem` is a presence probe and `mriStep_AttachLinsol` moves
//!    the ARKLS record into `ark_mem.ark_lmem`.
//!  * C `step_mem->jcur` is the shared [`ARKJcurPtr`] cell so that a
//!    preconditioner-setup routine reached re-entrantly through
//!    `arkLsSetup` writes through the same flag `mriStep_GetGammas` handed
//!    out (contract §"THE jcur SEAM").
//!  * The fused-op scratch arrays keep the C shape: `cvals` is a
//!    `Vec<sunrealtype>` sized `nfusedopvecs` and `Xvecs` a
//!    `Vec<Option<N_Vector>>` of the same length (C `calloc` leaves NULL
//!    slots, which map to `None`). [`mriStep_xvecs`] materialises the dense
//!    `&[N_Vector]` the fused vector kernels take.
//!  * Rust-forced renames: `step_mem->crate` is `crate_` and
//!    `MRIC->type` is `type_` (both are Rust keywords); `ark_mem->fn` is
//!    `ark_mem.fn_` (contract).

use std::any::Any;
use std::cell::{Cell, RefCell, RefMut};
use std::rc::Rc;

use crate::arkode::*;
use crate::arkode_impl::*;
use crate::arkode_io::*;
use crate::arkode_ls::*;
use crate::arkode_mri_tables::*;
use crate::arkode_mristep_io::*;
use crate::arkode_mristep_nls::*;
use sundials_core::sundials_adaptcontroller::*;
use sundials_core::sundials_context::SUNContext;
use sundials_core::sundials_errors::*;
use sundials_core::sundials_linearsolver::SUNLinearSolver_Type;
use sundials_core::sundials_math::*;
use sundials_core::sundials_nonlinearsolver::*;
use sundials_core::sundials_nvector::*;
use sundials_core::sundials_stepper::*;
use sundials_core::sundials_types::*;
use sundials_core::sundials_utils::*;
use sundials_core::sunnonlinsol_newton::SUNNonlinSol_Newton;

/*===============================================================
  MRIStep constants (arkode_mristep_impl.h)
  ===============================================================*/

/* Stage type identifiers */
pub const MRISTAGE_FIRST: i32 = -2;
pub const MRISTAGE_STIFF_ACC: i32 = -1;
pub const MRISTAGE_ERK_FAST: i32 = 0;
pub const MRISTAGE_ERK_NOFAST: i32 = 1;
pub const MRISTAGE_DIRK_NOFAST: i32 = 2;
pub const MRISTAGE_DIRK_FAST: i32 = 3;

/* The implicit-solver constants MAXCOR / CRDOWN / DGMAX / RDIV / MSBP /
   NLSCOEF are byte-identical duplicates of the ARKStep ones and live in
   `arkode_impl.rs` (contract §7). */

/*===============================================================
  Reusable MRIStep Error Messages (arkode_mristep_impl.h)
  ===============================================================*/

pub const MSG_MRISTEP_NO_MEM: &str = "Time step module memory is NULL.";
pub const MSG_NLS_INIT_FAIL: &str = "The nonlinear solver's init routine failed.";
pub const MSG_MRISTEP_NO_COUPLING: &str = "The MRIStepCoupling is NULL.";

/*===============================================================
  MRIStep user-supplied function types (arkode/arkode_mristep.h)
  ===============================================================*/

pub type MRIStepPreInnerFn = fn(
    t: sunrealtype,
    f_1d: &[N_Vector],
    nvecs: i32,
    user_data: &mut Option<Box<dyn Any>>,
) -> i32;

pub type MRIStepPostInnerFn =
    fn(t: sunrealtype, y: &N_Vector, user_data: &mut Option<Box<dyn Any>>) -> i32;

/*===============================================================
  MRIStep inner-stepper function types (arkode/arkode_mristep.h)
  ===============================================================*/

pub type MRIStepInnerEvolveFn = fn(
    stepper: &MRIStepInnerStepper,
    t0: sunrealtype,
    tout: sunrealtype,
    y: &N_Vector,
) -> i32;

pub type MRIStepInnerFullRhsFn = fn(
    stepper: &MRIStepInnerStepper,
    t: sunrealtype,
    y: &N_Vector,
    f: &N_Vector,
    mode: i32,
) -> i32;

pub type MRIStepInnerResetFn =
    fn(stepper: &MRIStepInnerStepper, tR: sunrealtype, yR: &N_Vector) -> i32;

pub type MRIStepInnerGetAccumulatedError =
    fn(stepper: &MRIStepInnerStepper, accum_error: &mut sunrealtype) -> i32;

pub type MRIStepInnerResetAccumulatedError = fn(stepper: &MRIStepInnerStepper) -> i32;

pub type MRIStepInnerSetRTol = fn(stepper: &MRIStepInnerStepper, rtol: sunrealtype) -> i32;

/*===============================================================
  MRI inner time stepper data structure (arkode_mristep_impl.h)
  ===============================================================*/

#[derive(Default, Clone)]
pub struct _MRIStepInnerStepper_Ops {
    pub evolve: Option<MRIStepInnerEvolveFn>,
    pub fullrhs: Option<MRIStepInnerFullRhsFn>,
    pub reset: Option<MRIStepInnerResetFn>,
    pub geterror: Option<MRIStepInnerGetAccumulatedError>,
    pub reseterror: Option<MRIStepInnerResetAccumulatedError>,
    pub setrtol: Option<MRIStepInnerSetRTol>,
}

pub type MRIStepInnerStepper_Ops = _MRIStepInnerStepper_Ops;

/// C `struct _MRIStepInnerStepper`. Handle model: `Rc` clone = C pointer
/// copy, `Rc::ptr_eq` = C pointer equality; every mutable field carries its
/// own `Cell`/`RefCell` so a holder of the handle can write through it
/// exactly as C writes through the pointer.
///
/// `content` is C `void* content`: every in-tree producer stores a SUNDIALS
/// handle there (`ARKodeMem` from `ARKodeCreateMRIStepInnerStepper`,
/// `SUNStepper` from `MRIStepInnerStepper_CreateFromSUNStepper`), so it is an
/// `Rc<dyn Any>` that `MRIStepInnerStepper_GetContent` hands out by clone
/// (a `Box` could not be handed out without moving it).
pub struct _MRIStepInnerStepper {
    /* stepper specific content and operations */
    pub content: RefCell<Option<Rc<dyn Any>>>,
    /// C `void* python` (Python bindings are out of scope); always `None`.
    pub python: RefCell<Option<Box<dyn Any>>>,
    pub ops: RefCell<MRIStepInnerStepper_Ops>,

    /* stepper context */
    pub sunctx: RefCell<SUNContext>,

    /* base class data */
    pub forcing: RefCell<Vec<N_Vector>>, /* array of forcing vectors            */
    pub nforcing: Cell<i32>,             /* number of forcing vectors active    */
    pub nforcing_allocated: Cell<i32>,   /* number of forcing vectors allocated */
    pub last_flag: Cell<i32>,            /* last stepper return flag            */
    pub tshift: Cell<sunrealtype>,       /* time normalization shift            */
    pub tscale: Cell<sunrealtype>,       /* time normalization scaling          */

    /* fused op workspace */
    pub vals: RefCell<Vec<sunrealtype>>,
    pub vecs: RefCell<Vec<Option<N_Vector>>>,

    /* Space requirements */
    pub lrw1: Cell<sunindextype>,
    pub liw1: Cell<sunindextype>,
    pub lrw: Cell<i64>,
    pub liw: Cell<i64>,
}

pub type MRIStepInnerStepper = Rc<_MRIStepInnerStepper>;

impl _MRIStepInnerStepper {
    /// C `malloc` + `memset(*stepper, 0, sizeof(**stepper))` in
    /// `MRIStepInnerStepper_Create` (which then assigns `last_flag`,
    /// `sunctx` and `python`).
    pub fn zeroed(sunctx: SUNContext) -> _MRIStepInnerStepper {
        _MRIStepInnerStepper {
            content: RefCell::new(None),
            python: RefCell::new(None),
            ops: RefCell::new(MRIStepInnerStepper_Ops::default()),
            sunctx: RefCell::new(sunctx),
            forcing: RefCell::new(Vec::new()),
            nforcing: Cell::new(0),
            nforcing_allocated: Cell::new(0),
            last_flag: Cell::new(0),
            tshift: Cell::new(ZERO),
            tscale: Cell::new(ZERO),
            vals: RefCell::new(Vec::new()),
            vecs: RefCell::new(Vec::new()),
            lrw1: Cell::new(0),
            liw1: Cell::new(0),
            lrw: Cell::new(0),
            liw: Cell::new(0),
        }
    }
}

/*===============================================================
  MRI time step module data structure (arkode_mristep_impl.h)
  ===============================================================*/

/// C `struct ARKodeMRIStepMemRec`, held BY VALUE in `ark_mem.step_mem`.
///
/// Deviations from the C layout, all forced and all documented at their use
/// sites: `crate` -> `crate_` (Rust keyword); `void* lmem` is gone (the ARKLS
/// record lives in `ark_mem.ark_lmem`, contract §4); `sunbooleantype jcur` is
/// the shared [`ARKJcurPtr`] cell; the `N_Vector*` arrays `Fse` / `Fsi` /
/// `forcing` are `Vec<N_Vector>` (empty == C `NULL`); the `int*` /
/// `sunrealtype*` arrays are `Vec` (empty == C `NULL`); `Xvecs` is
/// `Vec<Option<N_Vector>>` (C `calloc`'d NULL slots).
pub struct ARKodeMRIStepMemRec {
    /* MRI problem specification */
    pub fse: Option<ARKRhsFn>, /* y' = fse(t,y) + fsi(t,y) + ff(t,y) */
    pub fsi: Option<ARKRhsFn>,
    pub linear: sunbooleantype,         /* SUNTRUE if fi is linear        */
    pub linear_timedep: sunbooleantype, /* SUNTRUE if dfi/dy depends on t */
    pub explicit_rhs: sunbooleantype,   /* SUNTRUE if fse is provided     */
    pub implicit_rhs: sunbooleantype,   /* SUNTRUE if fsi is provided     */
    pub deduce_rhs: sunbooleantype,     /* SUNTRUE if fi is deduced after
                                        a nonlinear solve              */

    /* Outer RK method storage and parameters */
    pub Fse: Vec<N_Vector>,       /* explicit RHS at each stage               */
    pub Fsi: Vec<N_Vector>,       /* implicit RHS at each stage               */
    pub unify_Fs: sunbooleantype, /* Fse and Fsi point at the same memory     */
    pub fse_is_current: sunbooleantype,
    pub fsi_is_current: sunbooleantype,
    pub MRIC: Option<MRIStepCoupling>, /* slow->fast coupling table           */
    pub q: i32,                        /* method order                        */
    pub p: i32,                        /* embedding order                     */
    pub stages: i32,                   /* total number of stages              */
    pub nstages_active: i32,           /* number of active stage RHS vectors  */
    pub nstages_allocated: i32,        /* number of stage RHS vectors alloc'd */
    pub stage_map: Vec<i32>,           /* index map for stage RHS vectors     */
    pub stagetypes: Vec<i32>,          /* type flags for stages               */
    pub Ae_row: Vec<sunrealtype>,      /* equivalent explicit RK coeffs       */
    pub Ai_row: Vec<sunrealtype>,      /* equivalent implicit RK coeffs       */

    /* Algebraic solver data and parameters */
    pub sdata: Option<N_Vector>,          /* old stage data in residual        */
    pub zpred: Option<N_Vector>,          /* predicted stage solution          */
    pub zcor: Option<N_Vector>,           /* stage correction                  */
    pub NLS: Option<SUNNonlinearSolver>,  /* generic SUNNonlinearSolver object */
    pub ownNLS: sunbooleantype,           /* flag indicating ownership of NLS  */
    pub nls_fsi: Option<ARKRhsFn>,        /* fsi(t,y) used in the nonlin solver*/
    pub gamma: sunrealtype,               /* gamma = h * A(i,i)                */
    pub gammap: sunrealtype,              /* gamma at the last setup call      */
    pub gamrat: sunrealtype,              /* gamma / gammap                    */
    pub dgmax: sunrealtype,               /* call lsetup if |gamma/gammap-1| >= dgmax */
    pub predictor: i32,                   /* implicit prediction method to use */
    pub crdown: sunrealtype,              /* nonlin conv rate estimation const */
    pub rdiv: sunrealtype,                /* divergence if delnrm/delnrm_p > rdiv */
    /// C `step_mem->crate` (estimated nonlinear convergence rate); `crate`
    /// is a Rust keyword and cannot even be written raw.
    pub crate_: sunrealtype,
    pub delnrm_p: sunrealtype, /* norm of previous nonlinear solver update */
    pub delnrm: sunrealtype,   /* norm of current nonlinear solver update  */
    pub eRNrm: sunrealtype,    /* estimated residual norm, used in nonlin
                               and linear solver convergence tests      */
    pub nlscoef: sunrealtype,  /* coefficient in nonlin. convergence test  */
    pub msbp: i32,             /* positive => max # steps between lsetup
                               negative => call at each Newton iter     */
    pub nstlp: i64,            /* step number of last setup call           */
    pub maxcor: i32,           /* max num iterations for solving the
                               nonlinear equation                       */
    pub convfail: i32,         /* NLS fail flag (for interface routines)   */
    /// C `sunbooleantype jcur` — is the Jacobian info for the linear solver
    /// current? Shared cell so `step_getgammas` can hand out the same flag
    /// `arkLsSetup` / `arkLsPSetup` write through.
    pub jcur: ARKJcurPtr,
    pub stage_predict: Option<ARKStagePredictFn>, /* User-supplied stage predictor */
    pub istage: i32,                              /* stage index in nonlinear solve */

    /* Informational output for mriStep_GetStageIndex -- note that this
       may differ from istage, since that is used internally by the
       nonlinear solver, and it is manually modified during embedding
       stages to match the last internal stage index. */
    pub cur_stage: i32,

    /* Linear Solver Data (C `void* lmem` lives in `ark_mem.ark_lmem`) */
    pub linit: Option<ARKLinsolInitFn>,
    pub lsetup: Option<ARKLinsolSetupFn>,
    pub lsolve: Option<ARKLinsolSolveFn>,
    pub lfree: Option<ARKLinsolFreeFn>,

    /* Inner stepper */
    pub stepper: Option<MRIStepInnerStepper>,

    /* User-supplied pre and post inner evolve functions */
    pub pre_inner_evolve: Option<MRIStepPreInnerFn>,
    pub post_inner_evolve: Option<MRIStepPostInnerFn>,

    /* MRI adaptivity parameters */
    pub inner_rtol_factor: sunrealtype, /* prev control parameter               */
    pub inner_dsm: sunrealtype,         /* prev inner stepper accumulated error */
    pub inner_rtol_factor_new: sunrealtype, /* upcoming control parameter       */

    /* Counters */
    pub nfse: i64,        /* num fse calls                    */
    pub nfsi: i64,        /* num fsi calls                    */
    pub nsetups: i64,     /* num linear solver setup calls    */
    pub nls_iters: i64,   /* num nonlinear solver iters       */
    pub nls_fails: i64,   /* num nonlinear solver fails       */
    pub inner_fails: i64, /* num recov. inner solver fails    */
    pub nfusedopvecs: i32, /* length of cvals and Xvecs arrays */

    /* Data for using MRIStep with external polynomial forcing */
    pub expforcing: sunbooleantype, /* add forcing to explicit RHS */
    pub impforcing: sunbooleantype, /* add forcing to implicit RHS */
    pub tshift: sunrealtype,        /* time normalization shift    */
    pub tscale: sunrealtype,        /* time normalization scaling  */
    pub forcing: Vec<N_Vector>,     /* array of forcing vectors    */
    pub nforcing: i32,              /* number of forcing vectors   */

    /* Reusable arrays for fused vector operations */
    pub cvals: Vec<sunrealtype>,
    pub Xvecs: Vec<Option<N_Vector>>,
}

impl ARKodeMRIStepMemRec {
    /// C `calloc(1, sizeof(*step_mem))` in `MRIStepCreate`.
    pub fn zeroed() -> ARKodeMRIStepMemRec {
        ARKodeMRIStepMemRec {
            fse: None,
            fsi: None,
            linear: SUNFALSE,
            linear_timedep: SUNFALSE,
            explicit_rhs: SUNFALSE,
            implicit_rhs: SUNFALSE,
            deduce_rhs: SUNFALSE,
            Fse: Vec::new(),
            Fsi: Vec::new(),
            unify_Fs: SUNFALSE,
            fse_is_current: SUNFALSE,
            fsi_is_current: SUNFALSE,
            MRIC: None,
            q: 0,
            p: 0,
            stages: 0,
            nstages_active: 0,
            nstages_allocated: 0,
            stage_map: Vec::new(),
            stagetypes: Vec::new(),
            Ae_row: Vec::new(),
            Ai_row: Vec::new(),
            sdata: None,
            zpred: None,
            zcor: None,
            NLS: None,
            ownNLS: SUNFALSE,
            nls_fsi: None,
            gamma: ZERO,
            gammap: ZERO,
            gamrat: ZERO,
            dgmax: ZERO,
            predictor: 0,
            crdown: ZERO,
            rdiv: ZERO,
            crate_: ZERO,
            delnrm_p: ZERO,
            delnrm: ZERO,
            eRNrm: ZERO,
            nlscoef: ZERO,
            msbp: 0,
            nstlp: 0,
            maxcor: 0,
            convfail: 0,
            jcur: Rc::new(Cell::new(SUNFALSE)),
            stage_predict: None,
            istage: 0,
            cur_stage: 0,
            linit: None,
            lsetup: None,
            lsolve: None,
            lfree: None,
            stepper: None,
            pre_inner_evolve: None,
            post_inner_evolve: None,
            inner_rtol_factor: ZERO,
            inner_dsm: ZERO,
            inner_rtol_factor_new: ZERO,
            nfse: 0,
            nfsi: 0,
            nsetups: 0,
            nls_iters: 0,
            nls_fails: 0,
            inner_fails: 0,
            nfusedopvecs: 0,
            expforcing: SUNFALSE,
            impforcing: SUNFALSE,
            tshift: ZERO,
            tscale: ZERO,
            forcing: Vec::new(),
            nforcing: 0,
            cvals: Vec::new(),
            Xvecs: Vec::new(),
        }
    }
}

/// Downcast helper: view `ark_mem.step_mem` as the MRIStep memory record.
///
/// Panics if no stepper memory is attached or it is not an MRIStep record
/// (C would blindly cast the `void*` — UB maps to a deterministic panic).
/// NEVER hold the returned guard across `arkProcessError`, a user callback,
/// an `N_Vector` / matrix / linear-solver / nonlinear-solver operation, an
/// inner-stepper call, or another borrow of the same `ark_mem`.
pub fn mriStep_mem_mut(ark_mem: &ARKodeMem) -> RefMut<'_, ARKodeMRIStepMemRec> {
    RefMut::map(ark_mem.borrow_mut(), |m| {
        m.step_mem
            .as_mut()
            .expect("step_mem set")
            .downcast_mut::<ARKodeMRIStepMemRec>()
            .expect("MRIStep step memory")
    })
}

/// C `mriStep_AccessStepMem(ark_mem, fname, &step_mem)` reduced to its
/// presence check; the record itself is reached with [`mriStep_mem_mut`] at
/// each use site (contract §3).
fn mriStep_step_mem_ok(ark_mem: &ARKodeMem, fname: &str) -> i32 {
    let missing = ark_mem.borrow().step_mem.is_none();
    if missing {
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

/// Materialise the first `nvec` entries of the fused-op vector scratch as
/// the dense `&[N_Vector]` the `N_V*` kernels take (C hands the calloc'd
/// `step_mem->Xvecs` array over directly).
pub fn mriStep_xvecs(step_mem: &ARKodeMRIStepMemRec, nvec: i32) -> Vec<N_Vector> {
    step_mem.Xvecs[..nvec as usize]
        .iter()
        .map(|v| v.clone().expect("Xvecs entry set"))
        .collect()
}

/*===============================================================
  Callback invocation helpers

  C `void* user_data` is `Option<Box<dyn Any>>`: the box is taken out of
  the mem for the duration of the call and restored on EVERY path, and no
  borrow of `ark_mem` is held across the callback.
  ===============================================================*/

/// C `step_mem->fse(t, y, ydot, ark_mem->user_data)`.
pub fn mriStep_CallFse(
    ark_mem: &ARKodeMem,
    t: sunrealtype,
    y: &N_Vector,
    ydot: &N_Vector,
) -> i32 {
    let fse = { mriStep_mem_mut(ark_mem).fse }.expect("fse set");
    let mut user_data = ark_mem.borrow_mut().user_data.take();
    let retval = fse(t, y, ydot, &mut user_data);
    ark_mem.borrow_mut().user_data = user_data;
    retval
}

/// C `step_mem->fsi(t, y, ydot, ark_mem->user_data)`.
pub fn mriStep_CallFsi(
    ark_mem: &ARKodeMem,
    t: sunrealtype,
    y: &N_Vector,
    ydot: &N_Vector,
) -> i32 {
    let fsi = { mriStep_mem_mut(ark_mem).fsi }.expect("fsi set");
    let mut user_data = ark_mem.borrow_mut().user_data.take();
    let retval = fsi(t, y, ydot, &mut user_data);
    ark_mem.borrow_mut().user_data = user_data;
    retval
}

/// C `ark_mem->PreRhsFn(t, y, ark_mem->user_data)`.
pub fn mriStep_CallPreRhsFn(ark_mem: &ARKodeMem, t: sunrealtype, y: &N_Vector) -> i32 {
    let f = ark_mem.borrow().PreRhsFn.expect("PreRhsFn set");
    let mut user_data = ark_mem.borrow_mut().user_data.take();
    let retval = f(t, y, &mut user_data);
    ark_mem.borrow_mut().user_data = user_data;
    retval
}

/// C `ark_mem->PostProcessStageFn(t, y, ark_mem->user_data)`.
pub fn mriStep_CallPostProcessStageFn(
    ark_mem: &ARKodeMem,
    t: sunrealtype,
    y: &N_Vector,
) -> i32 {
    let f = ark_mem.borrow().PostProcessStageFn.expect("PostProcessStageFn set");
    let mut user_data = ark_mem.borrow_mut().user_data.take();
    let retval = f(t, y, &mut user_data);
    ark_mem.borrow_mut().user_data = user_data;
    retval
}

/// C `ark_mem->PostProcessStepFn(t, y, ark_mem->user_data)`.
pub fn mriStep_CallPostProcessStepFn(
    ark_mem: &ARKodeMem,
    t: sunrealtype,
    y: &N_Vector,
) -> i32 {
    let f = ark_mem.borrow().PostProcessStepFn.expect("PostProcessStepFn set");
    let mut user_data = ark_mem.borrow_mut().user_data.take();
    let retval = f(t, y, &mut user_data);
    ark_mem.borrow_mut().user_data = user_data;
    retval
}

/// C `step_mem->stage_predict(t, zpred, ark_mem->user_data)`.
pub fn mriStep_CallStagePredict(
    ark_mem: &ARKodeMem,
    t: sunrealtype,
    zpred: &N_Vector,
) -> i32 {
    let f = { mriStep_mem_mut(ark_mem).stage_predict }.expect("stage_predict set");
    let mut user_data = ark_mem.borrow_mut().user_data.take();
    let retval = f(t, zpred, &mut user_data);
    ark_mem.borrow_mut().user_data = user_data;
    retval
}

/*===============================================================
  Exported functions
  ===============================================================*/

pub fn MRIStepCreate(
    fse: Option<ARKRhsFn>,
    fsi: Option<ARKRhsFn>,
    t0: sunrealtype,
    y0: &N_Vector,
    stepper: &MRIStepInnerStepper,
    sunctx: &SUNContext,
) -> Option<ARKodeMem> {
    let mut retval: i32;

    /* Check that at least one of fse, fsi is supplied and is to be used*/
    if fse.is_none() && fsi.is_none() {
        arkProcessError(
            None,
            ARK_ILL_INPUT,
            line!() as i32,
            "MRIStepCreate",
            file!(),
            MSG_ARK_NULL_F,
        );
        return None;
    }

    /* Check that y0 is supplied: handled by the type system */
    /* Check that stepper is supplied: handled by the type system */
    /* Check that context is supplied: handled by the type system */

    /* Create ark_mem structure and set default values */
    let ark_mem = match arkCreate(sunctx) {
        Some(m) => m,
        None => {
            arkProcessError(
                None,
                ARK_MEM_NULL,
                line!() as i32,
                "MRIStepCreate",
                file!(),
                MSG_ARK_NO_MEM,
            );
            return None;
        }
    };

    /* Allocate ARKodeMRIStepMem structure, and initialize to zero */
    let step_mem = ARKodeMRIStepMemRec::zeroed();

    /* Attach step_mem structure and function pointers to ark_mem */
    {
        let mut m = ark_mem.borrow_mut();
        m.step_attachlinsol = Some(mriStep_AttachLinsol);
        m.step_disablelsetup = Some(mriStep_DisableLSetup);
        m.step_getlinmem = Some(mriStep_GetLmem);
        m.step_getimplicitrhs = Some(mriStep_GetImplicitRHS);
        m.step_getgammas = Some(mriStep_GetGammas);
        m.step_init = Some(mriStep_Init);
        m.step_fullrhs = Some(mriStep_FullRHS);
        m.step = Some(mriStep_TakeStepMRIGARK);
        m.step_setuserdata = Some(mriStep_SetUserData);
        m.step_printallstats = Some(mriStep_PrintAllStats);
        m.step_writeparameters = Some(mriStep_WriteParameters);
        m.step_setusecompensatedsums = None;
        m.step_resize = Some(mriStep_Resize);
        m.step_reset = Some(mriStep_Reset);
        m.step_free = Some(mriStep_Free);
        m.step_printmem = Some(mriStep_PrintMem);
        m.step_setdefaults = Some(mriStep_SetDefaults);
        m.step_computestate = Some(mriStep_ComputeState);
        m.step_setoptions = Some(mriStep_SetOptions);
        m.step_setorder = Some(mriStep_SetOrder);
        m.step_setnonlinearsolver = Some(mriStep_SetNonlinearSolver);
        m.step_setlinear = Some(mriStep_SetLinear);
        m.step_setnonlinear = Some(mriStep_SetNonlinear);
        m.step_setnlsrhsfn = Some(mriStep_SetNlsRhsFn);
        m.step_setdeduceimplicitrhs = Some(mriStep_SetDeduceImplicitRhs);
        m.step_setnonlincrdown = Some(mriStep_SetNonlinCRDown);
        m.step_setnonlinrdiv = Some(mriStep_SetNonlinRDiv);
        m.step_setdeltagammamax = Some(mriStep_SetDeltaGammaMax);
        m.step_setlsetupfrequency = Some(mriStep_SetLSetupFrequency);
        m.step_setpredictormethod = Some(mriStep_SetPredictorMethod);
        m.step_setmaxnonliniters = Some(mriStep_SetMaxNonlinIters);
        m.step_setnonlinconvcoef = Some(mriStep_SetNonlinConvCoef);
        m.step_setstagepredictfn = Some(mriStep_SetStagePredictFn);
        m.step_getnumrhsevals = Some(mriStep_GetNumRhsEvals);
        m.step_getnumlinsolvsetups = Some(mriStep_GetNumLinSolvSetups);
        m.step_getcurrentgamma = Some(mriStep_GetCurrentGamma);
        m.step_setadaptcontroller = Some(mriStep_SetAdaptController);
        m.step_getestlocalerrors = Some(mriStep_GetEstLocalErrors);
        m.step_getnonlinearsystemdata = Some(mriStep_GetNonlinearSystemData);
        m.step_getnumnonlinsolviters = Some(mriStep_GetNumNonlinSolvIters);
        m.step_getnumnonlinsolvconvfails = Some(mriStep_GetNumNonlinSolvConvFails);
        m.step_getnonlinsolvstats = Some(mriStep_GetNonlinSolvStats);
        m.step_setforcing = Some(mriStep_SetInnerForcing);
        m.step_getstageindex = Some(mriStep_GetStageIndex);
        m.step_supports_adaptive = SUNTRUE;
        m.step_supports_implicit = SUNTRUE;
        m.step_mem = Some(Box::new(step_mem));
    }

    /* Set default values for optional inputs */
    retval = mriStep_SetDefaults(&ark_mem);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(&ark_mem),
            retval,
            line!() as i32,
            "MRIStepCreate",
            file!(),
            "Error setting default solver options",
        );
        let mut mem = Some(ark_mem.clone());
        ARKodeFree(&mut mem);
        return None;
    }

    /* Allocate the general MRI stepper vectors using y0 as a template */
    /* NOTE: Fse, Fsi, inner_forcing, sdata, zpred and zcor will be allocated
       later on (based on the MRI method) */

    /* Copy the slow RHS functions into stepper memory */
    {
        let mut s = mriStep_mem_mut(&ark_mem);
        s.fse = fse;
        s.fsi = fsi;
        s.fse_is_current = SUNFALSE;
        s.fsi_is_current = SUNFALSE;

        /* Set implicit/explicit problem based on function pointers */
        s.explicit_rhs = fse.is_some();
        s.implicit_rhs = fsi.is_some();
    }

    /* Update the ARKODE workspace requirements */
    {
        let mut m = ark_mem.borrow_mut();
        m.liw += 49; /* fcn/data ptr, int, long int, sunindextype, sunbooleantype */
        m.lrw += 14;
    }

    /* Create a default Newton NLS object (just in case; will be deleted if
       the user attaches a nonlinear solver) */
    let implicit_rhs = {
        let mut s = mriStep_mem_mut(&ark_mem);
        s.NLS = None;
        s.ownNLS = SUNFALSE;
        s.implicit_rhs
    };

    if implicit_rhs {
        let ark_sunctx = ark_mem.borrow().sunctx.clone();
        let NLS = match SUNNonlinSol_Newton(y0, &ark_sunctx) {
            Some(nls) => nls,
            None => {
                arkProcessError(
                    Some(&ark_mem),
                    ARK_MEM_FAIL,
                    line!() as i32,
                    "MRIStepCreate",
                    file!(),
                    "Error creating default Newton solver",
                );
                let mut mem = Some(ark_mem.clone());
                ARKodeFree(&mut mem);
                return None;
            }
        };
        retval = ARKodeSetNonlinearSolver(&ark_mem, &NLS);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(&ark_mem),
                ARK_MEM_FAIL,
                line!() as i32,
                "MRIStepCreate",
                file!(),
                "Error attaching default Newton solver",
            );
            let mut mem = Some(ark_mem.clone());
            ARKodeFree(&mut mem);
            return None;
        }
        mriStep_mem_mut(&ark_mem).ownNLS = SUNTRUE;
    }

    /* Set the linear solver addresses to NULL (we check != NULL later) */
    {
        let mut s = mriStep_mem_mut(&ark_mem);
        s.linit = None;
        s.lsetup = None;
        s.lsolve = None;
        s.lfree = None;

        /* Initialize error norm  */
        s.eRNrm = ONE;

        /* Initialize all the counters */
        s.nfse = 0;
        s.nfsi = 0;
        s.nsetups = 0;
        s.nstlp = 0;
        s.nls_iters = 0;
        s.nls_fails = 0;
        s.inner_fails = 0;
    }
    /* C `step_mem->lmem = NULL`: the ARKLS record lives in ark_mem (§4) */
    ark_mem.borrow_mut().ark_lmem = None;

    /* Initialize fused op work space with sufficient storage for at least filling
       the full RHS on an ImEx problem -- must be allocate here as the full RHS
       is called before mriStep_Init when nesting MRI methods.
       The C calloc-failure branches are unreachable: a Rust allocation failure
       aborts rather than returning NULL. */
    let nfusedopvecs = {
        let mut s = mriStep_mem_mut(&ark_mem);
        s.nfusedopvecs = 3;
        s.cvals = vec![ZERO; s.nfusedopvecs as usize];
        s.Xvecs = vec![None; s.nfusedopvecs as usize];
        s.nfusedopvecs
    };
    {
        let mut m = ark_mem.borrow_mut();
        m.lrw += nfusedopvecs as i64;
        m.liw += nfusedopvecs as i64;
    }

    {
        let mut s = mriStep_mem_mut(&ark_mem);

        /* Initialize adaptivity parameters */
        s.inner_rtol_factor = ONE;
        s.inner_dsm = ONE;
        s.inner_rtol_factor_new = ONE;

        /* Initialize pre and post inner evolve functions */
        s.pre_inner_evolve = None;
        s.post_inner_evolve = None;

        /* Initialize external polynomial forcing data */
        s.expforcing = SUNFALSE;
        s.impforcing = SUNFALSE;
        s.forcing = Vec::new();
        s.nforcing = 0;
    }

    /* Initialize main ARKODE infrastructure (allocates vectors) */
    retval = arkInit(&ark_mem, t0, y0, FIRST_INIT);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(&ark_mem),
            retval,
            line!() as i32,
            "MRIStepCreate",
            file!(),
            "Unable to initialize main ARKODE infrastructure",
        );
        let mut mem = Some(ark_mem.clone());
        ARKodeFree(&mut mem);
        return None;
    }

    /* Attach the inner stepper memory */
    mriStep_mem_mut(&ark_mem).stepper = Some(stepper.clone());

    /* Check for required stepper functions */
    retval = mriStepInnerStepper_HasRequiredOps(stepper);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(&ark_mem),
            ARK_ILL_INPUT,
            line!() as i32,
            "MRIStepCreate",
            file!(),
            "A required inner stepper function is NULL",
        );
        let mut mem = Some(ark_mem.clone());
        ARKodeFree(&mut mem);
        return None;
    }

    /* return ARKODE memory */
    Some(ark_mem)
}

/*---------------------------------------------------------------
  MRIStepReInit:

  This routine re-initializes the MRIStep module to solve a new
  problem of the same size as was previously solved (all counter
  values are set to 0).

  NOTE: the inner stepper needs to be reinitialized before
  calling this function.
  ---------------------------------------------------------------*/
pub fn MRIStepReInit(
    arkode_mem: &ARKodeMem,
    fse: Option<ARKRhsFn>,
    fsi: Option<ARKRhsFn>,
    t0: sunrealtype,
    y0: &N_Vector,
) -> i32 {
    let mut retval: i32;

    /* access ARKodeMem and ARKodeMRIStepMem structures */
    retval = mriStep_step_mem_ok(arkode_mem, "MRIStepReInit");
    if retval != ARK_SUCCESS {
        return retval;
    }
    let ark_mem = arkode_mem;

    /* Check if ark_mem was allocated */
    if !ark_mem.borrow().MallocDone {
        arkProcessError(
            Some(ark_mem),
            ARK_NO_MALLOC,
            line!() as i32,
            "MRIStepReInit",
            file!(),
            MSG_ARK_NO_MALLOC,
        );
        return ARK_NO_MALLOC;
    }

    /* Check that at least one of fse, fsi is supplied and is to be used */
    if fse.is_none() && fsi.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!() as i32,
            "MRIStepReInit",
            file!(),
            MSG_ARK_NULL_F,
        );
        return ARK_ILL_INPUT;
    }

    /* Check that y0 is supplied: handled by the type system */

    /* Set implicit/explicit problem based on function pointers */
    let (implicit_rhs, have_nls) = {
        let mut s = mriStep_mem_mut(ark_mem);
        s.explicit_rhs = fse.is_some();
        s.implicit_rhs = fsi.is_some();
        (s.implicit_rhs, s.NLS.is_some())
    };

    /* Create a default Newton NLS object (just in case; will be deleted if
       the user attaches a nonlinear solver) */
    if implicit_rhs && !have_nls {
        let ark_sunctx = ark_mem.borrow().sunctx.clone();
        let NLS = match SUNNonlinSol_Newton(y0, &ark_sunctx) {
            Some(nls) => nls,
            None => {
                arkProcessError(
                    Some(ark_mem),
                    ARK_MEM_FAIL,
                    line!() as i32,
                    "MRIStepReInit",
                    file!(),
                    "Error creating default Newton solver",
                );
                let mut mem = Some(ark_mem.clone());
                ARKodeFree(&mut mem);
                return ARK_MEM_FAIL;
            }
        };
        retval = ARKodeSetNonlinearSolver(ark_mem, &NLS);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_FAIL,
                line!() as i32,
                "MRIStepReInit",
                file!(),
                "Error attaching default Newton solver",
            );
            let mut mem = Some(ark_mem.clone());
            ARKodeFree(&mut mem);
            return ARK_MEM_FAIL;
        }
        mriStep_mem_mut(ark_mem).ownNLS = SUNTRUE;
    }

    /* ReInitialize main ARKODE infrastructure */
    retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            retval,
            line!() as i32,
            "MRIStepReInit",
            file!(),
            "Unable to reinitialize main ARKODE infrastructure",
        );
        return retval;
    }

    /* Copy the input parameters into ARKODE state */
    {
        let mut s = mriStep_mem_mut(ark_mem);
        s.fse = fse;
        s.fsi = fsi;
        s.fse_is_current = SUNFALSE;
        s.fsi_is_current = SUNFALSE;

        /* Initialize all the counters */
        s.nfse = 0;
        s.nfsi = 0;
        s.nsetups = 0;
        s.nstlp = 0;
        s.nls_iters = 0;
        s.nls_fails = 0;
        s.inner_fails = 0;
    }

    /* C `if (step_mem->lmem) { arkLsInitializeCounters(step_mem->lmem); }` */
    let have_lmem = ark_mem.borrow().ark_lmem.is_some();
    if have_lmem {
        let mut arkls_mem = arkls_mem_mut(ark_mem);
        arkLsInitializeCounters(&mut arkls_mem);
    }

    ARK_SUCCESS
}

/*===============================================================
  Interface routines supplied to ARKODE
  ===============================================================*/

/*---------------------------------------------------------------
  mriStep_Resize:

  This routine resizes the memory within the MRIStep module.
  ---------------------------------------------------------------*/
pub fn mriStep_Resize(
    ark_mem: &ARKodeMem,
    y0: &N_Vector,
    _hscale: sunrealtype,
    _t0: sunrealtype,
    resize: Option<ARKVecResizeFn>,
    resize_data: &mut Option<Box<dyn Any>>,
) -> i32 {
    let mut retval: i32;

    /* access ARKodeMRIStepMem structure */
    retval = mriStep_step_mem_ok(ark_mem, "mriStep_Resize");
    if retval != ARK_SUCCESS {
        return retval;
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

    /* Resize Fse */
    let (have_Fse, nstages_allocated) = {
        let s = mriStep_mem_mut(ark_mem);
        (!s.Fse.is_empty(), s.nstages_allocated)
    };
    if have_Fse {
        let mut Fse = { std::mem::take(&mut mriStep_mem_mut(ark_mem).Fse) };
        let (mut lrw, mut liw) = {
            let m = ark_mem.borrow();
            (m.lrw, m.liw)
        };
        let ok = arkResizeVecArray(
            resize,
            resize_data,
            nstages_allocated,
            y0,
            &mut Fse,
            lrw_diff,
            &mut lrw,
            liw_diff,
            &mut liw,
        );
        {
            let mut m = ark_mem.borrow_mut();
            m.lrw = lrw;
            m.liw = liw;
        }
        {
            let mut s = mriStep_mem_mut(ark_mem);
            s.Fse = Fse;
        }
        if !ok {
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_FAIL,
                line!() as i32,
                "mriStep_Resize",
                file!(),
                "Unable to resize vector",
            );
            return ARK_MEM_FAIL;
        }
        let mut s = mriStep_mem_mut(ark_mem);
        if s.unify_Fs {
            s.Fsi = s.Fse.clone();
        }
    }

    /* Resize Fsi */
    let resize_Fsi = {
        let s = mriStep_mem_mut(ark_mem);
        !s.Fsi.is_empty() && !s.unify_Fs
    };
    if resize_Fsi {
        let mut Fsi = { std::mem::take(&mut mriStep_mem_mut(ark_mem).Fsi) };
        let (mut lrw, mut liw) = {
            let m = ark_mem.borrow();
            (m.lrw, m.liw)
        };
        let ok = arkResizeVecArray(
            resize,
            resize_data,
            nstages_allocated,
            y0,
            &mut Fsi,
            lrw_diff,
            &mut lrw,
            liw_diff,
            &mut liw,
        );
        {
            let mut m = ark_mem.borrow_mut();
            m.lrw = lrw;
            m.liw = liw;
        }
        {
            let mut s = mriStep_mem_mut(ark_mem);
            s.Fsi = Fsi;
        }
        if !ok {
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_FAIL,
                line!() as i32,
                "mriStep_Resize",
                file!(),
                "Unable to resize vector",
            );
            return ARK_MEM_FAIL;
        }
    }

    /* Resize the nonlinear solver interface vectors (if applicable) */
    {
        let mut sdata = { mriStep_mem_mut(ark_mem).sdata.clone() };
        if sdata.is_some() {
            let ok = arkResizeVec(
                ark_mem,
                resize,
                resize_data,
                lrw_diff,
                liw_diff,
                y0,
                &mut sdata,
            );
            mriStep_mem_mut(ark_mem).sdata = sdata;
            if !ok {
                arkProcessError(
                    Some(ark_mem),
                    ARK_MEM_FAIL,
                    line!() as i32,
                    "mriStep_Resize",
                    file!(),
                    "Unable to resize vector",
                );
                return ARK_MEM_FAIL;
            }
        }
    }
    {
        let mut zpred = { mriStep_mem_mut(ark_mem).zpred.clone() };
        if zpred.is_some() {
            let ok = arkResizeVec(
                ark_mem,
                resize,
                resize_data,
                lrw_diff,
                liw_diff,
                y0,
                &mut zpred,
            );
            mriStep_mem_mut(ark_mem).zpred = zpred;
            if !ok {
                arkProcessError(
                    Some(ark_mem),
                    ARK_MEM_FAIL,
                    line!() as i32,
                    "mriStep_Resize",
                    file!(),
                    "Unable to resize vector",
                );
                return ARK_MEM_FAIL;
            }
        }
    }
    {
        let mut zcor = { mriStep_mem_mut(ark_mem).zcor.clone() };
        if zcor.is_some() {
            let ok = arkResizeVec(
                ark_mem,
                resize,
                resize_data,
                lrw_diff,
                liw_diff,
                y0,
                &mut zcor,
            );
            mriStep_mem_mut(ark_mem).zcor = zcor;
            if !ok {
                arkProcessError(
                    Some(ark_mem),
                    ARK_MEM_FAIL,
                    line!() as i32,
                    "mriStep_Resize",
                    file!(),
                    "Unable to resize vector",
                );
                return ARK_MEM_FAIL;
            }
        }
    }

    /* If a NLS object was previously used, destroy and recreate default Newton
       NLS object (can be replaced by user-defined object if desired) */
    let recreate_nls = {
        let s = mriStep_mem_mut(ark_mem);
        s.NLS.is_some() && s.ownNLS
    };
    if recreate_nls {
        /* destroy existing NLS object */
        let old_nls = { mriStep_mem_mut(ark_mem).NLS.take() };
        retval = SUNNonlinSolFree(old_nls);
        if retval != ARK_SUCCESS {
            return retval;
        }
        {
            let mut s = mriStep_mem_mut(ark_mem);
            s.NLS = None;
            s.ownNLS = SUNFALSE;
        }

        /* create new Newton NLS object */
        let ark_sunctx = ark_mem.borrow().sunctx.clone();
        let NLS = match SUNNonlinSol_Newton(y0, &ark_sunctx) {
            Some(nls) => nls,
            None => {
                arkProcessError(
                    Some(ark_mem),
                    ARK_MEM_FAIL,
                    line!() as i32,
                    "mriStep_Resize",
                    file!(),
                    "Error creating default Newton solver",
                );
                return ARK_MEM_FAIL;
            }
        };

        /* attach new Newton NLS object */
        retval = ARKodeSetNonlinearSolver(ark_mem, &NLS);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_FAIL,
                line!() as i32,
                "mriStep_Resize",
                file!(),
                "Error attaching default Newton solver",
            );
            return ARK_MEM_FAIL;
        }
        mriStep_mem_mut(ark_mem).ownNLS = SUNTRUE;
    }

    /* Resize the inner stepper vectors */
    let stepper = { mriStep_mem_mut(ark_mem).stepper.clone() }.expect("stepper set");
    retval = mriStepInnerStepper_Resize(&stepper, resize, resize_data, lrw_diff, liw_diff, y0);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_FAIL,
            line!() as i32,
            "mriStep_Resize",
            file!(),
            "Unable to resize vector",
        );
        return ARK_MEM_FAIL;
    }

    /* reset nonlinear solver counters */
    {
        let mut s = mriStep_mem_mut(ark_mem);
        if s.NLS.is_some() {
            s.nsetups = 0;
        }
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_Reset:

  This routine resets the MRIStep module state to solve the same
  problem from the given time with the input state (all counter
  values are retained).  It is called after the main ARKODE
  infrastructure is reset.
  ---------------------------------------------------------------*/
pub fn mriStep_Reset(ark_mem: &ARKodeMem, tR: sunrealtype, yR: &N_Vector) -> i32 {
    /* access ARKodeMRIStepMem structure */
    let retval = mriStep_step_mem_ok(ark_mem, "mriStep_Reset");
    if retval != ARK_SUCCESS {
        return retval;
    }

    /* Reset the inner integrator with this same state */
    let stepper = { mriStep_mem_mut(ark_mem).stepper.clone() }.expect("stepper set");
    let retval = mriStepInnerStepper_Reset(&stepper, tR, yR);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            ARK_INNERSTEP_FAIL,
            line!() as i32,
            "mriStep_Reset",
            file!(),
            "Unable to reset the inner stepper",
        );
        return ARK_INNERSTEP_FAIL;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_ComputeState:

  Computes y based on the current prediction and given correction.
  ---------------------------------------------------------------*/
pub fn mriStep_ComputeState(ark_mem: &ARKodeMem, zcor: &N_Vector, z: &N_Vector) -> i32 {
    /* access ARKodeMRIStepMem structure */
    let retval = mriStep_step_mem_ok(ark_mem, "mriStep_ComputeState");
    if retval != ARK_SUCCESS {
        return retval;
    }

    let zpred = { mriStep_mem_mut(ark_mem).zpred.clone() }.expect("zpred set");
    N_VLinearSum(ONE, &zpred, ONE, zcor, z);

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_Free frees all MRIStep memory.
  ---------------------------------------------------------------*/
pub fn mriStep_Free(ark_mem: &ARKodeMem) {
    /* nothing to do if ark_mem is already NULL: handled by the type system */

    /* conditional frees on non-NULL MRIStep module */
    if ark_mem.borrow().step_mem.is_none() {
        return;
    }

    /* free the coupling structure and derived quantities */
    let MRIC = { mriStep_mem_mut(ark_mem).MRIC.take() };
    if let Some(MRIC) = MRIC {
        let mut Cliw: sunindextype = 0;
        let mut Clrw: sunindextype = 0;
        MRIStepCoupling_Space(&MRIC, &mut Cliw, &mut Clrw);
        MRIStepCoupling_Free(Some(MRIC));
        /* `step_mem->MRIC = NULL` performed by the `take()` above */
        {
            let mut m = ark_mem.borrow_mut();
            m.liw -= Cliw;
            m.lrw -= Clrw;
        }
        let stages = { mriStep_mem_mut(ark_mem).stages };
        let have_stagetypes = { !mriStep_mem_mut(ark_mem).stagetypes.is_empty() };
        if have_stagetypes {
            mriStep_mem_mut(ark_mem).stagetypes = Vec::new();
            ark_mem.borrow_mut().liw -= (stages + 1) as i64;
        }
        let have_stage_map = { !mriStep_mem_mut(ark_mem).stage_map.is_empty() };
        if have_stage_map {
            mriStep_mem_mut(ark_mem).stage_map = Vec::new();
            ark_mem.borrow_mut().liw -= stages as i64;
        }
        let have_Ae_row = { !mriStep_mem_mut(ark_mem).Ae_row.is_empty() };
        if have_Ae_row {
            mriStep_mem_mut(ark_mem).Ae_row = Vec::new();
            ark_mem.borrow_mut().lrw -= stages as i64;
        }
        let have_Ai_row = { !mriStep_mem_mut(ark_mem).Ai_row.is_empty() };
        if have_Ai_row {
            mriStep_mem_mut(ark_mem).Ai_row = Vec::new();
            ark_mem.borrow_mut().lrw -= stages as i64;
        }
    }

    /* free the nonlinear solver memory (if applicable) */
    let free_nls = {
        let s = mriStep_mem_mut(ark_mem);
        s.NLS.is_some() && s.ownNLS
    };
    if free_nls {
        let nls = { mriStep_mem_mut(ark_mem).NLS.clone() };
        let _ = SUNNonlinSolFree(nls);
        mriStep_mem_mut(ark_mem).ownNLS = SUNFALSE;
    }
    mriStep_mem_mut(ark_mem).NLS = None;

    /* free the linear solver memory */
    let lfree = { mriStep_mem_mut(ark_mem).lfree };
    if let Some(lfree) = lfree {
        let _ = lfree(ark_mem);
        /* C `step_mem->lmem = NULL` (the record lives in ark_mem, §4) */
        ark_mem.borrow_mut().ark_lmem = None;
    }

    /* free the sdata, zpred and zcor vectors */
    {
        let mut sdata = { mriStep_mem_mut(ark_mem).sdata.take() };
        if sdata.is_some() {
            arkFreeVec(ark_mem, &mut sdata);
            mriStep_mem_mut(ark_mem).sdata = None;
        }
    }
    {
        let mut zpred = { mriStep_mem_mut(ark_mem).zpred.take() };
        if zpred.is_some() {
            arkFreeVec(ark_mem, &mut zpred);
            mriStep_mem_mut(ark_mem).zpred = None;
        }
    }
    {
        let mut zcor = { mriStep_mem_mut(ark_mem).zcor.take() };
        if zcor.is_some() {
            arkFreeVec(ark_mem, &mut zcor);
            mriStep_mem_mut(ark_mem).zcor = None;
        }
    }

    /* free the RHS vectors */
    let (nstages_allocated, have_Fse) = {
        let s = mriStep_mem_mut(ark_mem);
        (s.nstages_allocated, !s.Fse.is_empty())
    };
    if have_Fse {
        let (lrw1, liw1, mut lrw, mut liw) = {
            let m = ark_mem.borrow();
            (m.lrw1, m.liw1, m.lrw, m.liw)
        };
        let mut Fse = { std::mem::take(&mut mriStep_mem_mut(ark_mem).Fse) };
        arkFreeVecArray(nstages_allocated, &mut Fse, lrw1, &mut lrw, liw1, &mut liw);
        {
            let mut m = ark_mem.borrow_mut();
            m.lrw = lrw;
            m.liw = liw;
        }
        let mut s = mriStep_mem_mut(ark_mem);
        s.Fse = Fse;
        if s.unify_Fs {
            s.Fsi = Vec::new();
        }
    }

    let have_Fsi = { !mriStep_mem_mut(ark_mem).Fsi.is_empty() };
    if have_Fsi {
        let (lrw1, liw1, mut lrw, mut liw) = {
            let m = ark_mem.borrow();
            (m.lrw1, m.liw1, m.lrw, m.liw)
        };
        let mut Fsi = { std::mem::take(&mut mriStep_mem_mut(ark_mem).Fsi) };
        arkFreeVecArray(nstages_allocated, &mut Fsi, lrw1, &mut lrw, liw1, &mut liw);
        {
            let mut m = ark_mem.borrow_mut();
            m.lrw = lrw;
            m.liw = liw;
        }
        mriStep_mem_mut(ark_mem).Fsi = Fsi;
    }

    /* free the reusable arrays for fused vector interface */
    let (have_cvals, have_Xvecs, nfusedopvecs) = {
        let s = mriStep_mem_mut(ark_mem);
        (!s.cvals.is_empty(), !s.Xvecs.is_empty(), s.nfusedopvecs)
    };
    if have_cvals {
        mriStep_mem_mut(ark_mem).cvals = Vec::new();
        ark_mem.borrow_mut().lrw -= nfusedopvecs as i64;
    }
    if have_Xvecs {
        mriStep_mem_mut(ark_mem).Xvecs = Vec::new();
        ark_mem.borrow_mut().liw -= nfusedopvecs as i64;
    }
    mriStep_mem_mut(ark_mem).nfusedopvecs = 0;

    /* free the time stepper module itself */
    ark_mem.borrow_mut().step_mem = None;
}

/*---------------------------------------------------------------
  mriStep_PrintMem:

  This routine outputs the memory from the MRIStep structure to
  a specified file pointer (useful when debugging).
  ---------------------------------------------------------------*/
pub fn mriStep_PrintMem(ark_mem: &ARKodeMem, outfile: &SUNFile) {
    /* access ARKodeMRIStepMem structure */
    let retval = mriStep_step_mem_ok(ark_mem, "mriStep_PrintMem");
    if retval != ARK_SUCCESS {
        return;
    }

    /* output integer quantities */
    let (q, p, istage, cur_stage, stages, maxcor, msbp, predictor, convfail, stagetypes) = {
        let s = mriStep_mem_mut(ark_mem);
        (
            s.q,
            s.p,
            s.istage,
            s.cur_stage,
            s.stages,
            s.maxcor,
            s.msbp,
            s.predictor,
            s.convfail,
            s.stagetypes.clone(),
        )
    };
    outfile.write_str(&format!("MRIStep: q = {q}\n"));
    outfile.write_str(&format!("MRIStep: p = {p}\n"));
    outfile.write_str(&format!("MRIStep: istage = {istage}\n"));
    outfile.write_str(&format!("MRIStep: cur_stage = {cur_stage}\n"));
    outfile.write_str(&format!("MRIStep: stages = {stages}\n"));
    outfile.write_str(&format!("MRIStep: maxcor = {maxcor}\n"));
    outfile.write_str(&format!("MRIStep: msbp = {msbp}\n"));
    outfile.write_str(&format!("MRIStep: predictor = {predictor}\n"));
    outfile.write_str(&format!("MRIStep: convfail = {convfail}\n"));
    outfile.write_str("MRIStep: stagetypes =");
    for i in 0..=stages {
        outfile.write_str(&format!(" {}", stagetypes[i as usize]));
    }
    outfile.write_str("\n");

    /* output long integer quantities */
    let (nfse, nfsi, nsetups, nstlp, nls_iters, nls_fails, inner_fails) = {
        let s = mriStep_mem_mut(ark_mem);
        (
            s.nfse,
            s.nfsi,
            s.nsetups,
            s.nstlp,
            s.nls_iters,
            s.nls_fails,
            s.inner_fails,
        )
    };
    outfile.write_str(&format!("MRIStep: nfse = {nfse}\n"));
    outfile.write_str(&format!("MRIStep: nfsi = {nfsi}\n"));
    outfile.write_str(&format!("MRIStep: nsetups = {nsetups}\n"));
    outfile.write_str(&format!("MRIStep: nstlp = {nstlp}\n"));
    outfile.write_str(&format!("MRIStep: nls_iters = {nls_iters}\n"));
    outfile.write_str(&format!("MRIStep: nls_fails = {nls_fails}\n"));
    outfile.write_str(&format!("MRIStep: inner_fails = {inner_fails}\n"));

    /* output boolean quantities */
    let (linear, linear_timedep, explicit_rhs, implicit_rhs, jcur, ownNLS) = {
        let s = mriStep_mem_mut(ark_mem);
        (
            s.linear,
            s.linear_timedep,
            s.explicit_rhs,
            s.implicit_rhs,
            s.jcur.get(),
            s.ownNLS,
        )
    };
    outfile.write_str(&format!("MRIStep: user_linear = {}\n", linear as i32));
    outfile.write_str(&format!(
        "MRIStep: user_linear_timedep = {}\n",
        linear_timedep as i32
    ));
    outfile.write_str(&format!("MRIStep: user_explicit = {}\n", explicit_rhs as i32));
    outfile.write_str(&format!("MRIStep: user_implicit = {}\n", implicit_rhs as i32));
    outfile.write_str(&format!("MRIStep: jcur = {}\n", jcur as i32));
    outfile.write_str(&format!("MRIStep: ownNLS = {}\n", ownNLS as i32));

    /* output sunrealtype quantities */
    outfile.write_str("MRIStep: Coupling structure:\n");
    let MRIC = { mriStep_mem_mut(ark_mem).MRIC.clone() };
    if let Some(MRIC) = MRIC {
        MRIStepCoupling_Write(&MRIC, outfile);
    }

    let (gamma, gammap, gamrat, crate_, delnrm_p, eRNrm, nlscoef, crdown, rdiv, dgmax) = {
        let s = mriStep_mem_mut(ark_mem);
        (
            s.gamma, s.gammap, s.gamrat, s.crate_, s.delnrm_p, s.eRNrm, s.nlscoef, s.crdown,
            s.rdiv, s.dgmax,
        )
    };
    outfile.write_str(&format!("MRIStep: gamma = {}\n", sun_format_g(gamma)));
    outfile.write_str(&format!("MRIStep: gammap = {}\n", sun_format_g(gammap)));
    outfile.write_str(&format!("MRIStep: gamrat = {}\n", sun_format_g(gamrat)));
    outfile.write_str(&format!("MRIStep: crate = {}\n", sun_format_g(crate_)));
    outfile.write_str(&format!("MRIStep: delnrm_p = {}\n", sun_format_g(delnrm_p)));
    outfile.write_str(&format!("MRIStep: eRNrm = {}\n", sun_format_g(eRNrm)));
    outfile.write_str(&format!("MRIStep: nlscoef = {}\n", sun_format_g(nlscoef)));
    outfile.write_str(&format!("MRIStep: crdown = {}\n", sun_format_g(crdown)));
    outfile.write_str(&format!("MRIStep: rdiv = {}\n", sun_format_g(rdiv)));
    outfile.write_str(&format!("MRIStep: dgmax = {}\n", sun_format_g(dgmax)));

    let (nstages_active, Ae_row, Ai_row) = {
        let s = mriStep_mem_mut(ark_mem);
        (s.nstages_active, s.Ae_row.clone(), s.Ai_row.clone())
    };
    outfile.write_str("MRIStep: Ae_row =");
    for i in 0..nstages_active {
        outfile.write_str(&format!(" {}", sun_format_g(Ae_row[i as usize])));
    }
    outfile.write_str("\n");
    outfile.write_str("MRIStep: Ai_row =");
    for i in 0..nstages_active {
        outfile.write_str(&format!(" {}", sun_format_g(Ai_row[i as usize])));
    }
    outfile.write_str("\n");

    /* SUNDIALS_DEBUG_PRINTVEC vector output is not enabled in this build */

    /* print the inner stepper memory */
    let stepper = { mriStep_mem_mut(ark_mem).stepper.clone() }.expect("stepper set");
    mriStepInnerStepper_PrintMem(&stepper, outfile);
}

/*---------------------------------------------------------------
  mriStep_AttachLinsol:

  This routine attaches the various set of system linear solver
  interface routines, data structure, and solver type to the
  MRIStep module.
  ---------------------------------------------------------------*/
pub fn mriStep_AttachLinsol(
    ark_mem: &ARKodeMem,
    linit: Option<ARKLinsolInitFn>,
    lsetup: Option<ARKLinsolSetupFn>,
    lsolve: Option<ARKLinsolSolveFn>,
    lfree: Option<ARKLinsolFreeFn>,
    _lsolve_type: SUNLinearSolver_Type,
    lmem: Option<Box<dyn Any>>,
) -> i32 {
    /* access ARKodeMRIStepMem structure */
    let retval = mriStep_step_mem_ok(ark_mem, "mriStep_AttachLinsol");
    if retval != ARK_SUCCESS {
        return retval;
    }

    /* free any existing system solver */
    let old_lfree = { mriStep_mem_mut(ark_mem).lfree };
    if let Some(old_lfree) = old_lfree {
        let _ = old_lfree(ark_mem);
    }

    /* Attach the provided routines, data structure and solve type */
    {
        let mut s = mriStep_mem_mut(ark_mem);
        s.linit = linit;
        s.lsetup = lsetup;
        s.lsolve = lsolve;
        s.lfree = lfree;

        /* Reset all linear solver counters */
        s.nsetups = 0;
        s.nstlp = 0;
    }
    /* C `step_mem->lmem = lmem`: the record is owned by ark_mem (§4) */
    ark_mem.borrow_mut().ark_lmem = lmem;

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_DisableLSetup:

  This routine NULLifies the lsetup function pointer in the
  MRIStep module.
  ---------------------------------------------------------------*/
pub fn mriStep_DisableLSetup(ark_mem: &ARKodeMem) {
    /* access ARKodeMRIStepMem structure */
    let retval = mriStep_step_mem_ok(ark_mem, "mriStep_DisableLSetup");
    if retval != ARK_SUCCESS {
        return;
    }

    /* nullify the lsetup function pointer */
    mriStep_mem_mut(ark_mem).lsetup = None;
}

/*---------------------------------------------------------------
  mriStep_GetLmem:

  This routine returns the system linear solver interface memory
  structure, lmem.

  Seam (§4): the ARKLS record is stored BY VALUE in `ark_mem.ark_lmem`, so
  this reports PRESENCE; `arkls_mem_mut(ark_mem)` reaches the record.
  ---------------------------------------------------------------*/
pub fn mriStep_GetLmem(ark_mem: &ARKodeMem) -> sunbooleantype {
    /* access ARKodeMRIStepMem structure, and return lmem */
    let retval = mriStep_step_mem_ok(ark_mem, "mriStep_GetLmem");
    if retval != ARK_SUCCESS {
        return SUNFALSE;
    }
    ark_mem.borrow().ark_lmem.is_some()
}

/*---------------------------------------------------------------
  mriStep_GetImplicitRHS:

  This routine returns the implicit RHS function pointer, fi.
  ---------------------------------------------------------------*/
pub fn mriStep_GetImplicitRHS(ark_mem: &ARKodeMem) -> Option<ARKRhsFn> {
    /* access ARKodeMRIStepMem structure, and return fi */
    let retval = mriStep_step_mem_ok(ark_mem, "mriStep_GetImplicitRHS");
    if retval != ARK_SUCCESS {
        return None;
    }
    let s = mriStep_mem_mut(ark_mem);
    if s.implicit_rhs {
        s.fsi
    } else {
        None
    }
}

/*---------------------------------------------------------------
  mriStep_GetGammas:

  This routine fills the current value of gamma, and states
  whether the gamma ratio fails the dgmax criteria.
  ---------------------------------------------------------------*/
pub fn mriStep_GetGammas(
    ark_mem: &ARKodeMem,
    gamma: &mut sunrealtype,
    gamrat: &mut sunrealtype,
    jcur: &mut Option<ARKJcurPtr>,
    dgamma_fail: &mut sunbooleantype,
) -> i32 {
    /* access ARKodeMRIStepMem structure */
    let retval = mriStep_step_mem_ok(ark_mem, "mriStep_GetGammas");
    if retval != ARK_SUCCESS {
        return retval;
    }

    /* set outputs */
    let (s_gamma, s_gamrat, s_jcur, s_dgmax) = {
        let s = mriStep_mem_mut(ark_mem);
        (s.gamma, s.gamrat, s.jcur.clone(), s.dgmax)
    };
    *gamma = s_gamma;
    *gamrat = s_gamrat;
    *jcur = Some(s_jcur);
    *dgamma_fail = SUNRabs(*gamrat - ONE) >= s_dgmax;

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_Init:

  This routine is called just prior to performing internal time
  steps (after all user "set" routines have been called) from
  within arkInitialSetup.

  With initialization type RESET_INIT, this routine does nothing.

  For other initialization types, this routine:
  - initializes and sets up the linear and nonlinear solvers
    (if applicable)
  - initializes and sets up the nonlinear solver (if applicable)
  - performs timestep adaptivity checks and initial setup,
    including setting the initial time step size if needed
  - sets the relevant TakeStep routine based on the current
    problem configuration
  - sets/checks the coefficient tables to be used
  - allocates any internal memory that depends on the MRI method
    structure or solver options

  With other initialization types, this routine does nothing.
  ---------------------------------------------------------------*/
pub fn mriStep_Init(ark_mem: &ARKodeMem, init_type: i32) -> i32 {
    let mut retval: i32;

    /* access ARKodeMRIStepMem structure */
    retval = mriStep_step_mem_ok(ark_mem, "mriStep_Init");
    if retval != ARK_SUCCESS {
        return retval;
    }

    /* immediately return if reset */
    if init_type == RESET_INIT {
        return ARK_SUCCESS;
    }

    /* initializations/checks for (re-)initialization call */
    if init_type == FIRST_INIT {
        /* enforce use of arkEwtSmallReal if using a fixed step size for
           an explicit method, an internal error weight function, and not performing
           accumulated temporal error estimation */
        let mut reset_efun: sunbooleantype = SUNTRUE;
        let implicit_rhs = { mriStep_mem_mut(ark_mem).implicit_rhs };
        if implicit_rhs {
            reset_efun = SUNFALSE;
        }
        let (fixedstep, user_efun, accum_type) = {
            let m = ark_mem.borrow();
            (m.fixedstep, m.user_efun, m.AccumErrorType)
        };
        if !fixedstep {
            reset_efun = SUNFALSE;
        }
        if user_efun {
            reset_efun = SUNFALSE;
        }
        if accum_type != ARK_ACCUMERROR_NONE {
            reset_efun = SUNFALSE;
        }
        if reset_efun {
            {
                let mut m = ark_mem.borrow_mut();
                m.user_efun = SUNFALSE;
                m.efun = Some(arkEwtSetSmallReal);
            }
            /* C `ark_mem->e_data = ark_mem`: a boxed handle clone playing the
               same role (the Rc cycle is broken in ARKodeFree) */
            let token: Box<dyn Any> = Box::new(ark_mem.clone());
            ark_mem.borrow_mut().e_data = Some(token);
        }

        /* Create coupling structure (if not already set) */
        retval = mriStep_SetCoupling(ark_mem);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!() as i32,
                "mriStep_Init",
                file!(),
                "Could not create coupling table",
            );
            return ARK_ILL_INPUT;
        }

        /* Check that coupling structure is OK */
        retval = mriStep_CheckCoupling(ark_mem);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!() as i32,
                "mriStep_Init",
                file!(),
                "Error in coupling table",
            );
            return ARK_ILL_INPUT;
        }

        let MRIC = { mriStep_mem_mut(ark_mem).MRIC.clone() }.expect("MRIC set");

        /* Attach correct TakeStep routine for this coupling table.
           (The C `default:` "Unknown method type" branch is unreachable:
           MRISTEP_METHOD_TYPE is a Rust enum.) */
        let mric_type = { MRIC.borrow().type_ };
        match mric_type {
            MRISTEP_EXPLICIT => ark_mem.borrow_mut().step = Some(mriStep_TakeStepMRIGARK),
            MRISTEP_IMPLICIT => ark_mem.borrow_mut().step = Some(mriStep_TakeStepMRIGARK),
            MRISTEP_IMEX => ark_mem.borrow_mut().step = Some(mriStep_TakeStepMRIGARK),
            MRISTEP_MERK => ark_mem.borrow_mut().step = Some(mriStep_TakeStepMERK),
            MRISTEP_SR => ark_mem.borrow_mut().step = Some(mriStep_TakeStepMRISR),
        }

        /* Request arkode ensure that ycur==yn upon entry to TakeStep function */
        ark_mem.borrow_mut().ensure_ycur = SUNTRUE;

        /* Retrieve/store method and embedding orders now that tables are finalized */
        let (mric_stages, mric_q, mric_p, mric_nmat) = {
            let c = MRIC.borrow();
            (c.stages, c.q, c.p, c.nmat)
        };
        {
            let mut s = mriStep_mem_mut(ark_mem);
            s.stages = mric_stages;
            s.q = mric_q;
            s.p = mric_p;
        }
        {
            let mut m = ark_mem.borrow_mut();
            let hadapt_mem = m.hadapt_mem.as_mut().expect("hadapt_mem set");
            hadapt_mem.q = mric_q;
            hadapt_mem.p = mric_p;
        }

        /* Ensure that if adaptivity or error accumulation is enabled, then
           method includes embedding coefficients */
        let (fixedstep, accum_type) = {
            let m = ark_mem.borrow();
            (m.fixedstep, m.AccumErrorType)
        };
        let p = { mriStep_mem_mut(ark_mem).p };
        if (!fixedstep || (accum_type != ARK_ACCUMERROR_NONE)) && (p <= 0) {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!() as i32,
                "mriStep_Init",
                file!(),
                "Temporal error estimation cannot be performed without embedding coefficients",
            );
            return ARK_ILL_INPUT;
        }

        /* allocate/fill derived quantities from MRIC structure */

        /* stage map */
        let (have_stage_map, stages) = {
            let s = mriStep_mem_mut(ark_mem);
            (!s.stage_map.is_empty(), s.stages)
        };
        if have_stage_map {
            mriStep_mem_mut(ark_mem).stage_map = Vec::new();
            ark_mem.borrow_mut().liw -= stages as i64;
        }
        mriStep_mem_mut(ark_mem).stage_map = vec![0i32; mric_stages as usize];
        ark_mem.borrow_mut().liw += mric_stages as i64;
        {
            let (mut stage_map, mut nstages_active) = {
                let mut s = mriStep_mem_mut(ark_mem);
                (std::mem::take(&mut s.stage_map), s.nstages_active)
            };
            retval = mriStepCoupling_GetStageMap(&MRIC, &mut stage_map, &mut nstages_active);
            let mut s = mriStep_mem_mut(ark_mem);
            s.stage_map = stage_map;
            s.nstages_active = nstages_active;
        }
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!() as i32,
                "mriStep_Init",
                file!(),
                "Error in coupling table",
            );
            return ARK_ILL_INPUT;
        }

        /* stage types */
        let have_stagetypes = { !mriStep_mem_mut(ark_mem).stagetypes.is_empty() };
        if have_stagetypes {
            mriStep_mem_mut(ark_mem).stagetypes = Vec::new();
            ark_mem.borrow_mut().liw -= stages as i64;
        }
        mriStep_mem_mut(ark_mem).stagetypes = vec![0i32; (mric_stages + 1) as usize];
        ark_mem.borrow_mut().liw += (mric_stages + 1) as i64;
        for j in 0..=mric_stages {
            let stagetype = mriStepCoupling_GetStageType(&MRIC, j);
            mriStep_mem_mut(ark_mem).stagetypes[j as usize] = stagetype;
        }

        /* explicit RK coefficient row */
        let have_Ae_row = { !mriStep_mem_mut(ark_mem).Ae_row.is_empty() };
        if have_Ae_row {
            mriStep_mem_mut(ark_mem).Ae_row = Vec::new();
            ark_mem.borrow_mut().lrw -= stages as i64;
        }
        mriStep_mem_mut(ark_mem).Ae_row = vec![ZERO; mric_stages as usize];
        ark_mem.borrow_mut().lrw += mric_stages as i64;

        /* implicit RK coefficient row */
        let have_Ai_row = { !mriStep_mem_mut(ark_mem).Ai_row.is_empty() };
        if have_Ai_row {
            mriStep_mem_mut(ark_mem).Ai_row = Vec::new();
            ark_mem.borrow_mut().lrw -= stages as i64;
        }
        mriStep_mem_mut(ark_mem).Ai_row = vec![ZERO; mric_stages as usize];
        ark_mem.borrow_mut().lrw += mric_stages as i64;

        /* Allocate reusable arrays for fused vector operations */
        let nforcing = { mriStep_mem_mut(ark_mem).nforcing };
        let fused_workspace_size: i32 = SUNMAX(3, 2 * mric_stages + 2 + nforcing);

        let nfusedopvecs = { mriStep_mem_mut(ark_mem).nfusedopvecs };
        if nfusedopvecs < fused_workspace_size {
            let (have_cvals, have_Xvecs) = {
                let s = mriStep_mem_mut(ark_mem);
                (!s.cvals.is_empty(), !s.Xvecs.is_empty())
            };
            if have_cvals {
                mriStep_mem_mut(ark_mem).cvals = Vec::new();
                ark_mem.borrow_mut().lrw -= nfusedopvecs as i64;
            }
            if have_Xvecs {
                mriStep_mem_mut(ark_mem).Xvecs = Vec::new();
                ark_mem.borrow_mut().liw -= nfusedopvecs as i64;
            }
            {
                let mut s = mriStep_mem_mut(ark_mem);
                s.nfusedopvecs = 0;

                /* The C calloc-failure branches are unreachable: a Rust
                   allocation failure aborts rather than returning NULL. */
                s.cvals = vec![ZERO; fused_workspace_size as usize];
                s.Xvecs = vec![None; fused_workspace_size as usize];
                s.nfusedopvecs = fused_workspace_size;
            }
            let mut m = ark_mem.borrow_mut();
            m.lrw += fused_workspace_size as i64;
            m.liw += fused_workspace_size as i64;
        }

        /* Retrieve/store method and embedding orders now that tables are finalized */
        {
            let mut s = mriStep_mem_mut(ark_mem);
            s.stages = mric_stages;
            s.q = mric_q;
            s.p = mric_p;

            /* If an MRISR method is applied to a non-ImEx problem, we "unify"
               the Fse and Fsi vectors to point at the same memory */
            s.unify_Fs = SUNFALSE;
            if (mric_type == MRISTEP_SR)
                && ((s.explicit_rhs && !s.implicit_rhs) || (!s.explicit_rhs && s.implicit_rhs))
            {
                s.unify_Fs = SUNTRUE;
            }
        }

        /* Allocate MRI RHS vector memory, update storage requirements */
        /*   Allocate Fse[0] ... Fse[nstages_active - 1] and           */
        /*   Fsi[0] ... Fsi[nstages_active - 1] if needed              */
        let (nstages_allocated, nstages_active, explicit_rhs, implicit_rhs, unify_Fs) = {
            let s = mriStep_mem_mut(ark_mem);
            (
                s.nstages_allocated,
                s.nstages_active,
                s.explicit_rhs,
                s.implicit_rhs,
                s.unify_Fs,
            )
        };
        if nstages_allocated < nstages_active {
            if nstages_allocated != 0 {
                if explicit_rhs {
                    let (lrw1, liw1, mut lrw, mut liw) = {
                        let m = ark_mem.borrow();
                        (m.lrw1, m.liw1, m.lrw, m.liw)
                    };
                    let mut Fse = { std::mem::take(&mut mriStep_mem_mut(ark_mem).Fse) };
                    arkFreeVecArray(nstages_allocated, &mut Fse, lrw1, &mut lrw, liw1, &mut liw);
                    {
                        let mut m = ark_mem.borrow_mut();
                        m.lrw = lrw;
                        m.liw = liw;
                    }
                    let mut s = mriStep_mem_mut(ark_mem);
                    s.Fse = Fse;
                    if s.unify_Fs {
                        s.Fsi = Vec::new();
                    }
                }
                if implicit_rhs {
                    let (lrw1, liw1, mut lrw, mut liw) = {
                        let m = ark_mem.borrow();
                        (m.lrw1, m.liw1, m.lrw, m.liw)
                    };
                    let mut Fsi = { std::mem::take(&mut mriStep_mem_mut(ark_mem).Fsi) };
                    arkFreeVecArray(nstages_allocated, &mut Fsi, lrw1, &mut lrw, liw1, &mut liw);
                    {
                        let mut m = ark_mem.borrow_mut();
                        m.lrw = lrw;
                        m.liw = liw;
                    }
                    let mut s = mriStep_mem_mut(ark_mem);
                    s.Fsi = Fsi;
                    if s.unify_Fs {
                        s.Fse = Vec::new();
                    }
                }
            }
            let ewt = ark_mem.borrow().ewt.clone().expect("ewt set");
            if explicit_rhs && !unify_Fs {
                let (lrw1, liw1, mut lrw, mut liw) = {
                    let m = ark_mem.borrow();
                    (m.lrw1, m.liw1, m.lrw, m.liw)
                };
                let mut Fse = { std::mem::take(&mut mriStep_mem_mut(ark_mem).Fse) };
                let ok = arkAllocVecArray(
                    nstages_active,
                    &ewt,
                    &mut Fse,
                    lrw1,
                    &mut lrw,
                    liw1,
                    &mut liw,
                );
                {
                    let mut m = ark_mem.borrow_mut();
                    m.lrw = lrw;
                    m.liw = liw;
                }
                mriStep_mem_mut(ark_mem).Fse = Fse;
                if !ok {
                    return ARK_MEM_FAIL;
                }
            }
            if implicit_rhs && !unify_Fs {
                let (lrw1, liw1, mut lrw, mut liw) = {
                    let m = ark_mem.borrow();
                    (m.lrw1, m.liw1, m.lrw, m.liw)
                };
                let mut Fsi = { std::mem::take(&mut mriStep_mem_mut(ark_mem).Fsi) };
                let ok = arkAllocVecArray(
                    nstages_active,
                    &ewt,
                    &mut Fsi,
                    lrw1,
                    &mut lrw,
                    liw1,
                    &mut liw,
                );
                {
                    let mut m = ark_mem.borrow_mut();
                    m.lrw = lrw;
                    m.liw = liw;
                }
                mriStep_mem_mut(ark_mem).Fsi = Fsi;
                if !ok {
                    return ARK_MEM_FAIL;
                }
            }
            if unify_Fs {
                let (lrw1, liw1, mut lrw, mut liw) = {
                    let m = ark_mem.borrow();
                    (m.lrw1, m.liw1, m.lrw, m.liw)
                };
                let mut Fse = { std::mem::take(&mut mriStep_mem_mut(ark_mem).Fse) };
                let ok = arkAllocVecArray(
                    nstages_active,
                    &ewt,
                    &mut Fse,
                    lrw1,
                    &mut lrw,
                    liw1,
                    &mut liw,
                );
                {
                    let mut m = ark_mem.borrow_mut();
                    m.lrw = lrw;
                    m.liw = liw;
                }
                mriStep_mem_mut(ark_mem).Fse = Fse;
                if !ok {
                    return ARK_MEM_FAIL;
                }
                let mut s = mriStep_mem_mut(ark_mem);
                s.Fsi = s.Fse.clone();
            }

            mriStep_mem_mut(ark_mem).nstages_allocated = nstages_active;
        }

        /* if any slow stage is implicit, allocate sdata, zpred, zcor vectors;
           if all stages explicit, free default NLS object, and detach all
           linear solver routines.  Note: step_mem->implicit_rhs will only equal
           SUNTRUE if an implicit table has been user-provided. */
        if implicit_rhs {
            let ewt = ark_mem.borrow().ewt.clone().expect("ewt set");
            {
                let mut sdata = { mriStep_mem_mut(ark_mem).sdata.clone() };
                let ok = arkAllocVec(ark_mem, &ewt, &mut sdata);
                mriStep_mem_mut(ark_mem).sdata = sdata;
                if !ok {
                    return ARK_MEM_FAIL;
                }
            }
            {
                let mut zpred = { mriStep_mem_mut(ark_mem).zpred.clone() };
                let ok = arkAllocVec(ark_mem, &ewt, &mut zpred);
                mriStep_mem_mut(ark_mem).zpred = zpred;
                if !ok {
                    return ARK_MEM_FAIL;
                }
            }
            {
                let mut zcor = { mriStep_mem_mut(ark_mem).zcor.clone() };
                let ok = arkAllocVec(ark_mem, &ewt, &mut zcor);
                mriStep_mem_mut(ark_mem).zcor = zcor;
                if !ok {
                    return ARK_MEM_FAIL;
                }
            }
        } else {
            let free_nls = {
                let s = mriStep_mem_mut(ark_mem);
                s.NLS.is_some() && s.ownNLS
            };
            if free_nls {
                let nls = { mriStep_mem_mut(ark_mem).NLS.take() };
                let _ = SUNNonlinSolFree(nls);
                let mut s = mriStep_mem_mut(ark_mem);
                s.NLS = None;
                s.ownNLS = SUNFALSE;
            }
            {
                let mut s = mriStep_mem_mut(ark_mem);
                s.linit = None;
                s.lsetup = None;
                s.lsolve = None;
                s.lfree = None;
            }
            /* C `step_mem->lmem = NULL` (the record lives in ark_mem, §4) */
            ark_mem.borrow_mut().ark_lmem = None;
        }

        /* Allocate inner stepper data */
        let stepper = { mriStep_mem_mut(ark_mem).stepper.clone() }.expect("stepper set");
        let ewt = ark_mem.borrow().ewt.clone().expect("ewt set");
        retval = mriStepInnerStepper_AllocVecs(&stepper, mric_nmat, &ewt);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!() as i32,
                "mriStep_Init",
                file!(),
                "Error allocating inner stepper memory",
            );
            return ARK_MEM_FAIL;
        }

        /* Override the interpolant degree (if needed), used in arkInitialSetup */
        let q = { mriStep_mem_mut(ark_mem).q };
        let interp_degree = ark_mem.borrow().interp_degree;
        if q > 1 && interp_degree > (q - 1) {
            /* Limit max degree to at most one less than the method global order */
            ark_mem.borrow_mut().interp_degree = q - 1;
        } else if q == 1 && interp_degree > 1 {
            /* Allow for linear interpolant with first order methods to ensure
               solution values are returned at the time interval end points */
            ark_mem.borrow_mut().interp_degree = 1;
        }

        /* Higher-order predictors require interpolation */
        let interp_type = ark_mem.borrow().interp_type;
        let predictor = { mriStep_mem_mut(ark_mem).predictor };
        if interp_type == ARK_INTERP_NONE && predictor != 0 {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!() as i32,
                "mriStep_Init",
                file!(),
                "Non-trival predictors require an interpolation module",
            );
            return ARK_ILL_INPUT;
        }
    }

    /* Call linit (if it exists) */
    let linit = { mriStep_mem_mut(ark_mem).linit };
    if let Some(linit) = linit {
        retval = linit(ark_mem);
        if retval != 0 {
            arkProcessError(
                Some(ark_mem),
                ARK_LINIT_FAIL,
                line!() as i32,
                "mriStep_Init",
                file!(),
                MSG_ARK_LINIT_FAIL,
            );
            return ARK_LINIT_FAIL;
        }
    }

    /* Initialize the nonlinear solver object (if it exists) */
    let have_nls = { mriStep_mem_mut(ark_mem).NLS.is_some() };
    if have_nls {
        retval = mriStep_NlsInit(ark_mem);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_NLS_INIT_FAIL,
                line!() as i32,
                "mriStep_Init",
                file!(),
                "Unable to initialize SUNNonlinearSolver object",
            );
            return ARK_NLS_INIT_FAIL;
        }
    }

    /* get timestep adaptivity type */
    let hcontroller = ark_mem
        .borrow()
        .hadapt_mem
        .as_ref()
        .expect("hadapt_mem set")
        .hcontroller
        .clone()
        .expect("hcontroller set");
    let adapt_type = SUNAdaptController_GetType(&hcontroller);

    let fixedstep = ark_mem.borrow().fixedstep;
    if fixedstep {
        /* Fixed step sizes: user must supply the initial step size */
        if ark_mem.borrow().hin == ZERO {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!() as i32,
                "mriStep_Init",
                file!(),
                "Timestep adaptivity disabled, but missing user-defined fixed stepsize",
            );
            return ARK_ILL_INPUT;
        }
    } else {
        /* ensure that a compatible adaptivity controller is provided */
        if (adapt_type != SUN_ADAPTCONTROLLER_MRI_H_TOL) && (adapt_type != SUN_ADAPTCONTROLLER_H)
        {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!() as i32,
                "mriStep_Init",
                file!(),
                "SUNAdaptController type is unsupported by MRIStep",
            );
            return ARK_ILL_INPUT;
        }

        /* Controller provides adaptivity (at least at the slow time scale):
           - verify that the MRI method includes an embedding, and
           - estimate initial slow step size (store in ark_mem->hin) */
        let mric_p = {
            let MRIC = { mriStep_mem_mut(ark_mem).MRIC.clone() }.expect("MRIC set");
            let p = MRIC.borrow().p;
            p
        };
        if mric_p <= 0 {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!() as i32,
                "mriStep_Init",
                file!(),
                "Timestep adaptivity enabled, but non-embedded MRI table specified",
            );
            return ARK_ILL_INPUT;
        }
    }

    /* Perform additional setup for (H,tol) controller */
    if adapt_type == SUN_ADAPTCONTROLLER_MRI_H_TOL {
        /* Verify that adaptivity type is supported by inner stepper */
        let stepper = { mriStep_mem_mut(ark_mem).stepper.clone() }.expect("stepper set");
        if !mriStepInnerStepper_SupportsRTolAdaptivity(&stepper) {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!() as i32,
                "mriStep_Init",
                file!(),
                "MRI H-TOL SUNAdaptController provided, but unsupported by inner stepper",
            );
            return ARK_ILL_INPUT;
        }

        /* initialize fast stepper to use the same relative tolerance as MRIStep */
        mriStep_mem_mut(ark_mem).inner_rtol_factor = ONE;
    }

    ARK_SUCCESS
}

/*------------------------------------------------------------------------------
  mriStep_ComputeH0:

  This utility routine computes the initial slow step size for MRI methods.

  It is assumed that the IVP is defined by multiple RHS functions,
     y'(t) = f(t,y) = fs(t,y)  + ff(t,y),
  where fs corresponds to dynamics that should be evolved directly by MRIStep,
  and ff corresponds to dynamics that will be evolved by an inner stepper.
  ----------------------------------------------------------------------------*/
pub fn mriStep_ComputeH0(ark_mem: &ARKodeMem, tout: sunrealtype, hin: &mut sunrealtype) -> i32 {
    let mut retval: i32;

    /*   tempv1 = fs(t0, y0) */
    let (tn, yn, tempv1) = {
        let m = ark_mem.borrow();
        (
            m.tn,
            m.yn.clone().expect("yn set"),
            m.tempv1.clone().expect("tempv1 set"),
        )
    };
    if mriStep_SlowRHS(ark_mem, tn, &yn, &tempv1, ARK_FULLRHS_START) != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            ARK_RHSFUNC_FAIL,
            line!() as i32,
            "mriStep_ComputeH0",
            file!(),
            "error calling slow RHS function(s)",
        );
        return ARK_RHSFUNC_FAIL;
    }
    retval = mriStep_Hin(ark_mem, tn, tout, &tempv1, hin);
    if retval != ARK_SUCCESS {
        retval = arkHandleFailure(ark_mem, retval);
        return retval;
    }

    ARK_SUCCESS
}

/*------------------------------------------------------------------------------
  mriStep_FullRHS:

  This is just a wrapper to call the user-supplied RHS functions,
  f(t,y) = fse(t,y) + fsi(t,y)  + ff(t,y).

  Note: this relies on the utility routine mriStep_UpdateF0 to update Fse[0]
  and Fsi[0] as appropriate (i.e., leveraging previous evaluations, etc.), and
  merely combines the resulting values together with ff to construct the output.

  However, in ARK_FULLRHS_OTHER mode, this routine must call all slow RHS
  functions directly, since that mode cannot reuse internally stored values.

   ARK_FULLRHS_OTHER -> called in the following circumstances:
                        (a) when estimating the initial time step size,
                        (b) for high-order dense output with the Hermite
                            interpolation module,
                        (c) by an "outer" stepper when MRIStep is used as an
                            inner solver), or
                        (d) when a high-order implicit predictor is requested from
                            the Hermite interpolation module within the time step
                            t_{n} \to t_{n+1}.

                        While instances (a)-(c) will occur in-between MRIStep time
                        steps, instance (d) can occur at the start of each internal
                        MRIStep stage.  Since the (t,y) input does not correspond
                        to an "official" time step, thus the RHS functions should
                        always be evaluated, and the values should *not* be stored
                        anywhere that will interfere with other reused MRIStep data
                        from one stage to the next (but it may use nonlinear solver
                        scratch space).

  Note that this routine always calls the fast RHS function, ff(t,y), in
  ARK_FULLRHS_OTHER mode.
  ----------------------------------------------------------------------------*/
pub fn mriStep_FullRHS(
    ark_mem: &ARKodeMem,
    t: sunrealtype,
    y: &N_Vector,
    f: &N_Vector,
    mode: i32,
) -> i32 {
    let mut nvec: i32;
    let mut retval: i32;

    /* access ARKodeMRIStepMem structure */
    retval = mriStep_step_mem_ok(ark_mem, "mriStep_FullRHS");
    if retval != ARK_SUCCESS {
        return retval;
    }

    let stepper = { mriStep_mem_mut(ark_mem).stepper.clone() }.expect("stepper set");

    /* ensure that inner stepper provides fullrhs function */
    let has_fullrhs = stepper.ops.borrow().fullrhs.is_some();
    if !has_fullrhs {
        arkProcessError(
            Some(ark_mem),
            ARK_RHSFUNC_FAIL,
            line!() as i32,
            "mriStep_FullRHS",
            file!(),
            MSG_ARK_MISSING_FULLRHS,
        );
        return ARK_RHSFUNC_FAIL;
    }

    /* perform RHS functions contingent on 'mode' argument */
    if mode == ARK_FULLRHS_START || mode == ARK_FULLRHS_END {
        /* update the internal storage for Fse[0] and Fsi[0] */
        retval = mriStep_UpdateF0(ark_mem, t, y, mode);
        if retval != 0 {
            arkProcessError(
                Some(ark_mem),
                ARK_RHSFUNC_FAIL,
                line!() as i32,
                "mriStep_FullRHS",
                file!(),
                &MSG_ARK_RHSFUNC_FAILED(t),
            );
            return ARK_RHSFUNC_FAIL;
        }

        /* evaluate fast component */
        retval = mriStepInnerStepper_FullRhs(&stepper, t, y, f, ARK_FULLRHS_OTHER);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_RHSFUNC_FAIL,
                line!() as i32,
                "mriStep_FullRHS",
                file!(),
                &MSG_ARK_RHSFUNC_FAILED(t),
            );
            return ARK_RHSFUNC_FAIL;
        }

        /* combine RHS vectors into output */
        let (explicit_rhs, implicit_rhs) = {
            let s = mriStep_mem_mut(ark_mem);
            (s.explicit_rhs, s.implicit_rhs)
        };
        if explicit_rhs && implicit_rhs {
            /* ImEx */
            let (cvals, Xvecs) = {
                let mut s = mriStep_mem_mut(ark_mem);
                s.cvals[0] = ONE;
                s.Xvecs[0] = Some(f.clone());
                s.cvals[1] = ONE;
                let v = s.Fse[0].clone();
                s.Xvecs[1] = Some(v);
                s.cvals[2] = ONE;
                let v = s.Fsi[0].clone();
                s.Xvecs[2] = Some(v);
                (s.cvals.clone(), mriStep_xvecs(&s, 3))
            };
            nvec = 3;
            let _ = N_VLinearCombination(nvec, &cvals, &Xvecs, f);
        } else if implicit_rhs {
            /* implicit */
            let v = { mriStep_mem_mut(ark_mem).Fsi[0].clone() };
            N_VLinearSum(ONE, &v, ONE, f, f);
        } else {
            /* explicit */
            let v = { mriStep_mem_mut(ark_mem).Fse[0].clone() };
            N_VLinearSum(ONE, &v, ONE, f, f);
        }
    } else if mode == ARK_FULLRHS_OTHER {
        /* compute the fast component (force new RHS computation) */
        nvec = 0;
        retval = mriStepInnerStepper_FullRhs(&stepper, t, y, f, ARK_FULLRHS_OTHER);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_RHSFUNC_FAIL,
                line!() as i32,
                "mriStep_FullRHS",
                file!(),
                &MSG_ARK_RHSFUNC_FAILED(t),
            );
            return ARK_RHSFUNC_FAIL;
        }
        {
            let mut s = mriStep_mem_mut(ark_mem);
            s.cvals[nvec as usize] = ONE;
            s.Xvecs[nvec as usize] = Some(f.clone());
        }
        nvec += 1;

        /* call the user-supplied pre-RHS function (if supplied) */
        if ark_mem.borrow().PreRhsFn.is_some() {
            retval = mriStep_CallPreRhsFn(ark_mem, t, y);
            if retval != 0 {
                return ARK_PRERHSFN_FAIL;
            }
        }

        let (explicit_rhs, implicit_rhs) = {
            let s = mriStep_mem_mut(ark_mem);
            (s.explicit_rhs, s.implicit_rhs)
        };

        /* compute the implicit component and store in sdata */
        if implicit_rhs {
            let sdata = { mriStep_mem_mut(ark_mem).sdata.clone() }.expect("sdata set");
            retval = mriStep_CallFsi(ark_mem, t, y, &sdata);
            mriStep_mem_mut(ark_mem).nfsi += 1;
            if retval != 0 {
                arkProcessError(
                    Some(ark_mem),
                    ARK_RHSFUNC_FAIL,
                    line!() as i32,
                    "mriStep_FullRHS",
                    file!(),
                    &MSG_ARK_RHSFUNC_FAILED(t),
                );
                return ARK_RHSFUNC_FAIL;
            }
            {
                let mut s = mriStep_mem_mut(ark_mem);
                s.cvals[nvec as usize] = ONE;
                s.Xvecs[nvec as usize] = Some(sdata);
            }
            nvec += 1;
        }

        /* compute the explicit component and store in ark_tempv2 */
        if explicit_rhs {
            let tempv2 = ark_mem.borrow().tempv2.clone().expect("tempv2 set");
            retval = mriStep_CallFse(ark_mem, t, y, &tempv2);
            mriStep_mem_mut(ark_mem).nfse += 1;
            if retval != 0 {
                arkProcessError(
                    Some(ark_mem),
                    ARK_RHSFUNC_FAIL,
                    line!() as i32,
                    "mriStep_FullRHS",
                    file!(),
                    &MSG_ARK_RHSFUNC_FAILED(t),
                );
                return ARK_RHSFUNC_FAIL;
            }
            {
                let mut s = mriStep_mem_mut(ark_mem);
                s.cvals[nvec as usize] = ONE;
                s.Xvecs[nvec as usize] = Some(tempv2);
            }
            nvec += 1;
        }

        /* Add external forcing components to linear combination */
        let (expforcing, impforcing) = {
            let s = mriStep_mem_mut(ark_mem);
            (s.expforcing, s.impforcing)
        };
        if expforcing || impforcing {
            let mut s = mriStep_mem_mut(ark_mem);
            mriStep_ApplyForcing(&mut s, t, ONE, &mut nvec);
        }

        /* combine RHS vectors into output */
        let (cvals, Xvecs) = {
            let s = mriStep_mem_mut(ark_mem);
            (s.cvals.clone(), mriStep_xvecs(&s, nvec))
        };
        let _ = N_VLinearCombination(nvec, &cvals, &Xvecs, f);
    } else {
        /* return with RHS failure if unknown mode is passed */
        arkProcessError(
            Some(ark_mem),
            ARK_RHSFUNC_FAIL,
            line!() as i32,
            "mriStep_FullRHS",
            file!(),
            "Unknown full RHS mode",
        );
        return ARK_RHSFUNC_FAIL;
    }

    ARK_SUCCESS
}

/*------------------------------------------------------------------------------
  mriStep_UpdateF0:

  This routine is called by mriStep_FullRHS to update the internal storage for
  Fse[0] and Fsi[0], incorporating forcing from a slower time scale as necessary.
  This supports the ARK_FULLRHS_START and ARK_FULLRHS_END "mode" values
  provided to mriStep_FullRHS, and contains all internal logic regarding whether
  RHS functions must be called, versus if the relevant data can just be copied.

  (See the C source for the full ARK_FULLRHS_START / ARK_FULLRHS_END commentary.)

  The C `ARKodeMRIStepMem step_mem` parameter is dropped: the record lives
  inside `ark_mem` and is reached through `mriStep_mem_mut` (an `&mut` to it
  could not coexist with the `&ARKodeMem` this routine also needs).
  ----------------------------------------------------------------------------*/
pub fn mriStep_UpdateF0(
    ark_mem: &ARKodeMem,
    t: sunrealtype,
    y: &N_Vector,
    mode: i32,
) -> i32 {
    let mut nvec: i32;
    let mut retval: i32;

    /* perform RHS functions contingent on 'mode' argument */
    if mode == ARK_FULLRHS_START {
        /* update the RHS components */

        let (fse_is_current, fsi_is_current, explicit_rhs, implicit_rhs, expforcing, impforcing) = {
            let s = mriStep_mem_mut(ark_mem);
            (
                s.fse_is_current,
                s.fsi_is_current,
                s.explicit_rhs,
                s.implicit_rhs,
                s.expforcing,
                s.impforcing,
            )
        };
        let fn_is_current = ark_mem.borrow().fn_is_current;

        /* call the user-supplied pre-RHS function (if supplied) */
        if ark_mem.borrow().PreRhsFn.is_some()
            && ((!fse_is_current || !fn_is_current) || (!fsi_is_current || !fn_is_current))
        {
            retval = mriStep_CallPreRhsFn(ark_mem, t, y);
            if retval != 0 {
                return ARK_PRERHSFN_FAIL;
            }
        }

        /*   implicit component */
        if implicit_rhs {
            /* if either ARKODE or MRIStep consider Fsi[0] stale, then recompute */
            if !fsi_is_current || !fn_is_current {
                let Fsi0 = { mriStep_mem_mut(ark_mem).Fsi[0].clone() };
                retval = mriStep_CallFsi(ark_mem, t, y, &Fsi0);
                mriStep_mem_mut(ark_mem).nfsi += 1;
                if retval != 0 {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_RHSFUNC_FAIL,
                        line!() as i32,
                        "mriStep_UpdateF0",
                        file!(),
                        &MSG_ARK_RHSFUNC_FAILED(t),
                    );
                    return ARK_RHSFUNC_FAIL;
                }
                mriStep_mem_mut(ark_mem).fsi_is_current = SUNTRUE;

                /* Add external forcing, if applicable */
                if impforcing {
                    let (cvals, Xvecs) = {
                        let mut s = mriStep_mem_mut(ark_mem);
                        s.cvals[0] = ONE;
                        let v = s.Fsi[0].clone();
                        s.Xvecs[0] = Some(v);
                        nvec = 1;
                        mriStep_ApplyForcing(&mut s, t, ONE, &mut nvec);
                        (s.cvals.clone(), mriStep_xvecs(&s, nvec))
                    };
                    let _ = N_VLinearCombination(Xvecs.len() as i32, &cvals, &Xvecs, &Fsi0);
                }
            }
        }

        /*   explicit component */
        if explicit_rhs {
            /* if either ARKODE or MRIStep consider Fse[0] stale, then recompute */
            if !fse_is_current || !fn_is_current {
                let Fse0 = { mriStep_mem_mut(ark_mem).Fse[0].clone() };
                retval = mriStep_CallFse(ark_mem, t, y, &Fse0);
                mriStep_mem_mut(ark_mem).nfse += 1;
                if retval != 0 {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_RHSFUNC_FAIL,
                        line!() as i32,
                        "mriStep_UpdateF0",
                        file!(),
                        &MSG_ARK_RHSFUNC_FAILED(t),
                    );
                    return ARK_RHSFUNC_FAIL;
                }
                mriStep_mem_mut(ark_mem).fse_is_current = SUNTRUE;

                /* Add external forcing, if applicable */
                if expforcing {
                    let (cvals, Xvecs) = {
                        let mut s = mriStep_mem_mut(ark_mem);
                        s.cvals[0] = ONE;
                        let v = s.Fse[0].clone();
                        s.Xvecs[0] = Some(v);
                        nvec = 1;
                        mriStep_ApplyForcing(&mut s, t, ONE, &mut nvec);
                        (s.cvals.clone(), mriStep_xvecs(&s, nvec))
                    };
                    let _ = N_VLinearCombination(Xvecs.len() as i32, &cvals, &Xvecs, &Fse0);
                }
            }
        }
    } else if mode == ARK_FULLRHS_END {
        /* compute the full RHS */
        if !ark_mem.borrow().fn_is_current {
            /* call the user-supplied pre-RHS function (if supplied) */
            if ark_mem.borrow().PreRhsFn.is_some() {
                retval = mriStep_CallPreRhsFn(ark_mem, t, y);
                if retval != 0 {
                    return ARK_PRERHSFN_FAIL;
                }
            }

            let (explicit_rhs, implicit_rhs, expforcing, impforcing) = {
                let s = mriStep_mem_mut(ark_mem);
                (s.explicit_rhs, s.implicit_rhs, s.expforcing, s.impforcing)
            };

            /* compute the implicit component */
            if implicit_rhs {
                let Fsi0 = { mriStep_mem_mut(ark_mem).Fsi[0].clone() };
                retval = mriStep_CallFsi(ark_mem, t, y, &Fsi0);
                mriStep_mem_mut(ark_mem).nfsi += 1;
                if retval != 0 {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_RHSFUNC_FAIL,
                        line!() as i32,
                        "mriStep_UpdateF0",
                        file!(),
                        &MSG_ARK_RHSFUNC_FAILED(t),
                    );
                    return ARK_RHSFUNC_FAIL;
                }
                mriStep_mem_mut(ark_mem).fsi_is_current = SUNTRUE;

                /* Add external forcing, as appropriate */
                if impforcing {
                    let (cvals, Xvecs) = {
                        let mut s = mriStep_mem_mut(ark_mem);
                        s.cvals[0] = ONE;
                        let v = s.Fsi[0].clone();
                        s.Xvecs[0] = Some(v);
                        nvec = 1;
                        mriStep_ApplyForcing(&mut s, t, ONE, &mut nvec);
                        (s.cvals.clone(), mriStep_xvecs(&s, nvec))
                    };
                    let _ = N_VLinearCombination(Xvecs.len() as i32, &cvals, &Xvecs, &Fsi0);
                }
            }

            /* compute the explicit component */
            if explicit_rhs {
                let Fse0 = { mriStep_mem_mut(ark_mem).Fse[0].clone() };
                retval = mriStep_CallFse(ark_mem, t, y, &Fse0);
                mriStep_mem_mut(ark_mem).nfse += 1;
                if retval != 0 {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_RHSFUNC_FAIL,
                        line!() as i32,
                        "mriStep_UpdateF0",
                        file!(),
                        &MSG_ARK_RHSFUNC_FAILED(t),
                    );
                    return ARK_RHSFUNC_FAIL;
                }
                mriStep_mem_mut(ark_mem).fse_is_current = SUNTRUE;

                /* Add external forcing, as appropriate */
                if expforcing {
                    let (cvals, Xvecs) = {
                        let mut s = mriStep_mem_mut(ark_mem);
                        s.cvals[0] = ONE;
                        let v = s.Fse[0].clone();
                        s.Xvecs[0] = Some(v);
                        nvec = 1;
                        mriStep_ApplyForcing(&mut s, t, ONE, &mut nvec);
                        (s.cvals.clone(), mriStep_xvecs(&s, nvec))
                    };
                    let _ = N_VLinearCombination(Xvecs.len() as i32, &cvals, &Xvecs, &Fse0);
                }
            }
        }
    } else {
        /* return with RHS failure if unknown mode is requested */
        arkProcessError(
            Some(ark_mem),
            ARK_RHSFUNC_FAIL,
            line!() as i32,
            "mriStep_UpdateF0",
            file!(),
            "Unknown full RHS mode",
        );
        return ARK_RHSFUNC_FAIL;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_TakeStepMRIGARK:

  This routine serves the primary purpose of the MRIStep module:
  it performs a single MRI step (with embedding, if possible).

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
pub fn mriStep_TakeStepMRIGARK(
    ark_mem: &ARKodeMem,
    dsmPtr: &mut sunrealtype,
    nflagPtr: &mut i32,
) -> i32 {
    let mut is: i32; /* current stage index        */
    /* Like C, `retval` is one reused local: the `_ => {}` match arms below
       (C `switch` cases with no `default:`) deliberately leave the previous
       value in place. */
    let mut retval: i32;
    let mut t0: sunrealtype;
    let mut tf: sunrealtype; /* start/end of each stage    */
    let mut calc_fslow: sunbooleantype;
    let mut need_inner_dsm: sunbooleantype;
    let do_embedding: sunbooleantype;
    let nested_mri: sunbooleantype;
    let mut nvec: i32;

    /* access the MRIStep mem structure */
    retval = mriStep_step_mem_ok(ark_mem, "mriStep_TakeStepMRIGARK");
    if retval != ARK_SUCCESS {
        return retval;
    }

    /* initialize algebraic solver convergence flag to success;
       error estimate to zero */
    *nflagPtr = ARK_SUCCESS;
    *dsmPtr = ZERO;

    /* determine whether embedding stage is needed */
    let (fixedstep, accum_type) = {
        let m = ark_mem.borrow();
        (m.fixedstep, m.AccumErrorType)
    };
    do_embedding = !fixedstep || (accum_type != ARK_ACCUMERROR_NONE);

    /* initialize the current stage index */
    {
        let mut s = mriStep_mem_mut(ark_mem);
        s.istage = 0;
        s.cur_stage = 0;
    }

    let stepper = { mriStep_mem_mut(ark_mem).stepper.clone() }.expect("stepper set");

    /* if MRI adaptivity is enabled: reset fast accumulated error,
       and send appropriate control parameter to the fast integrator */
    let hcontroller = ark_mem
        .borrow()
        .hadapt_mem
        .as_ref()
        .expect("hadapt_mem set")
        .hcontroller
        .clone()
        .expect("hcontroller set");
    let adapt_type = SUNAdaptController_GetType(&hcontroller);
    need_inner_dsm = SUNFALSE;
    if adapt_type == SUN_ADAPTCONTROLLER_MRI_H_TOL {
        need_inner_dsm = SUNTRUE;
        mriStep_mem_mut(ark_mem).inner_dsm = ZERO;
        retval = mriStepInnerStepper_ResetAccumulatedError(&stepper);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_INNERSTEP_FAIL,
                line!() as i32,
                "mriStep_TakeStepMRIGARK",
                file!(),
                "Unable to reset the inner stepper error estimate",
            );
            return ARK_INNERSTEP_FAIL;
        }
        let inner_rtol_factor = { mriStep_mem_mut(ark_mem).inner_rtol_factor };
        let reltol = ark_mem.borrow().reltol;
        retval = mriStepInnerStepper_SetRTol(&stepper, inner_rtol_factor * reltol);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_INNERSTEP_FAIL,
                line!() as i32,
                "mriStep_TakeStepMRIGARK",
                file!(),
                "Unable to set the inner stepper tolerance",
            );
            return ARK_INNERSTEP_FAIL;
        }
    }

    /* for adaptive computations, reset the inner integrator to the beginning of this step */
    if !fixedstep {
        let (tcur, ycur) = {
            let m = ark_mem.borrow();
            (m.tcur, m.ycur.clone().expect("ycur set"))
        };
        retval = mriStepInnerStepper_Reset(&stepper, tcur, &ycur);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_INNERSTEP_FAIL,
                line!() as i32,
                "mriStep_TakeStepMRIGARK",
                file!(),
                "Unable to reset the inner stepper",
            );
            return ARK_INNERSTEP_FAIL;
        }
    }

    /* call nonlinear solver setup if it exists */
    let NLS = { mriStep_mem_mut(ark_mem).NLS.clone() };
    if let Some(NLS) = NLS {
        if NLS.ops.borrow().setup.is_some() {
            let tempv3 = ark_mem.borrow().tempv3.clone().expect("tempv3 set");
            N_VConst(ZERO, &tempv3); /* set guess to 0 */
            let mut nls_mem: Option<Box<dyn Any>> = Some(Box::new(ark_mem.clone()));
            retval = SUNNonlinSolSetup(&NLS, &tempv3, &mut nls_mem);
            if retval < 0 {
                return ARK_NLS_SETUP_FAIL;
            }
            if retval > 0 {
                return ARK_NLS_SETUP_RECVR;
            }
        }
    }

    /* Evaluate the slow RHS functions if needed. NOTE: we decide between calling the
       full RHS function (if ark_mem->fn is non-NULL and MRIStep is not an inner
       integrator) versus just updating the stored values of Fse[0] and Fsi[0].  In
       either case, we use ARK_FULLRHS_START mode because MRIGARK methods do not
       evaluate the RHS functions at the end of the time step (so nothing can be
       leveraged). */
    let (expforcing, impforcing) = {
        let s = mriStep_mem_mut(ark_mem);
        (s.expforcing, s.impforcing)
    };
    nested_mri = expforcing || impforcing;
    let fn_is_null = ark_mem.borrow().fn_.is_none();
    if fn_is_null || nested_mri {
        let (tcur, ycur) = {
            let m = ark_mem.borrow();
            (m.tcur, m.ycur.clone().expect("ycur set"))
        };
        retval = mriStep_UpdateF0(ark_mem, tcur, &ycur, ARK_FULLRHS_START);
        if retval != 0 {
            return ARK_RHSFUNC_FAIL;
        }

        /* For a nested MRI configuration we might still need fn to create a predictor
           but it should be fn only for the current nesting level which is why we use
           UpdateF0 in this case rather than FullRHS */
        let fn_v = ark_mem.borrow().fn_.clone();
        let (explicit_rhs, implicit_rhs) = {
            let s = mriStep_mem_mut(ark_mem);
            (s.explicit_rhs, s.implicit_rhs)
        };
        if fn_v.is_some() && nested_mri && implicit_rhs {
            let fn_v = fn_v.expect("fn set");
            if implicit_rhs && explicit_rhs {
                let (Fsi0, Fse0) = {
                    let s = mriStep_mem_mut(ark_mem);
                    (s.Fsi[0].clone(), s.Fse[0].clone())
                };
                N_VLinearSum(ONE, &Fsi0, ONE, &Fse0, &fn_v);
            } else {
                let Fsi0 = { mriStep_mem_mut(ark_mem).Fsi[0].clone() };
                N_VScale(ONE, &Fsi0, &fn_v);
            }
        }
    } else if !fn_is_null && !ark_mem.borrow().fn_is_current {
        let (tcur, ycur, fn_v) = {
            let m = ark_mem.borrow();
            (
                m.tcur,
                m.ycur.clone().expect("ycur set"),
                m.fn_.clone().expect("fn set"),
            )
        };
        retval = mriStep_FullRHS(ark_mem, tcur, &ycur, &fn_v, ARK_FULLRHS_START);
        if retval != 0 {
            return ARK_RHSFUNC_FAIL;
        }
    }
    ark_mem.borrow_mut().fn_is_current = SUNTRUE;

    /* The first stage is the previous time-step solution, so its RHS
       is the [already-computed] slow RHS from the start of the step */

    let stages = { mriStep_mem_mut(ark_mem).stages };
    let MRIC = { mriStep_mem_mut(ark_mem).MRIC.clone() }.expect("MRIC set");

    /* Loop over remaining internal stages */
    is = 1;
    while is < stages - 1 {
        /* Set relevant stage times (including desired stage time for implicit solves)
           and stage index */
        let (tn, h) = {
            let m = ark_mem.borrow();
            (m.tn, m.h)
        };
        let (c_prev, c_cur) = {
            let c = MRIC.borrow();
            (c.c[(is - 1) as usize], c.c[is as usize])
        };
        t0 = tn + c_prev * h;
        tf = tn + c_cur * h;
        ark_mem.borrow_mut().tcur = tf;
        {
            let mut s = mriStep_mem_mut(ark_mem);
            s.istage = is;
            s.cur_stage = is;
        }

        /* Determine current stage type, and call corresponding routine; the
           vector ark_mem->ycur stores the previous stage solution on input, and
           should store the result of this stage solution on output. */
        let stagetype = { mriStep_mem_mut(ark_mem).stagetypes[is as usize] };
        match stagetype {
            MRISTAGE_ERK_FAST => {
                retval = mriStep_ComputeInnerForcing(ark_mem, is, t0, tf);
                if retval != ARK_SUCCESS {
                    return retval;
                }
                let (ycur, tempv2) = {
                    let m = ark_mem.borrow();
                    (
                        m.ycur.clone().expect("ycur set"),
                        m.tempv2.clone().expect("tempv2 set"),
                    )
                };
                retval = mriStep_StageERKFast(ark_mem, t0, tf, &ycur, &tempv2, need_inner_dsm);
                if retval != ARK_SUCCESS {
                    *nflagPtr = CONV_FAIL;
                }
            }
            MRISTAGE_ERK_NOFAST => {
                retval = mriStep_StageERKNoFast(ark_mem, is);
            }
            MRISTAGE_DIRK_NOFAST => {
                retval = mriStep_StageDIRKNoFast(ark_mem, is, nflagPtr);
            }
            MRISTAGE_DIRK_FAST => {
                retval = mriStep_StageDIRKFast(ark_mem, is, nflagPtr);
            }
            MRISTAGE_STIFF_ACC => {
                retval = ARK_SUCCESS;
            }
            _ => {}
        }
        if retval != ARK_SUCCESS {
            return retval;
        }

        /* apply user-supplied stage postprocessing function (if supplied) */
        if ark_mem.borrow().PostProcessStageFn.is_some() && (stagetype != MRISTAGE_STIFF_ACC) {
            let (tcur, ycur) = {
                let m = ark_mem.borrow();
                (m.tcur, m.ycur.clone().expect("ycur set"))
            };
            retval = mriStep_CallPostProcessStageFn(ark_mem, tcur, &ycur);
            if retval != 0 {
                return ARK_POSTPROCESS_STAGE_FAIL;
            }
        }

        /* conditionally reset the inner integrator with the modified stage solution */
        if stagetype != MRISTAGE_STIFF_ACC {
            let have_postprocess = ark_mem.borrow().PostProcessStageFn.is_some();
            if (stagetype != MRISTAGE_ERK_FAST) || have_postprocess {
                let ycur = ark_mem.borrow().ycur.clone().expect("ycur set");
                retval = mriStepInnerStepper_Reset(&stepper, tf, &ycur);
                if retval != ARK_SUCCESS {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_INNERSTEP_FAIL,
                        line!() as i32,
                        "mriStep_TakeStepMRIGARK",
                        file!(),
                        "Unable to reset the inner stepper",
                    );
                    return ARK_INNERSTEP_FAIL;
                }
            }
        }

        /* Compute updated slow RHS, except:
           1. if the stage is excluded from stage_map
           2. if the next stage has "STIFF_ACC" type, and temporal estimation is disabled */
        calc_fslow = SUNTRUE;
        let (stage_map_is, next_stagetype) = {
            let s = mriStep_mem_mut(ark_mem);
            (s.stage_map[is as usize], s.stagetypes[(is + 1) as usize])
        };
        if stage_map_is == -1 {
            calc_fslow = SUNFALSE;
        }
        if !do_embedding && (next_stagetype == MRISTAGE_STIFF_ACC) {
            calc_fslow = SUNFALSE;
        }
        if calc_fslow {
            let (explicit_rhs, implicit_rhs, deduce_rhs, expforcing, impforcing) = {
                let s = mriStep_mem_mut(ark_mem);
                (
                    s.explicit_rhs,
                    s.implicit_rhs,
                    s.deduce_rhs,
                    s.expforcing,
                    s.impforcing,
                )
            };

            /* call the user-supplied pre-RHS function (if supplied) */
            if ark_mem.borrow().PreRhsFn.is_some() {
                if explicit_rhs
                    || (implicit_rhs
                        && (!deduce_rhs || (stagetype != MRISTAGE_DIRK_NOFAST)))
                {
                    let (tcur, ycur) = {
                        let m = ark_mem.borrow();
                        (m.tcur, m.ycur.clone().expect("ycur set"))
                    };
                    retval = mriStep_CallPreRhsFn(ark_mem, tcur, &ycur);
                    if retval != 0 {
                        return ARK_PRERHSFN_FAIL;
                    }
                }
            }

            /* store implicit slow rhs  */
            if implicit_rhs {
                if !deduce_rhs || (stagetype != MRISTAGE_DIRK_NOFAST) {
                    let (tcur, ycur) = {
                        let m = ark_mem.borrow();
                        (m.tcur, m.ycur.clone().expect("ycur set"))
                    };
                    let Fsi_is =
                        { mriStep_mem_mut(ark_mem).Fsi[stage_map_is as usize].clone() };
                    retval = mriStep_CallFsi(ark_mem, tcur, &ycur, &Fsi_is);
                    mriStep_mem_mut(ark_mem).nfsi += 1;

                    if retval < 0 {
                        return ARK_RHSFUNC_FAIL;
                    }
                    if retval > 0 {
                        return ARK_UNREC_RHSFUNC_ERR;
                    }

                    /* Add external forcing to Fsi, if applicable */
                    if impforcing {
                        let (cvals, Xvecs) = {
                            let mut s = mriStep_mem_mut(ark_mem);
                            s.cvals[0] = ONE;
                            let v = s.Fsi[stage_map_is as usize].clone();
                            s.Xvecs[0] = Some(v);
                            nvec = 1;
                            mriStep_ApplyForcing(&mut s, tf, ONE, &mut nvec);
                            (s.cvals.clone(), mriStep_xvecs(&s, nvec))
                        };
                        let _ = N_VLinearCombination(
                            Xvecs.len() as i32,
                            &cvals,
                            &Xvecs,
                            &Fsi_is,
                        );
                    }
                } else {
                    let (gamma, zcor, sdata, Fsi_is) = {
                        let s = mriStep_mem_mut(ark_mem);
                        (
                            s.gamma,
                            s.zcor.clone().expect("zcor set"),
                            s.sdata.clone().expect("sdata set"),
                            s.Fsi[stage_map_is as usize].clone(),
                        )
                    };
                    N_VLinearSum(ONE / gamma, &zcor, -ONE / gamma, &sdata, &Fsi_is);
                }
            }

            /* store explicit slow rhs */
            if explicit_rhs {
                let (tcur, ycur) = {
                    let m = ark_mem.borrow();
                    (m.tcur, m.ycur.clone().expect("ycur set"))
                };
                let Fse_is = { mriStep_mem_mut(ark_mem).Fse[stage_map_is as usize].clone() };
                retval = mriStep_CallFse(ark_mem, tcur, &ycur, &Fse_is);
                mriStep_mem_mut(ark_mem).nfse += 1;

                if retval < 0 {
                    return ARK_RHSFUNC_FAIL;
                }
                if retval > 0 {
                    return ARK_UNREC_RHSFUNC_ERR;
                }

                /* Add external forcing to Fse, if applicable */
                if expforcing {
                    let (cvals, Xvecs) = {
                        let mut s = mriStep_mem_mut(ark_mem);
                        s.cvals[0] = ONE;
                        let v = s.Fse[stage_map_is as usize].clone();
                        s.Xvecs[0] = Some(v);
                        nvec = 1;
                        mriStep_ApplyForcing(&mut s, tf, ONE, &mut nvec);
                        (s.cvals.clone(), mriStep_xvecs(&s, nvec))
                    };
                    let _ =
                        N_VLinearCombination(Xvecs.len() as i32, &cvals, &Xvecs, &Fse_is);
                }
            }
        } /* compute slow RHS */

        is += 1;
    } /* loop over stages */

    /* perform embedded stage (if needed) */
    if do_embedding {
        is = stages;
        {
            let mut s = mriStep_mem_mut(ark_mem);
            s.istage = is;
            s.cur_stage = is;
        }

        /* Temporarily swap ark_mem->ycur and ark_mem->tempv4 pointers, copying
           data so that both hold the current ark_mem->ycur value.  This ensures
           that during this embedding "stage":
             - ark_mem->ycur will be the correct initial condition for the final stage.
             - ark_mem->tempv4 will hold the embedded solution vector. */
        {
            let (ycur, tempv4) = {
                let m = ark_mem.borrow();
                (
                    m.ycur.clone().expect("ycur set"),
                    m.tempv4.clone().expect("tempv4 set"),
                )
            };
            N_VScale(ONE, &ycur, &tempv4);
        }
        {
            let mut m = ark_mem.borrow_mut();
            let tmp = m.ycur.take();
            m.ycur = m.tempv4.take();
            m.tempv4 = tmp;
        }

        /* Reset ark_mem->tcur as the time value corresponding with the end of the step */
        /* Set relevant stage times (including desired stage time for implicit solves) */
        let (tn, h) = {
            let m = ark_mem.borrow();
            (m.tn, m.h)
        };
        let c_im2 = { MRIC.borrow().c[(is - 2) as usize] };
        t0 = tn + c_im2 * h;
        tf = tn + h;
        ark_mem.borrow_mut().tcur = tf;

        /* Determine embedding stage type, and call corresponding routine; the
           vector ark_mem->ycur stores the previous stage solution on input, and
           should store the result of this stage solution on output. */
        let stagetype = { mriStep_mem_mut(ark_mem).stagetypes[is as usize] };
        match stagetype {
            MRISTAGE_ERK_FAST => {
                retval = mriStep_ComputeInnerForcing(ark_mem, is, t0, tf);
                if retval != ARK_SUCCESS {
                    return retval;
                }
                let (ycur, tempv2) = {
                    let m = ark_mem.borrow();
                    (
                        m.ycur.clone().expect("ycur set"),
                        m.tempv2.clone().expect("tempv2 set"),
                    )
                };
                retval = mriStep_StageERKFast(ark_mem, t0, tf, &ycur, &tempv2, SUNFALSE);
                if retval != ARK_SUCCESS {
                    *nflagPtr = CONV_FAIL;
                }
            }
            MRISTAGE_ERK_NOFAST => {
                retval = mriStep_StageERKNoFast(ark_mem, is);
            }
            MRISTAGE_DIRK_NOFAST => {
                retval = mriStep_StageDIRKNoFast(ark_mem, is, nflagPtr);
            }
            MRISTAGE_DIRK_FAST => {
                retval = mriStep_StageDIRKFast(ark_mem, is, nflagPtr);
            }
            _ => {}
        }
        if retval != ARK_SUCCESS {
            return retval;
        }

        /* Swap back ark_mem->ycur with ark_mem->tempv4, and reset the inner integrator */
        {
            let mut m = ark_mem.borrow_mut();
            let tmp = m.ycur.take();
            m.ycur = m.tempv4.take();
            m.tempv4 = tmp;
        }
        let ycur = ark_mem.borrow().ycur.clone().expect("ycur set");
        retval = mriStepInnerStepper_Reset(&stepper, t0, &ycur);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_INNERSTEP_FAIL,
                line!() as i32,
                "mriStep_TakeStepMRIGARK",
                file!(),
                "Unable to reset the inner stepper",
            );
            return ARK_INNERSTEP_FAIL;
        }
    }

    /* Compute final stage (for evolved solution), along with error estimate */
    {
        is = stages - 1;
        {
            let mut s = mriStep_mem_mut(ark_mem);
            s.istage = is;
            s.cur_stage = is;
        }

        /* Set relevant stage times (including desired stage time for implicit solves) */
        let (tn, h) = {
            let m = ark_mem.borrow();
            (m.tn, m.h)
        };
        let c_im1 = { MRIC.borrow().c[(is - 1) as usize] };
        t0 = tn + c_im1 * h;
        tf = tn + h;
        ark_mem.borrow_mut().tcur = tf;

        /* Determine final stage type, and call corresponding routine; the
           vector ark_mem->ycur stores the previous stage solution on input, and
           should store the result of this stage solution on output. */
        let stagetype = { mriStep_mem_mut(ark_mem).stagetypes[is as usize] };
        match stagetype {
            MRISTAGE_ERK_FAST => {
                retval = mriStep_ComputeInnerForcing(ark_mem, is, t0, tf);
                if retval != ARK_SUCCESS {
                    return retval;
                }
                let (ycur, tempv2) = {
                    let m = ark_mem.borrow();
                    (
                        m.ycur.clone().expect("ycur set"),
                        m.tempv2.clone().expect("tempv2 set"),
                    )
                };
                retval = mriStep_StageERKFast(ark_mem, t0, tf, &ycur, &tempv2, need_inner_dsm);
                if retval != ARK_SUCCESS {
                    *nflagPtr = CONV_FAIL;
                }
            }
            MRISTAGE_ERK_NOFAST => {
                retval = mriStep_StageERKNoFast(ark_mem, is);
            }
            MRISTAGE_DIRK_NOFAST => {
                retval = mriStep_StageDIRKNoFast(ark_mem, is, nflagPtr);
            }
            MRISTAGE_DIRK_FAST => {
                retval = mriStep_StageDIRKFast(ark_mem, is, nflagPtr);
            }
            MRISTAGE_STIFF_ACC => {
                retval = ARK_SUCCESS;
            }
            _ => {}
        }
        if retval != ARK_SUCCESS {
            return retval;
        }

        /* apply user-supplied step postprocessing function (if supplied) */
        if ark_mem.borrow().PostProcessStepFn.is_some() && (stagetype != MRISTAGE_STIFF_ACC) {
            let (tcur, ycur) = {
                let m = ark_mem.borrow();
                (m.tcur, m.ycur.clone().expect("ycur set"))
            };
            retval = mriStep_CallPostProcessStepFn(ark_mem, tcur, &ycur);
            if retval != 0 {
                return ARK_POSTPROCESS_STEP_FAIL;
            }
        }

        /* conditionally reset the inner integrator with the modified stage solution */
        if stagetype != MRISTAGE_STIFF_ACC {
            let have_postprocess = ark_mem.borrow().PostProcessStepFn.is_some();
            if (stagetype != MRISTAGE_ERK_FAST) || have_postprocess {
                let ycur = ark_mem.borrow().ycur.clone().expect("ycur set");
                retval = mriStepInnerStepper_Reset(&stepper, tf, &ycur);
                if retval != ARK_SUCCESS {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_INNERSTEP_FAIL,
                        line!() as i32,
                        "mriStep_TakeStepMRIGARK",
                        file!(),
                        "Unable to reset the inner stepper",
                    );
                    return ARK_INNERSTEP_FAIL;
                }
            }
        }

        /* Compute temporal error estimate via difference between step
           solution and embedding, store in ark_mem->tempv1, and take norm. */
        if do_embedding {
            let (tempv4, ycur, tempv1, ewt) = {
                let m = ark_mem.borrow();
                (
                    m.tempv4.clone().expect("tempv4 set"),
                    m.ycur.clone().expect("ycur set"),
                    m.tempv1.clone().expect("tempv1 set"),
                    m.ewt.clone().expect("ewt set"),
                )
            };
            N_VLinearSum(ONE, &tempv4, -ONE, &ycur, &tempv1);
            *dsmPtr = N_VWrmsNorm(&tempv1, &ewt);
        }
    } /* loop over stages */

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_TakeStepMRISR:

  This routine performs a single MRISR step.

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
pub fn mriStep_TakeStepMRISR(
    ark_mem: &ARKodeMem,
    dsmPtr: &mut sunrealtype,
    nflagPtr: &mut i32,
) -> i32 {
    let mut stage: i32;
    let mut retval: i32; /* reusable return flag       */
    let mut embedding: sunbooleantype; /* flag indicating embedding  */
    let mut solution: sunbooleantype; /*   or solution stages       */
    let mut impl_corr: sunbooleantype; /* is slow correct. implicit? */
    let mut need_inner_dsm: sunbooleantype;
    let nested_mri: sunbooleantype;
    let mut nvec: i32;
    let tol: sunrealtype = 100.0 * SUN_UNIT_ROUNDOFF;

    /* access the MRIStep mem structure */
    retval = mriStep_step_mem_ok(ark_mem, "mriStep_TakeStepMRISR");
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
            m.tempv4.clone().expect("tempv4 set"),
            m.tempv2.clone().expect("tempv2 set"),
        )
    };

    /* initialize the current stage index */
    {
        let mut s = mriStep_mem_mut(ark_mem);
        s.istage = 0;
        s.cur_stage = 0;
    }

    let stepper = { mriStep_mem_mut(ark_mem).stepper.clone() }.expect("stepper set");

    /* if MRI adaptivity is enabled: reset fast accumulated error,
       and send appropriate control parameter to the fast integrator */
    let hcontroller = ark_mem
        .borrow()
        .hadapt_mem
        .as_ref()
        .expect("hadapt_mem set")
        .hcontroller
        .clone()
        .expect("hcontroller set");
    let adapt_type = SUNAdaptController_GetType(&hcontroller);
    need_inner_dsm = SUNFALSE;
    if adapt_type == SUN_ADAPTCONTROLLER_MRI_H_TOL {
        need_inner_dsm = SUNTRUE;
        mriStep_mem_mut(ark_mem).inner_dsm = ZERO;
        retval = mriStepInnerStepper_ResetAccumulatedError(&stepper);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_INNERSTEP_FAIL,
                line!() as i32,
                "mriStep_TakeStepMRISR",
                file!(),
                "Unable to reset the inner stepper error estimate",
            );
            return ARK_INNERSTEP_FAIL;
        }
        let inner_rtol_factor = { mriStep_mem_mut(ark_mem).inner_rtol_factor };
        let reltol = ark_mem.borrow().reltol;
        retval = mriStepInnerStepper_SetRTol(&stepper, inner_rtol_factor * reltol);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_INNERSTEP_FAIL,
                line!() as i32,
                "mriStep_TakeStepMRISR",
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
            (m.tcur, m.ycur.clone().expect("ycur set"))
        };
        retval = mriStepInnerStepper_Reset(&stepper, tcur, &ycur);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_INNERSTEP_FAIL,
                line!() as i32,
                "mriStep_TakeStepMRISR",
                file!(),
                "Unable to reset the inner stepper",
            );
            return ARK_INNERSTEP_FAIL;
        }
    }

    /* call nonlinear solver setup if it exists */
    let NLS = { mriStep_mem_mut(ark_mem).NLS.clone() };
    if let Some(NLS) = NLS {
        if NLS.ops.borrow().setup.is_some() {
            let tempv3 = ark_mem.borrow().tempv3.clone().expect("tempv3 set");
            N_VConst(ZERO, &tempv3); /* set guess to 0 */
            let mut nls_mem: Option<Box<dyn Any>> = Some(Box::new(ark_mem.clone()));
            retval = SUNNonlinSolSetup(&NLS, &tempv3, &mut nls_mem);
            if retval < 0 {
                return ARK_NLS_SETUP_FAIL;
            }
            if retval > 0 {
                return ARK_NLS_SETUP_RECVR;
            }
        }
    }

    /* Evaluate the slow RHS functions if needed. NOTE: we decide between calling the
       full RHS function (if ark_mem->fn is non-NULL and MRIStep is not an inner
       integrator) versus just updating the stored values of Fse[0] and Fsi[0].  In
       either case, we use ARK_FULLRHS_START mode because MRISR methods do not
       evaluate the RHS functions at the end of the time step (so nothing can be
       leveraged). */
    let (expforcing, impforcing) = {
        let s = mriStep_mem_mut(ark_mem);
        (s.expforcing, s.impforcing)
    };
    nested_mri = expforcing || impforcing;
    let fn_is_null = ark_mem.borrow().fn_.is_none();
    if fn_is_null || nested_mri {
        let (tcur, ycur) = {
            let m = ark_mem.borrow();
            (m.tcur, m.ycur.clone().expect("ycur set"))
        };
        retval = mriStep_UpdateF0(ark_mem, tcur, &ycur, ARK_FULLRHS_START);
        if retval != 0 {
            return ARK_RHSFUNC_FAIL;
        }

        /* For a nested MRI configuration we might still need fn to create a predictor
           but it should be fn only for the current nesting level which is why we use
           UpdateF0 in this case rather than FullRHS */
        let fn_v = ark_mem.borrow().fn_.clone();
        let (explicit_rhs, implicit_rhs) = {
            let s = mriStep_mem_mut(ark_mem);
            (s.explicit_rhs, s.implicit_rhs)
        };
        if fn_v.is_some() && nested_mri && implicit_rhs {
            let fn_v = fn_v.expect("fn set");
            if implicit_rhs && explicit_rhs {
                let (Fsi0, Fse0) = {
                    let s = mriStep_mem_mut(ark_mem);
                    (s.Fsi[0].clone(), s.Fse[0].clone())
                };
                N_VLinearSum(ONE, &Fsi0, ONE, &Fse0, &fn_v);
            } else {
                let Fsi0 = { mriStep_mem_mut(ark_mem).Fsi[0].clone() };
                N_VScale(ONE, &Fsi0, &fn_v);
            }
        }
    }
    if !fn_is_null && !ark_mem.borrow().fn_is_current {
        let (tcur, ycur, fn_v) = {
            let m = ark_mem.borrow();
            (
                m.tcur,
                m.ycur.clone().expect("ycur set"),
                m.fn_.clone().expect("fn set"),
            )
        };
        retval = mriStep_FullRHS(ark_mem, tcur, &ycur, &fn_v, ARK_FULLRHS_START);
        if retval != 0 {
            return ARK_RHSFUNC_FAIL;
        }
    }
    ark_mem.borrow_mut().fn_is_current = SUNTRUE;

    /* combine both RHS into FSE for ImEx problems, since MRISR fast forcing function
       only depends on Omega coefficients  */
    let (explicit_rhs, implicit_rhs) = {
        let s = mriStep_mem_mut(ark_mem);
        (s.explicit_rhs, s.implicit_rhs)
    };
    if implicit_rhs && explicit_rhs {
        let (Fse0, Fsi0) = {
            let s = mriStep_mem_mut(ark_mem);
            (s.Fse[0].clone(), s.Fsi[0].clone())
        };
        N_VLinearSum(ONE, &Fse0, ONE, &Fsi0, &Fse0);
    }

    /* The first stage is the previous time-step solution, so its RHS
       is the [already-computed] slow RHS from the start of the step */

    let stages = { mriStep_mem_mut(ark_mem).stages };
    let MRIC = { mriStep_mem_mut(ark_mem).MRIC.clone() }.expect("MRIC set");
    let accum_type = ark_mem.borrow().AccumErrorType;

    /* Determine how many stages will be needed */
    let max_stages: i32 = if fixedstep && (accum_type == ARK_ACCUMERROR_NONE) {
        stages
    } else {
        stages + 1
    };

    /* Loop over stages */
    stage = 1;
    while stage < max_stages {
        /* Determine if this is an "embedding" or "solution" stage */
        solution = stage == stages - 1;
        embedding = stage == stages;

        /* Set initial condition for this stage (all but first stage) */
        if stage > 1 {
            let (yn, ycur) = {
                let m = ark_mem.borrow();
                (
                    m.yn.clone().expect("yn set"),
                    m.ycur.clone().expect("ycur set"),
                )
            };
            N_VScale(ONE, &yn, &ycur);
        }

        /* Set current stage abscissa */
        let cstage: sunrealtype = if embedding {
            ONE
        } else {
            MRIC.borrow().c[stage as usize]
        };

        /* Set current stage time and index */
        let (tn, h) = {
            let m = ark_mem.borrow();
            (m.tn, m.h)
        };
        let tcur = tn + cstage * h;
        ark_mem.borrow_mut().tcur = tcur;
        {
            let mut s = mriStep_mem_mut(ark_mem);
            s.istage = stage;
            s.cur_stage = stage;
        }

        /* Compute forcing function for inner solver */
        retval = mriStep_ComputeInnerForcing(ark_mem, stage, tn, tcur);
        if retval != ARK_SUCCESS {
            return retval;
        }

        /* Reset the inner stepper on all but the first stage due to
           "stage-restart" structure */
        if stage > 1 {
            let ycur = ark_mem.borrow().ycur.clone().expect("ycur set");
            retval = mriStepInnerStepper_Reset(&stepper, tn, &ycur);
            if retval != ARK_SUCCESS {
                arkProcessError(
                    Some(ark_mem),
                    ARK_INNERSTEP_FAIL,
                    line!() as i32,
                    "mriStep_TakeStepMRISR",
                    file!(),
                    "Unable to reset the inner stepper",
                );
                return ARK_INNERSTEP_FAIL;
            }
        }

        /* Evolve fast IVP for this stage, potentially get inner dsm on
           all non-embedding stages */
        {
            let ycur = ark_mem.borrow().ycur.clone().expect("ycur set");
            retval = mriStep_StageERKFast(
                ark_mem,
                tn,
                tcur,
                &ycur,
                &ytemp,
                need_inner_dsm && !embedding,
            );
        }
        if retval != ARK_SUCCESS {
            *nflagPtr = CONV_FAIL;
            return retval;
        }

        /* perform MRISR slow/implicit correction */
        impl_corr = SUNFALSE;
        if implicit_rhs {
            /* determine whether implicit RHS correction will require an implicit solve */
            let g_ss = { MRIC.borrow().G[0][stage as usize][stage as usize] };
            impl_corr = SUNRabs(g_ss) > tol;

            /* perform implicit solve for correction */
            if impl_corr {
                /* update stage index for prediction and nonlinear solver if this is an "embedded" stage */
                if embedding {
                    mriStep_mem_mut(ark_mem).istage = stage - 1;
                }

                /* Call predictor for current stage solution (result placed in zpred) */
                let (istage, zpred) = {
                    let s = mriStep_mem_mut(ark_mem);
                    (s.istage, s.zpred.clone().expect("zpred set"))
                };
                retval = mriStep_Predict(ark_mem, istage, &zpred);
                if retval != ARK_SUCCESS {
                    return retval;
                }

                /* If a user-supplied predictor routine is provided, call that here
                   Note that mriStep_Predict is *still* called, so this user-supplied
                   routine can just "clean up" the built-in prediction, if desired. */
                let have_stage_predict = { mriStep_mem_mut(ark_mem).stage_predict.is_some() };
                if have_stage_predict {
                    retval = mriStep_CallStagePredict(ark_mem, tcur, &zpred);
                    if retval < 0 {
                        return ARK_USER_PREDICT_FAIL;
                    }
                    if retval > 0 {
                        return TRY_AGAIN;
                    }
                }

                /* fill sdata with explicit contributions to correction */
                let ycur = ark_mem.borrow().ycur.clone().expect("ycur set");
                let (cvals, Xvecs, sdata) = {
                    let mut s = mriStep_mem_mut(ark_mem);
                    s.cvals[0] = ONE;
                    s.Xvecs[0] = Some(ycur);
                    s.cvals[1] = -ONE;
                    let v = s.zpred.clone().expect("zpred set");
                    s.Xvecs[1] = Some(v);
                    for j in 0..stage {
                        let g = MRIC.borrow().G[0][stage as usize][j as usize];
                        s.cvals[(j + 2) as usize] = h * g;
                        let v = s.Fsi[j as usize].clone();
                        s.Xvecs[(j + 2) as usize] = Some(v);
                    }
                    let sdata = s.sdata.clone().expect("sdata set");
                    (s.cvals.clone(), mriStep_xvecs(&s, stage + 2), sdata)
                };
                retval = N_VLinearCombination(stage + 2, &cvals, &Xvecs, &sdata);
                if retval != 0 {
                    return ARK_VECTOROP_ERR;
                }

                /* Update gamma for implicit solver */
                let firststage = ark_mem.borrow().firststage;
                {
                    let mut s = mriStep_mem_mut(ark_mem);
                    s.gamma = h * g_ss;
                    if firststage {
                        s.gammap = s.gamma;
                    }
                    s.gamrat = if firststage { ONE } else { s.gamma / s.gammap };
                }

                /* perform implicit solve (result is stored in ark_mem->ycur); return
                   with positive value on anything but success */
                *nflagPtr = mriStep_Nls(ark_mem, *nflagPtr);
                if *nflagPtr != ARK_SUCCESS {
                    return TRY_AGAIN;
                }
            }
            /* perform explicit update for correction */
            else {
                let ycur = ark_mem.borrow().ycur.clone().expect("ycur set");
                let (cvals, Xvecs) = {
                    let mut s = mriStep_mem_mut(ark_mem);
                    s.cvals[0] = ONE;
                    s.Xvecs[0] = Some(ycur.clone());
                    for j in 0..stage {
                        let g = MRIC.borrow().G[0][stage as usize][j as usize];
                        s.cvals[(j + 1) as usize] = h * g;
                        let v = s.Fsi[j as usize].clone();
                        s.Xvecs[(j + 1) as usize] = Some(v);
                    }
                    (s.cvals.clone(), mriStep_xvecs(&s, stage + 1))
                };
                retval = N_VLinearCombination(stage + 1, &cvals, &Xvecs, &ycur);
                if retval != 0 {
                    return ARK_VECTOROP_ERR;
                }
            }
        }

        /* apply user-supplied stage or step postprocessing function (if supplied),
           and reset the inner integrator with the modified stage solution */
        let (have_post_stage, have_post_step) = {
            let m = ark_mem.borrow();
            (m.PostProcessStageFn.is_some(), m.PostProcessStepFn.is_some())
        };
        if !solution && !embedding && have_post_stage {
            let (tcur_now, ycur) = {
                let m = ark_mem.borrow();
                (m.tcur, m.ycur.clone().expect("ycur set"))
            };
            retval = mriStep_CallPostProcessStageFn(ark_mem, tcur_now, &ycur);
            if retval != 0 {
                return ARK_POSTPROCESS_STAGE_FAIL;
            }
            let (tcur_now, ycur) = {
                let m = ark_mem.borrow();
                (m.tcur, m.ycur.clone().expect("ycur set"))
            };
            retval = mriStepInnerStepper_Reset(&stepper, tcur_now, &ycur);
            if retval != ARK_SUCCESS {
                arkProcessError(
                    Some(ark_mem),
                    ARK_INNERSTEP_FAIL,
                    line!() as i32,
                    "mriStep_TakeStepMRISR",
                    file!(),
                    "Unable to reset the inner stepper",
                );
                return ARK_INNERSTEP_FAIL;
            }
        } else if solution && have_post_step {
            let (tcur_now, ycur) = {
                let m = ark_mem.borrow();
                (m.tcur, m.ycur.clone().expect("ycur set"))
            };
            retval = mriStep_CallPostProcessStepFn(ark_mem, tcur_now, &ycur);
            if retval != 0 {
                return ARK_POSTPROCESS_STEP_FAIL;
            }
            let (tcur_now, ycur) = {
                let m = ark_mem.borrow();
                (m.tcur, m.ycur.clone().expect("ycur set"))
            };
            retval = mriStepInnerStepper_Reset(&stepper, tcur_now, &ycur);
            if retval != ARK_SUCCESS {
                arkProcessError(
                    Some(ark_mem),
                    ARK_INNERSTEP_FAIL,
                    line!() as i32,
                    "mriStep_TakeStepMRISR",
                    file!(),
                    "Unable to reset the inner stepper",
                );
                return ARK_INNERSTEP_FAIL;
            }
        }

        /* Compute updated slow RHS (except for final solution or embedding) */
        if !solution && !embedding {
            let deduce_rhs = { mriStep_mem_mut(ark_mem).deduce_rhs };

            /* call the user-supplied pre-RHS function (if supplied) */
            if ark_mem.borrow().PreRhsFn.is_some() {
                if explicit_rhs || (implicit_rhs && (!deduce_rhs || !impl_corr)) {
                    let (tcur_now, ycur) = {
                        let m = ark_mem.borrow();
                        (m.tcur, m.ycur.clone().expect("ycur set"))
                    };
                    retval = mriStep_CallPreRhsFn(ark_mem, tcur_now, &ycur);
                    if retval != 0 {
                        return ARK_PRERHSFN_FAIL;
                    }
                }
            }

            /* store implicit slow rhs */
            if implicit_rhs {
                if !deduce_rhs || !impl_corr {
                    let (tcur_now, ycur) = {
                        let m = ark_mem.borrow();
                        (m.tcur, m.ycur.clone().expect("ycur set"))
                    };
                    let Fsi_stage = { mriStep_mem_mut(ark_mem).Fsi[stage as usize].clone() };
                    retval = mriStep_CallFsi(ark_mem, tcur_now, &ycur, &Fsi_stage);
                    mriStep_mem_mut(ark_mem).nfsi += 1;

                    if retval < 0 {
                        return ARK_RHSFUNC_FAIL;
                    }
                    if retval > 0 {
                        return ARK_UNREC_RHSFUNC_ERR;
                    }

                    /* Add external forcing to Fsi[stage], if applicable */
                    if impforcing {
                        let (cvals, Xvecs) = {
                            let mut s = mriStep_mem_mut(ark_mem);
                            s.cvals[0] = ONE;
                            let v = s.Fsi[stage as usize].clone();
                            s.Xvecs[0] = Some(v);
                            nvec = 1;
                            mriStep_ApplyForcing(&mut s, tcur_now, ONE, &mut nvec);
                            (s.cvals.clone(), mriStep_xvecs(&s, nvec))
                        };
                        let _ = N_VLinearCombination(
                            Xvecs.len() as i32,
                            &cvals,
                            &Xvecs,
                            &Fsi_stage,
                        );
                    }
                } else {
                    let (gamma, zcor, sdata, Fsi_stage) = {
                        let s = mriStep_mem_mut(ark_mem);
                        (
                            s.gamma,
                            s.zcor.clone().expect("zcor set"),
                            s.sdata.clone().expect("sdata set"),
                            s.Fsi[stage as usize].clone(),
                        )
                    };
                    N_VLinearSum(ONE / gamma, &zcor, -ONE / gamma, &sdata, &Fsi_stage);
                }
            }

            /* store explicit slow rhs */
            if explicit_rhs {
                let (tcur_now, ycur) = {
                    let m = ark_mem.borrow();
                    (m.tcur, m.ycur.clone().expect("ycur set"))
                };
                let Fse_stage = { mriStep_mem_mut(ark_mem).Fse[stage as usize].clone() };
                retval = mriStep_CallFse(ark_mem, tcur_now, &ycur, &Fse_stage);
                mriStep_mem_mut(ark_mem).nfse += 1;

                if retval < 0 {
                    return ARK_RHSFUNC_FAIL;
                }
                if retval > 0 {
                    return ARK_UNREC_RHSFUNC_ERR;
                }

                /* Add external forcing to Fse[stage], if applicable */
                if expforcing {
                    let (cvals, Xvecs) = {
                        let mut s = mriStep_mem_mut(ark_mem);
                        s.cvals[0] = ONE;
                        let v = s.Fse[stage as usize].clone();
                        s.Xvecs[0] = Some(v);
                        nvec = 1;
                        mriStep_ApplyForcing(&mut s, tcur_now, ONE, &mut nvec);
                        (s.cvals.clone(), mriStep_xvecs(&s, nvec))
                    };
                    let _ =
                        N_VLinearCombination(Xvecs.len() as i32, &cvals, &Xvecs, &Fse_stage);
                }
            }

            /* combine both RHS into Fse for ImEx problems since
               fast forcing function only depends on Omega coefficients */
            if implicit_rhs && explicit_rhs {
                let (Fse_stage, Fsi_stage) = {
                    let s = mriStep_mem_mut(ark_mem);
                    (s.Fse[stage as usize].clone(), s.Fsi[stage as usize].clone())
                };
                N_VLinearSum(ONE, &Fse_stage, ONE, &Fsi_stage, &Fse_stage);
            }
        }

        /* If this is the solution stage, archive for error estimation */
        if solution {
            let ycur = ark_mem.borrow().ycur.clone().expect("ycur set");
            N_VScale(ONE, &ycur, &ytilde);
        }

        stage += 1;
    } /* loop over stages */

    /* if temporal error estimation is enabled: compute estimate via difference between
       step solution and embedding, store in ark_mem->tempv1, store norm in dsmPtr, and
       copy solution back to ycur */
    if !fixedstep || (accum_type != ARK_ACCUMERROR_NONE) {
        let (ycur, tempv1, ewt) = {
            let m = ark_mem.borrow();
            (
                m.ycur.clone().expect("ycur set"),
                m.tempv1.clone().expect("tempv1 set"),
                m.ewt.clone().expect("ewt set"),
            )
        };
        N_VLinearSum(ONE, &ytilde, -ONE, &ycur, &tempv1);
        *dsmPtr = N_VWrmsNorm(&tempv1, &ewt);
        N_VScale(ONE, &ytilde, &ycur);
    }

    ARK_SUCCESS
}
