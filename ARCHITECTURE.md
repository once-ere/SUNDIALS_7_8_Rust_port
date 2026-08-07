# ARCHITECTURE — cross-module contracts

Pinned decisions. Read before modifying shared types; extend (append) when
a new contract is fixed, never silently change an existing one.

## Crate graph

`sundials_core` ← {`cvode_rs`, `cvodes_rs`, `kinsol_rs`, `ida_rs`,
`idas_rs`, `arkode_rs`}. Solver crates never depend on each other.

## Module naming

One Rust module per upstream C file, named after its base name. Impl
headers (`*_impl.h`) and public `include/<pkg>/*.h` content fold into the
matching module. `.def` X-macro tables become `const` table data in the
including module.

## Core type model (fixed in Phase 1)

- `sunrealtype` = `f64`, `sunindextype` = `i64`, `suncountertype` = `i64`,
  `sunbooleantype` = `bool`, `SUNErrCode` = `i32`, `SUNComm` = `i32`.
- **Handle model.** Every C heap object reached through a pointer
  (`N_Vector`, `SUNMatrix`, `SUNLinearSolver`, `SUNNonlinearSolver`,
  `SUNAdaptController`, `SUNContext`, solver mems, …) becomes
  `pub type X = Rc<X_>` where `X_` holds `content: RefCell<Box<dyn Any>>`
  (C `void* content`) plus an ops struct of plain `fn` pointers (where the
  C API has one) and the `sunctx` handle. Cloning the `Rc` is the C
  pointer copy; `Rc::ptr_eq` is C pointer equality.
- Ops are plain `fn` pointers taking `&Handle` arguments — identical call
  shape to C — and mutate through the `RefCell`. User-supplied override
  implementations keep working exactly as in C.
- `SUNContext` owns the error handler stack, logger, and profiler exactly
  as in C.
- `user_data`: `Option<Box<dyn Any>>`, passed to callbacks as
  `Option<&mut dyn Any>` (C: `void*`). Solver internals `Option::take`
  the box out of the mem record around each callback invocation, so the
  callback gets exclusive access without re-borrowing the mem.
- Solver mems: the public handle is `Rc<RefCell<CVodeMemRec>>`-style;
  public API functions borrow once at entry and pass `&mut CVodeMemRec`
  internally — matching C's `cv_mem->` style with zero borrow churn.
- User callbacks are plain `fn` pointer types matching the C signature
  argument-for-argument (same names in the same order in Rust signatures).
  Do not change a callback signature without updating every example.

## Aliasing / copy-back rule

Where C aliases internal state with user buffers (CVODE `cv_y`/`yout`,
IDA `ida_yy`/`yret` etc., ARKODE `ycur`/`yout`), the Rust port copies the
internal buffer to the user buffer at every return path — success, early
error, and root-return alike.

Vector ops: free functions (`N_VLinearSum(a, &x, b, &y, &z)`) mirror the
C call shape for **all** call sites and are alias-safe by construction —
implementations detect operand aliasing (`Rc::ptr_eq`) and take a single
mutable borrow for the aliased case. This satisfies the in-place-method
contract trivially: the free function *is* safe under aliasing, so C call
sites translate 1:1 whether or not they alias. In-place convenience
methods may exist additionally but are never required for safety.

## Formatting

`sundials_core::sundials_utils::{fmt_e, fmt_f, fmt_g}` implement C
`printf` `%e/%f/%g` exactly (default precision 6; `%g` strips trailing
zeros; exponent `e±dd`, at least two digits; `inf`/`nan` lowercase).
Width variants `fmt_ew/fmt_fw/fmt_gw(x, width, prec)` right-justify with
spaces (C `%W.Pe`). Never use Rust `{:e}`.

## OS mapping

`sundials_profiler.c` timers → `std::time::Instant`;
`sundials_futils.c` → `std::fs`; env access in logger/CLI → `std::env`.
Public API and observable behavior identical to C.

## Error/return conventions

Exact C flag names and values (`CV_SUCCESS`, `CV_MEM_NULL`, `IDA_SUCCESS`,
`KIN_SUCCESS`, `ARK_SUCCESS`, `SUN_SUCCESS`, …). Negative = fatal,
positive = recoverable. Functions that return flags in C return the same
integer type in Rust; output parameters in C (`T *out`) become `&mut T`
in the same argument position with the same name.

## Established porting patterns (locked during Phase 1)

- **Content downcast**: every implementation module defines
  `fn content_mut(X) -> RefMut<'_, ContentStruct>` via
  `RefMut::map(x.content.borrow_mut(), downcast_mut)`. Public accessor
  macros (`NV_DATA_S`, `SM_DATA_D`, …) are functions returning `RefMut`
  guards; drop the guard before any other op on the same object.
- **Granular borrow rule**: never hold a RefCell borrow (of a mem, a
  solver content, or vector data) across a call that can re-enter it —
  callbacks (RHS, ATimes, Psolve, Jacobian) reach integrator state
  through their own handle. Iterative-solver `solve` ops move ALL content
  state into locals at entry (`Option::take` for `Box<dyn Any>` callback
  data, `mem::take` for arrays, `Rc` clones for vectors), run the C
  algorithm inside a closure returning the flag, and restore + write back
  (numiters/resnorm/zeroguess/last_flag) at one exit point. Final flag
  values are identical to C's multi-return-path writes because
  logging-level-2 builds have no observable effects in between.
- **`SUNCheck*`/`SUNAssert`**: release no-ops; call sites evaluate the
  call and continue (`let _ = f(...)` where C had `(void)`).
  `SUNLogInfo`/`SUNLogDebug`/`SUNLogExtraDebug*` compile away entirely at
  logging level 2 and are omitted at translation time.
- **CLI parsing**: `argv: &[String]` with `argv[0]` = program name; C
  `atoi` maps to `s.trim().parse().unwrap_or(0)`; prefix matching is
  literal (`<id>.` with no leading dashes).
- **C output params**: `T *out` → `&mut T` same position/name;
  functions returning object pointers return `Option<Handle>`
  (NULL = `None`). Constructors that C would fail with NULL return
  `None`.
- **fmt helpers**: `fmt_e/f/g(x, prec)`, width variants
  `fmt_ew/fmt_fw/fmt_gw(x, width, prec)`, and `sun_format_e/g/sg(x)` for
  the `SUN_FORMAT_E/G/SG` macros ("% .15e" / "%.15g" / "%+.15g").
- **Vector arrays**: C `N_Vector*` → `&[N_Vector]` (handles are Rc
  clones); `N_Vector**` → `&[Vec<N_Vector>]`. Row-wise Hessenberg
  `sunrealtype**` → `&mut [Vec<f64>]`; column-pointer arrays
  (`SUNDlsMat cols`) → `dls_cols()` chunks_mut views.

## Accepted deviation classes (adversarially verified, Phase 1)

These divergences from the release-C reference build are deliberate,
verified unobservable on any path a valid serial example takes, and must
be applied CONSISTENTLY in later phases:

1. **Kept failure-path checks.** Where release C compiles out
   `SUNAssert`/`SUNCheckCall` (silently proceeding on misuse), ported
   modules may keep the check as a plain `if` returning the error code
   (fn-pointer-presence checks, `file_name` guards, propagated sub-call
   errors in Set*/Initialize forwarding). Only observable on invalid
   usage or malformed CLI input.
2. **Ownership snapshots.** `SUNMemoryHelper_Wrap`/`_Alias` take/clone
   owned `Vec<u8>` buffers instead of aliasing raw pointers; no
   write-through. Any future consumer ported from C that mutates a
   wrapped buffer afterward must mutate through the SUNMemory handle.
3. **C-locale ASCII whitespace.** `atoi`/`atol`/`SUNStrToReal` skip only
   ASCII whitespace (matching C-locale `isspace`/`strtod`), implemented
   via `trim_start_matches` — never Unicode `trim`.
4. **Unsigned wrap.** C `size_t` counter arithmetic that can underflow
   maps to `wrapping_sub` (never a panicking `-=`).
5. **C UB → deterministic panic.** NULL deref / out-of-bounds /
   double-free in C map to Rust panics at the same site.
6. **`user_data` pointer-snapshot sites (Phase 2).** C code that
   snapshots the raw `user_data` pointer and reuses it later
   (CVODE `cv_e_data = cv_user_data` in `cvInitialSetup`; CVLS
   `P_data = cv_user_data` at `CVodeSetLinearSolver`) cannot alias a
   `Box`: the port passes the CURRENT `cv_user_data` box at call time
   instead. Divergent only when `CVodeSetUserData` is called
   mid-integration after the snapshot point — no reference example
   does this, and the Rust behavior matches the documented SUNDIALS
   semantics. Same class: `void*`-returning getters
   (`CVodeGetUserData`, `CVodeGetNonlinearSystemData`) SWAP the box
   with the caller's out-param; the caller must hand it back (via
   `CVodeSetUserData` or a second swap) before the integrator next
   invokes a user callback.
7. **Hoisted callback fn-pointers within one evaluation.** DQ loops
   (`cvLsDQJtimes` retries, dense/band DQJac column loops) copy the
   RHS/jt fn pointer to a local before the loop where C re-reads the
   field each iteration; a callback re-entrantly swapping the fn
   mid-evaluation would take effect one call later than in C. This is
   the locked move-state-into-locals pattern; observable by no valid
   example.
