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

- `sunrealtype` = `f64`, `sunindextype` = `i64`, `sunbooleantype` = `bool`.
- C "object with ops table" structs (N_Vector, SUNMatrix, SUNLinearSolver,
  SUNNonlinearSolver, SUNAdaptController, …) become Rust structs holding
  their content directly; polymorphism over implementations is expressed
  the same way the C code does — an ops table of plain `fn` pointers where
  the C API exposes one, so user-supplied overrides keep working.
- `SUNContext` owns the error handler stack, logger, and profiler exactly
  as in C.
- `user_data`: `Option<Box<dyn Any>>`, passed to callbacks as
  `Option<&mut dyn Any>` (C: `void*`).
- User callbacks are plain `fn` pointer types matching the C signature
  argument-for-argument (same names in the same order in Rust signatures).
  Do not change a callback signature without updating every example.

## Aliasing / copy-back rule

Where C aliases internal state with user buffers (CVODE `cv_y`/`yout`,
IDA `ida_yy`/`yret` etc., ARKODE `ycur`/`yout`), the Rust port copies the
internal buffer to the user buffer at every return path — success, early
error, and root-return alike.

Vector ops: in-place methods handle aliased operands
(`linear_sum_with`-style); free functions (`N_VLinearSum`) mirror the C
call shape for distinct operands. Which one an internal call site uses is
decided by what the C source actually aliases at that site.

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
