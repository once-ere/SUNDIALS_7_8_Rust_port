# sundials_7_8_0__rs — workspace rules

Pure-Rust port of SUNDIALS 7.8.0. The upstream C tree is the parent
directory (`../src/`, `../include/`, `../examples/`) and is **read-only**.
This workspace is its own git repo; git is the undo mechanism.

## Hard rules

1. **Fidelity first.** Line-by-line faithful translation: control flow,
   constants, tolerances, heuristics, error/return codes, and argument
   lists (names, order, meaning) match the parent C function exactly.
   Preserve arithmetic order — acceptance is byte-identical printed output.
2. Zero `unsafe`, zero FFI, zero external crates (std only), zero warnings
   from `cargo build --workspace`.
3. Never stub a missing symbol — its definition is under `../src/` or
   `../include/`; port it into `sundials_core`.
4. Public API keeps exact C names and return-flag conventions
   (`CV_SUCCESS = 0`; negative = fatal, positive = recoverable). Crate
   roots carry `#![allow(non_snake_case, non_camel_case_types,
   non_upper_case_globals)]`.
5. All float output goes through
   `sundials_core::sundials_utils::{fmt_e, fmt_f, fmt_g}` — never `{:e}`.
6. C buffer aliasing (e.g. CVODE `cv_y` / user `yout`): copy back at
   **every** return path, including early-error and rootfinding exits.
   All of CVODE(S), IDA(S), ARKODE do this.
7. Once a crate's examples verify green they stay green — the cumulative
   regression gate runs `tools/verify_examples.sh` for all crates ported
   so far at every phase gate.

## Module layout

- Module = C file base name + `.rs` (`cvodes_nls_stg1.c` →
  `crates/cvodes_rs/src/cvodes_nls_stg1.rs`; `arkode_impl.h` →
  `arkode_impl.rs`). Public `include/` headers fold into the matching
  module.
- Solver crates re-export every shared `sundials_core` module at root and
  provide a flat prelude so examples can `use cvode_rs::*;`.
- One `[[example]]` entry per translated example; example name = C base
  name.
- `user_data` is `Option<Box<dyn Any>>`; callbacks are plain `fn`
  pointers. Aliasing vector ops get in-place methods; free functions
  (`N_VLinearSum`) serve distinct operands.

## Workflow

- Commit after every ported file (or small coherent group); tag phase
  gates (`phase2-cvode-green`, …).
- After EVERY build/test/run: `… 2>&1 | tee <log>` then **Read the log**
  before the next edit. Never re-run a command that produced no output.
- Max two attempts per failing command, then switch strategy.
- Read each in-scope C file exactly once, at translation time, completely.
  Never read excluded paths (GPU/MPI/KLU/LAPACK/Fortran/xbraid trees).
- Update `PROGRESS.md` (per-file status: todo | ported | building |
  committed) and `VERIFICATION.md` (per-variant status) as units land.
- Resume after context loss from this file + `PROGRESS.md` + `git log` —
  do not re-explore the tree.

## Verification

`tools/verify_examples.sh [crate|all|list]` parses the upstream
CMakeLists tuples (199 variants), builds release examples, runs each
variant with exact argv, diffs against `../examples/...` references
(noise-filtered symmetrically), and writes `logs/summary.txt`. Read only
the summary; open individual diffs only for non-IDENTICAL lines.
Documented exceptions: `*L` examples (LAPACK→native) and iterative-solver
last-digit drift — one-line justification each in `VERIFICATION.md`.
CLI-option variants use bare `<solverid>.<key>` tokens (no leading
dashes); the parser prefix-matches literally.
