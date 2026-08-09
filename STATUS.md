# STATUS — where the port stands

Snapshot of the finished state. For the public guide read `sundials.md`; for
per-variant evidence `VERIFICATION.md`; for per-file status `PROGRESS.md`.
Resume after an interruption from this file + `PROGRESS.md` + `git log`.

## Library ports — complete

All seven crates are ported and build warning-free.

| crate | modules | library | example gate |
|---|---:|---|---|
| `sundials_core` | 51 | committed | n/a (23 lib unit tests) |
| `cvode_rs` | 12 | committed | swept — 12/24 identical, 9 documented, 3 excluded |
| `cvodes_rs` | 16 | committed | swept — 23/39 identical, 10 documented, 6 excluded |
| `kinsol_rs` | 8 | committed | swept — 19/22 identical, 1 documented, 2 excluded |
| `ida_rs` | 8 | committed | swept — 10/14 identical, 1 documented, 3 excluded |
| `idas_rs` | 12 | committed | swept — 12/22 identical, 4 documented, 6 excluded |
| `arkode_rs` | 34 | committed | **not swept** — 78 variants outstanding |

Verified first-hand on 2026-08-09:

* `cargo clean -p` for all seven crates, then `cargo build --workspace` →
  all seven compile, **zero warnings**.
* `cargo test --workspace --lib` → **25 passed, 0 failed**
  (23 `sundials_core`, 1 `cvodes_rs`, 1 `idas_rs`).
* `tools/verify_examples.sh all` → **76 IDENTICAL / 25 documented exceptions /
  20 excluded** across the 121 variants of the five swept crates.
  Reproduces the statuses recorded in `VERIFICATION.md` exactly; no
  regressions.
* `grep unsafe` over `crates/*/src` and `crates/*/examples` → 0.
  `Cargo.lock` holds exactly the 7 workspace packages — no external crates.

## What remains

1. **ARKODE example gate.** The 78 ARKODE variants have never been run against
   their references. Outstanding work, in order:
   a. finish translating the remaining `examples/arkode/C_serial` programs;
   b. add an `[[example]]` entry per program to `crates/arkode_rs/Cargo.toml`
      (without it `tools/verify_examples.sh` skips the crate);
   c. fix the examples that do not compile — at last check
      `ark_KrylovDemo_prec.rs` fails with `E0384` (immutable `nrmfac`
      reassigned, ~line 276);
   d. run `tools/verify_examples.sh arkode_rs`, then `all` as the cumulative
      regression gate;
   e. record every non-IDENTICAL variant in `VERIFICATION.md` with evidence,
      using the established bar: byte-identity against a locally built
      pristine upstream C binary, never by tuning the example.
   Example translation is in flight in a parallel workstream; re-check
   `crates/arkode_rs/examples/` and `git log` before starting.
2. **ERKStep discrete adjoint** (ARCHITECTURE §12). `ERKStepCreateAdjointStepper`,
   `erkStep_TakeStep_Adjoint`, `erkStep_fe_Adj` and `erkStep_SUNStepperReInit`
   are untranslated, so `erkStep_Init` installs `erkStep_TakeStep`
   unconditionally. The original blocker (`nvector_manyvector`) is gone and the
   ARKStep adjoint cluster is ported, so this is translation work only. Adding
   it must restore the `do_adjoint` branch in `erkStep_Init` in the same change.
3. **Own `LICENSE` file.** The workspace is a BSD-3-Clause derivative of
   SUNDIALS but carries no licence file of its own — it inherits the one at the
   root of the containing upstream tree. Add a copy (plus the MIT notice for
   the ARM/musl `pow` port in `sundials_math.rs`) before the crates are
   distributed independently.
4. **`PROGRESS.md` bookkeeping.** Two stale entries, both cosmetic:
   the four `arkode_*.def` files are marked `todo` although their tables are
   folded into the corresponding modules per the ARCHITECTURE rule (verified:
   34/34 ERK, 32/32 DIRK, 29/29 MRI, 8/8 splitting coefficient sets present);
   and the 18 `cvode_rs` example lines are marked `todo` although that crate
   passed its gate at `phase2-cvode-green`.
5. **Unexercised code.** Every `*_bbdpre` module is compile-only — BBD
   preconditioning is MPI-only upstream, so no serial reference example can
   regression-test it. Same for the excluded KLU/SuperLU paths.
6. **kinsol duplicate QR-add.** `kinsol_orth::kinQRAdd` and an inlined copy in
   `kinsol.rs` implement the same ICWY update. Verified semantically
   equivalent; worth collapsing, but only once `kinAnalytic_fp --m_aa 2
   --orth_aa 1` can prove it byte-for-byte (it currently verifies IDENTICAL,
   so the proof is available).

## Acceptance rules that still bind

* Byte-identical stdout against the upstream `.out`, noise-filtered
  symmetrically, is the only pass condition (`tools/verify_examples.sh`).
* When a shipped `.out` cannot be reproduced, the fallback bar is byte-identity
  against a locally built pristine upstream C binary — CMake Release, clang,
  `-O3 -DNDEBUG -ffp-contract=off`, logging level 2, error checks off,
  profiling off, serial. Three recurring reference-side classes: `ref-libm`,
  stale upstream reference, LAPACK→native.
* Zero `unsafe`, zero FFI, zero external crates, zero build warnings.
* Once a crate's examples verify green they stay green; the cumulative gate is
  `tools/verify_examples.sh all` at every phase gate.

## Remaining sequence

arkode examples → `tools/verify_examples.sh arkode_rs` → cumulative
`tools/verify_examples.sh all` → tag `v1-complete` → push to `origin`
(https://github.com/once-ere/SUNDIALS_7_8_Rust_port).
