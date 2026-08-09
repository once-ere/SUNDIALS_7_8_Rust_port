# STATUS — complete

The port is finished. All eight phases are done, every crate is verified
against the upstream reference outputs, and the cumulative gate passes.
For the public guide read `sundials.md`; for per-variant evidence
`VERIFICATION.md`; for per-file status `PROGRESS.md`.

## Final state

| crate | modules | examples | byte-identical | documented | excluded |
|---|---:|---:|---:|---:|---:|
| `sundials_core` | 51 | — | — | — | — |
| `cvode_rs` | 12 | 18 | 12 | 9 | 3 |
| `cvodes_rs` | 16 | 26 | 23 | 10 | 6 |
| `kinsol_rs` | 8 | 9 | 19 | 1 | 2 |
| `ida_rs` | 8 | 8 | 10 | 1 | 3 |
| `idas_rs` | 12 | 13 | 12 | 4 | 6 |
| `arkode_rs` | 34 | 34 | 51 | 27 | 0 |
| **total** | **141** | **108** | **127** | **52** | **20** |

Verified first-hand on 2026-08-09, on the final tree:

* `cargo build --workspace` → all seven crates, **zero warnings**.
* `cargo test --workspace --lib` → **25 passed, 0 failed**.
* `tools/verify_examples.sh all` → **127 IDENTICAL / 52 documented /
  20 excluded** over all **199** variants. Zero FAIL, zero NO-REF, zero
  NO-BINARY, zero regressions.
* Zero `unsafe`, zero FFI, zero external crates (`Cargo.lock` holds exactly
  the seven workspace packages).

## What the 52 documented divergences are

None is a port defect. Each is a case where the shipped `.out` cannot be
reproduced by its own C source on this machine, established by building the
pristine upstream C locally (CMake Release, clang, `-O3 -DNDEBUG
-ffp-contract=off`, logging level 2, error checks off, profiling off,
monitoring on, serial) and comparing against that. Three classes:

* **`ref-libm`** — the reference embeds the generating host's glibc
  `sin`/`exp`/`pow` rounding inside the integration feedback loop. A one-ulp
  difference forks the step-size trajectory. Several `.out` files in one
  family require *mutually incompatible* libm versions, so no single build
  can match them all.
* **stale upstream reference** — the `.out` predates its own source: the
  `SUN_TABLE_WIDTH` 28→29 change (whitespace-only diffs across a statistics
  block), a changed format string, or regeneration on a new CI host.
* **LAPACK→native** — the `*L` examples run the native dense/band solvers, so
  factorisation arithmetic differs in the last digit.

## Notable defects the byte-identity gate caught

The gate paid for itself; these were all found by output comparison, not by
review:

1. **Deterministic `pow`, unfused FMA** (`sundials_math.rs`) — the final
   `scale + scale * tmp` must be a single fused multiply-add, as in the glibc
   build the references came from. Unfused it rounds the wrong way near a
   midpoint, forking `ark_robertson` into a 228-line diff. Fusing it turned
   eight ARKODE variants identical.
2. **Platform `pow` at all** — the original `f64::powf` differs from glibc by
   one ulp on rare inputs inside the step-size heuristics; replacing it with
   the ported ARM/musl algorithm fixed three examples at once.
3. **Newton solver retry loop** — an initial-residual failure fell into the
   jbad-retry block instead of breaking out, spinning forever on a recoverable
   RHS flag (`cvRoberts_dns_negsol` hung).
4. **IDAS sensitivity `user_data`** — the "`None` means pass the integrator's
   `user_data`" fallback was never implemented at six call sites, so
   user-supplied sensitivity residuals got `None` and panicked.
5. **CVODES `cv_p`** — the sensitivity parameter array was an owned copy, so
   internal difference-quotient perturbations never reached the user's RHS.
   Now shared as `Rc<RefCell<Vec<sunrealtype>>>`, with a regression test whose
   negative control reproduces the original defect.

## Known limitations

1. **Unexercised code.** Every `*_bbdpre` module is compile-only — BBD
   preconditioning is MPI-only upstream, so no serial reference example can
   regression-test it. Same for the excluded KLU/SuperLU paths.
2. **Adjoint steppers are compiler-checked only.** `ERKStepCreateAdjointStepper`
   and its cluster are translated and build clean, but no serial reference
   example exercises them, so line-by-line reading is the only check they have
   had. Their `user_data` does not alias the forward memory (deviation class 6)
   — a caller porting from C must call `SUNAdjointStepper_SetUserData` itself.
3. **Accepted deviations.** Thirteen numbered classes in `ARCHITECTURE.md`,
   each verified unobservable on any path a valid serial example takes.

## Rules that still bind future work

* Byte-identical stdout against the upstream `.out`, noise-filtered
  symmetrically, is the only pass condition (`tools/verify_examples.sh`).
* When a shipped `.out` cannot be reproduced, the fallback bar is byte-identity
  against a locally built pristine upstream C binary — never tune an example to
  match a reference.
* Zero `unsafe`, zero FFI, zero external crates, zero build warnings.
* Once a crate's examples verify green they stay green; the cumulative gate is
  `tools/verify_examples.sh all`.
