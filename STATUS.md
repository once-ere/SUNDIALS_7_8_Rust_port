# STATUS — resume point

Snapshot for resuming after an interruption. Read this plus `PROGRESS.md`,
`VERIFICATION.md` and `git log` before restarting any work.

## Library ports (all six crates build warning-free)

| crate | modules | state |
|---|---|---|
| `sundials_core` | 51 (incl. the 9 Phase-1 deferrals) | committed, 20/20 tests |
| `cvode_rs` | 12 | committed + examples verified (Phase 2 gate) |
| `cvodes_rs` | 17 | committed; verification pass incomplete |
| `kinsol_rs` | 8 | committed; verification pass incomplete |
| `ida_rs` | 9 | committed; verification pass incomplete |
| `idas_rs` | 13 | committed; verification pass incomplete |
| `arkode_rs` | — | port in flight (33 C files, ~54k lines) |

`cargo build --workspace` → zero warnings. `cargo test --workspace --lib` → green.

## Verification state

Only `cvode_rs` has been through the full example gate:
12/21 variants byte-IDENTICAL, 9 documented exceptions (see
`VERIFICATION.md`). Examples for the other crates are written or in
flight but not yet swept.

## Known open items

1. **cvodes `cv_p` aliasing** — C shares the user's parameter pointer so
   internal difference-quotient perturbations reach the user RHS. Fix in
   flight: share the array as `Rc<RefCell<Vec<sunrealtype>>>`. Until it
   lands, FSA examples that pass a real parameter array compute wrong
   sensitivities.
2. **Adversarial verification incomplete** for cvodes / kinsol / ida /
   idas — the fan-outs were cut short by a session limit. Re-run per
   crate (the Phase-2 pattern: one verifier per module group, instructed
   to refute, then a fixer that re-checks each finding against the C).
3. **Unexercised code**: every `*_bbdpre` module is compile-only — no
   serial reference example uses BBD preconditioning (it is MPI-only
   upstream), so the example gate can never regression-test it.
4. **kinsol duplicate QR-add** — `kinsol_orth::kinQRAdd` and an inlined
   copy inside `kinsol.rs` implement the same ICWY update; verified
   semantically equivalent, worth collapsing once
   `kinAnalytic_fp` with `orth_aa 1` can prove it byte-for-byte.

## Reference-environment exception classes (established in Phase 2)

Recorded in `VERIFICATION.md`; expect them to recur in every crate:
`ref-libm` (shipped `.out` embeds the generating host's glibc
sin/exp/pow rounding inside the integration feedback loop), stale
upstream references, and LAPACK→native last-digit drift. The bar for
accepting one is byte-identity against a locally built pristine
upstream C binary (`-O3 -ffp-contract=off`, reference config).

## Remaining sequence

arkode port → arkode examples → per-crate verification sweeps → full
cumulative `tools/verify_examples.sh all` → `sundials.md` public guide →
tag `v1-complete` → push to the GitHub remote (`origin`, already wired).
