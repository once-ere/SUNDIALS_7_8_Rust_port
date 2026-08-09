# SUNDIALS 7.8.0 — pure-Rust port

A line-by-line translation of [SUNDIALS](https://github.com/LLNL/sundials)
7.8.0 into safe Rust. **No `unsafe`, no FFI, no external crates, no build
warnings.** Acceptance is byte-identical printed output against the upstream C
reference examples.

**→ Read [`sundials.md`](sundials.md) — the full guide: crate map, worked
example, C-to-Rust API conventions, verification results and known
deviations.**

## Headline facts

* 7 crates: `sundials_core` plus `cvode_rs`, `cvodes_rs`, `kinsol_rs`,
  `ida_rs`, `idas_rs`, `arkode_rs`. Solver crates depend on the core, never on
  each other.
* 141 modules, one per upstream C file, keeping the exact C function names,
  constants and return-flag conventions (`CV_SUCCESS = 0`; negative fatal,
  positive recoverable).
* Serial only. No MPI, GPU, KLU, SuperLU, LAPACK, Fortran or XBraid backends.
* Example gate: of the 199 reference `(example, argv)` variants,
  **76 are byte-identical, 25 are documented reference-side exceptions
  (0 port defects), 20 are excluded as KLU/SuperLU**, and the **78 ARKODE
  variants are not yet swept** — the ARKODE library is ported and warning-free
  but unproven against reference output.
* `cargo build --workspace` → zero warnings. `cargo test --workspace --lib` →
  25 passed.

## Quick start

```sh
cargo build --workspace
cargo run -p cvode_rs --example cvRoberts_dns
tools/verify_examples.sh cvode_rs        # or `all`, or `list`
```

```rust
use cvode_rs::prelude::*;
```

## Other documents

| file | contents |
|---|---|
| [`sundials.md`](sundials.md) | public guide — start here |
| [`ARCHITECTURE.md`](ARCHITECTURE.md) | handle model, locked porting patterns, numbered deviation classes |
| [`VERIFICATION.md`](VERIFICATION.md) | per-variant matrix and the evidence behind every exception |
| [`PROGRESS.md`](PROGRESS.md) | per-file port status |
| [`STATUS.md`](STATUS.md) | what is done, what remains, how to resume |

## Licence

Derivative work of SUNDIALS, **BSD-3-Clause**, Copyright © 2002–2026 Lawrence
Livermore National Security, Southern Methodist University, University of
Maryland Baltimore County and the SUNDIALS contributors. The deterministic
`pow` in `crates/sundials_core/src/sundials_math.rs` is a port of the ARM
optimized-routines / musl implementation, **MIT**, Copyright © 2018 Arm
Limited. Not an LLNL product; not endorsed by the SUNDIALS project. See
`sundials.md` §8.
