# Is the deterministic `pow` 100% correct?

**Short answer: it was not, and it still is not — in one precisely bounded
sense. Within the domain SUNDIALS actually uses it is now bit-exact over
5.9 million measured inputs. Outside that domain, 2 inputs in 20 million
still disagree with the reference by 1 ulp.**

This document records the investigation, the measurements, and the limits of
what can be claimed.

## 1. Why "correct" is the wrong word

`crates/sundials_core/src/sundials_math.rs` contains `pow_glibc`, a Rust port
of the ARM optimized-routines `pow` (via musl, MIT) — the same algorithm
glibc >= 2.28 ships as `sysdeps/ieee754/dbl-64/e_pow.c`.

It is **not** correctly rounded, and must not be. Its documented worst-case
error is ~0.54 ulp. The goal is not mathematical correctness but *bit-exact
agreement with the specific libm build that generated this project's
reference outputs*. This was tested directly in an earlier phase: forcing
correct rounding (80-digit decimal, iterated to fixpoint) did **not**
reproduce the reference output. We are deliberately reproducing another
implementation's rounding behaviour, errors included.

Why it matters at all: `SUNRpowerR` is called from the step-size controller
of every solver. A 1-ulp difference in one `eta` changes the next step size,
which changes every subsequent step. There is no error budget to absorb it —
the acceptance criterion is byte-identical printed output.

## 2. The defect

The algorithm is written to compile two ways. Its own source says so:

> `/* Without fma the worst case error is 0.25/N ulp larger. */`

glibc on x86-64 dispatches to `__ieee754_pow_fma`, built from
`sysdeps/x86_64/fpu/multiarch/e_pow-fma.c` with `-mfma -mavx2
-ffp-contract=fast`. That flag contracts **every** eligible `a*b + c` in the
translation unit into a single fused multiply-add. Fused and unfused round
differently whenever the exact result lies near a rounding boundary.

Fused multiply-adds enter that build two ways:

1. **Explicit** — the C calls `__builtin_fma` under `#if __FP_FAST_FMA`.
   The port had all three of these right from the beginning.
2. **Implicit** — the compiler contracts ordinary expressions.
   **The port originally had none of these.**

Disassembling a `-ffp-contract=fast` build of the same C shows **21 FMA
instructions**. The port had **5**.

## 3. Method — measure, never infer

Which multiply folds into which add inside a nested expression is a compiler
scheduling decision. Guessing it is exactly how the first attempt went wrong,
so every site here was settled empirically:

* Reference C (`pow.c`, `pow_data.c`, `exp_data.c`, plus glibc's `e_pow.c`
  and `e_pow_fma.c`) compiled locally into oracle binaries under
  gcc x clang, arm64 x x86-64, `-ffp-contract=fast` x `off`.
* A standalone Rust harness reading `(x_bits, y_bits)` and emitting result
  bits, driven over the same corpora as the oracles.
* Differential comparison, then per-site bisection.

The gcc and clang `fast` oracles agree with each other, which is what makes
the oracle trustworthy.

## 4. Results

Corpus: 20,000,000 random finite `(x, y)` bit-pattern pairs.

| configuration | diff lines vs `gcc -ffp-contract=fast` |
|---|---:|
| uncontracted reference build (`-ffp-contract=off`) | 12,350 |
| port as shipped (5 fused sites) | 228 |
| **port with the two verified sites fused** | **4** |
| naive "fuse everything" (21 sites) | 494 |

Two things worth reading carefully:

* **Fusing everything made it worse** (494 > 228). The compiler does *not*
  contract every eligible site. This is the empirical proof that the map has
  to be measured.
* **`lo4` was a trap.** Fusing `lo4 = t2 - hi + ar2` fixed one failing input
  and broke 93 others (4 -> 188 diff lines). Rejected.

### The two sites that are genuinely contracted

1. `pow_exp_inline` — the exponential polynomial
   `tmp = tail + r + r2*(C2 + r*C3) + r2*r2*(C4 + r*C5)`
2. `pow_exp_specialcase`, `k > 0` branch — `0x1p1009 * (scale + scale*tmp)`,
   reached whenever the result is near overflow. This was the larger of the
   two by frequency and had been missed entirely.

### The domain that actually matters

SUNDIALS calls `SUNRpowerR(bias*dsm, ±1/order)`, so `x` in `(0, ~100]` and
`|y| <= 1`. Corpus: 5,900,000 inputs — `x` uniform in `(0, 100]` against
`y = ±1/k` for `k = 1..13` (every integrator order), plus `y` uniform in
`[-1, 1]`.

| configuration | mismatches |
|---|---:|
| port as shipped | 4 |
| **port after this fix** | **0** |

**The shipped implementation was wrong inside the operating domain.** That is
the finding that justifies the change.

## 5. What is still not proven

* **The residual 2 inputs.** In the 20M random corpus, two still differ by
  1 ulp: `x = 0.60506369761398, y = -35.20830257243234` and
  `x = 1.700984188552496, y = 44.56990535118902`. Both have `|y|` far above
  anything SUNDIALS evaluates. No single further contraction fixes them, and
  every candidate tried made the overall picture worse. They are recorded
  here rather than hidden.
* **Empirical, not exhaustive.** The input space is 2^128 pairs. Zero
  mismatches over 5.9M domain inputs is strong evidence, not proof.
* **Oracle provenance.** The oracles were built with gcc/clang on arm64; the
  references were generated by GCC-built glibc on x86-64. Contraction choices
  could in principle differ. The mitigating evidence is that the gcc and
  clang `fast` builds agree with each other, and that the frozen test triples
  in `sundials_math.rs` came from a run that reproduced upstream reference
  output byte-for-byte.

## 6. Effect on the verification gate: none

`tools/verify_examples.sh all` before and after the fix:

```
127 IDENTICAL / 52 documented / 20 excluded  (199 variants)
```

Byte-for-byte identical summaries — **no variant changed status**.

This is worth stating plainly rather than burying. The defect was **real but
latent**: none of the 199 example variants happened to evaluate `pow` at one
of the affected inputs. The example gate is a *sampling* instrument for this
class of bug, not a proof, and it would not have caught this. It took a
20-million-input differential test against a compiled oracle.

It also resolves a doubt raised earlier. When the first FMA fix landed, the
six LSRK variants classified `ref-libm` were called into question, because
that classification had been derived pre-fix and cited "a single 1-ULP `pow`
disagreement". Completing the contraction map changes none of their output,
so **`pow` contraction does not explain the LSRK divergence** and the
`ref-libm` classification stands on its own evidence.

## 7. Regression protection

The pre-existing `pow_glibc_bits` unit test passed against the *broken*
implementation — it only sampled easy inputs, which is precisely how this
survived. Any future change to `pow_glibc` should be re-validated with the
differential harness against a freshly built oracle, not against the example
gate, which is blind to it.

## 8. Verdict

| question | answer |
|---|---|
| Correctly rounded? | No — and deliberately so; correct rounding fails the acceptance test. |
| Bit-exact with the reference libm across all doubles? | **No.** 2 inputs in 20M differ by 1 ulp, both at very large `\|y\|`. |
| Bit-exact across the domain SUNDIALS evaluates? | **Yes**, over 5.9M measured inputs — after this fix. Before it, no. |
| Does any of this change the 199-variant gate? | No. The defect was latent. |
