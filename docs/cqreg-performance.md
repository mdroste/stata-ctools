# Bounded cqreg performance checks

The changes target allocations and order statistics that become costly as the
observation count grows. The estimator, quantile definitions, convergence
tolerance, and variance formulas are unchanged.

- The default Frisch–Newton solver now uses its existing lightweight workspace.
  Legacy IPM directions and per-thread observation buffers remain available for
  experimental preprocessing, which uses the other solver.
- Both workspaces allocate only `N` doubles for the weighted right-hand side,
  replacing an unused `N × K` allocation.
- Predictor and corrector steps reuse the same Cholesky factor. The numerical
  QR retry skips cross-products and normal-equation right-hand sides it never uses.
- Response quantiles use selection instead of a full sort. Residual quantiles
  use three-way partitioning so repeated values cannot cause quadratic scanning.
  A work budget falls back to sorting the remaining partition on difficult inputs.
- The pure C matrix-vector fallback now honors BLAS's `beta=0` convention:
  it overwrites the output without reading its uninitialized contents.

## Measurements

Apple Silicon, Accelerate and four OpenMP threads; medians of three native runs
against a saved copy of the previous source. These are individual components,
not complete Stata command timings.

| Component | Size | Before | After |
|---|---|---:|---:|
| Response quantile | 200,000 observations | 15.887 ms | 1.759 ms |
| Main FN solve, including crossover | 200,000 observations, 12 coefficients | 102.921 ms | 103.264 ms |
| Residual sparsity with ties | 12,000 observations | 3.774 ms | 0.013 ms |
| Main plus auxiliary aligned workspace allocations | 200,000 observations, 12 coefficients, 4 threads | 120.007 MB | 59.203 MB |

Coefficients, iteration counts, objectives, quantiles and sparsity estimates
were identical in the before/after comparisons. Ordinary solver timing was
essentially unchanged; the substantial measured gains are in selection and
workspace size. A cache-blocked weighted cross-product experiment was slower
and was not retained.

For large samples, the leading allocation saving during the auxiliary fit is
`8*N*(10 + 2*K + T)` bytes, where `K` includes the constant and `T` is the old
workspace's thread count. At five million observations, 12 coefficients and
four threads, this removes about **1.52 GB of requested workspace allocation**.
This is an allocation calculation, not a measured reduction in resident memory
or a timing extrapolation. No million-observation fits were run.

## Reproducing the checks

From the repository root:

```sh
python3 validation/test_cqreg_native.py
CTOOLS_SANITIZERS=address,undefined python3 validation/test_cqreg_native.py
CTOOLS_TEST_OPENMP_PREFIX=/opt/homebrew/opt/libomp \
  python3 validation/test_cqreg_native.py --benchmark --rows 200000 --columns 12
```

Use `--baseline-dir /path/to/saved/cqreg` for an automated comparison against
an earlier copy of `src/cqreg`. The benchmark defaults to 20,000 observations
and three repetitions, with caps of 200,000 observations and five repetitions.
The default correctness checks use at most 8,193 values and 303-observation
regressions, including forced QR solves, auxiliary solves, and nonconvergence.

`validation/validate_cqreg_performance.do` adds a 603-observation Stata smoke
test against `qreg`, covering sample filtering, three quantiles, all density
methods, robust and clustered VCE, ties, and iteration-limit errors. Its optional
argument selects an isolated ado/plugin build directory. Invoke Stata through
the `stata` shell alias and exit cleanly with `exit, clear`.

Both native BLAS and pure C paths passed AddressSanitizer and
UndefinedBehaviorSanitizer. Native OpenMP checks passed as well. The Stata smoke
test passed using a temporary serial plugin build: the installed static OpenMP
runtime targets macOS 26, so the distribution build's macOS 11 compatibility
check correctly rejected that runtime. Tracked plugin binaries were not replaced.
