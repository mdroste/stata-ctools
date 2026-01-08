# cqreg VCE Debugging Plan

## Executive Summary

The `cqreg` command produces coefficients that match Stata's `qreg`, but standard errors differ. Based on analysis of the codebase and reference implementations (R's `quantreg` package, Julia's `QuantileRegressions.jl`), the issue is likely in the **sparsity estimation** and potentially in how the VCE formula is applied.

## Key References

- [R quantreg package source](https://github.com/cran/quantreg/blob/master/R/quantreg.R)
- [R bandwidth.rq documentation](https://rdrr.io/cran/quantreg/man/bandwidth.rq.html)
- [Stata qreg manual](https://www.stata.com/manuals/rqreg.pdf)
- Hall & Sheather (1988), JRSS(B)

---

## Phase 1: Verify Coefficient Correctness

**Goal**: Confirm coefficients match `qreg` exactly before debugging VCE.

**Test Script**:
```stata
sysuse auto, clear
qreg price mpg weight
matrix b_qreg = e(b)

cqreg price mpg weight
matrix b_cqreg = e(b)

* Compare
matrix diff = b_qreg - b_cqreg
matrix list diff
```

**Expected**: Coefficient differences should be < 1e-6.

---

## Phase 2: Fix Hall-Sheather Bandwidth Formula

**Issue Found**: In `cqreg_sparsity.c:187`, there's a spurious `h *= 2.0` scaling factor that is NOT in the original Hall-Sheather formula.

**Current Code** (`cqreg_bandwidth_hsheather`):
```c
ST_double h = pow(z_alpha, 2.0/3.0) * pow(numer / denom, 1.0/3.0) * pow((ST_double)N, -1.0/3.0);
h *= 2.0;  // <-- THIS IS WRONG
return h;
```

**Correct Formula** (from R's `quantreg`):
```r
n^(-1/3) * qnorm(1 - alpha/2)^(2/3) * ((1.5 * f0^2)/(2 * x0^2 + 1))^(1/3)
```

Where:
- `x0 = qnorm(p)` (p = quantile, e.g., 0.5)
- `f0 = dnorm(x0)`
- `alpha = 0.05` (default)

**Fix**: Remove the `h *= 2.0` line entirely.

**Verification**:
```stata
* In R, compute bandwidth for comparison:
* bandwidth.rq(p=0.5, n=74, hs=TRUE, alpha=0.05)
* Should get approximately 0.156 for n=74, p=0.5
```

---

## Phase 3: Verify Sparsity Estimation

**Background**: Sparsity `s = 1/f(0)` where `f(0)` is the density of residuals at zero. Stata uses the "difference quotient" method:

```
s = (r_{[n*(q+h)]} - r_{[n*(q-h)]}) / (2h)
```

Where `r_{[k]}` is the k-th order statistic of residuals.

**Current Implementation** (`cqreg_estimate_sparsity` in `cqreg_sparsity.c`):
```c
ST_double x_lo = cqreg_sample_quantile(sp->sorted_resid, N, q_lo);
ST_double x_hi = cqreg_sample_quantile(sp->sorted_resid, N, q_hi);
sp->sparsity = (x_hi - x_lo) / dq;
```

**Potential Issues**:
1. The `cqreg_sample_quantile` uses R's "type 7" interpolation - verify Stata uses the same.
2. The bandwidth `h` from Phase 2 is used as a probability/quantile width, but may need adjustment.

**Debugging Steps**:
1. Print intermediate values: `h`, `q_lo`, `q_hi`, `x_lo`, `x_hi`, `sparsity`
2. Compare with `qreg, vce(iid)` output and extract Stata's implied sparsity from SE

**Implied Sparsity Calculation**:
```stata
* From qreg output, we can back-calculate sparsity:
* V = s^2 * q*(1-q) * (X'X)^{-1}
* SE = sqrt(V[k,k])
* s = SE / sqrt(q*(1-q) * (X'X)^{-1}[k,k])
```

---

## Phase 4: Create Diagnostic Comparison Script

**Create `benchmark/debug_vce.do`**:

```stata
* debug_vce.do - Compare cqreg vs qreg VCE
clear all
sysuse auto

* Run qreg
quietly qreg price mpg weight, quantile(0.5)
local se_mpg_qreg = _se[mpg]
local se_wt_qreg = _se[weight]
local se_cons_qreg = _se[_cons]

* Store qreg V matrix
matrix V_qreg = e(V)

* Run cqreg
quietly cqreg price mpg weight, quantile(0.5) verbose

* Compare
display "=== Standard Error Comparison ==="
display "          qreg        cqreg       ratio"
display "mpg:    " %9.4f `se_mpg_qreg' "  " %9.4f _se[mpg] "  " %6.3f _se[mpg]/`se_mpg_qreg'
display "weight: " %9.4f `se_wt_qreg' "  " %9.4f _se[weight] "  " %6.3f _se[weight]/`se_wt_qreg'
display "_cons:  " %9.4f `se_cons_qreg' "  " %9.4f _se[_cons] "  " %6.3f _se[_cons]/`se_cons_qreg'

* Display intermediate values from cqreg
display ""
display "=== cqreg Diagnostics ==="
display "N:        " %8.0f __cqreg_N
display "Sparsity: " %12.6f __cqreg_sparsity
display "Bandwidth:" %12.6f __cqreg_bandwidth

* Compute implied sparsity from qreg
* V = s^2 * q*(1-q) * (X'X)^{-1}
* For median: q*(1-q) = 0.25
* s^2 = V / (0.25 * (X'X)^{-1})
matrix XtX_inv = syminv(X'*X)  // Need to compute this
```

---

## Phase 5: Fix IID VCE Formula

**Current Formula** (`cqreg_vce_iid` in `cqreg_vce.c`):
```c
ST_double scale = sparsity * sparsity * q * (1.0 - q);
V[i * K + j] = scale * XtX_inv[i * K + j];
```

**Correct Formula**:
```
V = s^2 * tau*(1-tau) * (X'X)^{-1}
```

This looks correct IF the sparsity is correct. The issue is likely upstream in sparsity estimation.

**Alternative VCE Formulas to Consider**:

1. **Koenker-Bassett (1978) IID**: Exactly what we have.

2. **Powell Sandwich (Robust)**:
```
V = (X'X)^{-1} * X' * D * X * (X'X)^{-1}
```
Where D is diagonal with D[i,i] = tau^2 if r_i > 0, (1-tau)^2 if r_i < 0.

3. **Hendricks-Koenker (Cluster)**:
```
V = (X'X)^{-1} * Omega * (X'X)^{-1}
Omega = sum_g (X_g' * psi_g * psi_g' * X_g)
psi_i = tau - I(r_i < 0)
```

---

## Phase 6: Implement Robust VCE

**Issue**: Current `cqreg_vce_robust` mixes concepts incorrectly.

**Correct Robust VCE** (from quantreg's "nid" method):

1. Re-estimate at `tau + h` and `tau - h`
2. Compute fitted values difference: `dyhat = X * (beta_high - beta_low)`
3. Estimate local density: `f_i = 2h / max(dyhat_i, epsilon)`
4. Construct sandwich: `V = tau*(1-tau) * (X' * diag(f_i) * X)^{-1} * X'X * (X' * diag(f_i) * X)^{-1}`

**Simplified Alternative** (Powell-type):
```
V = (X'X)^{-1} * (sum_i w_i * X_i * X_i') * (X'X)^{-1}
```
Where `w_i = tau^2` if `r_i > 0`, else `(1-tau)^2`.

---

## Phase 7: Implement Cluster-Robust VCE

**Current Implementation Looks Correct** but verify:

1. Influence function: `psi_i = tau - I(r_i < 0)`
2. Cluster scores: `S_g = sum_{i in g} psi_i * X_i`
3. Meat: `M = sum_g S_g * S_g'`
4. Small-sample adjustment: `G/(G-1) * (N-1)/(N-K)`
5. Sandwich: `V = (X'X)^{-1} * M * (X'X)^{-1}`

**Test with**:
```stata
qreg price mpg weight, vce(cluster rep78)
cqreg price mpg weight, vce(cluster rep78)
```

---

## Phase 8: Multi-Dataset Testing

**Test Cases**:

| Dataset | Model | Quantiles | VCE Types |
|---------|-------|-----------|-----------|
| auto | price ~ mpg | 0.25, 0.5, 0.75 | IID, robust |
| auto | price ~ mpg weight | 0.5 | IID, robust, cluster(rep78) |
| nlsw88 | wage ~ ttl_exp | 0.1, 0.5, 0.9 | IID, robust |
| Generated (N=1000) | y ~ x1 x2 | 0.5 | IID, robust |
| Generated (N=10000) | y ~ x1 x2 x3 | 0.5 | IID, robust |

**Pass Criteria**:
- Coefficients within 1e-6 of qreg
- Standard errors within 1% of qreg

---

## Debugging Checklist

- [ ] Remove `h *= 2.0` from bandwidth calculation
- [ ] Verify bandwidth matches R's `bandwidth.rq(p, n, hs=TRUE)`
- [ ] Verify sample quantile computation matches Stata's convention
- [ ] Print and compare sparsity values with implied sparsity from qreg
- [ ] Verify (X'X)^{-1} computation is numerically stable
- [ ] Test edge cases: extreme quantiles (0.1, 0.9), small N
- [ ] Document any intentional differences from qreg

---

## Key Code Locations

| Component | File | Function |
|-----------|------|----------|
| Hall-Sheather bandwidth | `cqreg_sparsity.c:164-190` | `cqreg_bandwidth_hsheather` |
| Sparsity estimation | `cqreg_sparsity.c:319-373` | `cqreg_estimate_sparsity` |
| Sample quantile | `cqreg_sparsity.c:264-280` | `cqreg_sample_quantile` |
| IID VCE | `cqreg_vce.c:239-295` | `cqreg_vce_iid` |
| Robust VCE | `cqreg_vce.c:301-425` | `cqreg_vce_robust` |
| Cluster VCE | `cqreg_vce.c:431-542` | `cqreg_vce_cluster` |
| VCE dispatcher | `cqreg_vce.c:548-584` | `cqreg_compute_vce` |

---

## Expected Timeline

1. **Phase 1-2**: Fix bandwidth - immediate win expected
2. **Phase 3-4**: Verify sparsity - may require iteration
3. **Phase 5**: IID VCE should work once sparsity is correct
4. **Phase 6-7**: Robust/Cluster VCE - more complex debugging
5. **Phase 8**: Comprehensive validation
