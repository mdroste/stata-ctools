# HAC VCE Debugging Plan for civreghdfe

## Problem Summary

civreghdfe's HAC VCE differs from ivreghdfe by ~13% (ratio 1.134 for V[1,1]). The element-wise ratios vary (1.13, 0.93, 0.86), indicating the issue is in the structure of the meat matrix, not a simple scaling factor.

**Key findings:**
- Coefficients match perfectly
- Robust (HC) VCE matches perfectly
- HAC VCE differs
- Discrepancy exists even without fixed effects

## Root Cause Analysis

### What's Correct
1. Kernel weight formula: `w(j) = 1 - j/bw` for Bartlett - matches ivreg2
2. HAC formula structure: `S = Γ(0) + Σ w_j * (Γ(j) + Γ(j)')`
3. Score vectors: `u_i = PzX_i * e_i`
4. Observation ordering is correctly handled
5. DOF adjustment formula: `N / df_r`

### Likely Sources of Discrepancy

#### 1. Score Vector Definition
- **civreghdfe**: Uses `u = PzX * e` where `PzX = Z(Z'Z)^{-1}Z'X`
- **ivreghdfe**: Uses `s = Z * e` for the S matrix, then transforms

While mathematically equivalent, the order of operations differs:
- civreghdfe: `meat = (PzX)' * Ω * PzX`, then `V = XkX^{-1} * meat * XkX^{-1}`
- ivreghdfe: `S = Z' * Ω * Z`, then `V = 1/N * transform(S)` with normalized matrices

#### 2. Matrix Normalization
ivreghdfe uses normalized matrices:
- `QZZ = Z'Z / N` (normalized)
- Final VCE has `1/N` factor

civreghdfe uses unnormalized matrices:
- `XtPzX = X'PzX` (unnormalized)
- No `1/N` in final VCE

This could cause a scaling difference if not balanced correctly.

#### 3. Centering Option
ivreghdfe has a `center` option for HAC (centers residuals before outer product).
civreghdfe has `center` parameter but it's marked TODO and not implemented.

## Debugging Steps

### Step 1: Add Verbose HAC Output to civreghdfe

In `civreghdfe_vce.c`, add debug output (when verbose > 1):

```c
if (verbose > 1) {
    /* Print meat matrix sum before DOF adjustment */
    SF_display("HAC meat matrix (before sandwich):\n");
    for (j = 0; j < K_total; j++) {
        for (k = 0; k < K_total; k++) {
            char buf[50];
            sprintf(buf, "  meat[%d,%d] = %e\n", j, k, meat[j * K_total + k]);
            SF_display(buf);
        }
    }

    /* Print XkX_inv diagonal */
    SF_display("XkX_inv diagonal:\n");
    for (j = 0; j < K_total; j++) {
        char buf[50];
        sprintf(buf, "  XkX_inv[%d,%d] = %e\n", j, j, XkX_inv[j * K_total + j]);
        SF_display(buf);
    }
}
```

### Step 2: Compare with ivreghdfe Matrices

Run ivreghdfe and capture:
```stata
ivreghdfe y (x = z), absorb(fe) kernel(bartlett) bw(2)
matrix S = e(S)    // HAC meat matrix (K_iv x K_iv)
matrix list S
```

### Step 3: Verify Transformation Formula

The key equation to verify:
```
civreghdfe_meat = X'Z * (Z'Z)^{-1} * S * (Z'Z)^{-1} * Z'X
```

If this relationship holds, the meat matrices are consistent.

### Step 4: Check Normalization

Test if adding/removing 1/N factor fixes the discrepancy:
```c
/* Test: divide meat by N */
for (i = 0; i < K_total * K_total; i++) {
    meat[i] /= (ST_double)N;
}
```

## Proposed Fix Approaches

### Approach A: Match ivreghdfe's Score Vector Approach

Instead of computing `meat = (PzX)' * Ω * PzX`, compute:
1. `S = Z' * Ω * Z` (K_iv x K_iv)
2. `meat = X'Z * (Z'Z)^{-1} * S * (Z'Z)^{-1} * Z'X` (K_total x K_total)

This matches ivreghdfe's computation order exactly.

### Approach B: Add Normalization Factor

If the issue is pure normalization:
1. Track the relationship between ivreghdfe and civreghdfe meat matrices
2. Apply correction factor (e.g., 1/N if needed)

### Approach C: Implement Centering

Implement the TODO center option:
```c
if (center) {
    /* Center residuals: e_centered = e - mean(e) */
    ST_double e_mean = 0.0;
    for (i = 0; i < N; i++) e_mean += resid[i];
    e_mean /= N;
    for (i = 0; i < N; i++) resid[i] -= e_mean;
}
```

## Testing Plan

1. Add debug output and compare intermediate matrices
2. Implement fix
3. Re-run validation tests
4. Verify all 9 failing HAC tests pass
5. Verify robust VCE still matches (regression test)

## Files to Modify

1. `src/civreghdfe/civreghdfe_vce.c` - HAC meat computation
2. `src/civreghdfe/civreghdfe_tests.c` - Diagnostic statistics (idstat, widstat, sargan)
3. `build/civreghdfe.ado` - If any Stata-side changes needed

## Key Numerical Finding

From testing on the auto dataset with bw=2, absorb(foreign):

| Metric | Value |
|--------|-------|
| V_civ / V_iv | 1.134 |
| V_civ_no_dof / V_iv | 1.073 |
| DOF adjustment | 1.057 |

**Interpretation:**
- ~7% of the discrepancy comes from the meat matrix computation
- ~6% comes from DOF adjustment differences
- Total: ~13% discrepancy

This rules out a simple missing factor. The meat matrix itself differs by 7%.

## Recommended Next Step

Add verbose debug output to `civreghdfe_vce.c` to print:
1. Sum of meat matrix diagonal elements
2. Frobenius norm of meat matrix
3. First few elements of the u (score) vectors

Compare these with ivreghdfe's e(S) matrix to identify where the 7% difference originates.

## Priority

High - 7 of 9 failing tests are HAC-related
