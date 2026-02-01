# ctools New Features

This document describes features available **only** in ctools commands—functionality that does not exist in the standard Stata commands or popular user-written alternatives they replace.

---

## cencode / cdecode

Stata's native `encode` and `decode` commands only process one variable at a time and require creating new variables.

**ctools adds:**

- **Varlist support**: Encode or decode multiple variables in a single command
- **Replace option**: Modify variables in-place instead of creating new ones

```stata
* Encode three variables at once, replacing in-place
cencode var1 var2 var3, replace
```

---

## cbinscatter

The original `binscatter` package uses the "classic" method for covariate control, which residualizes both x and y before creating bins.

**ctools adds:**

- **Binsreg methodology**: `method(binsreg)` implements the Cattaneo et al. (2024) "On Binscatter" approach, which creates bins on raw x values and then residualizes within bins. This avoids distortions that can occur with the classic method when the conditional expectation function is nonlinear.

```stata
* Use the binsreg covariate control method
cbinscatter y x, controls(z) method(binsreg)
```

---

## csort

Stata's native `sort` command provides no user control over memory usage or algorithm selection.

**ctools adds:**

- **Algorithm selection**: `algorithm()` option to choose from `auto`, `counting`, `lsd`, `msd`, `ips4o`, `timsort`, `sample`, or `merge`
- **Streaming mode**: `stream()` option reduces memory usage for wide datasets by permuting non-key variables sequentially

```stata
* Use counting sort for integer data
csort year, algorithm(counting)

* Streaming mode for memory-constrained sorting
csort id, stream(4)
```

---

## cqreg

Stata's native `qreg` command does not support fixed effects absorption.

**ctools adds:**

- **HDFE support**: `absorb()` option to absorb high-dimensional fixed effects in quantile regression

```stata
* Quantile regression with two-way fixed effects
cqreg ln_wage age tenure, absorb(idcode year) quantile(0.5)
```
