# Project Instructions

## Environment

- **Stata CLI**: Available via the bash alias `stata`. Use directly.
  - Example: `stata -b do myscript.do`
- **Source Directory**: `./src/` - C source files
- **Build Directory**: `./build/` - Compiled plugins (`.plugin`) and Stata files (`.ado`, `.sthlp`)
- **Benchmark Directory**: `./benchmark/` - Stata `.do` files for speed/correctness benchmarks

## Building the Plugin

```bash
make              # Build for current platform
make all          # Build all platform plugins
make native       # Build optimized for current architecture
make check        # Check build dependencies
make clean        # Remove compiled files
```

Output plugins by platform:
- `ctools_mac_arm.plugin` - macOS Apple Silicon
- `ctools_mac_x86.plugin` - macOS Intel
- `ctools_windows.plugin` - Windows x64
- `ctools_linux.plugin` - Linux x64
- `ctools.plugin` - Generic/native

Requirements: Xcode CLT (macOS), optionally `brew install libomp` for OpenMP support.

## Project Overview

**ctools** is a suite of C-accelerated Stata programs with multithreading support. These are high-performance drop-in replacements for native Stata commands.

### Commands

| Command | Replaces | Description |
|---------|----------|-------------|
| `csort` | `sort` | Parallel LSD radix sort |
| `cmerge` | `merge` | C-accelerated merge (1:1, m:1, 1:m, m:m) |
| `cimport` | `import delimited` | Multi-threaded CSV import |
| `cexport` | `export delimited` | Parallel CSV export |
| `creghdfe` | `reghdfe` | HDFE regression with CG solver |
| `cqreg` | `qreg` | Quantile regression with IPM solver |

### Design Principles

1. **Speed is critical** - Every operation is optimized for performance
2. **Modular architecture** - Core routines (load, store, sort) are reusable across commands
3. **Parallel by default** - Uses pthreads for parallel I/O and radix sort

## Source Code Organization

### Core Infrastructure (`./src/`)

| File | Purpose |
|------|---------|
| `stplugin.c/h` | Stata Plugin Interface (SPI) - **DO NOT MODIFY** |
| `ctools_plugin.c` | Main dispatcher - routes commands to handlers |
| `ctools_config.h` | Performance tuning parameters (threads, thresholds) |
| `ctools_types.h/c` | Core data structures (`stata_data`, `stata_variable`) |
| `ctools_error.h/c` | Unified error reporting utilities |
| `ctools_threads.h/c` | Thread pool utilities for parallel execution |
| `ctools_timer.h/c` | Cross-platform high-resolution timing |
| `ctools_data_load.c` | Stata → C data transfer (parallel, 8x unrolled) |
| `ctools_data_store.c` | C → Stata data transfer (parallel, 8x unrolled) |
| `ctools_sort_radix_lsd.c` | Parallel LSD radix sort (numeric + string) |

### Command Modules

Each command has its own subdirectory with `*_impl.c` and `*_impl.h`:

- `./src/csort/` - Sort implementation
- `./src/cmerge/` - Merge implementation
- `./src/cimport/` - CSV import implementation
- `./src/cexport/` - CSV export implementation
- `./src/creghdfe/` - HDFE regression (multiple files for OLS, VCE, solver, etc.)
- `./src/cqreg/` - Quantile regression (IPM solver, VCE, sparsity estimation)

### Build Output (`./build/`)

| File | Type |
|------|------|
| `*.plugin` | Compiled C plugins |
| `*.ado` | Stata command wrappers |
| `*.sthlp` | Stata help files |

## Key Data Structures

### `stata_data` (ctools_types.h)
Main container for dataset in C memory:
```c
typedef struct {
    size_t nobs;            // Number of observations
    size_t nvars;           // Number of variables
    stata_variable *vars;   // Array of variables [nvars]
    size_t *sort_order;     // Permutation array [nobs]
} stata_data;
```

### `stata_variable` (ctools_types.h)
Single variable's data:
```c
typedef struct {
    stata_vartype type;     // STATA_TYPE_DOUBLE or STATA_TYPE_STRING
    size_t nobs;            // Number of observations
    union {
        double *dbl;        // Numeric: contiguous double[nobs]
        char **str;         // String: char*[nobs], heap-allocated
    } data;
    size_t str_maxlen;      // Max string length (string vars only)
} stata_variable;
```

## Typical Data Flow

1. **Load**: `ctools_data_load()` - Reads from Stata via SPI (parallel by variable)
2. **Process**: Command-specific logic (sort, merge, etc.)
3. **Store**: `ctools_data_store()` - Writes back to Stata via SPI (parallel by variable)
4. **Free**: `stata_data_free()` - Releases all allocated memory

## Performance Optimizations

- **8x loop unrolling** in data load/store for better pipelining
- **Parallel variable I/O** - one thread per variable when nvars >= 2
- **LSD radix sort** with parallel histogram and scatter phases
- **Pointer swapping** instead of memcpy between radix passes
- **Early exit** for uniform byte distributions (skip no-op passes)
- **Pre-cached string lengths** for string sorting
- **Aligned memory allocation** (64-byte cache line alignment)
- **SD_FASTMODE** compile flag disables SPI bounds checking

## Configuration (`ctools_config.h`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `MIN_PARALLEL_BYTES` | 1MB | Threshold for parallel I/O |
| `MIN_OBS_PER_THREAD` | 100,000 | Min observations per thread for sort |
| `NUM_THREADS` | 8 | Max threads for parallel operations |
| `CACHE_LINE_SIZE` | 64 | Alignment for allocations |
| `RADIX_BITS` | 8 | Radix size (256 buckets per pass) |

## Debugging Tips

- Add `verbose` option to any command for timing breakdown
- Use `timeit` option where available for detailed timing
- Check Stata scalars after commands (e.g., `_csort_time_load`)
- Plugin errors return Stata error codes (198 = syntax, 920 = memory, etc.)

## Adding a New Command

1. Create `./src/newcmd/newcmd_impl.c` and `newcmd_impl.h`
2. Implement `ST_retcode newcmd_main(const char *args)`
3. Add dispatch case in `ctools_plugin.c`
4. Add sources to Makefile (`NEWCMD_SRCS`, `NEWCMD_HEADERS`)
5. Create `./build/newcmd.ado` wrapper
6. Create `./build/newcmd.sthlp` help file

## cqreg: Quantile Regression

### Overview

`cqreg` is a C-accelerated quantile regression command that aims to replicate Stata's `qreg` output exactly while providing better performance for large datasets.

### Source Files (`./src/cqreg/`)

| File | Purpose |
|------|---------|
| `cqreg_main.c` | Entry point, argument parsing, orchestration |
| `cqreg_fn.c` | Frisch-Newton interior point method (IPM) solver |
| `cqreg_ipm.c` | IPM state management and utilities |
| `cqreg_vce.c` | Variance-covariance estimation (IID, robust, cluster) |
| `cqreg_sparsity.c` | Sparsity/density estimation for VCE |
| `cqreg_linalg.c` | Linear algebra (Cholesky, matrix ops) |
| `cqreg_hdfe.c` | High-dimensional fixed effects support |
| `cqreg_types.c` | Data structures and state management |

### Algorithm

The solver uses the Frisch-Newton primal-dual interior point method:
1. Initialize from OLS solution
2. Iterate with Mehrotra predictor-corrector steps
3. Converge when duality gap < tolerance (1e-12)

Reference: QuantileRegressions.jl (Julia implementation)

### VCE Computation

The variance-covariance matrix depends on sparsity estimation:
```
V = τ(1-τ) * s² * (X'X)⁻¹  (IID case)
```
where `s = 1/f(0)` is the sparsity (inverse density at the quantile).

### VCE Types and Density Methods

| VCE Type | denmethod | Description | Matches |
|----------|-----------|-------------|---------|
| IID | residual | Difference quotient on residuals | `qreg vce(iid, residual)` exactly |
| IID | fitted | Kernel density on residuals | Similar to qreg default |
| Robust | residual | Score-based sandwich with scalar sparsity | Equals IID for median |
| Robust | fitted | Powell sandwich with per-obs densities | `qreg vce(robust)` (within 1%) |
| Cluster | any | Cluster-robust with influence functions | Unique to cqreg |

**Syntax:**
```stata
cqreg y x, vce(iid) denmethod(residual)   * Matches qreg vce(iid, residual) exactly
cqreg y x, vce(robust) denmethod(fitted)  * Matches qreg vce(robust) closely
cqreg y x, vce(cluster grp)               * Cluster VCE (qreg doesn't support this)
```

**Important Notes:**

1. **denmethod(residual)** uses difference quotient sparsity estimation, excluding LP basis observations (|residual| < 1e-8). This matches `qreg vce(iid, residual)` exactly.

2. **denmethod(fitted)** uses kernel density estimation. For robust VCE, it implements the Powell sandwich with per-observation density weights, closely matching `qreg vce(robust)`.

3. **Cluster VCE** is NOT supported by Stata's `qreg` command. cqreg provides this as additional functionality with proper influence-function-based standard errors.

### Key Functions

- `cqreg_fn_solve()` - Main IPM solver loop
- `cqreg_estimate_sparsity()` - Difference quotient sparsity estimator (residual method)
- `cqreg_estimate_fitted_density()` - Per-observation density estimation (fitted method)
- `cqreg_vce_iid()` - IID variance computation with scalar sparsity
- `cqreg_vce_robust()` - Robust sandwich variance with scalar sparsity
- `cqreg_vce_robust_fitted()` - Powell sandwich with per-observation densities
- `cqreg_vce_cluster()` - Cluster-robust variance

### Testing

```stata
* Compare cqreg vs qreg
sysuse auto, clear

* Exact match with qreg residual method
qreg price mpg, vce(iid, residual)
cqreg price mpg, vce(iid) denmethod(residual)
* Coefficients and standard errors match exactly

* Close match with qreg robust (within 1%)
qreg price mpg, vce(robust)
cqreg price mpg, vce(robust) denmethod(fitted)

* Cluster VCE (unique to cqreg)
cqreg price mpg, vce(cluster rep78)
```

### Implementation Details

**Sparsity Estimation** (`cqreg_sparsity.c`):
- Uses Hall-Sheather bandwidth formula
- Excludes LP basis observations (|residual| < 1e-8)
- Uses adjusted N for bandwidth computation
- Employs Stata-compatible ceiling(N*q) quantile method

**Reference Implementations:**
- R's `quantreg` package: [bandwidth.rq](https://rdrr.io/cran/quantreg/man/bandwidth.rq.html)
- Julia's `QuantileRegressions.jl`

**Hall-Sheather Bandwidth Formula:**
```
h = n^(-1/3) * z_{1-α/2}^(2/3) * ((1.5 * φ(x_p)²) / (2*x_p² + 1))^(1/3)
```
Where `x_p = Φ⁻¹(p)`, `φ` = standard normal PDF, `α = 0.05`.
