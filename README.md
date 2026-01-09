# ctools

High-performance C-accelerated drop-in replacements for Stata commands.

## Overview

**ctools** provides drop-in replacements for common Stata commands that run substantially faster on large datasets. All commands are implemented in optimized C with OpenMP parallelization and compile to native plugins for Windows, macOS (Intel and Apple Silicon), and Linux.

| Stata command | ctools command | Purpose | Typical speedup |
| --- | --- | --- | ---: |
| `sort` | `csort` | Sort dataset | **2-4x** |
| `import delimited` | `cimport delimited` | Import CSV/TSV files | **50x** |
| `export delimited` | `cexport delimited` | Export CSV/TSV files | **25x** |
| `merge` | `cmerge` | Merge (join) datasets | **2-5x** |
| `reghdfe` | `creghdfe` | High-dimensional FE regression | **10x** |
| `qreg` | `cqreg` | Quantile regression | **4x** |

## Installation

### From GitHub (recommended)

```stata
net install ctools, from("https://raw.githubusercontent.com/your-username/stata-ctools/main/build")
```

### Manual Installation

1. Download or clone this repository
2. Copy the contents of the `build/` directory to a location on Stata's adopath (e.g., your personal ado directory)

The package automatically detects your operating system and architecture, loading the appropriate precompiled plugin.

## Commands

### csort - High-performance sorting

Drop-in replacement for Stata's `sort` command using a parallel least significant digit (LSD) radix sort.

```stata
* Basic usage
csort varname

* Sort by multiple variables
csort var1 var2 var3

* With timing output
csort myvar, verbose timeit
```

### cimport delimited - Fast CSV import

High-performance replacement for `import delimited` with multi-threaded, super efficient data reading.

```stata
* Import CSV file
cimport delimited using "data.csv", clear

* Import tab-delimited file
cimport delimited using "data.tsv", clear delimiters(tab)

* Import specific row range
cimport delimited using "bigdata.csv", clear rowrange(1000:2000)

* Fast mode (writes temp DTA, then loads)
cimport delimited using "bigdata.csv", clear fast verbose
```

**Options:**
- `delimiters(chars)` - Field delimiter (default: comma)
- `varnames(1|nonames)` - How to read variable names
- `case(preserve|lower|upper)` - Variable name case
- `rowrange([start][:end])` - Import specific rows
- `fast` - Use fast DTA mode for large files
- `verbose` - Show progress and throughput

### cexport delimited - Fast CSV export

High-performance replacement for `export delimited` with multi-threaded, super-efficient data writing. `cexport delimited` replaces `export delimited`.

```stata
* Export all variables
cexport delimited using "output.csv", replace

* Export selected variables
cexport delimited id name value using "output.csv", replace

* Tab-delimited output
cexport delimited using "output.tsv", delimiter(tab) replace

* Export subset
cexport delimited using "subset.csv" if year > 2020, replace
```

**Options:**
- `delimiter(char)` - Field delimiter (default: comma)
- `replace` - Overwrite existing file
- `novarnames` - Omit header row
- `quote` - Quote all string fields
- `verbose`, `timeit` - Show progress/timing

### cmerge - Fast dataset merging

Drop-in replacement for `merge`, leveraging ctools for internal sorting. `cmerge` replaces `merge`.

```stata
* One-to-one merge
cmerge 1:1 id using "data2.dta"

* Many-to-one merge
cmerge m:1 state year using "statedata.dta", keep(match)

* Keep only specific variables from using dataset
cmerge 1:1 id using "fulldata.dta", keepusing(var1 var2)

```

**Options:**
- `keep(match|master|using)` - Which observations to keep
- `assert(results)` - Verify merge results
- `generate(varname)` - Name of merge indicator (default: `_merge`)
- `nogenerate` - Don't create merge indicator
- `keepusing(varlist)` - Variables to keep from using data
- `sorted` - Skip sorting (data already sorted)
- `verbose` - Show timing information

### creghdfe - High-dimensional fixed effects regression

Drop-in replacement for `reghdfe` written in pure C with a variety of speed improvements. `creghdfe` replaces `reghdfe`.


```stata
* Basic regression with fixed effects
creghdfe y x1 x2, absorb(firm year)

* Clustered standard errors
creghdfe y x1 x2, absorb(firm year) vce(cluster firm)

* Robust standard errors
creghdfe y x1 x2, absorb(firm year) vce(robust)

* With weights
creghdfe y x1 x2 [aw=weight], absorb(firm year)
creghdfe y x1 x2 [fw=freq], absorb(firm year)
creghdfe y x1 x2 [pw=sampwt], absorb(firm year)  // forces robust VCE

* Verbose timing
creghdfe y x1 x2, absorb(firm year) verbose
```

**Options:**
- `absorb(varlist)` - Fixed effects to absorb (required)
- `vce(unadjusted|robust|cluster varlist)` - Variance estimator
- Supports `aweight`, `fweight`, and `pweight`
- `verbose`, `timeit` - Show progress/timing

### cqreg - Quantile regression

Drop-in replacement for `qreg` using an Interior Point Method (IPM) solver. Supports high-dimensional fixed effects via the `absorb()` option.

```stata
* Median regression (default)
cqreg price mpg weight

* 25th percentile
cqreg price mpg weight, quantile(0.25)

* 75th percentile with robust SEs
cqreg price mpg weight, quantile(0.75) vce(robust)

* Quantile regression with fixed effects
webuse nlswork, clear
cqreg ln_wage age ttl_exp tenure, absorb(idcode)

* Two-way FE with clustering
cqreg ln_wage age ttl_exp tenure, absorb(idcode year) vce(cluster idcode)
```

**Options:**
- `quantile(#)` - Quantile to estimate (0-1, default: 0.5)
- `absorb(varlist)` - Fixed effects to absorb
- `vce(iid|robust|cluster varname)` - Variance estimator
- `bwmethod(hsheather|bofinger|chamberlain)` - Bandwidth for sparsity
- `tolerance(#)` - Convergence tolerance (default: 1e-8)
- `maxiter(#)` - Maximum IPM iterations (default: 50)
- `verbose`, `timeit` - Show progress/timing

## Building from Source

You almost certainly do not need to compile this program from the source files, as pre-compiled plugins are automatically built by GitHub and installeed when you follow the installation instructions above. If you want to, then it is highly recommended that your compiler work with OpenMP. The Makefile should automatically detect your platforms and compile with OpenMP when available.

### Requirements

| Platform | Requirements |
| --- | --- |
| macOS | Xcode Command Line Tools, Homebrew with `libomp` |
| Windows | MSVC/MinGW, or cross-compile with `llvm-mingw` |
| Linux | GCC with OpenMP (`libgomp`) |

### Build Commands

```bash
make              # Build for current platform
make all          # Build all platform plugins
make check        # Check build dependencies
make clean        # Remove compiled files
```

### Output Files

| File | Platform |
| --- | --- |
| `build/ctools_mac_arm.plugin` | macOS Apple Silicon |
| `build/ctools_mac_x86.plugin` | macOS Intel |
| `build/ctools_windows.plugin` | Windows x64 |
| `build/ctools_linux.plugin` | Linux x64 |

## Project Structure

```
ctools/
├── src/                        # C source files
│   ├── stplugin.c/h            # Stata Plugin Interface
│   ├── ctools_plugin.c         # Main dispatcher
│   ├── ctools_data_load.c      # Parallel data loading
│   ├── ctools_sort_radix_lsd.c # Parallel LSD radix sort
│   ├── csort/                  # Sort command
│   ├── cmerge/                 # Merge command
│   ├── cimport/                # CSV import
│   ├── cexport/                # CSV export
│   ├── creghdfe/               # HDFE regression
│   └── cqreg/                  # Quantile regression (IPM solver)
├── build/                      # Compiled plugins, .ado, .sthlp files
├── benchmark/                  # Performance test .do files
├── Makefile                    # Build system
└── README.md
```

## Performance Notes

- All commands use OpenMP for parallel computation when available
- Data loading uses 8-way loop unrolling for cache efficiency
- Sorting uses a parallel LSD radix sort optimized for Stata's data types
- `creghdfe` uses a conjugate gradient solver with Kaczmarz-style sweeping
- `cqreg` uses a primal-dual Interior Point Method with Mehrotra predictor-corrector

## Compatibility

- Stata 14.0 or later
- Windows (x64), macOS (Intel and Apple Silicon), Linux (x64)
- OpenMP support recommended but not required

## License

MIT License

## Contributing

Contributions are welcome. Please open an issue or submit a pull request.

## Citation

If you use ctools in your research, please cite:

```bibtex
@software{ctools,
  title = {ctools: High-performance C-accelerated commands for Stata},
  url = {https://github.com/mdroste/stata-ctools}
}
```
