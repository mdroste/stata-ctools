# stata-ctools
C-accelerated drop-in replacements for Stata programs

## About

ctools provides drop-in replacements for a variety of Stata programs that run substantially faster for large datasets.
.

| Stata command         | ctools command       | Purpose                        | Typical speedup  |
| --------------------- | -------------------- | ------------------------------ | ---------------: |
| `sort`                | `csort`              | Sort dataset                   | **2-4×**         |
| `import delimited`    | `cimport delimited`  | Import text-delimited          | **50×**          |
| `export delimited`    | `cexport delimited`  | Export text-delimited          | **25×**          |
| `merge`               | `cmerge`             | Merge (join) datsets           | **2-5×**         |
| `binscatter`          | `cbinscatter`        | Binned scatter plots           | **10×**          |
| `reghdfe`             | `creghdfe`           | High-dimensional FE regression | **10×**          |
| `qreg`                | `cqreg`              | Quantile regression            | **4×**          |

## Quick start

You can begin using these programs by installing directly from this Git repository with the Stata command:
```stata
net install ctools
```

Running this command places all of the relevant Stata files (.ado and .sthlp), along with precompiled C plugins for Windows, Mac OSX, and Unix directly in your ADO path. Stata will automatically detect your operating system and use the correct precompiled plugin. 

If for whatever reason you would like to compile these plugins from scratch, please refer to the BUILD readme. 

Alternatively, you can install by downloading this repository and moving the contents of the ./build directory somewhere on Stata's adopath.


