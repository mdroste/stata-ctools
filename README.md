# stata-ctools
C-accelerated drop-in replacements for Stata programs

## About

ctools provides drop-in replacements for a variety of Stata programs that run substantially faster for large datasets.

| Replaces (Stata command) | Use instead (stata-ctools) | Typical speedup |
| ------------------------ | -------------------------- | --------------: |
| `sort`                   | `csort`                    |         **2×**  |
| `import delimited`       | `cimport delimited`        |         **40×** |
| `export delimited`       | `cexport delimited`        |         **20×** |
| `merge`                  | `cmerge`                   |         **4×**  |
| `binscatter`             | `cbinscatter`              |         **10×** |
| `reghdfe`                | `creghdfe`                 |         **15×** |

## Quick start

You can begin using these programs by installing directly from this Git repository with the Stata command:
```stata
net install ctools
```

Running this command places all of the relevant Stata files (.ado and .sthlp), along with precompiled C plugins for Windows, Mac OSX, and Unix directly in your ADO path. Stata will automatically detect your operating system and use the correct precompiled plugin. If for whatever reason you would like to compile these plugins from scratch, please refer to the BUILD readme. 

Alternatively, you can install by downloading this repository and movin the contents of the ./build directory somewhere on Stata's adopath.

