# stata-ctools
C-accelerated drop-in replacements for Stata programs

| Replaces (Stata command) | Use instead (stata-ctools) | Typical speedup |
| ------------------------ | -------------------------- | --------------: |
| `sort`                   | `csort`                    |         **10×** |
| `import delimited`       | `cimport delimited`        |         **40×** |
| `export delimited`       | `cexport delimited`        |         **20×** |
| `merge`                  | `cmerge`                   |         **10×** |
| `binscatter`             | `cbinscatter`              |         **10×** |
| `reghdfe`                | `creghdfe`                 |         **10×** |

Behind the scenes, ctools has a variety of utility programs that may be useful for the development of future Stata programs.
