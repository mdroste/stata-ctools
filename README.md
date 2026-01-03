# stata-ctools
C-accelerated drop-in replacements for Stata programs

| Replaces (Stata command) | Use instead (stata-ctools) | Typical speedup |
| ------------------------ | -------------------------- | --------------: |
| `sort`                   | `csort`                    |         **10×** |
| `import delimited`       | `cimport delimited`        |         **10×** |
| `export delimited`       | `cexport delimited`        |         **10×** |
| `merge`                  | `cmerge`                   |         **10×** |
| `binscatter`             | `cbinscatter`              |         **10×** |
| `reghdfe`                | `creghdfe`                 |         **10×** |
