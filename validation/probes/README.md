# Diagnostic probes

These Stata scripts help investigate specific command behavior. Run them manually;
the automated correctness suite is listed in `validation/suites.json`.

| Script | Purpose |
| --- | --- |
| `merge.do` | Compares `merge` and `cmerge` on datasets with 25 million observations |
| `nobatch.do` | Exercises `csort` with batch transfers enabled and disabled |
| `which_cppml.do` | Shows which `cppmlhdfe` ado-file Stata loads and prints its opening lines |
| `cppml_pathology_debug.do` | Compares PPML coefficients and covariance matrices across convergence settings |
| `cppml_pathology_minchange.do` | Examines PPML covariance differences as convergence tolerances change |
| `cwinsor_macro_rounding.do` | Checks percentile cutoff rounding against `winsor2` |

Run probes from the repository root in a fresh Stata session, with `build/` on
the ado-path and the required reference packages installed. Review each script's
inputs and settings before running it. The merge and sorting probes create large
datasets; the sorting probe also changes the `ctools_nobatch` global.

Save logs and generated output in `temp/<task>/`. For batch runs, use
`validation/support/stata_runner.py` with a driver that closes its log and ends
with `exit, clear`. See [Stata setup](../README.md#run-the-tests).

Add repeatable regression tests to the appropriate test directory and register
them in `validation/suites.json`.
