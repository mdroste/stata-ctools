# Diagnostic scripts

These Stata scripts investigate specific command behavior. Run them manually.
The automated tests are listed in `validation/suites.json`.

| Script | What it checks |
| --- | --- |
| `merge.do` | Compares `merge` and `cmerge` on datasets with 25 million observations |
| `nobatch.do` | Runs `csort` with batch data transfers enabled and disabled |
| `which_cppml.do` | Shows which `cppmlhdfe` ado-file Stata loads and prints its opening lines |
| `cppml_pathology_debug.do` | Compares PPML coefficients and covariance matrices across convergence settings |
| `cppml_pathology_minchange.do` | Examines PPML covariance differences as convergence tolerances change |
| `cwinsor_macro_rounding.do` | Checks percentile cutoff rounding against `winsor2` |

Run these scripts from the repository root in a fresh Stata session, with
`build/` on the ado-path and the commands used for comparison installed. Check
each script's inputs and settings first. The merge and sorting scripts create
large datasets; the sorting script also changes the `ctools_nobatch` global.

Save logs and generated output in `temp/<task>/`. For batch runs, use
`validation/support/stata_runner.py` with a script that closes its log and ends
with `exit, clear`. See [Stata setup](../README.md#run-the-tests).

When a script identifies a bug, add a test for it to the appropriate test
directory and list it in `validation/suites.json`.
