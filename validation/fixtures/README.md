# Offline validation data

`manifest.json` pins 24 official Stata 18 example datasets by source URL and
SHA-256. Prepare the cache once from the repository root:

```sh
python3 scripts/fetch_validation_data.py
python3 scripts/fetch_validation_data.py --check-only
```

Downloads go to the ignored `data/` directory. A changed upstream file fails
checksum verification; a Stata test never downloads a replacement. Review
upstream changes before using the maintainer-only `--record` option to change
the manifest. The binary datasets are cached locally rather than redistributed
as part of the ctools installation package.

`ctools_fixture` in `validation/validate_setup.do` loads these files. Missing
fixtures fail the component. Built-in `sysuse` examples remain local to Stata.
The import suite uses available official examples for generic CSV round trips;
obsolete, unavailable example names are not treated as skipped checks.

Run the complete gate with `python3 validation/run_stata_audit.py`. On this
machine, run that Python command through a scoped shell outside the sandbox;
it invokes Stata only through `/bin/zsh -lic 'oldstata ...'` and always exits
Stata cleanly. The reference packages listed in the driver's preflight must
already be installed. The gate requires each of 20 components to reach its
completion marker, rejects failed checks and unexpected skips, and verifies a
unique driver marker so an old successful log cannot satisfy a new run.

The 48 documented matching standard-error comparisons use different methods
and are recorded as exclusions, not passing checks. `test_release_gate.py`
verifies that aborted components, stale logs, missing summaries, failed checks,
and unexplained skips cannot pass the gate. Publication in CI requires a
successful licensed Stata job configured with `CTOOLS_STATA_RUNNER`; skipping
that job cannot publish or commit release binaries.
