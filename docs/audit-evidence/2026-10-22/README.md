# October 22 audit evidence

These files support the itemized repository audit. All data are synthetic or existing repository fixtures. Production files were not edited by the audit.

- `estimation/`: IV, PPML and quantile regression probes and logs. `cqreg_single` records an excluded candidate: native Stata also permits one-cluster missing covariance.
- `data/`: conversion/sampling/merge and binscatter probes and logs. `bins_probes` intercepts graph construction to inspect its syntax.
- `io/`: fixture generator, native parser harness and compiler driver, and licensed import/export probe/log.
- `shared/`: range statistics, matching option/input tests, matching allocation-failure harness, and baseline validation drivers/logs.

Some source files intentionally refer to `/tmp/ctools-*` working directories and the original repository path. Adjust those paths and the candidate plugin location before rerunning on another checkout. The initial shared root probe includes captured unavailable-reference errors; these are not treated as reference comparisons. `range64.txt` is the end-to-end optimized-path reproducer; the four-row case in `root_probes.txt` is a passing direct-path control.

From the repository root on macOS, the two standalone native probes can be reproduced with:

```sh
clang -std=c11 -O2 -D_DARWIN_C_SOURCE -DSYSTEM=APPLEMAC -Isrc -Wl,-dead_strip docs/audit-evidence/2026-10-22/shared/range_cancel.c -o /tmp/ctools-range-probe
/tmp/ctools-range-probe
clang -std=c11 -O1 -g -fsanitize=address,undefined -D_DARWIN_C_SOURCE -DSYSTEM=APPLEMAC -Isrc -Wl,-dead_strip docs/audit-evidence/2026-10-22/shared/match_oom.c -o /tmp/ctools-match-oom-probe
/tmp/ctools-match-oom-probe
```

Expected matching output before a fix:

```text
fail_growth=0 rc=0 matched=1 treated_support=1
fail_growth=1 rc=0 matched=0 treated_support=0
```

The second line is the defect: allocation failure is reported as successful no-match output. The mock only replaces transport and selected allocation calls; matching logic is the actual production implementation.

The native baseline log records 37 attempts, 26 successful suites and 11 failures caused by the absent OpenMP runtime/header. The two component logs record 2,219 successful assertions in 12 completed components, plus rangejoin blocked by its missing reference. Do not interpret setup-blocked drivers as successful tests.
