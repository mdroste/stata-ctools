# Manual diagnostic probes

These retained programs were previously named `validation/temp_*.do`. They are
historical debugging aids, not release tests, and may rely on interactive state
or obsolete paths. Useful cases have separate named regressions in the suite
registry. Run probes only after checking their paths and inputs, and route new
logs and temporary files to `temp/<task>/`.

All Stata execution uses the `stata` alias through `/bin/zsh -lic`; batch drivers
must end with `exit, clear`. Do not invoke an application binary directly.
