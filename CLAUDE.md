# Project Instructions

## Important Conventions

**Always call Stata from bash using the stata alias.** You do not need to find the Stata installation path. stata is defined as the correct alias. 

**Always use relative paths from the project root.** Never use `cd` to change directories.

```bash
# Good: automatically logs output to benchmark/test_per.log
stata -b do benchmark/test_perf.do

# Bad: no need to change directory
cd benchmark && stata -b do test_perf.do
```

## Building

```bash
make              # Build
make clean        # Remove compiled files
make check        # Check dependencies
```

Output: `build/ctools_mac_arm.plugin`, `build/ctools_mac_x86.plugin`, etc.

The Stata ctools program calls should automatically detect a user's OS/architecture and load the appropriate compiled plugin.

## Project Structure

```
src/                    # C source files
  stplugin.c/h          # Stata Plugin Interface (DO NOT MODIFY)
  ctools_plugin.c       # Main dispatcher
  ctools_data_load.c    # Parallel data loading (8x unrolled)
  ctools_sort_radix_lsd.c  # Parallel LSD radix sort
  csort/                # Sort command
  cmerge/               # Merge command
  cimport/              # CSV import
  cexport/              # CSV export
  creghdfe/             # HDFE regression
  cqreg/                # Quantile regression
build/                  # Compiled plugins, .ado, .sthlp files
benchmark/              # Stata .do files for testing
```

## Included Commands

| Command | Replaces | Description |
|---------|----------|-------------|
| `csort` | `sort` | Parallel LSD radix sort |
| `cmerge` | `merge` | C-accelerated merge |
| `cimport` | `import delimited` | Multi-threaded CSV import |
| `cexport` | `export delimited` | Parallel CSV export |
| `creghdfe` | `reghdfe` | HDFE regression with CG solver |
| `cqreg` | `qreg` | Quantile regression with IPM solver |

### Debug Flags for cqreg

In `cqreg_fn.c`:
- `FN_TIMING 1` - Prints detailed timing breakdown per IPM phase
- `IPM_DEBUG 0` - Verbose iteration logging (generates large output)

In `cqreg_main.c`:
- `MAIN_DEBUG 0` - Entry/exit logging for main functions

In `cqreg_vce.c`:
- `VCE_DEBUG 0` - VCE computation logging

Set to 1 to enable, 0 to disable. Rebuild after changing.

### Performance Notes

- Focus on highly performant, parallelized, and modular commands.
- Use 'common' ctools C programs, like ctools_data_load and ctools_data_store, to work with data wherever possible for individual ctools commands.
- Refer to ctools src files for the common data structures you should work with.
- Never modify the stata function inteface ./src/stplugin.c/h files.
- Optimize memory allocation and access patterns to speed up the programs.
- Use OpenMP and the ctools thread management programs wherever possible.
- ctools commands should support a verbose option reporting timing of the plugin components.

## Adding a New Command

1. Create `src/newcmd/newcmd_impl.c` and `newcmd_impl.h`
2. Implement `ST_retcode newcmd_main(const char *args)`
3. Add dispatch case in `src/ctools_plugin.c`
4. Add sources to `Makefile` (`NEWCMD_SRCS`, `NEWCMD_HEADERS`)
5. Create `build/newcmd.ado` and `build/newcmd.sthlp`

## Debugging

- Add `verbose` option to commands for timing breakdown
- Plugin errors: 198 = syntax, 920 = memory
- Run benchmarks: `stata -b do benchmark/test_breakdown.do`
