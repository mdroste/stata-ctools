# Row tiles in the shared loader and bulk store

Status: implemented in `ctools_data_load_ex()` (`src/ctools_data_io.c`) and
`ctools_io_data_store_impl()` (`src/ctools_data_store.c`). Every command that
reads or writes through `ctools_data_load()`, `ctools_data_load_ex()` or
`ctools_data_store*()` uses them; no command passes a flag for it. All
measurements below come from StataNow 19.5 MP on a Ryzen 9 9950X (16 cores,
32 threads, 94 GB RAM, two channels of DDR5-4800) under Windows 11 on
2026-11-27, with the plugin built by MSYS2 GCC 16. macOS and Linux have not
been measured; [Other machines](#other-machines) says what to expect and how
to measure.

## Why

Stata keeps the values of one observation next to each other in memory. A
plugin reads them through the SPI one value at a time (`SF_vdata()`), so the
order of the reads decides how often each observation's memory is fetched:

- column by column (the loader's chunked scheduler, one variable per task),
  every observation's cache lines are fetched once per variable;
- in row tiles (`ctools_io_load_tiled_variables()`: 512 observations, all
  variables), they are fetched once for all the variables of the load.

Writes behave the same way, with a read of each line before it is written.
The values loaded or stored do not depend on the order; only the time does.
The persistent thread pool that the chunked scheduler uses is also rebuilt in
every plugin call, which the tiles' OpenMP loop does not need.

## Rules

| Path | Rows tiled when | Before 2026-11-27 |
|---|---|---|
| Load, numeric variables only | 2 or more variables (`IO_AUTO_TILE_MIN_VARS`) and 200,000 or more selected rows (`IO_AUTO_TILE_MIN_ROWS`) | only with `CTOOLS_LOAD_ROW_TILES`, 128+ variables, or 256+ bytes of string columns |
| Load with a string variable | as before: `CTOOLS_LOAD_ROW_TILES`, or 256+ bytes of string columns | the same |
| Bulk store (`ctools_data_store*`) | 2 to 4 variables with enough threads, or 4+ tileable variables (no strL) and 200,000+ rows | 2 to 4 variables, or 256+ bytes of string columns |

One variable is loaded or stored row-parallel in either case. A load the
caller marks with `CTOOLS_LOAD_ROW_TILES` is tiled from two variables and
200,000 selected rows on, strings included (not strL).

`IO_AUTO_TILES` (`src/ctools_transport_internal.h`, default 1) switches the
two defaults off when it is defined as 0. It exists for A/B measurements on
other machines (see below); the native test
`test_transport_adaptive_native` fails when it is 0.

## Measurements on Windows

### The loader and store alone

A benchmark-only subcommand loaded the plugin's variables once per plugin
call, as a command does, alternating the two routes (four calls per
repetition, ABBA). Every loaded value was hashed and compared. Median
milliseconds:

| Load | Column route | Row tiles | Ratio |
|---|---:|---:|---:|
| 10M rows, 2 doubles of 64 in memory | 48.7 | 34.6 | 0.71 |
| 10M rows, 4 doubles | 95.0 | 50.5 | 0.53 |
| 10M rows, 8 doubles | 162.9 | 81.8 | 0.50 |
| 10M rows, 12 doubles | 250.0 | 113.5 | 0.45 |
| 10M rows, 24 doubles | 450.3 | 208.9 | 0.46 |
| 10M rows, 12 doubles, nothing else in memory | 244.1 | 97.7 | 0.40 |
| 10M rows, 6 of byte/int/long/float/double | 129.6 | 52.9 | 0.41 |
| 10M rows, 3 str12/str24 | 131.2 | 130.2 | 0.99 |
| 1M rows, 2 to 24 doubles | 4.4 to 30.9 | 3.1 to 22.0 | 0.37 to 0.81 |
| 300K rows, 2 to 24 doubles | 1.1 to 9.8 | 0.9 to 8.2 | 0.48 to 0.98 |

The rows above load every row or select them with a sample marker. A marker
selecting half of the rows gave ratios of 0.52 to 0.93. With an `if`
condition both routes spend most of their time in Stata's per-observation
`if` evaluation (340 ms for 10M rows), so the ratio there is 0.76 to 1.03.
String loads, which keep their old rule, were at 0.94 to 1.09. With 1, 2 or
4 threads the ratio was 0.53 to 0.87.

| Bulk store, 10M rows | Column chunks | Row tiles | Ratio |
|---|---:|---:|---:|
| 5 doubles | 156.6 | 53.1 | 0.34 |
| 8 doubles | 235.6 | 70.9 | 0.30 |
| 12 doubles | 344.5 | 93.1 | 0.27 |
| 24 doubles | 616.1 | 164.1 | 0.27 |
| 6 of mixed storage types | 179.9 | 77.6 | 0.43 |
| 4 doubles and 2 str12 | 192.2 | 91.7 | 0.48 |

At 1M rows the store ratios were 0.35 to 0.56. Two to four variables already
used row tiles. The dataset signature was the same before and after every
store.

### Commands

Wall time of the command, 10M rows, the median of 6 runs in 3 interleaved
rounds, build of 2026-11-27 without and with the two rules. Every result
matched bit for bit (e(b), e(V), generated variables, sort order).

| Command | Without (s) | With (s) | Change |
|---|---:|---:|---:|
| `civreghdfe y x1 (xe = z1 z2), absorb(firm grp)` | 0.893 | 0.776 | -13% |
| `civreghdfe y x1-x6 (xe = z1 z2), absorb(id year) vce(cluster id)` | 4.561 | 4.147 | -9% |
| `cppmlhdfe yc x1 x2, absorb(firm grp)` | 1.973 | 1.923 | -3% |
| `cbinscatter y x1, controls(x2 x3 x4) absorb(firm grp)` | 1.026 | 0.917 | -11% |
| `chdfe x1-x5, absorb(id year) generate(hd_)` | 0.569 | 0.508 | -11% |
| `cegen double sx = total(x2), by(grp year)` | 0.160 | 0.138 | -14% |
| `cegen double mx1 = mean(x1), by(firm)` | 0.112 | 0.101 | -10% |
| `cwinsor x1-x6, cuts(1 99) prefix(w_)` | 0.839 | 0.768 | -9% |
| `cdistinct id year grp` | 0.904 | 0.814 | -10% |
| `ccollapse (mean) x1-x8, by(id)` (1M rows out) | 0.204 | 0.186 | -9% |
| `cmerge m:1 firm using ...` | 0.317 | 0.307 | -3% |
| `csort rid1 x1`, `csample 50`, `cbsample`, `creghdfe` | | | 0% |

csort, creghdfe and most of cmerge already read and wrote through row tiles
or their own fused paths; csample and cbsample wrote only their sample flag
in these runs.

## Other machines

The cause, Stata's observation-major layout read one value per SPI call, is
the same on every platform, so the direction of the result should hold on
macOS and Linux. The size need not:

- Apple silicon has 128-byte cache lines, more memory bandwidth per core
  and larger shared caches than this PC. Fetching an observation's memory
  again for every variable costs less there, so the ratio may be closer
  to 1.
- When the time goes to the SPI calls themselves rather than to memory, the
  order matters less. This is the case with an `if` condition, and may be
  the case on a machine with fast memory.
- macOS builds use libomp, and the chunked scheduler's pthread pool behaves
  differently from Windows' threads.

Nothing here has been measured on a Mac. `validation/bench_transport.py`,
written for macOS, compiles the loader and store of two source trees into
benchmark plugins and times `ctools_data_load()` and `ctools_data_store()` in
Stata. To compare the two rules on a Mac, from the repository root:

```sh
# a copy of the sources with the old rules
rm -rf /tmp/rowtiles-old && cp -R src /tmp/rowtiles-old
printf '#define IO_AUTO_TILES 0\n' | cat - src/ctools_transport_internal.h \
    > /tmp/rowtiles-old/ctools_transport_internal.h

# build both and write the Stata do-file (Homebrew llvm and libomp, or set
# CC and LIBOMP_PREFIX)
python3 -m validation.bench transport \
    --variant baseline=/tmp/rowtiles-old --variant candidate=src \
    --output /tmp/rowtiles/run1 --rows 10000000 --columns 12 \
    --shapes numeric mixed --widths 8 32 --modes identity filtered \
    --hints 1 --repetitions 5 --threads 8

# run it (the path to Stata varies)
/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -b do /tmp/rowtiles/run1/run.do

# medians and paired speedups (baseline / candidate)
python3 -m validation.bench summarize-transport /tmp/rowtiles /tmp/rowtiles/rowtiles
```

`/tmp/rowtiles/rowtiles_summary.csv` then holds `load_speedup` and
`store_speedup` for each shape and mode; above 1 means the row tiles are
faster. The `identity` mode loads every row and stores through the bulk
store; `filtered` loads with an `if` condition and stores one variable at a
time, so its store is the same in both builds. The `numeric` shape cycles
through byte, int, long, float and double with missing values; `mixed`
alternates them with string variables, which keep the old load rule.

## Tests

`validation/native/transport.suite.c` (`test_transport_adaptive_native`) loads
three numeric variables with exactly 200,000 selected rows, filtered and not,
and checks that the tiled allocation is reached (fault injection makes it fail
and the column route take over without changing a value), and two variables
one row short of the gate, which must not reach it. Its other cases cover the
older gates, allocation retries and errors. On Windows the native transport
suites do not build through `validation/run.py` (it compiles them as Unix
code); they were built by hand with `-DSYSTEM=STWIN32` for these checks.
