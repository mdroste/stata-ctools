# Interpolation, interval joins, and string splitting: native/SSC comparisons

These measurements compare `cipolate` with native `ipolate`, `csplit` with native
`split`, and `crangejoin` with SSC `rangejoin`. They support the approximate
speedups in the [README](../README.md). The separate
[optimization report](PERFORMANCE_NEWCOMMANDS.md) compares earlier and later C
implementations and is not the baseline for these speedups.

## Results

Times are end-to-end seconds, reported as the median of four measured calls.
Speedup is the reference median divided by the ctools median. Each scale has
one full-size warmup per implementation before the measured calls.

| Command | Input rows per dataset | Reference (s) | ctools (s) | Speedup |
|---|---:|---:|---:|---:|
| `ipolate` | 1,000,000 | 0.860 | 0.053 | 16.2× |
| `ipolate` | 5,000,000 | 4.604 | 0.246 | 18.7× |
| `ipolate` | 20,000,000 | 20.590 | 1.127 | 18.3× |
| `split` | 1,000,000 | 3.559 | 0.169 | 21.0× |
| `split` | 5,000,000 | 14.223 | 0.613 | 23.2× |
| `split` | 20,000,000 | 62.777 | 3.116 | 20.1× |
| `rangejoin` | 1,000,000 | 2.857 | 0.186 | 15.4× |
| `rangejoin` | 5,000,000 | 13.949 | 0.764 | 18.3× |
| `rangejoin` | 20,000,000 | 57.027 | 2.849 | 20.0× |

For interpolation and splitting, input and output row counts are equal. The
join has the stated number of observations **in each input**, with 1,998,000,
9,990,000, and 39,960,000 output rows respectively. The README rounds the range
of speedups across these three scales, not a confidence interval.

These are synthetic workloads modeled on firm panels, transaction records, and
daily offers. They use relatively narrow datasets and fixed-width strings.
The join inputs are already ordered by firm and date; the interpolation input
is interleaved by firm. Different sorting requirements, group structures,
delimiters, string widths, join expansion, memory pressure, and hardware can
change the results substantially. In particular, these results do not measure
`csplit`'s native `strL` fallback, optional `destring`, wide joins, or dense
many-to-many matches. Measurements came from a shared workstation; all measured
calls, including slower calls, are retained.

## Workloads and exact commands

The [generator and benchmark loop](../validation/benchmark_readme_commands.do)
define the data completely. Random and sort seeds are 95183 at every scale.

| Workload | Input structure | Options and output |
|---|---|---|
| Firm-month sales | N rows, 200 months per firm; 5,000 / 25,000 / 100,000 firms. Seven variables: long row ID and firm, int month and employment, byte sector, str2 region, double sales. Time-major monthly extracts interleave firms. Sales have a firm-size component, trend, lognormal noise, and independently missing values with probability 35%. | Separate linear interpolation by firm with endpoint extrapolation; one new double variable. |
| Transaction records | N rows, four variables: long row ID, str64 composite code, int quantity, double amount. Codes vary by row and contain region, sales channel, SKU, account, currency, and date. | Six fields separated by a literal pipe; six new strings. No numeric conversion. |
| Inquiries and offers | N master rows and N using rows, each with five variables. There are 500 dates per firm (2,000 / 10,000 / 40,000 firms). Both inputs are ordered by firm and date. Master has long row ID and firm, double date, int quantity, str2 region. Using has long offer ID and firm, int date, double price, str3 venue. Master dates lie halfway between successive daily offers. | Inclusive ±1-day join within firm. Two offers per inquiry except the final inquiry per firm, which has one. Nine output variables. |

The timed commands are:

```stata
ipolate sales month, generate(sales_i) by(firm) epolate
cipolate sales month, generate(sales_i) by(firm) epolate threads(12)

split code, generate(part) parse("|")
csplit code, generate(part) parse("|") threads(12)

rangejoin date -1 1 using offers, by(firm) keepusing(offer_id price venue)
crangejoin date -1 1 using offers, by(firm) keepusing(offer_id price venue) threads(12)
```

## Machine, software, and measurement

- Measured September 25, 2026 on an Apple M4 Pro, 12 logical processors, 48 GiB
  RAM, macOS 26.7 (25G229).
- StataNow/MP 18.5 for Apple Silicon, revision 26 February 2025, with
  `set processors 12`. ctools also uses an explicit `threads(12)` cap.
- Native `ipolate` 1.3.5 (23 April 2020), native `split` 2.1.0
  (2 November 2021), SSC `rangejoin` 1.1.3 (13 April 2021) and its dependency
  `rangestat` 1.1.1 (9 May 2017). These dependencies are needed only for the
  reference benchmark, not for `crangejoin` itself.
- ctools 1.0.2, local frozen source identity `readme-bench-20260925`, built with
  Homebrew clang 21.1.8, `-O3`, LTO, strict floating-point flags, and a macOS
  11-compatible static OpenMP runtime. The actual plugin SHA-256 is
  `d57175ad005b20bfd0a10d9a99d9e5f41dd51be6bad11ca0a7439f89a2c38aed`.
  This is a frozen working-tree build, not an assertion that a release tag or
  the moving checkout contains exactly the same code. The
  [manifest](benchmarks/readme_commands_manifest.json) hashes the sources,
  wrappers, plugin, reference commands, runtime archive, and benchmark harness.
- One Stata session per scale, run sequentially. Within each workload,
  reference and ctools run in alternating AB, BA, AB, BA pairs, after warmups.
  No observations are dropped from the timing summary.
- Timing surrounds the complete command: ado work, C loading/computation/
  storage, allocations, cleanup, and the join's using-file read. Generating
  data, reloading the starting dataset between calls, sorting solely for
  comparisons, and checking results are outside the timed region. Inputs and
  using files have been read previously; this is a warm-cache comparison.
- A [benchmark-only clock plugin](../validation/benchmark_clock.c) uses
  `CLOCK_MONOTONIC`. Earlier runs used a clock-changing wrapper, so the
  monotonic timer keeps the reported timings comparable across those runs.

## Correctness and evidence

At every scale, warmup results are compared across **all output cells** using
Stata's `cf _all`. Interpolation and join outputs are first put in common row
order outside the timed region. The benchmark also asserts the expected row
count after every call. All nine full-output comparisons passed exactly.
These checks establish equality for the measured
workloads; they do not replace the commands' broader correctness suites.

- [Every measured call (CSV)](benchmarks/readme_commands_timings.csv)
- [Exact-comparison and completion markers, plus raw timings](benchmarks/readme_commands_evidence.txt)
- [Build, source, reference, and harness identities](benchmarks/readme_commands_manifest.json)

## Reproduction

Build an isolated copy of the source and ado files using the
[platform build instructions](PLATFORMS.md), keeping the compiled plugin in its
`build/` directory. Install SSC `rangejoin` and `rangestat` for the reference
side, and provide the directory containing `rangejoin.ado` to the preparer.
On this machine the reference copy and installed dependency were supplied as:

```sh
python3 validation/prepare_readme_benchmarks.py \
  --snapshot /private/tmp/ctools-readme-bench/frozen \
  --reference /private/tmp/ctools-rangejoin-reference \
  --reference-dependency "/Users/Mike/Library/Application Support/Stata/ado/plus/r/rangestat.ado" \
  --reference-dependency "/Users/Mike/Library/Application Support/Stata/ado/plus/r/rangestat_run.ado" \
  --output /private/tmp/ctools-readme-reproduction
```

Preparation compiles the clock and writes three driver files; it never launches
Stata. Defaults are `--sizes 1000000 5000000 20000000 --threads 12 --reps 4`.
Run the printed `stata` commands sequentially. Each driver captures errors and exits cleanly.
Then collect the results:

```sh
python3 validation/summarize_readme_benchmarks.py /private/tmp/ctools-readme-reproduction \
  --csv /private/tmp/ctools-readme-reproduction/timings.csv \
  --evidence /private/tmp/ctools-readme-reproduction/evidence.txt
```

The summarizer rejects incomplete or failed drivers, missing exact comparisons,
missing or duplicate trials, wrong row counts, and invalid elapsed times.
The frozen build and original local logs for this run remain under
`/private/tmp/ctools-readme-bench`; the compact public evidence omits Stata's
license details. Temporary benchmark datasets are erased after use.
