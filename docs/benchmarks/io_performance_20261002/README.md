# cimport / cexport performance, 2 October 2026

The main gains are in default CSV import of numeric ASCII files and new XLSX
workbook export. CSV export is a control: its implementation was left unchanged.
XLSX import improvements are smaller because decompression and serial metadata
and data transfer still contribute substantially to wall time.

## Changes and evidence

* CSV import retains the mapped input when the detector chooses an ASCII-compatible
  ISO-8859 or KOI8-R charset and the **entire file** is ASCII. Numeric CSV files in
  this campaign are detected as ISO-8859-2. Previously that name missed the fast
  path, spending most of the command in `iconv` and allocating a conversion buffer.
  The detected charset returned to Stata is preserved; explicit encodings and
  files with non-ASCII bytes still use the converter. A regression places a
  non-ASCII byte beyond the detector's 8,000-byte sample and compares decoded
  bytes against the original converter.
* CSV verbose timing now includes charset recognition/conversion in the plugin
  total and displays it separately. Previously it was reported as Stata overhead.
* XLSX export uses libdeflate 1.25 CRC-32 in place of miniz's byte-at-a-time checksum.
  The existing libdeflate subset contained only decompression, despite the old
  checksum comment claiming acceleration. Unmodified upstream CRC sources, tables,
  runtime dispatch, and portable fallback were added at the **same version**.
  Sampling showed 4,323 checksum leaf samples in the baseline XLSX export profile.
  Compression itself remains miniz DEFLATE; XML, values, and compression level are
  unchanged. `libdeflate_crc32` matches the old CRC, including incremental chunks.
* XLSX compression workers are capped by the number of data chunks, excluding
  the tiny preamble/footer. A one-chunk worksheet no longer starts a thread pool
  just to compress its XML serially in one worker. This benefits small workbooks.
* XLSX import recognizes complete compact `<v>…</v></c>` and
  `<is><t>…</t></is></c>` bodies. It performs the same value conversion after
  recognizing the complete envelope. Namespaces, formulas, rich text, attributes,
  extra elements, and all other layouts retain the general XML walk. Every prefix
  of representative compact and general cells is checked under ASan/UBSan.

## Method

Apple M4 Pro, 12 logical CPUs, 48 GiB RAM, macOS ARM64; Stata/MP 19.5 with four
Stata processors. The ctools thread settings are specified separately below. Both frozen
plugins use the same shared-source snapshot, Clang `-O3`, LTO, `-mcpu=apple-m1`,
strict floating point, and static LLVM OpenMP 21.1.8 built for macOS 11.
The repository was also undergoing a structural refactor, so the experiment
froze shared sources before applying only the I/O changes. The complete source
hash inventory is in `source-manifest.json`; only the three I/O C files and the
added CRC implementation differ. Both benchmark builds use identical ado files;
the corrected verbose display is validated separately in the working package.

A separate Stata process loads each frozen build in **B–C–C–B** order. Each
case/thread/operation gets one warmup and three measured iterations per session,
for six observations per build in the main table. A benchmark-only monotonic
clock measures the complete command, including ado work and file commit. Dataset
creation and restoration are outside the timer. These are **warm filesystem-cache**
measurements, not cold-storage throughput or peak-memory measurements. No other
benchmark or test from this task ran alongside the measured sessions.

The reported times are medians. Raw measurements, min/max, logs, and plugin/ado
SHA-256 hashes accompany each table. Small changes in these short runs can be
noise: even the unchanged CSV exporter varies. The large import/export gains
persist at one million observations and in two full campaigns.

Workloads use a fixed seed, Gaussian numeric doubles, approximately 1% numeric
missing values, and where present text containing commas, quotes, ampersands,
and angle brackets:

| Case | Rows | Numeric columns | String columns |
|---|---:|---:|---:|
| numeric | 200,000 | 10 | 0 |
| mixed | 200,000 | 10 | 5 |
| wide | 30,000 | 100 | 0 |
| strings | 50,000 | 10 | 5, doubled text |
| large | 1,000,000 | 10 | 0 |
| small | 1,000 | 10 | 0 |

## Main command timings

Eight threads for delimited import/export and XLSX export. XLSX import uses its
automatic worker count (`threads=0` in the CSV); its ado does not expose a
`threads()` option. The accompanying CSV also contains one-thread timings.

| Case | Operation | Baseline seconds | Candidate seconds | Speedup |
|---|---|---:|---:|---:|
| numeric | export_csv | 0.023330 | 0.024344 | 0.96× |
| numeric | import_csv | 0.444663 | 0.013293 | 33.45× |
| numeric | export_xlsx | 0.281135 | 0.094156 | 2.99× |
| numeric | import_xlsx | 0.122690 | 0.122618 | 1.00× |
| mixed | export_csv | 0.028834 | 0.033185 | 0.87× |
| mixed | import_csv | 0.034356 | 0.035408 | 0.97× |
| mixed | export_xlsx | 0.510707 | 0.148385 | 3.44× |
| mixed | import_xlsx | 0.251043 | 0.244411 | 1.03× |
| wide | export_csv | 0.045692 | 0.044059 | 1.04× |
| wide | import_csv | 0.668891 | 0.021261 | 31.46× |
| wide | export_xlsx | 0.409335 | 0.141459 | 2.89× |
| wide | import_xlsx | 0.185104 | 0.175543 | 1.05× |
| strings | export_csv | 0.013029 | 0.013596 | 0.96× |
| strings | import_csv | 0.014352 | 0.014056 | 1.02× |
| strings | export_xlsx | 0.163838 | 0.044690 | 3.67× |
| strings | import_xlsx | 0.070148 | 0.069512 | 1.01× |

The numeric CSV gains depend on the charset selected by automatic detection.
Mixed/string files already taking the ASCII fast path show no corresponding
large improvement. XLSX import is approximately unchanged on numeric data and
modestly faster on the wide/mixed cases; the result does not support a broad
claim of large Excel-import gains.

## Small files

Forty measured repetitions per build (20 per session), after the worker cap:

| Operation | Threads | Baseline seconds | Candidate seconds | Speedup |
|---|---:|---:|---:|---:|
| export_csv | 1 | 0.001418 | 0.001504 | 0.94× |
| import_csv | 1 | 0.003182 | 0.001143 | 2.78× |
| export_xlsx | 1 | 0.005188 | 0.003970 | 1.31× |
| export_csv | 8 | 0.001380 | 0.001324 | 1.04× |
| import_csv | 8 | 0.003215 | 0.001038 | 3.10× |
| export_xlsx | 8 | 0.006525 | 0.003901 | 1.67× |
| import_xlsx | 0 | 0.003337 | 0.003531 | 0.94× |

## One million rows

Four measured repetitions per build, using the final candidate. These results
retain all default file-commit behavior and use the same 10 numeric columns.

| Operation | Threads | Baseline seconds | Candidate seconds | Speedup |
|---|---:|---:|---:|---:|
| export_csv | 1 | 0.408679 | 0.411869 | 0.99× |
| import_csv | 1 | 2.389927 | 0.254092 | 9.41× |
| export_xlsx | 1 | 3.529942 | 2.557223 | 1.38× |
| export_csv | 8 | 0.079072 | 0.081740 | 0.97× |
| import_csv | 8 | 2.204018 | 0.057554 | 38.30× |
| export_xlsx | 8 | 1.322974 | 0.405027 | 3.27× |
| import_xlsx | 0 | 0.599274 | 0.578195 | 1.04× |

## Native Stata reference and externally written XLSX

Native `import`/`export` commands were also timed on the numeric and mixed
200,000-row data, using their default processor configuration. These are
contextual medians of two measured runs after one warmup, not an ABBA comparison.
The CSV input and ctools-written XLSX input match those in the main campaign.

| Case | Native operation | Seconds |
|---|---|---:|
| numeric | export_csv | 0.785643 |
| numeric | import_csv | 0.330234 |
| numeric | export_xlsx | 5.995009 |
| numeric | import_xlsx | 3.182980 |
| mixed | export_csv | 0.861598 |
| mixed | import_csv | 0.631895 |
| mixed | export_xlsx | 8.067674 |
| mixed | import_xlsx | 6.457406 |

Stata's native exporter supplied an additional XLSX input, so the parser was
also checked on an independent writer. Five measurements after warmup, baseline
then candidate (less controlled than ABBA); all candidate imports match the
baseline with `cf _all`:

| Native-written XLSX | Baseline seconds | Candidate seconds | Speedup |
|---|---:|---:|---:|
| numeric | 0.136773 | 0.134592 | 1.02× |
| mixed | 0.283413 | 0.245866 | 1.15× |

## Reproduction

Keep two frozen build directories. Build both with the same compiler flags and
OpenMP runtime. The public benchmark runner prepares inputs, compiles the clock,
checks build hashes before each session, validates fresh Stata completion nonces,
and compares exported artifacts:

```sh
python3 validation/prepare_io_performance.py \
  --baseline /absolute/path/to/baseline \
  --candidate /absolute/path/to/candidate \
  --output temp/io_benchmark --run
```

Use `--cases large --reps 2` for the large check and `--cases small --reps 20`
for stable small-file timings. Stata is launched through
`/bin/zsh -lic 'stata -q -b do ...'`. Each session exits with `exit, clear`.
The benchmark does not install or alter the user's Stata configuration.

## Correctness and integration

Both frozen builds produce identical CSV bytes and decompressed XLSX parts;
ZIP checksums are verified. Candidate imports match baseline values with
`cf _all`. The main, large, and small campaigns check 48, 12, and 12 export
comparisons respectively, in addition to the Stata value comparisons.

The final working package passes:

* 338 `validate_cimport.do` tests and 213 `validate_cexport.do` tests, with no
  failures or skips; the corrected CSV verbose timing also runs successfully.
* 196 exact Excel option comparisons and 15 formula/calculation-chain structural
  comparisons from `run_io_excel_options.py`.
* New ASan/UBSan regressions for CRC equivalence, bounded compact/general XML
  parsing, and ASCII/non-ASCII charset handling.
* Existing native I/O option, binary-format, parity, XLSX edit, and XLSX long-text
  tests, including their sanitizer checks.
* macOS ARM64 plugin dependency checks, package-source and release-metadata
  synchronization checks. The CRC source also compiles for macOS x86-64;
  runtime performance and accelerated CRC execution were tested on ARM64 only.

The final integration run exposed two protocol strings accidentally renamed by
the concurrent refactor (`scan`/`load` became `cio_scan`/`cio_load`). Restoring the
existing dispatcher protocol fixed the import failures; all the counts above are
from the rebuilt package after that correction. The frozen benchmark builds
already used the correct protocol. Existing validation expectations were not
changed. Logs and hashes of the final plugin, package ado files, and relevant
sources are under `validation/` in this evidence directory.

Legacy `.xls` is covered by correctness tests but was not separately profiled or
benchmarked. Its implementation was not optimized in this pass.

## Remaining costs

CSV export's sampled CPU time is mostly numeric formatting and streamed SPI
reads. This pass preserves its formatting and I/O behavior. XLSX export still
spends substantial CPU time in DEFLATE compression; changing compressors or
compression policy needs separate size/runtime and compatibility measurements.
XLSX import still includes whole-part decompression and serial transfer/type
inference. Broader changes to those stages were not needed for the measured
wins retained here.
