# Command compatibility

ctools implements selected interfaces of the listed commands. Similar syntax does not imply identical options, numerical methods, or inference. The installed help files describe accepted syntax. `aw`, `fw`, `pw`, and `iw` below denote Stata weight types; a dash means no Stata weight expression is accepted.

| Command | Reference | Weights | Supported scope and material differences |
|---|---|---|---|
| `cimport` | `import` file-format subcommands | — | C readers for delimited, XLS/XLSX, dBase, SAS/catalogs, XPORT5/8, SAV/ZSAV and SHP. Long delimited/Excel/statistical text and binary SHP headers use strL. Unicode/code-page encodings, delimiter sets and numeric locales are supported. FRED/Haver are excluded; see [the I/O parity audit](IO_PARITY.md). |
| `cexport` | `export` file-format subcommands | — | C writers for delimited, XLS/XLSX, dBase, SAS, XPORT5/8, SAV and SHP. Staged file publication and tested format-specific metadata/date/missing handling. Excel headings are opt-in; XLSX/BIFF8 XLS workbook editing and style preservation are implemented. Native SAS writer parity is unverified here. See [the I/O parity audit](IO_PARITY.md); full native parity is not established. |
| `csort` | `sort` | — | Stable ascending sort; rejects `if`/`in`. `nosortedby` skips the wrapper's final native sort used to establish Stata's sorted metadata. `stream(#)` limits variables loaded together. Fixed strings only. |
| `cmerge` | `merge` | — | Declared unique key sides are checked before modification; failures restore the master. Shared storage types are widened; key string/numeric types must match even with `force`. See help for merge types and conflict options; frames are used when available, with a preserve/restore path on older Stata. |
| `csample` | `sample` | — | Percentage or `count(#)` sampling without replacement; both commands support `by()`. Random selections need not coincide under a common seed. Native count syntax differs. |
| `cbsample` | `bsample` | `weight()` output | Bootstrap observations/clusters, optionally within strata. The positional count specifies draws per stratum, in observations or clusters. `weight()` overwrites an existing numeric variable with frequencies. |
| `cencode` | `encode` | — | Literal value labels; existing label sets use native `encode` to preserve its extension/noextend semantics. Fixed strings and text strL up to 2045 bytes are readable. |
| `cdecode` | `decode` | — | Native Stata decoding engine, including literal/multiline/long labels and labeled missing codes; multiple-variable and transactional replace extensions. `threads()` is a compatibility no-op. |
| `cdestring` | `destring` | — | Numeric conversion of fixed strings and bounded text strL; use the documented `ignore()`, `force`, and output options. |
| `cipolate` | `ipolate` | — | Linear interpolation, duplicate-x averaging, by()/by-prefix, epolate, if/in. Preserves input row order. ctools cipolate is distinct from the SSC cubic interpolator of the same name. |
| `csplit` | `split` | — | Fixed strings use C; strL input uses native split. Supports parse(), notrim, limit(), if/in, and native destring conversion options. Outputs are staged until success. |
| `crangejoin` | `rangejoin` 1.1.3 | — | Inclusive interval joins, by(), keepusing(), prefix()/suffix()/all. Preserves unmatched master rows and stable using-key order; restores data on error. Numeric/fixed string data only; output limited to 2^31−1 observations. |
| `cwinsor` | `winsor2`, `gstats winsor` | — | Winsorization/trimming, groups, and named outputs. Generated values are missing outside if/in; replace leaves those rows unchanged. Failed operations roll back. |
| `cbinscatter` | `binscatter`, selected `binsreg` behavior | aw/fw/pw/iw | Numeric plotted variables; factor controls and expression weights. Redundant controls are omitted; iterative projection sweeps absorb effects. No claim of identical `binsreg` inference. |
| `crangestat` | `rangestat` | — | Implemented range statistics with interval/by/excludeself; not arbitrary user Mata functions. |
| `creghdfe` | `reghdfe` | aw/fw/pw | Linear regression with absorbed effects. Residual/FE predictions require residuals saved at estimation; `xb` excludes absorbed effects. Solver and degrees-of-freedom options are not a complete `reghdfe` interface. |
| `civreghdfe` | `ivreghdfe` | aw/fw/pw | Implemented IV estimators and VCE options. `center` is rejected. Residual predictions require saved estimation residuals; see help for unsupported estimator/VCE combinations. |
| `cqreg` | `qreg` | — | Conditional quantiles with IID, robust, or clustered VCE; no weights or `absorb()`. Explicit factor indicators are allowed. `predict` supports xb/residuals/stdp. The IPM solver and density estimates can differ from native `qreg`; nonconvergence returns an error. |
| `cpsmatch` | `psmatch2` | — | ATT from one outcome; no ATE/ATU or observation-weight syntax. No-replacement matching requires nearest neighbor with one neighbor. `r(att_se)` is an approximate matched-sample SE, not the `psmatch2` estimator. |
| `cpplmhdfe` | `ppmlhdfe` | aw/fw/pw | Native C PPML with optional FE, categorical interactions, heterogeneous slopes, native simplex/ReLU/mu separation, saved individual FE, solver/initialization controls, exposure/offset, and robust or multiway clustered VCE. Reports the constant, z/Wald inference, absorbed DF, and final sample. `d()` enables fitted means and residual predictions. Hidden reference debugging options and every numerical iteration path are not identical; see [PPML details](README_cpplmhdfe.md). The public command is spelled `cpplmhdfe`. |

For `cpsmatch`, the reported SE is `sqrt(s1²/n1 + s0²/n0)`, using matching-weighted sample variances and observation counts in the matched treated and control groups. It does not account for estimated propensity scores or repeated use of a control in the same way as matching-specific variance estimators. Treat comparisons of ATT and comparisons of its SE as separate validation questions.

Shared plugin data transfer supports numeric and fixed-width string variables up to `str2045`. It also reads textual `strL` values up to 2045 bytes through Stata's length-aware API; its general binary/oversized strL reads and strL writes are rejected. The dedicated Excel/statistical/SHP engines have separate length-aware transfer for longer text and binary headers. Delimited import reports failed writes as incomplete data, and shared write failures report that data may be partially modified. Do not continue analysis after such an error without restoring or reloading the data.

## Performance evidence

Runtime depends on rows, columns, variable types, hardware, thread count, options, and the reference command. No universal speedup is promised. Published benchmark results should include the date, CPU/OS, Stata edition/version, reference package versions, ctools source revision, dimensions/types, seed or input data, options, thread counts, and timing method. Keep the correctness comparison and the benchmark log with each result. Existing benchmark scripts in `validation/` provide starting points; a speed claim without this provenance is not a release guarantee.

OLS, IV, and PPML evaluate numeric weight expressions once. Numeric cluster labels
are grouped as full doubles, so fractional, negative, and large labels preserve
their partitions. OLS/PPML return the final retained `e(sample)`; with fweights,
`e(N)` reports the frequency-weight sum.

PPML screens all-zero FE groups and singleton chains, then uses native simplex
and ReLU for general regressor/FE separation; fitted-mean screening is available
with `separation(mu)`. The reduced sample is rebuilt before
estimation and prediction. No runtime fallback to `ppmlhdfe` is used. See the
PPML details for the corrected pathological separation fixture, supported
advanced controls, and remaining hidden-option limits.

Output names for `creghdfe` (residuals, group IDs, and saved effects), `crangestat`,
and generated `cencode` variables must be distinct new names. These commands
restore the dataset on failure. Empty-data encoding follows the same validation
contract. OLS/IV/PPML projection failures and exhausted PPML IRLS fits return 430.

Excel daily dates use their workbook's 1900 or 1904 epoch; serial 60 in the 1900
epoch maps to 28 February 1900 in the public cimport path. Datetime cells import as milliseconds.
Excel does not represent leap seconds, so `%tC` export uses `cofC()` conversion.

A one-cluster covariance is undefined. The shared OLS/IV covariance path returns 498; PPML rejects fewer than two retained clusters with 459. Neither posts a successful zero covariance for this case.
