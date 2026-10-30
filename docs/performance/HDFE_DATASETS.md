# creghdfe and reghdfe on the HDFE benchmark collection

Sergio Correia's [HDFE benchmark collection](https://www.scorreia.com/data/hdfe/index.html) is a set of 15 datasets for testing estimators that absorb high-dimensional fixed effects. Each has an outcome `y`, regressors `x1` and `x2`, and two identifiers, `id1` and `id2`, that define a bipartite graph. Several carry further identifiers. The graphs differ in the property that makes absorption easy or hard, how well the two sets of units are connected: they range from a complete bipartite design to a single path through 10,001 nodes.

[`validation/performance/datasets.do`](../../validation/performance/datasets.do) runs `creghdfe` and `reghdfe` on every dataset with one, two, three and four fixed effects, with regrouped identifiers, and with unadjusted, robust and clustered covariances, and compares what the two commands store. This page walks through those tests and what they show. The tables are generated from the run described under [Environment](#environment). The underlying rows are in [`benchmarks/hdfe_datasets_results.csv`](../benchmarks/hdfe_datasets_results.csv) and [`benchmarks/hdfe_datasets_oracle.csv`](../benchmarks/hdfe_datasets_oracle.csv).

The suite is opt-in. It needs about 227 MB of downloaded data and is not part of the offline release gate.

## Results in brief

<!-- BEGIN GENERATED: totals -->
- Checks: 131 passed, 0 failed.
- Regression comparisons: 113 specifications, 113 pass, 0 fail, 0 errors.
- Smallest agreement: coefficients 9.3 digits over every default-settings run; among the passing comparisons, covariance 8.0 digits and fit statistics 8.2 digits.
- Default-settings wall clock over those 113 specifications, one run each: reghdfe 225.6 s, creghdfe 46.5 s.
<!-- END GENERATED: totals -->

- **Samples and singletons.** The estimation sample and the number of singleton observations dropped are identical in every specification, including `directors`, where iterated singleton removal discards 493,782 of 525,012 rows.
- **Coefficients.** `e(b)` agrees in every specification at each command's default settings.
- **Mobility groups.** On all 14 two-way graphs the two commands count the same number of mobility groups, and their `groupvar()` identifiers define the same partition of the rows.
- **Covariances and fit statistics.** These agree in every specification. Robust and clustered covariances are compared on converged solutions, for the reason given in [section 3](#3-robust-and-clustered-covariances).
- **Absorbed degrees of freedom.** `e(df_a)` matches in every specification. That includes the two cases a connected graph cannot show: clustering on an absorbed effect when there are several mobility groups, and a third or later effect that forms several mobility groups with an earlier one. See [Absorbed degrees of freedom](#absorbed-degrees-of-freedom).
- **Exact solutions.** On `credit2` and `soccer` both commands reproduce a direct dummy-variable regression. On the path graph neither command reaches seven digits at default tolerance. `creghdfe` does with `tolerance(1e-10)`.

## Reproducing

```sh
python3 scripts/fetch_validation_data.py --collection hdfe   # once: 15 files, checksum-pinned
python3 validation/bench.py hdfe                      # run the comparisons
python3 validation/bench.py hdfe --write-doc          # also regenerate the tables below
```

The first command downloads the Stata files into the ignored `validation/fixtures/data/hdfe/` and verifies each against the SHA-256 in [`validation/fixtures/hdfe_manifest.json`](../../validation/fixtures/hdfe_manifest.json). The suite never downloads anything itself. It needs `reghdfe` and `ftools` installed, uses the plugin in `build/` unless `--build-dir` names another build, and took about eleven minutes on the machine below. The runner exits with an error when any comparison fails.

## Environment

<!-- BEGIN GENERATED: environment -->
- Run date: 5 Oct 2026
- Machine: Apple M4 Pro, 12 logical cores, Darwin 25.6.0 (arm64)
- Stata: 19.5 MP (born 12 Aug 2026), licensed for 4 processors
- Reference: reghdfe `version 6.13.1 10Jan2026`, ftools `version 2.50.0 09jan2026`
- ctools: creghdfe `version 1.0.2 9feb2026 github.com/mdroste/stata-ctools`, source revision a6eb6d3 with uncommitted changes
- Plugin: `build/ctools_mac_arm.plugin`, SHA-256 `516583390bd32025505e03dc964f1a661ab428f0ceff0f30d6df7f5c85d158f1`, default thread count
<!-- END GENERATED: environment -->

## The datasets

Figures are those published with the collection. Components are connected components of the full `(id1, id2)` graph. Both commands drop singleton observations before estimating, so the mobility-group counts in the tables below are smaller.

| Dataset | Kind | Rows | Unique edges | `id1` levels | `id2` levels | Components | Other identifiers |
|---|---|---:|---:|---:|---:|---:|---|
| `credit2` | borrower by month, subsample of `credit` | 19,094 | 19,094 | 601 | 140 | 1 | `ruc`, `t` |
| `credit` | borrower by month | 516,810 | 516,810 | 13,946 | 140 | 1 | `ruc`, `t` |
| `soccer` | home by away player | 73,487 | 73,355 | 528 | 521 | 1 | |
| `synthetic-complete` | complete bipartite | 500,000 | 500,000 | 1,000 | 500 | 1 | |
| `synthetic-uniform-easy` | uniform random matching | 500,000 | 475,566 | 216,007 | 20 | 1 | |
| `synthetic-uniform-hard` | uniform random matching | 500,000 | 499,992 | 122,702 | 122,655 | 194 | |
| `synthetic-uniform-harder` | uniform random matching | 500,000 | 499,997 | 216,007 | 216,046 | 13,213 | |
| `synthetic-assortative` | CEO by firm, assortative | 499,155 | 196,245 | 126,101 | 94,933 | 33,042 | `year`, `ceo_type`, `firm_type` |
| `synthetic-zigzag` | path graph | 10,002 | 10,000 | 5,000 | 5,001 | 1 | |
| `enron` | email sender by recipient | 367,662 | 367,662 | 36,692 | 36,692 | 1,887 | |
| `github` | open-source contributions | 548,843 | 41,751 | 27,970 | 24,358 | 13,232 | `id3`, `t` |
| `patents` | patent citation pairs | 500,008 | 500,008 | 101,837 | 362,995 | 25,229 | |
| `workers` | firm by worker | 504,315 | 253,929 | 28,864 | 218,390 | 19,994 | `id3` (job title), `id4` (year) |
| `schools` | student by teacher | 413,444 | 406,998 | 206,722 | 11,960 | 15 | `id3` (year) |
| `directors` | firm-director matches | 525,012 | 524,819 | 395,925 | 400,850 | 281,627 | |

In the credit files `ruc` is the raw 11-digit borrower number behind `id1`, and `t` is the month behind `id2`, stored as a string. In `github`, `t` is a commit timestamp string, and the file also has two more regressors, `len_commit_msg` and `num_files`.

## How a comparison is scored

Each specification is one test. Both commands run on the same data in memory, and the test passes only if all of the following hold.

- **Exact.** `e(N)`, `e(num_singletons)`, `e(df_a)`, `e(df_r)`, `e(df_m)`, `e(rank)`, `e(N_hdfe)` and `e(N_clust)` are equal. For a two-way model without clustering, the number of mobility groups is equal too.
- **Seven significant digits**, the project default, where digits of agreement are `-log10(|a - b| / max(|a|, |b|))`. This applies to `e(r2)`, `e(r2_a)`, `e(r2_within)`, `e(rss)`, `e(tss)`, `e(mss)`, `e(F)`, `e(rmse)`, and every element of `e(b)` and `e(V)`.

The "Digits" columns report the smallest agreement over the elements of `e(b)`, of `e(V)`, and over the fit statistics, capped at 15. A covariance column therefore reflects its worst element.

Unless stated otherwise both commands use their defaults: `tolerance(1e-8)`, `reghdfe`'s accelerated alternating projections and `creghdfe`'s conjugate-gradient solver. Timings are one wall-clock run of each full command, `reghdfe` first, on data already in memory. They are recorded for context and never asserted. See [Timing notes](#timing-notes).

## Walkthrough

### 1. One fixed effect

`absorb(id1)` and `absorb(id2)`, on all 15 datasets. With one effect the projection is a group mean and nothing iterates. These rows test the parts around the solver: a single pass of singleton removal, the level counts that become absorbed degrees of freedom, and the regression arithmetic. `patents` with `absorb(id2)` drops more than half its rows, because most `id2` levels occur once.

<!-- BEGIN GENERATED: oneway -->
| Dataset | Specification | N used | Singletons dropped | Absorbed df | Digits, b | Digits, V | Digits, fit | reghdfe (s) | creghdfe (s) | Result |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| credit2 | `absorb(id1)` | 19,075 | 19 | 582 | 13.7 | 14.3 | 10.8 | 0.13 | 0.02 | pass |
| credit2 | `absorb(id2)` | 19,094 | 0 | 140 | 14.2 | 13.8 | 10.7 | 0.02 | 0.02 | pass |
| credit | `absorb(id1)` | 516,154 | 656 | 13,290 | 13.6 | 13.5 | 10.2 | 0.08 | 0.04 | pass |
| credit | `absorb(id2)` | 516,810 | 0 | 140 | 13.6 | 13.2 | 9.8 | 0.12 | 0.03 | pass |
| soccer | `absorb(id1)` | 73,487 | 0 | 528 | 13.8 | 14.1 | 9.4 | 0.02 | 0.01 | pass |
| soccer | `absorb(id2)` | 73,487 | 0 | 521 | 14.1 | 13.5 | 9.4 | 0.02 | 0.01 | pass |
| synthetic-complete | `absorb(id1)` | 500,000 | 0 | 1,000 | 13.8 | 14.3 | 10.4 | 0.07 | 0.03 | pass |
| synthetic-complete | `absorb(id2)` | 500,000 | 0 | 500 | 13.2 | 13.8 | 10.9 | 0.11 | 0.03 | pass |
| synthetic-uniform-easy | `absorb(id1)` | 432,350 | 67,650 | 148,357 | 13.3 | 13.6 | 9.4 | 0.08 | 0.04 | pass |
| synthetic-uniform-easy | `absorb(id2)` | 500,000 | 0 | 20 | 13.5 | 13.6 | 9.6 | 0.11 | 0.03 | pass |
| synthetic-uniform-hard | `absorb(id1)` | 490,795 | 9,205 | 113,497 | 13.6 | 13.5 | 9.8 | 0.08 | 0.04 | pass |
| synthetic-uniform-hard | `absorb(id2)` | 490,833 | 9,167 | 113,488 | 13.1 | 13.5 | 9.9 | 0.12 | 0.03 | pass |
| synthetic-uniform-harder | `absorb(id1)` | 432,350 | 67,650 | 148,357 | 13.1 | 13.5 | 10.9 | 0.08 | 0.04 | pass |
| synthetic-uniform-harder | `absorb(id2)` | 432,884 | 67,116 | 148,930 | 13.6 | 13.6 | 10.6 | 0.13 | 0.04 | pass |
| synthetic-assortative | `absorb(id1)` | 474,059 | 25,096 | 101,005 | 13.7 | 15.0 | 12.9 | 0.12 | 0.03 | pass |
| synthetic-assortative | `absorb(id2)` | 485,167 | 13,988 | 80,945 | 13.5 | 15.0 | 12.4 | 0.08 | 0.04 | pass |
| synthetic-zigzag | `absorb(id1)` | 10,002 | 0 | 5,000 | 13.0 | 15.0 | 10.1 | 0.01 | 0.00 | pass |
| synthetic-zigzag | `absorb(id2)` | 10,002 | 0 | 5,001 | 13.8 | 13.9 | 10.5 | 0.01 | 0.00 | pass |
| enron | `absorb(id1)` | 356,451 | 11,211 | 25,481 | 13.6 | 13.6 | 10.1 | 0.06 | 0.03 | pass |
| enron | `absorb(id2)` | 356,451 | 11,211 | 25,481 | 13.8 | 13.7 | 9.1 | 0.09 | 0.02 | pass |
| github | `absorb(id1)` | 542,516 | 6,327 | 21,643 | 13.3 | 13.7 | 9.5 | 0.09 | 0.04 | pass |
| github | `absorb(id2)` | 546,204 | 2,639 | 21,719 | 13.9 | 13.9 | 10.4 | 0.13 | 0.04 | pass |
| patents | `absorb(id1)` | 491,628 | 8,380 | 93,457 | 13.5 | 13.6 | 9.6 | 0.10 | 0.04 | pass |
| patents | `absorb(id2)` | 223,472 | 276,536 | 86,459 | 13.4 | 12.8 | 9.9 | 0.12 | 0.04 | pass |
| workers | `absorb(id1)` | 499,145 | 5,170 | 23,694 | 13.3 | 13.4 | 10.1 | 0.08 | 0.04 | pass |
| workers | `absorb(id2)` | 457,079 | 47,236 | 171,154 | 13.6 | 12.0 | 9.4 | 0.12 | 0.04 | pass |
| schools | `absorb(id1)` | 413,444 | 0 | 206,722 | 13.3 | 13.4 | 9.8 | 0.07 | 0.03 | pass |
| schools | `absorb(id2)` | 413,347 | 97 | 11,863 | 13.3 | 13.4 | 10.1 | 0.12 | 0.03 | pass |
| directors | `absorb(id1)` | 215,513 | 309,499 | 86,426 | 13.6 | 13.7 | 10.4 | 0.07 | 0.03 | pass |
| directors | `absorb(id2)` | 203,784 | 321,228 | 79,622 | 13.4 | 13.8 | 9.8 | 0.09 | 0.04 | pass |
<!-- END GENERATED: oneway -->

### 2. Two fixed effects: the benchmark graph

`absorb(id1 id2)`, the specification the collection is built around. Three things happen here that do not with one effect.

1. **Singleton removal iterates.** Dropping a row that is alone in its `id1` level can leave another row alone in its `id2` level. The counts match on every graph.
2. **Mobility groups are counted.** After singletons are removed, each connected component of the graph makes one coefficient of the second effect redundant. Both commands report this in `e(dof_table)`. The test also requests `groupvar()` from both commands and requires the two identifiers to define the same partition of the rows. Group labels are not required to match.
3. **The projection iterates.** The `reghdfe` iteration count is a rough guide to how hard the graph is, from 2 on the complete design to several hundred on the sparse real networks.

`synthetic-zigzag` is not in this table. Its two-way model is saturated with two regressors, and it is treated in [section 6](#6-exact-solutions).

<!-- BEGIN GENERATED: twoway -->
| Dataset | Rows | Singletons dropped | N used | Absorbed df | Mobility groups | reghdfe iterations | Digits, b | Digits, V | Digits, fit | reghdfe (s) | creghdfe (s) | Result |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| credit2 | 19,094 | 19 | 19,075 | 721 | 1 | 9 | 14.2 | 13.9 | 10.3 | 0.05 | 0.01 | pass |
| credit | 516,810 | 656 | 516,154 | 13,429 | 1 | 8 | 13.6 | 13.2 | 9.8 | 0.75 | 0.09 | pass |
| soccer | 73,487 | 0 | 73,487 | 1,048 | 1 | 6 | 13.7 | 13.7 | 9.5 | 0.11 | 0.01 | pass |
| synthetic-complete | 500,000 | 0 | 500,000 | 1,499 | 1 | 2 | 13.7 | 13.8 | 10.2 | 0.60 | 0.06 | pass |
| synthetic-uniform-easy | 500,000 | 67,650 | 432,350 | 148,376 | 1 | 5 | 13.3 | 13.3 | 9.5 | 0.79 | 0.06 | pass |
| synthetic-uniform-hard | 500,000 | 19,713 | 480,287 | 225,450 | 1 | 19 | 11.6 | 13.8 | 10.3 | 1.17 | 0.12 | pass |
| synthetic-uniform-harder | 500,000 | 182,255 | 317,745 | 236,585 | 1 | 48 | 12.7 | 13.3 | 9.2 | 1.35 | 0.15 | pass |
| synthetic-assortative | 499,155 | 33,934 | 465,221 | 154,058 | 24,320 | 167 | 13.1 | 11.1 | 11.4 | 3.87 | 0.67 | pass |
| enron | 367,662 | 21,168 | 346,494 | 50,329 | 243 | 60 | 12.5 | 11.3 | 9.8 | 1.23 | 0.21 | pass |
| github | 548,843 | 7,906 | 540,937 | 31,190 | 11,699 | 205 | 10.3 | 11.2 | 9.3 | 3.99 | 1.30 | pass |
| patents | 500,008 | 323,893 | 176,115 | 115,710 | 1,162 | 200 | 12.2 | 10.7 | 10.1 | 1.98 | 0.33 | pass |
| workers | 504,315 | 52,264 | 452,051 | 174,996 | 14,679 | 293 | 11.4 | 9.6 | 9.6 | 5.12 | 1.28 | pass |
| schools | 413,444 | 196 | 413,248 | 218,471 | 14 | 124 | 12.3 | 11.5 | 9.2 | 2.53 | 0.19 | pass |
| directors | 525,012 | 493,782 | 31,230 | 21,366 | 3,645 | 233 | 10.2 | 10.7 | 10.3 | 0.42 | 0.07 | pass |
<!-- END GENERATED: twoway -->

### 3. Robust and clustered covariances

The same two-way model with `vce(robust)`, `vce(cluster id1)` and `vce(cluster id2)`. Clustering on an absorbed identifier makes that effect nested within the clusters, so it is removed from the absorbed degrees of freedom.

These comparisons are tested on converged solutions. The reason is numerical. A sandwich covariance sums products of projected regressors and residuals, so an error in the projection enters it at first order. It enters the coefficients and the unadjusted covariance only at second order. At `tolerance(1e-8)` both commands deliver coefficients to nine or more digits on these graphs, but on the weakly connected ones their robust covariances are several digits less accurate. The column "Digits, V (default settings)" shows the resulting agreement between the two default runs, which for the robust rows falls as low as five digits.

Tightening the tolerance of both commands does not settle the comparison, because `reghdfe`'s default projection stops improving. In runs made while writing this suite on `schools`, `synthetic-assortative`, `directors`, `github` and `synthetic-uniform-harder`, each command's default robust or clustered covariance agreed with `reghdfe`'s own `technique(lsmr)` solution to between 4.6 and 6.6 digits. At `tolerance(1e-12)` `creghdfe` agreed with it to between 8.2 and 9.8 digits, while `reghdfe`'s default projection stayed between 6.1 and 8.0 and barely changed its iteration count from `tolerance(1e-10)`. So the reference for the test is `reghdfe` with `technique(lsmr) tolerance(1e-13)`, compared with `creghdfe` with `tolerance(1e-12)`, under the same seven-digit rule. The remaining columns report that comparison. Timings are from the default runs.

Every comparison passes. On the seven graphs with several mobility groups the clustered rows depend on how the non-nested effect is counted, the first rule under [Absorbed degrees of freedom](#absorbed-degrees-of-freedom).

<!-- BEGIN GENERATED: vce -->
| Dataset | Specification | Absorbed df, reghdfe | Absorbed df, creghdfe | Digits, V (default settings) | Digits, b | Digits, V | reghdfe (s) | creghdfe (s) | Result |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| credit2 | `absorb(id1 id2) vce(robust)` | 721 | 721 | 14.1 | 13.4 | 12.7 | 0.04 | 0.01 | pass |
| credit2 | `absorb(id1 id2) vce(cluster id1)` | 139 | 139 | 14.1 | 13.4 | 12.6 | 0.02 | 0.01 | pass |
| credit2 | `absorb(id1 id2) vce(cluster id2)` | 581 | 581 | 13.8 | 13.4 | 12.1 | 0.02 | 0.01 | pass |
| credit | `absorb(id1 id2) vce(robust)` | 13,429 | 13,429 | 13.3 | 13.1 | 11.9 | 0.76 | 0.09 | pass |
| credit | `absorb(id1 id2) vce(cluster id1)` | 139 | 139 | 13.2 | 13.1 | 12.5 | 0.32 | 0.08 | pass |
| credit | `absorb(id1 id2) vce(cluster id2)` | 13,289 | 13,289 | 13.2 | 13.1 | 12.4 | 0.31 | 0.08 | pass |
| soccer | `absorb(id1 id2) vce(robust)` | 1,048 | 1,048 | 13.1 | 13.5 | 11.8 | 0.11 | 0.01 | pass |
| soccer | `absorb(id1 id2) vce(cluster id1)` | 520 | 520 | 14.0 | 13.5 | 10.9 | 0.04 | 0.01 | pass |
| soccer | `absorb(id1 id2) vce(cluster id2)` | 527 | 527 | 14.0 | 13.5 | 12.1 | 0.05 | 0.01 | pass |
| synthetic-complete | `absorb(id1 id2) vce(robust)` | 1,499 | 1,499 | 13.2 | 13.7 | 13.2 | 0.60 | 0.05 | pass |
| synthetic-complete | `absorb(id1 id2) vce(cluster id1)` | 499 | 499 | 15.0 | 13.7 | 11.5 | 0.19 | 0.05 | pass |
| synthetic-complete | `absorb(id1 id2) vce(cluster id2)` | 999 | 999 | 11.9 | 13.7 | 11.4 | 0.19 | 0.05 | pass |
| synthetic-uniform-easy | `absorb(id1 id2) vce(robust)` | 148,376 | 148,376 | 13.1 | 13.1 | 12.1 | 0.79 | 0.06 | pass |
| synthetic-uniform-easy | `absorb(id1 id2) vce(cluster id1)` | 19 | 19 | 13.3 | 13.1 | 11.8 | 0.35 | 0.06 | pass |
| synthetic-uniform-easy | `absorb(id1 id2) vce(cluster id2)` | 148,356 | 148,356 | 12.5 | 13.1 | 11.1 | 0.23 | 0.06 | pass |
| synthetic-uniform-hard | `absorb(id1 id2) vce(robust)` | 225,450 | 225,450 | 11.7 | 11.3 | 10.5 | 1.17 | 0.11 | pass |
| synthetic-uniform-hard | `absorb(id1 id2) vce(cluster id1)` | 112,707 | 112,707 | 13.8 | 11.3 | 11.2 | 0.65 | 0.12 | pass |
| synthetic-uniform-hard | `absorb(id1 id2) vce(cluster id2)` | 112,742 | 112,742 | 14.1 | 11.3 | 11.7 | 0.67 | 0.12 | pass |
| synthetic-uniform-harder | `absorb(id1 id2) vce(robust)` | 236,585 | 236,585 | 7.2 | 13.4 | 10.2 | 1.32 | 0.14 | pass |
| synthetic-uniform-harder | `absorb(id1 id2) vce(cluster id1)` | 118,465 | 118,465 | 5.7 | 13.4 | 8.7 | 0.90 | 0.14 | pass |
| synthetic-uniform-harder | `absorb(id1 id2) vce(cluster id2)` | 118,119 | 118,119 | 7.2 | 13.4 | 10.4 | 0.86 | 0.14 | pass |
| synthetic-assortative | `absorb(id1 id2) vce(robust)` | 154,058 | 154,058 | 5.4 | 13.4 | 8.2 | 3.80 | 0.66 | pass |
| synthetic-assortative | `absorb(id1 id2) vce(cluster id1)` | 78,874 | 78,874 | 5.9 | 13.4 | 9.1 | 3.76 | 0.67 | pass |
| synthetic-assortative | `absorb(id1 id2) vce(cluster id2)` | 99,502 | 99,502 | 6.1 | 13.4 | 8.9 | 3.49 | 0.66 | pass |
| enron | `absorb(id1 id2) vce(robust)` | 50,329 | 50,329 | 7.4 | 12.7 | 11.1 | 1.23 | 0.20 | pass |
| enron | `absorb(id1 id2) vce(cluster id1)` | 25,285 | 25,285 | 7.7 | 12.7 | 11.8 | 0.91 | 0.20 | pass |
| enron | `absorb(id1 id2) vce(cluster id2)` | 25,285 | 25,285 | 7.8 | 12.7 | 11.5 | 1.00 | 0.23 | pass |
| github | `absorb(id1 id2) vce(robust)` | 31,190 | 31,190 | 6.6 | 12.3 | 9.3 | 3.99 | 1.29 | pass |
| github | `absorb(id1 id2) vce(cluster id1)` | 21,353 | 21,353 | 8.0 | 12.3 | 10.4 | 3.94 | 1.29 | pass |
| github | `absorb(id1 id2) vce(cluster id2)` | 21,534 | 21,534 | 6.5 | 12.3 | 10.8 | 3.87 | 1.29 | pass |
| patents | `absorb(id1 id2) vce(robust)` | 115,710 | 115,710 | 8.4 | 11.3 | 11.3 | 1.76 | 0.28 | pass |
| patents | `absorb(id1 id2) vce(cluster id1)` | 67,523 | 67,523 | 7.8 | 11.3 | 10.0 | 1.47 | 0.27 | pass |
| patents | `absorb(id1 id2) vce(cluster id2)` | 49,347 | 49,347 | 8.5 | 11.3 | 11.0 | 1.47 | 0.27 | pass |
| workers | `absorb(id1 id2) vce(robust)` | 174,996 | 174,996 | 7.9 | 12.0 | 10.1 | 5.09 | 1.27 | pass |
| workers | `absorb(id1 id2) vce(cluster id1)` | 169,253 | 169,253 | 9.0 | 12.0 | 9.9 | 4.76 | 1.27 | pass |
| workers | `absorb(id1 id2) vce(cluster id2)` | 20,420 | 20,420 | 7.6 | 12.0 | 10.5 | 5.32 | 1.27 | pass |
| schools | `absorb(id1 id2) vce(robust)` | 218,471 | 218,471 | 6.2 | 12.7 | 9.8 | 2.53 | 0.18 | pass |
| schools | `absorb(id1 id2) vce(cluster id1)` | 11,860 | 11,860 | 6.2 | 12.7 | 9.8 | 2.06 | 0.19 | pass |
| schools | `absorb(id1 id2) vce(cluster id2)` | 206,623 | 206,623 | 7.3 | 12.7 | 9.1 | 1.97 | 0.19 | pass |
| directors | `absorb(id1 id2) vce(robust)` | 21,366 | 21,366 | 6.7 | 10.7 | 9.6 | 0.42 | 0.07 | pass |
| directors | `absorb(id1 id2) vce(cluster id1)` | 12,084 | 12,084 | 5.6 | 10.7 | 8.6 | 0.38 | 0.07 | pass |
| directors | `absorb(id1 id2) vce(cluster id2)` | 12,925 | 12,925 | 6.9 | 10.7 | 10.9 | 0.38 | 0.07 | pass |
<!-- END GENERATED: vce -->

### 4. Other pairings and constructed groups

Two-way models on identifiers other than the benchmark pair, and on groups built from them.

- **`credit2`, `credit`.** `absorb(ruc id2)` replaces `id1` with the borrower number it was built from. The partition is the same, so the answer must be, but the values exceed 2^31 and are stored as doubles. `absorb(id1 year)` uses a year taken from the month string. `absorb(id1_year id2)` absorbs borrower-year cells and months, which splits the graph into one mobility group per year.
- **`synthetic-assortative`.** `absorb(ceo_type firm_type)` is a two-by-two design. `absorb(id1 year)` and `absorb(id2 year)` pair each side of the match with time.
- **`github`.** The benchmark pair with four regressors, and `absorb(id3 month)` with a month taken from the timestamp string.
- **`workers`.** `absorb(id1 id3)` is firm by job title. `absorb(id2 firm_year)` is worker by firm-year cell.
- **`schools`.** `absorb(id1 teacher_year)` is student by teacher-year cell.

All of these pass.

<!-- BEGIN GENERATED: regroup -->
| Dataset | Specification | N used | Singletons dropped | Absorbed df | Mobility groups | Digits, b | Digits, V | Digits, fit | reghdfe (s) | creghdfe (s) | Result |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| credit2 | `absorb(ruc id2)` | 19,075 | 19 | 721 | 1 | 14.2 | 13.9 | 10.3 | 0.04 | 0.01 | pass |
| credit2 | `absorb(id1 year)` | 19,075 | 19 | 593 | 1 | 14.0 | 14.6 | 10.5 | 0.03 | 0.01 | pass |
| credit2 | `absorb(id1_year id2)` | 18,966 | 128 | 2,231 | 12 | 14.5 | 14.1 | 10.3 | 0.04 | 0.01 | pass |
| credit | `absorb(ruc id2)` | 516,154 | 656 | 13,429 | 1 | 13.6 | 13.2 | 9.8 | 0.74 | 0.09 | pass |
| credit | `absorb(id1 year)` | 516,154 | 656 | 13,301 | 1 | 13.6 | 13.4 | 11.2 | 0.46 | 0.09 | pass |
| credit | `absorb(id1_year id2)` | 513,600 | 3,210 | 55,302 | 12 | 13.2 | 13.4 | 10.0 | 0.77 | 0.06 | pass |
| synthetic-assortative | `absorb(ceo_type firm_type)` | 499,155 | 0 | 3 | 1 | 13.3 | 15.0 | 12.8 | 0.37 | 0.06 | pass |
| synthetic-assortative | `absorb(id1 year)` | 474,059 | 25,096 | 101,014 | 1 | 13.5 | 15.0 | 12.6 | 0.92 | 0.06 | pass |
| synthetic-assortative | `absorb(id2 year)` | 485,167 | 13,988 | 80,954 | 1 | 13.6 | 15.0 | 12.4 | 0.79 | 0.07 | pass |
| github | `y x1 x2 len_commit_msg num_files, absorb(id1 id2)` | 540,937 | 7,906 | 31,190 | 11,699 | 9.5 | 9.6 | 9.3 | 6.30 | 1.30 | pass |
| github | `absorb(id3 month)` | 548,227 | 616 | 1,880 | 1 | 13.4 | 14.4 | 9.6 | 0.46 | 0.07 | pass |
| workers | `absorb(id1 id3)` | 498,835 | 5,480 | 27,663 | 154 | 10.2 | 9.8 | 9.1 | 4.98 | 1.35 | pass |
| workers | `absorb(id2 firm_year)` | 442,505 | 61,810 | 192,181 | 12,158 | 10.9 | 8.5 | 9.3 | 6.76 | 0.82 | pass |
| schools | `absorb(id1 teacher_year)` | 412,718 | 726 | 229,155 | 99 | 10.0 | 9.7 | 8.5 | 5.66 | 0.43 | pass |
<!-- END GENERATED: regroup -->

### 5. Three and four fixed effects

The datasets with further identifiers add them as a third and fourth effect. All of these pass. A time effect is usually connected to the earlier effects in a single component and then costs one redundant coefficient. Job title in `workers`, `id3` in `github`, and year in the credit files, where year is nested within month, are not, and their rows depend on the second rule under [Absorbed degrees of freedom](#absorbed-degrees-of-freedom). The clustered rows are tested on converged solutions as in section 3.

<!-- BEGIN GENERATED: multiway -->
| Dataset | Specification | Tested settings | N used | Absorbed df, reghdfe | Absorbed df, creghdfe | Digits, b | Digits, V | reghdfe (s) | creghdfe (s) | Result |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| credit2 | `absorb(id1 id2 year)` | default | 19,075 | 721 | 721 | 14.7 | 13.9 | 0.05 | 0.01 | pass |
| credit | `absorb(id1 id2 year)` | default | 516,154 | 13,429 | 13,429 | 13.6 | 13.2 | 1.16 | 0.11 | pass |
| synthetic-assortative | `absorb(id1 id2 year)` | default | 465,221 | 154,067 | 154,067 | 13.0 | 11.2 | 6.21 | 0.99 | pass |
| synthetic-assortative | `absorb(id1 id2 year) vce(cluster id2)` | converged | 465,221 | 99,511 | 99,511 | 13.1 | 9.0 | 5.33 | 0.99 | pass |
| github | `absorb(id1 id2 id3)` | default | 540,674 | 32,004 | 32,004 | 11.0 | 11.4 | 8.58 | 2.73 | pass |
| github | `absorb(id1 id2 month)` | default | 540,937 | 31,196 | 31,196 | 10.2 | 11.1 | 7.03 | 2.07 | pass |
| github | `absorb(id1 id2 id3 month)` | default | 540,674 | 32,010 | 32,010 | 11.0 | 11.1 | 12.96 | 3.78 | pass |
| workers | `absorb(id1 id2 id4)` | default | 452,051 | 174,998 | 174,998 | 11.6 | 9.7 | 8.77 | 1.80 | pass |
| workers | `absorb(id1 id2 id3)` | default | 451,626 | 178,005 | 178,005 | 9.4 | 8.1 | 13.00 | 2.82 | pass |
| workers | `absorb(id1 id2 id3 id4)` | default | 451,626 | 178,007 | 178,007 | 9.3 | 8.0 | 18.35 | 3.66 | pass |
| workers | `absorb(id1 id2 id3 id4) vce(cluster id1)` | converged | 451,626 | 172,296 | 172,296 | 12.0 | 9.3 | 18.11 | 3.71 | pass |
| schools | `absorb(id1 id2 id3)` | default | 413,248 | 218,474 | 218,474 | 11.9 | 11.3 | 4.41 | 0.36 | pass |
| schools | `absorb(id1 id2 id3) vce(cluster id2)` | converged | 413,248 | 206,626 | 206,626 | 13.2 | 9.1 | 3.69 | 0.35 | pass |
<!-- END GENERATED: multiway -->

### 6. Exact solutions

Agreement between two iterative solvers does not show that either is right. Three datasets allow a check against a regression that does not iterate. Each row compares one run with that exact solution: the slopes, their covariance block, the residual sum of squares and the residual degrees of freedom.

- **`credit2` and `soccer`** are small enough for `areg y x1 x2 i.id2, absorb(id1)`, which solves the second effect as explicit dummies. Both commands match it.
- **`synthetic-zigzag`** is a path: 10,000 distinct `(id1, id2)` cells linking 10,001 nodes. On a tree, any vector of cell means can be written as an `id1` effect plus an `id2` effect, so projecting on the two effects is the same as projecting on the cells, and `areg y x1, absorb(cell)` is exact. Two cells hold two rows each. All the within variation is in those four rows.

  At default tolerance both commands stop about one percent from the exact slope, which is two digits. This is the stopping rule meeting a graph whose projections converge extremely slowly. It is not a difference between the commands, and it is recorded, not tested. With `tolerance(1e-10)`, `creghdfe` reproduces the exact solution, and that is the test. In runs made while writing this suite, `reghdfe`'s default projection stayed 0.6 percent from the exact slope at every tolerance down to `1e-14`, and its `technique(lsmr)` reproduced it.

- **Saturated model.** With both `x1` and `x2` the zigzag model has 10,002 rows, 10,000 absorbed coefficients and two regressors: no residual degrees of freedom. `reghdfe` posts coefficients with a missing covariance. `creghdfe` returns `r(498)`. The suite asserts both behaviours, so a change to either is noticed.

<!-- BEGIN GENERATED: oracle -->
| Dataset | Specification | Run | Residual df | Exact residual df | Digits, b | Digits, V | Digits, RSS | Result |
|---|---|---|---:|---:|---:|---:|---:|---|
| credit2 | `y x1 x2, absorb(id1 id2)` | `reghdfe` | 18,352 | 18,352 | 14.3 | 15.0 | 15.0 | recorded |
| credit2 | `y x1 x2, absorb(id1 id2)` | `creghdfe` | 18,352 | 18,352 | 14.8 | 14.2 | 15.0 | pass |
| soccer | `y x1 x2, absorb(id1 id2)` | `reghdfe` | 72,437 | 72,437 | 13.9 | 13.5 | 15.0 | recorded |
| soccer | `y x1 x2, absorb(id1 id2)` | `creghdfe` | 72,437 | 72,437 | 14.0 | 13.5 | 13.9 | pass |
| synthetic-zigzag | `y x1, absorb(id1 id2)` | `reghdfe` | 1 | 1 | 2.0 | 2.5 | 2.5 | recorded |
| synthetic-zigzag | `y x1, absorb(id1 id2)` | `creghdfe` | 1 | 1 | 1.9 | 2.5 | 2.4 | recorded |
| synthetic-zigzag | `y x1, absorb(id1 id2)` | `creghdfe tolerance(1e-10)` | 1 | 1 | 10.3 | 8.5 | 8.5 | pass |
<!-- END GENERATED: oracle -->

## Absorbed degrees of freedom

Both commands report, for each absorbed effect, its levels minus its redundant coefficients, and `e(df_a)` is the total. Two of `reghdfe`'s rules for the redundant count only matter on graphs like these, and the suite is what brought `creghdfe` into line with them.

1. **Effects nested within the clusters are set aside first.** With `absorb(id1 id2) vce(cluster id1)`, `id1` costs no degrees of freedom, and mobility groups are then counted only among the effects that remain. `id2` is the only one, so it loses a single redundant coefficient: `e(df_a)` is the number of `id2` levels minus 1, not minus the number of mobility groups between `id1` and `id2`.
2. **Later effects are counted pairwise.** A third or later effect is compared with every earlier effect that is not nested, and its redundant count is the largest number of mobility groups it forms with any of them. In the credit files each year contains whole months, so all 12 year levels are redundant. `reghdfe` marks such counts as possibly too low in its table.

`dofadjustments()` selects these adjustments as in `reghdfe`: `all`, `none`, or any of `firstpair` or `pairwise`, `clusters` and `continuous`.

Builds of `creghdfe` from before these rules subtracted the first pair's mobility groups even when one of the pair was nested, and one coefficient for each later effect. Coefficients were the same. `reghdfe`'s clustered standard errors were up to 9.5 percent larger on `directors`, with 3,645 mobility groups in 31,230 observations, and 0.003 percent larger on `schools`, with 14 in 413,248.

## Differences found

### 1. Saturated models

With no residual degrees of freedom `reghdfe` posts estimates and a missing covariance, and `creghdfe` returns `r(498)`, the error it also returns for a one-cluster covariance (see [COMPATIBILITY.md](../COMPATIBILITY.md)).

### 2. Default tolerance on path-like graphs

On `synthetic-zigzag` both commands report convergence at `tolerance(1e-8)` while about one percent from the exact slope. Neither default is adequate for a graph of that shape. `creghdfe` reaches the exact solution when the tolerance is tightened.

### 3. Stored results the suite does not compare

The fourth column of `e(dof_table)` differs: `reghdfe` flags redundant counts that are only bounds, while `creghdfe` does not compute the flag and leaves it missing. Since 27 October 2026, `e(N_hdfe_extended)` and `e(df_a_initial)` follow `reghdfe` (the absorbed effects; the levels before any redundancy is removed). Earlier versions of `creghdfe` stored the mobility groups and the count after redundancy there.

## Timing notes

The totals above and the per-row seconds are single wall-clock runs taken during a correctness run, with the desktop in use. They include each command's parsing, singleton detection, degrees-of-freedom computation and estimation. `reghdfe` runs under Stata/MP with the licensed processor count, and `creghdfe`'s plugin uses its default thread count on the same machine. They show how the two commands behaved on these graphs on this machine on this day. They are not a controlled benchmark, and, as [COMPATIBILITY.md](../COMPATIBILITY.md) says, no general speedup is promised.

## Files

| File | Purpose |
|---|---|
| [`validation/performance/datasets.do`](../../validation/performance/datasets.do) | The comparisons |
| [`validation/bench_profiles.py`](../../validation/bench_profiles.py) | Checks the data cache, runs Stata, regenerates the tables here |
| [`scripts/fetch_validation_data.py`](../../scripts/fetch_validation_data.py) | `--collection hdfe` downloads and verifies the datasets |
| [`validation/fixtures/hdfe_manifest.json`](../../validation/fixtures/hdfe_manifest.json) | Source URL and SHA-256 of each dataset |
| [`benchmarks/hdfe_datasets_results.csv`](../benchmarks/hdfe_datasets_results.csv) | One row per run: settings, sample, degrees of freedom of both commands, digits, seconds, differences |
| [`benchmarks/hdfe_datasets_oracle.csv`](../benchmarks/hdfe_datasets_oracle.csv) | The exact-solution comparisons |

The datasets are Sergio Correia's. Their sources and construction are described on the [collection's pages](https://www.scorreia.com/data/hdfe/index.html). They are downloaded to a local cache and are not redistributed with ctools.
