# crangestat

`crangestat` computes statistics over inclusive numeric-key windows, optionally
within groups. It adds result variables and preserves observation order.

```stata
crangestat (mean) mean_price=price (sd) sd_price=price, ///
    interval(date -30 0) by(firm)
```

For a row with key `date`, this example uses observations in the same firm whose
keys lie between `date - 30` and `date`. A bound is a number (added to the key),
`.` for an unbounded endpoint, or a numeric variable holding the bound itself for
each observation, as in `rangestat`.
Without an explicit output name, `(mean) price` creates `price_mean`.

Supported statistics include `count`, `obs`, `missing`, `mean`, `sum`, `min`,
`max`, `sd`, `variance`, `median`, `iqr`, `first`, `last`, `firstnm`, `lastnm`,
and the percentiles listed in [the help file](../../build/crangestat.sthlp).
Percentiles, `median` and `iqr` use Stata's `_pctile` definition. `if` and `in`
restrict the observations that enter windows and receive results, as in
`rangestat`. The flexible statistics `(cov)`, `(corr)`, `(reg)` and user-written
Mata functions are computed in Mata with `rangestat`'s output names. `casewise`
drops, in each window, the observations with a missing value of any variable of
the same statistic, as `rangestat` does, for every statistic. `describe`,
`local(name)` and `nocheck` are accepted as in `rangestat`, and `r(newvars)` and
`r(rangestat_version)` are returned. Variable ranges may have spaces around the
dash (`x1 - x3`). `excludeself` omits the current
row from its window; `threads(#)` limits workers and `verbose` reports timings.
`by()` variables may be strL, with text of any length and binary values.

## Validation and failure behavior

Output targets must be distinct new variable names, including when a target aliases a source or interval key. Failed or interrupted operations restore the dataset without writing a snapshot to disk (the result variables are dropped again); option `snapshot` uses `preserve` instead and option `fast` undoes nothing. Missing interval keys are excluded from the population used to compute statistics, including unbounded windows and excludeself, on both sides of the parallel execution threshold.
