# crangestat

`crangestat` computes statistics over inclusive numeric-key windows, optionally
within groups. It adds result variables and preserves observation order.

```stata
crangestat (mean) mean_price=price (sd) sd_price=price, ///
    interval(date -30 0) by(firm)
```

For a row with key `date`, this example uses observations in the same firm whose
keys lie between `date - 30` and `date`. Bounds must be numeric constants or `.`
for an unbounded endpoint; observation-specific bound variables are unsupported.
Without an explicit output name, `(mean) price` creates `price_mean`.

Supported statistics include `count`, `mean`, `sum`, `min`, `max`, `sd`,
`variance`, `median`, `iqr`, `first`, `last`, `firstnm`, `lastnm`, and the
percentiles listed in [the help file](../../build/crangestat.sthlp).
User-written Mata functions are unsupported. `excludeself` omits the current
row from its window; `threads(#)` limits workers and `verbose` reports timings.

## Validation and failure behavior

Output targets must be distinct new variable names, including when a target aliases a source or interval key. Failed or interrupted operations restore the dataset without writing a snapshot to disk (the result variables are dropped again); option `snapshot` uses `preserve` instead and option `fast` undoes nothing. Missing interval keys are excluded from the population used to compute statistics, including unbounded windows and excludeself, on both sides of the parallel execution threshold.
