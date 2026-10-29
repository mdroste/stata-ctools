# crangestat

## Validation and failure behavior

Output targets must be distinct new variable names, including when a target aliases a source or interval key. Failed or interrupted operations restore the dataset without writing a snapshot to disk (the result variables are dropped again); option `snapshot` uses `preserve` instead and option `fast` undoes nothing. Missing interval keys are excluded from the population used to compute statistics, including unbounded windows and excludeself, on both sides of the parallel execution threshold.
