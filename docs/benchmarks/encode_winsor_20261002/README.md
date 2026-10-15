# Encoding and winsor performance evidence

See [the report](../../PERFORMANCE_ENCODE_WINSOR.md) for interpretation, methods,
limitations, and the outstanding integrated winsor validation failures.

- `integrated_encode` and `integrated_winsor` compare frozen builds containing the
  same concurrent compatibility changes, with only this optimization reversed in
  the baseline. These are the primary checkout measurements.
- `encode_1m`, `encode_3m`, `winsor_1m`, and small-input campaigns isolate the original
  performance changes before concurrent compatibility edits.
- `pilot` is exploratory; its supposed existing-label control did not retain the
  label through save/use. That fixture was corrected for every later campaign.
- `*references*` compare ctools directly with winsor2 and gstats on identical inputs.
- `profile*` retain macOS stack samples and their drivers. Counts include waiting
  threads and differing numbers of completed commands; use them diagnostically.
- `integrated_test*` and `attribution_baseline*` retain the same winsor test failures
  with and without the optimization. Assertions and tolerances were not relaxed.
- `optimizations.patch` isolates this task's five production-file changes from the
  already dirty starting checkout. `concurrent_changes.patch` records the unrelated
  compatibility changes present in the integrated pair.
- Frozen source/executable hashes are retained. Large datasets and actual binaries
  remain under `temp/encode_performance_20261002/` on the workstation.

Logs use `.txt` here. License serial numbers and licensee names are redacted.
`SHA256SUMS` covers the retained evidence bytes after that redaction.
