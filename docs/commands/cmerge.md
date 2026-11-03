# cmerge

Merge Stata datasets.

## Overview

`cmerge` joins datasets by matching key variables. It uses C code for data loading, sorting, and matching, and calls Stata's `merge` for datasets containing `strL` variables.

## Syntax

```stata
cmerge merge_type varlist using filename [, options]
```

Where `merge_type` is one of:
- `1:1` - One-to-one merge
- `m:1` - Many-to-one merge
- `1:m` - One-to-many merge
- `m:m` - Pair repeated keys sequentially, repeating the shorter group's last row; this does not form all combinations

As with `merge`, `1:n`, `n:1`, `n:n`, `n:m` and `m:n` are accepted for these types.
Use `cmerge 1:1 _n using filename` to merge by observation number.

Rows come out in `merge`'s order: the master observations in key order (the first
output row of each), then every other row, that is the using-only rows and the
later `1:m` or `m:m` pairings of a master row, in key order. A `1:1` result
without using-only rows is sorted by the keys. Where `merge` sorts a master whose
keys have ties, `merge`'s order within the ties depends on `set sortseed`;
`cmerge` keeps their existing order.

## Options

| Option | Description |
|--------|-------------|
| `keep(results)` | Keep `master`, `using`, `match`, `match_update`, `match_conflict`, or codes 1–5 |
| `assert(results)` | Verify merge results match expected values |
| `generate(varname)` | Name for merge indicator variable (default: `_merge`) |
| `nogenerate` | Do not create merge indicator variable |
| `keepusing(varlist)` | Variables to keep from using dataset |
| `sorted` | Assert data is already sorted (skip internal sorting) |
| `force` | Retain master type for incompatible non-key variables; key types must match |
| `noreport` | Suppress merge result table |
| `verbose` | Display timing and progress information |
| `nolabel` | Do not copy value labels from using data (`nolabels` is accepted, as in `merge`) |
| `nonotes` | Do not copy notes from using data |
| `update` | Fill missing master values from same-named using variables |
| `replace` | Replace conflicting nonmissing values; requires `update` |
| `preserve_order(#)` | Preserve master observation order when set to 1 |
| `threads(#)` | Maximum threads; zero uses the ctools thread setting |
| `snapshot` | Restore the master after a failure from a `preserve` snapshot (written to disk above `c(max_preservemem)`) |
| `fast` | Keep no copy of the master; a failure or Break can leave the data partly merged |

## Examples

```stata
* One-to-one merge
cmerge 1:1 id using data2.dta

* Many-to-one merge keeping only matched observations
cmerge m:1 state year using statedata.dta, keep(match)

* Keep only specific variables from using dataset
cmerge 1:1 id using fulldata.dta, keepusing(var1 var2)

* Merge with verbose output and no merge variable
cmerge m:1 id using lookup.dta, nogenerate verbose

* Assert all observations should match
cmerge 1:1 id using data2.dta, assert(match)

* Retain master types when non-key variables have incompatible types
cmerge m:1 code using lookup.dta, force
```

## Performance

`cmerge` achieves its speed advantage through a three-phase approach:

### Phase 1: Load Using Dataset
- Loads the using dataset into C memory
- Uses parallel data loading with 8-way loop unrolling

### Phase 2: Plan and Execute Merge
- Plans the join before the master is changed: loads the master keys alone, sorts them with a high-performance radix sort (if not pre-sorted), and joins them in a single pass; key, sortedness and uniqueness errors stop here
- If the result is the master in its existing order, no other master data are read
- Otherwise stages whole master observations and writes each output observation once, then fills the using variables in parallel

### Phase 3: Apply Options
- Applies `keep()` and `assert()` filters
- Renames placeholder variables to match using dataset

### Implementation details

- **Parallel 4-pass LSD radix sort**: Sorts key variables using 8-bit radix (256 buckets) with OpenMP-parallelized histogram and scatter phases; pointer swapping between passes avoids copying
- **Cache-line aligned histograms**: 64-byte aligned per-thread histograms prevent false sharing during parallel counting
- **Hybrid group search**: Uses linear scan for the first 16 elements, then switches to binary search for larger groups—exploits branch prediction for small groups
- **Specialized numeric key comparison**: Fast path skips string comparison logic when all key variables are numeric
- **Pre-sorted detection**: Automatically detects if data is already sorted and skips the sort phase entirely
- **Parallel data loading**: SIMD-accelerated (AVX2/SSE2/NEON) bulk load with 16x loop unrolling from Stata's Plugin Interface
- **Overflow-safe allocation**: `ctools_safe_malloc2` checks allocation sizes for integer overflow before requesting memory
- **Arena allocator for strings**: Bulk-allocated string storage with O(1) allocation and single bulk free, reducing `malloc` overhead
- **Persistent thread pool**: Worker threads are reused across sort and merge phases

## Merge Indicator Values

The `_merge` variable (or variable specified by `generate()`) indicates:

| Value | Label | Meaning |
|-------|-------|---------|
| 1 | master only | Observation from master, no match in using |
| 2 | using only | Observation from using, no match in master |
| 3 | matched | Observation matched without an update or conflict |
| 4 | match update | Missing master values updated, without a conflict |
| 5 | match conflict | Conflicting nonmissing values found |

## Stored Results

`cmerge` stores the following in `r()`:

### Scalars
| Result | Description |
|--------|-------------|
| `r(N)` | Number of observations after merge |
| `r(N_1)` | Number of observations only in master |
| `r(N_2)` | Number of observations only in using |
| `r(N_3)` | Matches without updates or conflicts |
| `r(N_4)` | Updated matches (with `update`) |
| `r(N_5)` | Conflicting matches (with `update`) |
| `r(time)` | Elapsed time in seconds |

### Macros
| Result | Description |
|--------|-------------|
| `r(using)` | Name of using file |
| `r(keyvars)` | Key variable names |

## Technical Notes

- Key variables must both be numeric or both be string, even with `force`
- The using file must be a Stata `.dta` file (extension can be omitted)
- Variable name conflicts: master values are kept by default
- Missing values in key variables are treated as valid match keys

## Comparison with Native `merge`

| Feature | Stata `merge` | `cmerge` |
|---------|---------------|----------|
| Implementation | Stata + Mata | Pure C |
| Parallelization | No | Yes (OpenMP) |
| Sort Algorithm | Stata default | Radix sort |

Memory use depends on the dataset, options, and number of threads.

## See Also

- [ctools Overview](../../README.md)
- [csort](csort.md) - High-performance sorting

## Validation and failure behavior

Shared numeric and fixed-string variables are widened before updates or using-only rows are written. A long/float combination uses double to preserve both inputs. Key variables must have matching string/numeric types even with force. For incompatible non-key types, force retains the master type and treats the using-side values as missing. String writes retain the shared plugin limit of str2045; strL output is unsupported. If `assert()` fails, the merged result stays in memory, as with `merge`; any other failed or interrupted merge (including Break) restores the master dataset.

By default, `cmerge` keeps an in-memory copy of the master observations until the merge succeeds. If the command fails or is interrupted, it restores the observations and reverses changes to variables and dataset attributes. Merges involving `strL`, and certain merges with empty datasets, keep a separate copy in a Stata frame. If that copy does not fit in memory, the merge fails before changing the dataset. On Stata 15 and earlier, those merges use `preserve`.

Two options change this:

- `snapshot` preserves the whole dataset for every merge, as earlier versions did (written to disk above `c(max_preservemem)`).
- `fast` keeps no copy at all, like `collapse, fast`. A failure or Break after the merge has started changing the data is reported, but the data are left partly merged. Key errors are still found before any change. It saves memory, and time only where a copy would otherwise be taken (`keep()` or `update` on data already in key order, strL variables, empty data).
