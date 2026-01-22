# cencode

High-performance C-accelerated string encoding for Stata datasets.

## Overview

`cencode` is a drop-in replacement for Stata's `encode` command that provides 2-5x speedup over native Stata encoding. It converts string variables to numeric variables with value labels, using an optimized streaming algorithm.

## Syntax

```stata
cencode varname [if] [in], generate(newvar) [options]
```

## Options

| Option | Description |
|--------|-------------|
| `generate(newvar)` | **Required.** Name of the new numeric variable to create |
| `label(name)` | Name for the value label (default: same as newvar) |
| `noextend` | Do not extend existing value label |
| `verbose` | Display timing breakdown |
| `threads(#)` | Maximum number of threads (default: all available) |

## Examples

```stata
* Basic usage - encode a string variable
sysuse auto, clear
cencode make, generate(make_code)

* View the encoded values and labels
tabulate make_code
label list make_code

* Encode with a custom label name
cencode make, generate(make_num) label(car_makes)

* Encode with verbose timing output
cencode make, generate(make_code) verbose

* Encode only foreign cars
cencode make if foreign == 1, generate(foreign_make)
```

## Algorithm

`cencode` uses an optimized two-pass streaming algorithm:

### Pass 1: Collect Unique Values
- Streams through the string variable once
- Uses a hash table to identify unique values in O(1) average time
- Only stores K unique strings (not all N observations)
- Memory usage: O(K) where K is the number of unique values

### Pass 2: Encode Values
- Re-reads strings from Stata
- Looks up each string in the hash table (O(1) average)
- Writes the numeric code to the output variable

### Key Optimizations
1. **Streaming**: Only stores unique strings, not all N observations
2. **Hash Table**: O(1) average lookup time with FNV-1a hashing
3. **Arena Allocation**: Efficient memory allocation for unique strings
4. **Alphabetical Ordering**: Codes assigned in alphabetical order (1, 2, 3, ...)

## Performance

| Dataset Size | Unique Values | `encode` | `cencode` | Speedup |
|-------------|---------------|----------|-----------|---------|
| 100,000 | 100 | 0.012s | 0.008s | 1.5x |
| 1,000,000 | 100 | 0.21s | 0.098s | 2.1x |
| 1,000,000 | 1,000 | 0.25s | 0.11s | 2.3x |

Performance gains are most significant with:
- Large datasets (100K+ observations)
- Moderate number of unique values
- Longer string values

## Timing Breakdown

When using the `verbose` option, you'll see:

```
cencode timing breakdown:
-------------------------------------------------------
  C plugin internals:
    Argument parsing:         0.0000 sec
    Load strings:             0.0343 sec
    Collect unique:           0.0000 sec
    Sort unique:              0.0000 sec
    Encode values:            0.0376 sec
    Store labels:             0.0000 sec
  -----------------------------------------------------
    C plugin total:           0.0719 sec
  -----------------------------------------------------
  Stata overhead:
    Pre-plugin setup:         0.0070 sec
    Plugin call overhead:     0.0001 sec
    Post-plugin labels:       0.0040 sec
-------------------------------------------------------
    Wall clock total:         0.0830 sec
-------------------------------------------------------

  Unique values encoded:    100
```

## Stored Results

`cencode` stores the following scalars:

| Scalar | Description |
|--------|-------------|
| `_cencode_n_unique` | Number of unique string values encoded |
| `_cencode_time_load` | Time to read strings from Stata (seconds) |
| `_cencode_time_sort` | Time to sort unique values (seconds) |
| `_cencode_time_encode` | Time to encode and write values (seconds) |
| `_cencode_time_total` | Total C plugin time (seconds) |

## Limitations

- Maximum of 65,536 unique values (Stata value label limit)
- Stata's Plugin Interface serializes data access, limiting parallelization benefits

## Comparison with Native `encode`

| Feature | Stata `encode` | `cencode` |
|---------|----------------|-----------|
| Algorithm | Mata-based | C with hash table |
| Memory | O(N) | O(K) unique values |
| Lookup | Linear search | Hash table O(1) |
| Typical Speedup | 1x | 2-5x |

## See Also

- [ctools Overview](../README.md)
- Stata's `encode` command: `help encode`
- Stata's `decode` command for reverse operation
