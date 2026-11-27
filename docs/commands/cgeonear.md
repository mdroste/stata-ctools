# cgeonear

Nearest neighbors and distances by geodetic distance, with the syntax,
options and results of Robert Picard's `geonear` and `geodist`, in one
command.

`cgeonear` is pinned to `geonear` 2.0.3 (September 4, 2019) and `geodist`
1.1.0 (June 18, 2019). See [Attribution](../ATTRIBUTION.md#cgeonear).

## Overview

`cgeonear` has two forms. Given a using file, it finds for each base location
in the data in memory the nearest neighbor locations in that file: the
nearest `nearcount()` of them and, in long form, every one within a distance
(`geonear`). Given four coordinates, it computes the distance between two
points in every observation, or between two given points (`geodist`). A
`using` file (or a `neighbors()` option) selects the first form.

Distances are geodetic, the length of the shortest path along the surface of
a model of the earth: great-circle distances on a sphere (haversine formula)
or distances on a reference ellipsoid (Vincenty's formula). `geonear`
defaults to the sphere and `geodist` to the WGS 1984 ellipsoid, and so do the
two forms of `cgeonear`.

## Syntax

```stata
* nearest neighbors (geonear)
cgeonear baseid baselat baselon using nborfile, neighbors(nborid nborlat nborlon) ///
    [wide long nearcount(#) ignoreself genstub(str) within(#) limit(#) ///
    miles radius(#) ellipsoid a(#) f(#) ops(#) report(#) maxbases(str) ///
    threads(#) verbose]

* distances (geodist)
cgeonear lat1 lon1 lat2 lon2 [if] [in], generate(newvar) [miles sphere radius(#) ///
    ellipsoid(#1,#2) maxiter(#) threads(#) verbose]
cgeonear lat1 lon1 lat2 lon2 [, miles sphere radius(#) ellipsoid(#1,#2) maxiter(#)]
```

Coordinates are signed decimal degrees, positive north and east. In the
distance form each coordinate is a numeric variable, a numeric scalar or a
number; with no variable, `cgeonear` displays the distance and returns
`r(distance)` (and `r(iterations)` on the ellipsoid).

## Nearest-neighbor options

| Option | Description |
|--------|-------------|
| `wide` | One pair of variables per neighbor, added to the data in memory (the default) |
| `long` | Replace the data with one observation per base and neighbor |
| `nearcount(#)` | The `#` nearest neighbors (default 1; `0` in long form with `within()`) |
| `ignoreself` | Ignore the neighbor with the base's identifier |
| `genstub(str)` | Wide form: name the variables `str`, `km_to_str` (default `nid`) |
| `within(#)` | Long form: also every neighbor within `#` km (miles with `miles`) |
| `limit(#)` | Long form: at most `#` neighbors per base, the nearest |
| `miles` | Distances in miles |
| `radius(#)` | Sphere radius in km (default 6371) |
| `ellipsoid` | WGS 1984 ellipsoid |
| `a(#) f(#)` | Another ellipsoid: semi-major axis (m), reciprocal of the flattening |
| `ops(#)`, `report(#)`, `maxbases()` | `geonear`'s region size and progress reports; accepted, no effect |
| `threads(#)` | Maximum threads; zero uses the ctools thread setting |
| `verbose` | Display where the time went |

## Distance options

| Option | Description |
|--------|-------------|
| `generate(newvar)` | The new variable of distances (required when a coordinate is a variable) |
| `miles` | Distances in miles |
| `sphere` | Great-circle distances (the default is the WGS 1984 ellipsoid) |
| `radius(#)` | Great-circle distances on a sphere of radius `#` km |
| `ellipsoid(#1,#2)` | Another ellipsoid: semi-major axis (m), reciprocal of the flattening |
| `maxiter(#)` | Bound on Vincenty's iterations (default 25; undocumented in `geodist`) |
| `threads(#)`, `verbose` | As above |

## Examples

```stata
* The nearest other cell tower of each tower
use cell, clear
cgeonear ctid ctlat ctlon using cell.dta, neighbors(ctid ctlat ctlon) ignoreself

* The nearest tower and every tower within 5 km of each block group
use blockgroups, clear
cgeonear bgid bglat bglon using cell.dta, neighbors(ctid ctlat ctlon) long within(5)

* The three nearest block groups of each tower, on the ellipsoid, in miles
use cell, clear
cgeonear ctid ctlat ctlon using blockgroups.dta, n(bgid bglat bglon) near(3) ///
    genstub(bg) ellipsoid miles

* Distances to a fixed point, and between two points
cgeonear 42.265837 -83.748696 lat lon, generate(d)
cgeonear 42.265837 -83.748696 42.207667 -83.356022, miles
```

## What the nearest-neighbor form gives

**Wide form.** For one neighbor the new variables are `nid` and `km_to_nid`;
for several, `nid1 km_to_nid1 nid2 km_to_nid2 ...`, nearest first. The
identifiers have `nborid`'s storage type. Every observation of a base
identifier gets its neighbors, also observations of that identifier whose
coordinates are missing; observations with a missing identifier, or one
without coordinates, get missing values. The data are then ordered by
`baseid`.

**Long form.** The data become one observation per base location and
neighbor, with `baseid`, `nborid` and `km_to_nborid`, sorted by `baseid`,
distance and `nborid`. As in `geonear`, a neighbor at distance 0 (the same
place) is always included, `within()` and `nearcount()` combine (a neighbor
that meets either is kept), and `limit()` then keeps the nearest.

**Ties.** Neighbors at the same distance are ranked by ascending `nborid`,
the order `geonear` documents.

**Metadata.** The long-form result keeps the data label, characteristics and
notes of the base data, everything that describes `baseid` and `nborid`
(formats, variable and value labels, notes, characteristics), and the value
labels attached to the six location variables. The using file's dataset notes
are appended (a note already present is not repeated) and its other dataset
characteristics are added where the data have none. In wide form only the
dataset notes and characteristics are merged that way, as in `geonear`.

**Errors.** Base and neighbor locations are checked as `geonear` checks them,
with its messages and return codes: identifiers whose coordinates differ
(r(459)), no location (r(459)), string coordinates (r(7)), coordinates out of
range (r(198)), fewer neighbors than `nearcount()` needs (r(2001)). A failure
leaves the data as they were.

## Differences from geonear and geodist

- Wide form: a variable that already has a result variable's name is an error
  (r(110)); `geonear` keeps its values and adds the other variables. A
  `_merge` variable in the data is allowed.
- Wide form: observations of one `baseid` keep their order. Data already in
  `baseid` order are left as they are (reading the base locations finds
  out); others are reordered by csort and marked as sorted by `baseid`.
  `geonear` leaves ties in no particular order and clears the sort marker.
- Ties are always ranked by `nborid`. In wide form with `nearcount()` above 1
  or `ignoreself`, `geonear` can rank neighbors at the same distance but at
  different places by their distance from the middle of its region instead.
- With `ellipsoid` in wide form, `geonear` measures some pairs from the
  neighbor to the base, and Vincenty's formula is not symmetric in its last
  digits: such distances differ by up to about 2e-13 relative.
- A neighbor whose Vincenty distance does not converge (nearly antipodal
  points) ranks after every other; `geonear` can then return one neighbor
  twice in wide form.
- `strL` identifiers are refused (r(109)), and `f()` must exceed 1.
- `cgeonear` shows the threads used where `geonear` shows its number of
  regions, and no progress lines.
- `r()` is cleared after the nearest-neighbor form; `geonear` leaves the
  results of the last command it ran internally.
- Distance form: none known. The undocumented `maxiter()` behaves as in
  `geodist`, including its scalar-mode report of a failure when the iteration
  converges on its last pass.

## Accuracy

The haversine and Vincenty formulas follow the Stata expressions of `geonear`
and `geodist` operation by operation (Stata evaluates `a*b/c` as `a*(b/c)`
and `x^2` as `x*x`, and so does the plugin). Stata's `tan()` is `sin()/cos()`
and its `asin(y)` is `atan(y/sqrt(1-y*y))`; on Windows its `sin()`, `cos()`
and `atan2()` are those of the Universal C Runtime, which the plugin calls
too. On StataNow 19.5 for Windows, 100,000 random arguments gave the same
`sin`, `cos`, `tan`, `atan`, `asin` and `sqrt` in every case and the same
`atan2` in all but 7, where the last bit differs. The validation suite
therefore finds `cgeonear`'s results identical to `geonear`'s (except the
wide-form ellipsoid case above) and its distance form identical to `geodist`'s in
all but a few distances in 10,000, which differ by one or two units in the
last place.

## Performance

`geonear` splits the base locations into ever smaller regions, keeps for each
the neighbors that may be near it, and computes the distances within each
region in Stata. `cgeonear` makes three plugin calls:

1. **Base locations.** The identifier and coordinates of every observation
   are read in one pass and put in identifier order (the radix sorts of
   `csort`, skipped when the data are already in order); one location per
   identifier, checked as `geonear` checks them.
2. **Neighbor locations.** Only the three neighbor variables are loaded from
   the file, in a frame of their own (Stata 16 or later), and read the same
   way.
3. **Search.** The neighbors are placed in a balanced k-d tree of their
   positions in three dimensions (on the unit sphere, or earth-centred
   coordinates on the ellipsoid). The straight line between two points is
   never longer than the path along the surface, so a branch of the tree
   whose box is farther away than the current k-th nearest distance cannot
   hold a nearer neighbor. Each base location is searched in parallel, in the
   order of a Morton curve so that consecutive searches share the tree's
   nodes in cache; the exact distance is computed only for the few neighbors
   that pass the straight-line bound, and the ranking uses that distance.
   Wide-form results are written in one parallel pass over the
   observations; long-form pairs are laid out in identifier order and
   written into the new dataset.

The distance form reads the coordinate variables and writes the distances
in one parallel pass each.

Memory: about 40 bytes per base location and 95 per neighbor location while
the search runs (plus the text of string identifiers), 4 bytes per
observation of the master data, and the results (12 bytes per neighbor of
each base in wide form, 16 per pair in long form).

Where the time goes, from `verbose`, for the README workload at 10M
observations (10M base locations, 50,000 neighbors, the nearest of each; 16
threads on a Ryzen 9 9950X shared with other work): reading the base
locations 0.23 s, the search 0.43 s (4 geodetic distances per base on
average), writing the two variables 0.04 s. Data that are not already in
`baseid` order are then reordered by `csort`, about 0.3 s more for the
profile's 20 variables.

The README speedups come from the benchmark profile, which runs both commands
on its 20-variable datasets (see the
[published results](../benchmarks/readme_profile/results.md) and
[how the comparison is run](../performance/README.md#compare-nearest-neighbors-and-distances)):
the nearest of `N`/200 places for every observation, and one ellipsoid
distance per observation. Smallest time of three passes, StataNow 19.5 MP with
4 processors for `geonear` and `geodist`, 16 threads for ctools, on a machine
shared with other work:

| Stata command | Observations | Stata | `cgeonear` | Speedup |
|---------------|-------------:|------:|-----------:|--------:|
| `geonear` | 1M | 37.3 s | 0.086 s | 433× |
| `geonear` | 10M | 357 s | 0.78 s | 455× |
| `geodist` | 1M | 1.73 s | 0.038 s | 45× |
| `geodist` | 10M | 17.1 s | 0.27 s | 63× |
| `geodist` | 25M | 43.7 s | 0.72 s | 61× |

`geonear` is not run at 25M (it takes over ten minutes there).
