# cfuzzyjoin

Join records whose keys are similar rather than equal, with the syntax,
options, scores and results of three community commands:

| Subcommand | Replaces | What it does |
|------------|----------|--------------|
| `cfuzzyjoin link` | `reclink` (Michael Blasnik) | The best fuzzy match of each master record, scored from several variables |
| `cfuzzyjoin match` | `matchit` (Julio D. Raffo) | Every pair of similar strings of two datasets, or the similarity of two columns |
| `cfuzzyjoin group` | `strgroup` (Julian Reif) | Groups of similar values of a string variable (Levenshtein distance) |

`cfuzzyjoin` is pinned to `reclink` 1.7, `matchit` 1.5.2 and `strgroup`
1.0.4, all from SSC. See [Attribution](../ATTRIBUTION.md#cfuzzyjoin).

## Syntax

```stata
* record linkage (reclink)
cfuzzyjoin link varlist using filename, idmaster(varname) idusing(varname) gen(newvar) ///
    [wmatch(numlist) wnomatch(numlist) orblock(varlist|none) required(varlist) ///
     exactstr(varlist) exclude(filename) minscore(#) minbigram(#) uvarlist(varlist) ///
     uprefix(string) _merge(newvar) debug(filename) verbose threads(#)]

* similarity join of two datasets (matchit)
cfuzzyjoin match idmaster txtmaster using filename, idusing(varname) txtusing(varname) ///
    [similmethod(simfunc) score(jaccard|simple|minsimple) weights(simple|root|log) ///
     wgtfile(filename) threshold(#) stopwordsauto swthreshold(#) diagnose ///
     generate(newvar) override time verbose threads(#)]

* similarity of two columns (matchit without using)
cfuzzyjoin match varname1 varname2 [, similmethod() score() weights() wgtfile() generate(newvar)]

* groups of similar strings (strgroup)
[by varlist:] cfuzzyjoin group varname [if] [in], generate(newvar) threshold(#) ///
    [normalize(shorter|longer|none) first noclean verbose threads(#)]
```

## link: record linkage

For each observation in memory (the master data), `cfuzzyjoin link` finds the
records of the using file most similar to it on the variables of *varlist*.
Each variable is scored from 0 to 1: strings with reclink's bigram comparator
(the share of the master string's bigrams found in the using string, with
Winkler's bonus for an agreeing start), numbers and the strings of
`exactstr()` as equal or not. The scores become agreement and disagreement
weights (`wmatch()`, `wnomatch()`) with a sharp change at `minbigram()`, and
the score of the pair is agreement / (agreement + disagreement). Records equal
on every variable are perfect matches (score 1) and take their using record
out of the search. The using records with the highest score, when it exceeds
`minscore()`, are the match, all of them when they tie.

The result keeps every master observation and adds:

| Variable | Content |
|----------|---------|
| `gen()` | The score (`float`, `%5.4f`); missing without a match |
| `idusing()` | The identifier of the matched using record |
| `U`*name* | Each matching variable of the using file, after its master counterpart (`uprefix()` changes the prefix) |
| other using variables | With their labels, formats, notes and characteristics, unless the master has a variable of that name |
| `_merge` | 3 matched, 1 not matched (`_merge()` renames it) |

The result is sorted by `idusing()`, with the unmatched master observations
last, as `reclink` leaves it. A master record that ties between several using
records appears once for each.

**Blocking** saves comparisons: `orblock()` compares a master record only with
using records that equal it on at least one of the listed variables (the
default when *varlist* has four or more variables; `orblock(none)` turns it
off), and `required()` only with those that equal it on all of the listed
variables. `exclude()` leaves out the records of earlier matches, so that a
fast, strict pass can be followed by a slower one on what is left.

## match: similarity joins

Each string becomes a vector of grams, chosen by `similmethod()`:

| `similmethod()` | Grams of "John Smith" |
|-----------------|-----------------------|
| `bigram` (default) | `Jo oh hn n␣ ␣S Sm mi it th` |
| `ngram, 3` | `Joh ohn hn␣ n␣S ␣Sm Smi mit ith` |
| `ngram_circ, 3` | the same plus `th␣ h␣J ␣Jo` (the end wraps to the start) |
| `token` | `John Smith` |
| `cotoken`, `scotoken` | pairs of consecutive tokens (`John Smith`), and with `scotoken` the tokens too |
| `firstgram, 3` | `Joh Smi` |
| `tkngram, 3` | the tokens and the 3-grams |
| `soundex`, `soundex_nara` | `J525` |
| `soundex_fk`, `soundex_ext` | the full soundex key `J5253`; its prefixes `J5 J52 J525 J5253` |
| `nysiis_fk` | `JANSNAT` |
| `token_soundex`, `tokenwrap, "f"` | a phonetic key of each token (`J500 S530`) |

The similarity of two vectors with dot product *N* and squared norms *D1*,
*D2* is `score(jaccard)` *N*/√(*D1·D2*) (default), `score(simple)`
2*N*/(*D1*+*D2*) or `score(minsimple)` *N*/min(*D1*,*D2*). `weights()` lets
rare grams count more: each gram is weighted by 1/*f*, 1/√*f* or
1/(log *f* + 1) of its frequency *f* in both string variables (or in
`wgtfile()`). `stopwordsauto` drops grams more frequent than `swthreshold()`
per observation.

With `using`, every pair sharing a gram and scoring at least `threshold()`
(default 0.5) is kept, and the pairs replace the data in memory:
*idmaster*, *txtmaster*, `idusing()`, `txtusing()` and the score
(`similscore`, or `generate()`), in the storage types `compress` leaves, as
`matchit` makes them. The pairs come in master order and, within a master
string, in using order. Without `using`, the score of the two variables of
each observation is added as a new `double` variable.

## group: groups of similar strings

Two strings match when their Levenshtein edit distance (bytes, unit costs),
divided by the length of the shorter string (`normalize(shorter)`, the
default), of the longer one (`longer`) or by nothing (`none`: the threshold
is a number of edits), is at most `threshold()`, which is used rounded to seven
decimal places, as `strgroup`'s plugin uses it. Groups are the sets of strings joined by
chains of matches, numbered from 1 by their first observation. With `by`, the
groups stay within by-groups and are numbered by-group after by-group; with
`first`, strings are also split by their first letter (other characters, then
a to z), as `strgroup` does. Leading and trailing blanks are removed first
unless `noclean` is given; missing (and blank) strings get a missing group.

## Examples

```stata
* Link two lists of firms on name, address and ZIP code
use firms_a, clear
cfuzzyjoin link name address zip using firms_b, idmaster(id) idusing(rid) gen(score) ///
    wmatch(10 8 2) wnomatch(8 4 1)

* Every pair of names at least 0.7 similar (bigrams), weighting rare grams more
use names_a, clear
cfuzzyjoin match id name using names_b, idusing(rid) txtusing(rname) threshold(0.7) weights(log)

* Similarity of two columns with trigrams
cfuzzyjoin match name1 name2, similmethod(ngram, 3) generate(sim3)

* Group spellings of firm names within states
bysort state: cfuzzyjoin group name, generate(firm) threshold(0.2)
```

## Results compared with the original commands

The validation suite (`validation/stata/cfuzzyjoin.do`) runs each original
command and its subcommand on the same data and compares everything they
leave.

- **link**: reclink's scores to the last bit (including reclink's rounding of
  each step to `float` and Stata's reading of `a*b/c` as `a*(b/c)`), the same
  matches, and the same variables in the same order, with the same types,
  formats, labels, notes and characteristics. Master strings reach the bigram
  comparator as reclink's `syntax` passes them: trimmed, with runs of blanks
  and tabs made one blank.
- **match**: the same pairs and scores as `matchit`, to the last bit without
  weights (every sum is an exact integer) and within 1e-13 (relative) with
  them (sums run over grams in another order), the same variable names (with
  matchit's renaming of clashing names) and storage types.
- **group**: strgroup's group numbers, with `by`, `if`/`in`, `first`, the
  three normalizations, long strings and boundary thresholds.

Results do not depend on the number of threads.

## Differences from the original commands

- **link**: when a using value combination shared by several using records is
  the best match of several master records, every master record gets every
  such record; `reclink`'s many-to-many `merge` pairs them off in an order
  that depends on the sort. Master values are compared as values, where
  `reclink` writes them into macros: master strings with double quotes work
  (`reclink` stops with r(132)), and equal numbers that are not integers agree
  (`reclink` compares a `float` 0.1 with the text of its value and scores the
  pair as disagreeing). A master record whose or-blocking
  strings are all empty has no candidates, where `reclink` stops with a syntax
  error. `wnomatch()` weights must be numbers. `cfuzzyjoin link` stores
  `r()` results.
- **match**: the using file may omit `.dta`. Gram frequencies missing from
  `wgtfile()` count as 1, also with `using`, where `matchit` stops with a Mata
  error. `diagnose` lists the grams both datasets have. Only the built-in
  similarity, score and weight functions exist (`matchit` also calls Mata
  functions a user adds). `flag()` and `keepmata` are accepted and do
  nothing. When `txtusing()` is named like the score variable, the
  score is named `similscore1`; `matchit` never returns (its renaming loop
  compares the requested name instead of the candidate).
- **group**: blanks are trimmed unless `noclean`, as `strgroup`'s help file
  describes; `strgroup` 1.0.4 never trims (its `noclean` has no effect). Blank
  strings get a missing group. There is no 10,000-observation limit (`force`
  is accepted) and no warning about duplicate values. When no observation is
  selected, the new variable is all missing (`strgroup` stops with r(198)).
  `strgroup`'s plugin reads only the first 244 bytes of each string;
  `cfuzzyjoin group` compares whole strings, so strings longer than that can
  be grouped differently.

Strings are compared up to 2,045 bytes; longer `strL` values are compared on
their first 2,045 bytes.

## Performance

`validation/performance/cfuzzyjoin.do` times each original command and its
subcommand on the same generated firm names with typos (two or three words
from a vocabulary of 25 syllables; 40% with a deleted, doubled or swapped
byte): `strgroup` with `threshold(0.2)`, `matchit` (two datasets) with
`threshold(0.6)`, and `reclink` on name, city and ZIP code with
`wmatch(10 2 4) wnomatch(8 1 2)` and no blocking. The original commands ran
once and `cfuzzyjoin` best of ten, on a Ryzen 9 9950X (StataNow 19.5 MP,
16 threads, the machine shared with other work). Every result was checked
against the original command's outside the timed region. Times below 10 ms are
close to the resolution of Stata's timer (1 ms).

| Original | Observations | Original | cfuzzyjoin | Speedup | Same result |
|---|---:|---:|---:|---:|:---:|
| `reclink` | 1,000 | 7.19 s | 0.015 s | 479× | yes |
| `reclink` | 2,000 | 27.33 s | 0.035 s | 781× | yes |
| `reclink` | 5,000 | 86.22 s | 0.160 s | 539× | yes |
| `matchit` | 1,000 | 19.77 s | 0.007 s | 2,825× | yes |
| `matchit` | 2,000 | 85.12 s | 0.009 s | 9,458× | yes |
| `matchit` | 5,000 | 469.59 s | 0.029 s | 16,193× | yes |
| `strgroup` | 1,000 | 0.08 s | 0.001 s | 84× | yes |
| `strgroup` | 2,000 | 0.40 s | 0.002 s | 200× | yes |
| `strgroup` | 5,000 | 2.67 s | 0.006 s | 446× | yes |

Observations are per dataset (master and using alike). The original commands
compare every pair in Stata or Mata, so their time grows with the square of
the data: at 50,000 observations each, `reclink` and `matchit` would take
hours by these timings, and `strgroup` would need 10 GB for its matrix of all
pairs. There `cfuzzyjoin group` takes 0.48 s, `cfuzzyjoin match` 2.2 s to
score 2.1 billion candidate pairs and keep the 11 million at or above 0.6,
and `cfuzzyjoin link` 1.2 s with `required(city)`.

Where the time goes:

- **link** indexes the distinct value combinations of the using variables and
  the blocking variables once, scores the candidates of each master record in
  parallel with the bigrams of every using string precomputed, and lays out
  the result in the plugin; Stata adds the variables (an empty copy of the
  using data is appended for their labels and formats) and sorts the rows.
- **match** turns every string into gram ids once (in parallel, numbered in
  order of first appearance so that sums do not depend on the thread count),
  indexes the using grams, and accumulates the dot products of each master
  string over the using strings that share a gram, in a dense array per
  thread; only the pairs kept are sorted.
- **group** collapses equal strings, sorts the rest by length, and compares a
  string only with longer ones until the length difference alone exceeds the
  allowed distance, skipping pairs already in one group and pairs ruled out by
  the characters only one of them has; the remaining distances use Myers's
  bit-parallel algorithm, and a lock-free union-find joins the groups.

To rerun the comparison:

```sh
python3 validation/bench.py do cfuzzyjoin_performance --output temp/timings-fuzzy "1000 2000 5000" "50000"
```
