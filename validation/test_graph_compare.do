* Test script to visually compare cbinscatter vs binscatter graphs
* Run this to generate side-by-side comparisons

clear all
set more off

* Add ctools to adopath
adopath + build

* Load test data
sysuse auto, clear

* Test 1: Basic comparison (single series with linear fit)
* This should show navy markers and maroon line for both
di _n "{hline 60}"
di "Test 1: Basic scatter - comparing graph appearance"
di "{hline 60}"

* binscatter
binscatter price mpg, nquantiles(10) name(bs_basic, replace) title("binscatter: price vs mpg")

* cbinscatter
cbinscatter price mpg, nquantiles(10) name(cbs_basic, replace) title("cbinscatter: price vs mpg")

graph combine bs_basic cbs_basic, rows(1) name(compare_basic, replace)

* Test 2: Multiple groups (by foreign)
di _n "{hline 60}"
di "Test 2: By-groups scatter - comparing graph appearance"
di "{hline 60}"

* binscatter
binscatter price mpg, nquantiles(10) by(foreign) name(bs_by, replace) title("binscatter: by(foreign)")

* cbinscatter
cbinscatter price mpg, nquantiles(10) by(foreign) name(cbs_by, replace) title("cbinscatter: by(foreign)")

graph combine bs_by cbs_by, rows(1) name(compare_by, replace)

* Test 3: No fit line
di _n "{hline 60}"
di "Test 3: No fit line - comparing graph appearance"
di "{hline 60}"

* binscatter
binscatter price mpg, nquantiles(10) linetype(none) name(bs_nofit, replace) title("binscatter: linetype(none)")

* cbinscatter
cbinscatter price mpg, nquantiles(10) linetype(none) name(cbs_nofit, replace) title("cbinscatter: linetype(none)")

graph combine bs_nofit cbs_nofit, rows(1) name(compare_nofit, replace)

* Test 4: Quadratic fit
di _n "{hline 60}"
di "Test 4: Quadratic fit - comparing graph appearance"
di "{hline 60}"

* binscatter
binscatter price mpg, nquantiles(10) linetype(qfit) name(bs_qfit, replace) title("binscatter: linetype(qfit)")

* cbinscatter
cbinscatter price mpg, nquantiles(10) linetype(qfit) name(cbs_qfit, replace) title("cbinscatter: linetype(qfit)")

graph combine bs_qfit cbs_qfit, rows(1) name(compare_qfit, replace)

* Test 5: Controls
di _n "{hline 60}"
di "Test 5: With controls - comparing graph appearance"
di "{hline 60}"

* binscatter
binscatter price mpg, nquantiles(10) controls(weight) name(bs_ctrl, replace) title("binscatter: controls(weight)")

* cbinscatter
cbinscatter price mpg, nquantiles(10) controls(weight) name(cbs_ctrl, replace) title("cbinscatter: controls(weight)")

graph combine bs_ctrl cbs_ctrl, rows(1) name(compare_ctrl, replace)

di _n "{hline 60}"
di "Visual comparison complete. Check the combined graphs."
di "{hline 60}"
