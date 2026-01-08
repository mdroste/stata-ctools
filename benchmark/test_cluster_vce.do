clear all
set more off
adopath + "../build"
sysuse auto, clear
drop if missing(price, mpg, weight, rep78)

display _n "=== TESTING CLUSTER VCE ===" _n

* Check if qreg supports cluster VCE
capture qreg price mpg weight, quantile(0.5) vce(cluster rep78)
if _rc != 0 {
    display "qreg does NOT support vce(cluster), error code: " _rc
    display "This is expected - Stata's qreg doesn't have cluster VCE"
}
else {
    display "qreg cluster VCE succeeded (unexpected)"
    matrix list e(V)
}

* Try bsqreg (bootstrap quantile regression) as alternative
display _n "=== Testing bootstrap QR as reference ===" _n
display "(skipping bootstrap - takes too long for quick test)"

* Test cqreg cluster VCE
display _n "=== Testing cqreg vce(cluster rep78) ===" _n
capture cqreg price mpg weight, quantile(0.5) vce(cluster rep78)
if _rc == 0 {
    display "cqreg cluster VCE succeeded"
    display "SE[mpg]    = " %12.6f _se[mpg]
    display "SE[weight] = " %12.6f _se[weight]
    display "Number of clusters: " e(N_clust)
    matrix list e(V)
}
else {
    display "cqreg cluster VCE failed, error code: " _rc
}

* For reference, get IID standard errors
display _n "=== Reference: cqreg IID VCE ===" _n
quietly cqreg price mpg weight, quantile(0.5) vce(iid)
display "SE[mpg]    = " %12.6f _se[mpg]
display "SE[weight] = " %12.6f _se[weight]
