* Test cqreg vs qreg
* Simple benchmark on sysuse auto

clear all
set more off

* Add build directory to adopath
adopath + "../build"

* Load test data
sysuse auto, clear

* Display data info
di "Number of observations: " _N
summarize price mpg weight

* Run native qreg first
di ""
di "============================================"
di "Running native qreg..."
di "============================================"
qreg price mpg weight
matrix b_qreg = e(b)
matrix V_qreg = e(V)
di "qreg coefficients:"
matrix list b_qreg
di "qreg VCE diagonal:"
matrix list V_qreg

* Now try cqreg
di ""
di "============================================"
di "Running cqreg..."
di "============================================"

cqreg price mpg weight, verbose

di "cqreg rc = " _rc

if _rc == 0 {
    di ""
    di "cqreg coefficients:"
    matrix list e(b)

    di ""
    di "cqreg VCE:"
    matrix list e(V)

    di ""
    di "Comparison of coefficients:"
    di "Variable     qreg          cqreg         diff"
    di "------------------------------------------------"
    forval i = 1/3 {
        local diff = b_qreg[1,`i'] - e(b)[1,`i']
        di "coef `i':     " %12.4f b_qreg[1,`i'] "  " %12.4f e(b)[1,`i'] "  " %12.6f `diff'
    }
}
else {
    di as error "cqreg failed with error code: " _rc
}

di ""
di "Test complete."
