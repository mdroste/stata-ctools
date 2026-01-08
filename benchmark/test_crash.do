* Test FN solver crash at various sizes
clear all
set more off
adopath + "../build"

di "{hline 50}"
di "Testing FN Solver Stability"
di "{hline 50}"

foreach n in 1000 2000 3000 4000 5000 {
    clear
    set seed 12345
    set obs `n'
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen y = 1 + 2*x1 - x2 + rnormal()*2

    di _newline "Testing N=`n'..."

    * Try running cqreg
    capture noisily cqreg y x1 x2

    if _rc == 0 {
        di "  SUCCESS at N=`n'"
        matrix list e(b)
    }
    else {
        di "  FAILED at N=`n' with rc=" _rc
    }
}

di _newline "All tests complete"
