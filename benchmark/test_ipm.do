* Test IPM implementation
clear all
set more off
adopath + "../build"

sysuse auto, clear

di _newline
di as text "{hline 70}"
di as text "Testing IPM implementation: cqreg vs qreg (median regression)"
di as text "{hline 70}"

* Run qreg
di _newline
di as result "=== Stata's qreg ===" _newline
qreg price mpg weight, nolog
matrix b_qreg = e(b)

* Run cqreg
di _newline
di as result "=== ctools cqreg (IPM) ===" _newline
cqreg price mpg weight
matrix b_cqreg = e(b)

* Compare
di _newline
di as text "{hline 70}"
di as text "Coefficient Comparison:"
di as text "{hline 70}"
di as text "Variable     " _col(20) "qreg" _col(35) "cqreg" _col(50) "Diff"

local names : colnames b_qreg
local i = 1
local max_diff = 0
foreach v of local names {
    local diff = b_cqreg[1,`i'] - b_qreg[1,`i']
    local abs_diff = abs(`diff')
    if `abs_diff' > `max_diff' {
        local max_diff = `abs_diff'
    }
    di as text "`v'" _col(20) as result %12.6f b_qreg[1,`i'] _col(35) as result %12.6f b_cqreg[1,`i'] _col(50) as result %12.6f `diff'
    local i = `i' + 1
}

di _newline
di as text "{hline 70}"
di as text "Maximum absolute difference: " as result %g `max_diff'
if `max_diff' < 1e-4 {
    di as result "SUCCESS: Difference < 1e-4"
}
else if `max_diff' < 1e-2 {
    di as text "CLOSE: Difference < 1e-2 but >= 1e-4"
}
else {
    di as error "FAIL: Difference >= 1e-2"
}
di as text "{hline 70}"
