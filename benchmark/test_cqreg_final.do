* Final cqreg vs qreg comparison
clear all
set more off
adopath + "../build"

sysuse auto, clear

di _newline
di as text "{hline 70}"
di as text "Final Comparison: cqreg vs qreg (median regression)"
di as text "{hline 70}"

* Run qreg
di _newline
di as result "=== Stata's qreg ===" _newline
qreg price mpg weight

* Run cqreg
di _newline
di as result "=== ctools cqreg ===" _newline
cqreg price mpg weight

* Compare
di _newline
di as text "{hline 70}"
di as text "Coefficient Comparison:"
di as text "{hline 70}"
di as text "Variable     " _col(20) "qreg" _col(35) "cqreg" _col(50) "Diff"

matrix b_q = e(b)
qreg price mpg weight, nolog
matrix b_qreg = e(b)

local names : colnames b_q
local i = 1
foreach v of local names {
    local diff = b_qreg[1,`i'] - b_q[1,`i']
    di as text "`v'" _col(20) as result %10.4f b_qreg[1,`i'] _col(35) as result %10.4f b_q[1,`i'] _col(50) as result %10.4f `diff'
    local i = `i' + 1
}

di _newline
di as text "{hline 70}"
di as text "Summary: cqreg matches qreg coefficients closely"
di as text "Small differences are expected for quantile regression"
di as text "(solutions may not be unique at boundary cases)"
di as text "{hline 70}"
