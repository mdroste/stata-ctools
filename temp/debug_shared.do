adopath ++ "./build"
clear
set obs 10
gen str10 var_a = "cat" + string(mod(_n, 3))
gen str10 var_b = "cat" + string(mod(_n, 3))
capture label drop shared
capture cencode var_a, generate(a_code) label(shared)
di "cencode var_a rc = `=_rc'"
capture cencode var_b, generate(b_code) label(shared)
di "cencode var_b rc = `=_rc'"
capture drop a_code b_code
capture label drop shared
capture encode var_a, generate(a_code) label(shared)
di "encode var_a rc = `=_rc'"
capture encode var_b, generate(b_code) label(shared)
di "encode var_b rc = `=_rc'"
