adopath ++ build
webuse nlswork, clear
xtset idcode year
civreghdfe ln_wage tenure (ttl_exp = age), absorb(idcode) vce(cluster idcode) verbose
