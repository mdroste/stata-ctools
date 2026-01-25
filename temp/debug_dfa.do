adopath ++ build
webuse nlswork, clear
xtset idcode year

* Check if nested detection works
civreghdfe ln_wage tenure (ttl_exp = age), absorb(idcode) vce(cluster idcode) verbose noid | head -30
