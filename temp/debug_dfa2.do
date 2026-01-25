adopath ++ build
webuse nlswork, clear
xtset idcode year
civreghdfe ln_wage tenure (ttl_exp = age), absorb(idcode) vce(cluster idcode) verbose noid
di ""
di "After civreghdfe:"
di "e(df_a) = " e(df_a)
di "e(df_r) = " e(df_r)
di "e(N) = " e(N)
di "e(N_clust) = " e(N_clust)
di "e(N_hdfe) = " e(N_hdfe)
