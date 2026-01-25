adopath ++ build
sysuse auto, clear

* ivreghdfe with cluster
ivreghdfe price (mpg = weight), absorb(foreign) vce(cluster foreign)
di "ivreghdfe: r2=" e(r2) ", r2_a=" e(r2_a) ", N=" e(N) ", df_a=" e(df_a)
di "ivreghdfe: sdofminus=" e(sdofminus) ", dofminus=" e(dofminus)

* civreghdfe with cluster
sysuse auto, clear
civreghdfe price (mpg = weight), absorb(foreign) vce(cluster foreign) noid
di ""
di "civreghdfe: r2=" e(r2) ", r2_a=" e(r2_a) ", N=" e(N) ", df_a=" e(df_a)

* For nested FE in cluster, ivreghdfe uses df_a=0 but sdofminus should still be 1
* Let me check what sdofminus ivreghdfe is using
di ""
di "For nested FE (absorb=cluster), ivreghdfe uses:"
di "  HDFE.df_a = 0 (nested), absorb_ct = max(1,0) = 1"
di "  sdofminus = absorb_ct = 1"
di "  Formula: r2_a = 1 - (1-r2)*N/(N-K-sdofminus) = 1 - (1-0.1386)*74/(74-1-1) = " ///
   1 - (1-0.1386)*74/(74-1-1)
