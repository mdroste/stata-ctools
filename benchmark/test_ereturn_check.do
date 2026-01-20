* Check what e() values ivreghdfe returns

clear all
set more off
adopath + "build"

sysuse auto, clear

di _n "{hline 78}"
di "Checking ivreghdfe e() returns with robust VCE"
di "{hline 78}" _n

ivreghdfe price (mpg = weight length), absorb(foreign) vce(robust)

di _n "All stored e() scalars:"
ereturn list

di _n _n "Specific Sargan/Hansen related:"
di "e(sargan):   " e(sargan)
di "e(sargan_p): " e(sargan_p)
di "e(j):        " e(j)
di "e(jp):       " e(jp)
di "e(hansen):   " e(hansen)

di _n _n "{hline 78}"
di "Now checking civreghdfe e() returns"
di "{hline 78}" _n

civreghdfe price (mpg = weight length), absorb(foreign) vce(robust)

di _n "All stored e() scalars:"
ereturn list

di _n _n "Specific Sargan/Hansen related:"
di "e(sargan):   " e(sargan)
di "e(sargan_p): " e(sargan_p)
di "e(sargan_df):" e(sargan_df)
