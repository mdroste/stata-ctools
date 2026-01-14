* Test multi-endogenous case with verbose
clear all
set more off
adopath ++ "build"

sysuse auto, clear
civreghdfe price (mpg weight = length turn displacement), absorb(foreign) verbose
di as text "underid_stat = " e(idstat)
di as text "cd_f = " e(cd_f)
