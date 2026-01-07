* Quick test of csort
clear all
set more off
adopath + "../build"

sysuse auto, clear
list price in 1/5
csort price
list price in 1/5
