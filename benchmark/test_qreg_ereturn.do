clear all
set more off
sysuse auto, clear

* Run qreg with default (fitted) method
qreg price mpg

* Show all stored results
ereturn list
