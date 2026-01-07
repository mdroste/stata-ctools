clear all
adopath ++ "./build"
sysuse auto, clear

di "Testing csort..."
csort price
di "csort completed!"
