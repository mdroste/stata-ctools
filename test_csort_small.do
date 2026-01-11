adopath + "./build"

* Use small dataset
webuse autosize, clear
desc, short

di "Calling csort..."
csort make
di "csort completed!"
list
