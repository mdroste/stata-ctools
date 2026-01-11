adopath + "./build"
webuse autosize, clear
di "Master data loaded:"
desc, short
list make in 1/6
di "Now calling cmerge..."
cmerge 1:1 make using https://www.stata-press.com/data/r18/autoexpense, verbose
di "cmerge completed!"
