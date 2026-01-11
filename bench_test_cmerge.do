* Add ctools to adopath
adopath + "./build"

cap program drop bench_merge
program bench_merge
	syntax anything using/, [opts(string asis)]
	di "Benchmarking merge and cmerge with command: `anything'"
	preserve
	merge `anything' using `using', `opts'
	tempfile results1
	qui save `results1'
	restore
	cmerge `anything' using `using', `opts'
	tempfile results2
	qui save `results2'
end

webuse autosize, clear
bench_merge 1:1 make using https://www.stata-press.com/data/r18/autoexpense

/*
* 1:1 match merge
webuse autosize, clear
merge 1:1 make using https://www.stata-press.com/data/r18/autoexpense

* Perform 1:1 match merge, requiring there to be only matches
webuse autosize, clear
merge 1:1 make using https://www.stata-press.com/data/r18/autoexpense

* Perform 1:1 match merge, keeping only matches and squelching the _merge variable
webuse autosize, clear
merge 1:1 make using https://www.stata-press.com/data/r18/autoexpense, keep(match) nogen

* Perform m:1 match merge with sforce in memory
webuse sforce, clear
merge m:1 region using https://www.stata-press.com/data/r18/dollars

* Perform m:1 match merge, illustrating update option
webuse overlap1, clear
merge m:1 id using https://www.stata-press.com/data/r18/overlap2, update

* Perform m:1 match merge, illustrating update replace option
webuse overlap1, clear
merge m:1 id using https://www.stata-press.com/data/r18/overlap2, update replace

* Perform 1:m match merge, illustrating update replace option
webuse overlap2, clear
merge 1:m id using https://www.stata-press.com/data/r18/overlap1, update replace

* Sequential merge by observation
webuse sforce, clear
merge 1:1 _n using https://www.stata-press.com/data/r18/dollars
*/
