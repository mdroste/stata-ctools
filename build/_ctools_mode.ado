*! version 1.2.0 20261029 github.com/mdroste/stata-ctools
*! Failure-handling mode in a ctools wrapper's options.
* Sets local mode in the caller to snapshot, fast or empty (the default: an
* exact undo without a snapshot of the dataset). Options are parsed after the
* first comma outside parentheses and quotes, so if-expressions such as
* inlist(g,1,2) are skipped. A command line syntax cannot parse leaves the
* default; the wrapper's own syntax then reports the error.
program define _ctools_mode
    version 14.1
    gettoken lhs 0 : 0, parse(",") bind
    capture syntax [, SNAPshot FAST *]
    if _rc {
        c_local mode
        exit
    }
    if "`snapshot'" != "" & "`fast'" != "" {
        di as error "options snapshot and fast may not be combined"
        exit 198
    }
    c_local mode `snapshot'`fast'
end
