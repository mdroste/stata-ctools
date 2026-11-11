*! version 1.2.6 20261111 github.com/mdroste/stata-ctools
* estat after creghdfe, named in e(estat_cmd), with the subcommands of
* reghdfe's estat: summarize (without its header, as reghdfe), vce, and ic
* with e(df_m) + 1 degrees of freedom (the regressors and the constant);
* also bootstrap after vce(bootstrap).
program define _creghdfe_estat, rclass
    version 14.1
    if "`e(cmd)'" != "creghdfe" error 301
    gettoken key 0 : 0, parse(", ")
    local lkey = length(`"`key'"')
    if `"`key'"' == substr("summarize", 1, max(2, `lkey')) {
        syntax [anything] [, noheader *]
        * noheader works around estat_summ's header; it is always on
        estat_summ `anything', `options' noheader
    }
    else if `"`key'"' == "vce" {
        vce `0'
    }
    else if `"`key'"' == substr("bootstrap", 1, max(4, `lkey')) & "`e(vce)'" == "bootstrap" {
        * after vce(bootstrap) or the bootstrap prefix: Stata's normal,
        * percentile and bias-corrected intervals from the e() results
        estat_default bootstrap `0'
    }
    else if `"`key'"' == "ic" {
        syntax [, *]
        estat_default ic, df(`=e(df_m) + 1') `options'
    }
    else {
        di as error `"invalid subcommand `key'"'
        exit 321
    }
    return add
end
