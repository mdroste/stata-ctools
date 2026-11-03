*! version 1.2.4 20261103 github.com/mdroste/stata-ctools
* estat after cppmlhdfe, named in e(estat_cmd), with the subcommands and
* conventions of ppmlhdfe's (reghdfe's) estat: summarize, vce, and ic with
* e(df_m) + 1 degrees of freedom.
program define _cppmlhdfe_estat, rclass
    version 14.1
    if "`e(cmd)'" != "cppmlhdfe" error 301
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
