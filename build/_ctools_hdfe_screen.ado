*! version 1.3.0 20261124 github.com/mdroste/stata-ctools
*! Native sample preparation before factor expansion or weight staging
program define _ctools_hdfe_screen
    version 14.1
    gettoken action 0 : 0
    if "`action'" == "plan" {
        syntax anything(name=model)
        * list strips factor operators without choosing levels or rejecting
        * invalid categories on observations that will be excluded later.
        * Keep time-series operators on the original marksample path. Remove
        * parentheses first so ib(first), ib(last), ib(freq) remain eligible.
        local ops = ustrregexra(`" `model'"', "\([^)]*\)", "")
        local ok = !regexm(`"`ops'"', "[ #][^ .#]*[LlFfDdSs][^ .#]*[.]")
        local inputs
        if `ok' {
            tempname rhold
            _return hold `rhold'
            capture fvrevar `model', list
            local ok = (_rc == 0)
            if `ok' {
                local inputs `r(varlist)'
                capture confirm numeric variable `inputs'
                local ok = (_rc == 0)
            }
            _return restore `rhold'
        }
        c_local __screen_ok = `ok'
        c_local __screen_inputs `inputs'
        exit
    }
    if "`action'" != "mark" error 198
    syntax varlist(numeric), TOUSE(varname numeric) PLUGIN(string) [THREADS(integer 0)]
    local inputs : list uniq varlist
    capture program ctools_plugin, plugin using("`plugin'")
    if _rc != 0 & _rc != 110 exit 601
    local threadopt
    if `threads' > 0 local threadopt threads(`threads')
    plugin call ctools_plugin `inputs' `touse', "hdfe_markout `threadopt'"
end
