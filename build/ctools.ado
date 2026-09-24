*! version 1.0.2 20260920 github.com/mdroste/stata-ctools

program define ctools, rclass
    version 14.1

    syntax [, Version ENVironment_check Verbose UPDATE]

    if "`update'" != "" {
        di as text ""
        di as text "Checking for ctools updates..."
        di as text ""

        capture noisily net install ctools, from("https://raw.githubusercontent.com/mdroste/stata-ctools/main/build") replace

        if _rc == 0 {
            di as text ""
            di as text "{bf:ctools has been updated successfully.}"
            di as text "You may need to restart Stata for all changes to take effect."
        }
        else if _rc == 631 {
            di as text ""
            di as text "{bf:ctools is already up to date.}"
        }
        else {
            di as error ""
            di as error "Update failed with error code `=_rc'."
            di as error "Please check your internet connection and try again."
        }
        exit
    }

    local __ado_version "1.0.2"
    if "`environment_check'" != "" | "`version'" != "" {
        * Load the plugin if not already loaded
        _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0

        capture noisily plugin call ctools_plugin, "version"
        local plugin_rc = _rc
        if `plugin_rc' == 0 {
            di as text "ctools ado `__ado_version'; plugin `ctools_plugin_version' (`ctools_plugin_revision')"
            return local version "`__ado_version'"
            return local plugin_version "`ctools_plugin_version'"
            return local plugin_revision "`ctools_plugin_revision'"
        }
        else {
            di as error "ctools plugin identity could not be read."
            exit `plugin_rc'
        }
        exit
    }

    di as text ""
    di as result " ██████╗████████╗ ██████╗  ██████╗ ██╗     ███████╗"
    di as result "██╔════╝╚══██╔══╝██╔═══██╗██╔═══██╗██║     ██╔════╝"
    di as result "██║        ██║   ██║   ██║██║   ██║██║     ███████╗"
    di as result "██║        ██║   ██║   ██║██║   ██║██║     ╚════██║"
    di as result "╚██████╗   ██║   ╚██████╔╝╚██████╔╝███████╗███████║"
    di as result " ╚═════╝   ╚═╝    ╚═════╝  ╚═════╝ ╚══════╝╚══════╝"
    di as text ""
    di as text "  {it:C-accelerated tools for Stata}{col 56}v1.0.2"
    di as text "{hline 60}"
    di as text ""
    di as text "{ul:Data Management}"
    di as text "  {help cimport:cimport}      Import delimited/Excel data"
    di as text "  {help cexport:cexport}      Export delimited/Excel data"
    di as text "  {help csort:csort}        Sort dataset"
    di as text "  {help cmerge:cmerge}       Merge datasets"
    di as text "  {help csample:csample}      Sample without replacement"
    di as text "  {help cbsample:cbsample}     Sample with replacement"
    di as text "  {help cencode:cencode}      String to labeled numeric"
    di as text "  {help cdecode:cdecode}      Labeled numeric to string"
    di as text "  {help cdestring:cdestring}    String to numeric"
    di as text "  {help cwinsor:cwinsor}      Winsorize variables"
    di as text "  {help crangestat:crangestat}   Range statistics"
    di as text ""
    di as text "{ul:Estimation}"
    di as text "  {help creghdfe:creghdfe}     OLS with multi-way FE"
    di as text "  {help civreghdfe:civreghdfe}   2SLS/GMM with multi-way FE"
    di as text "  {help cqreg:cqreg}        Quantile regression"
    di as text "  {help cpplmhdfe:cpplmhdfe}    PPML with multi-way FE"
    di as text "  {help cpsmatch:cpsmatch}     Propensity score matching"
    di as text ""
    di as text "{ul:Visualization}"
    di as text "  {help cbinscatter:cbinscatter}  Binned scatter plots"
    di as text ""
    di as text "{ul:Utilities}"
    di as text "  {cmd:ctools, environment_check} - Check plugin status"
    di as text "  {cmd:ctools, update}            - Update to latest version"
    di as text ""
    di as text "For help: {cmd:help} {it:command}    Source: {browse github.com/mdroste/stata-ctools}"
    di as text "{hline 60}"
end
