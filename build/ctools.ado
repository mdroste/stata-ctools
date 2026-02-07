*! version 0.9.1 06Feb2026
*! ctools: C-accelerated tools for Stata
*! Main info program for the ctools suite

program define ctools
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

    if "`environment_check'" != "" {
        * Load the plugin if not already loaded
        capture program list ctools_plugin
        if _rc != 0 {
            local __os = c(os)
            local __machine = c(machine_type)
            local __is_mac = 0
            if "`__os'" == "MacOSX" {
                local __is_mac = 1
            }
            else if strpos(lower("`__machine'"), "mac") > 0 {
                local __is_mac = 1
            }

            local __plugin = ""
            if "`__os'" == "Windows" {
                local __plugin "ctools_windows.plugin"
            }
            else if `__is_mac' {
                local __is_arm = 0
                if strpos(lower("`__machine'"), "apple") > 0 | strpos(lower("`__machine'"), "arm") > 0 | strpos(lower("`__machine'"), "silicon") > 0 {
                    local __is_arm = 1
                }
                if `__is_arm' == 0 {
                    tempfile __archfile
                    quietly shell uname -m > "`__archfile'" 2>&1
                    tempname __fh
                    file open `__fh' using "`__archfile'", read text
                    file read `__fh' __archline
                    file close `__fh'
                    capture erase "`__archfile'"
                    if strpos("`__archline'", "arm64") > 0 {
                        local __is_arm = 1
                    }
                }
                if `__is_arm' {
                    local __plugin "ctools_mac_arm.plugin"
                }
                else {
                    local __plugin "ctools_mac_x86.plugin"
                }
            }
            else if "`__os'" == "Unix" {
                local __plugin "ctools_linux.plugin"
            }
            else {
                local __plugin "ctools.plugin"
            }

            capture program ctools_plugin, plugin using("`__plugin'")
            if _rc != 0 & _rc != 110 & "`__plugin'" != "ctools.plugin" {
                capture program ctools_plugin, plugin using("ctools.plugin")
            }
        }

        capture program list ctools_plugin
        if _rc == 0 {
            di as text "ctools plugin is loaded and ready."
        }
        else {
            di as error "ctools plugin could NOT be loaded."
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
    di as text "  {it:C-accelerated tools for Stata}{col 56}v0.9.1"
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
