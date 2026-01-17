*! version 1.0.0 17Jan2026
*! ctools: C-accelerated tools for Stata
*! Main info program for the ctools suite

program define ctools
    version 14.0

    syntax [, Version ENVironment_check Verbose]

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
    di as text "{bf:ctools} - C-accelerated tools for Stata"
    di as text "{hline 50}"
    di as text "Version 1.1.0"
    di as text ""
    di as text "Available commands:"
    di as text "  {help csort:csort}      - Replaces -sort-"
    di as text "  {help cmerge:cmerge}     - Replaces -merge-"
    di as text "  {help cimport:cimport}    - Replaces -import delimited-"
    di as text "  {help cexport:cexport}    - Replaces -export delimited-"
    di as text "  {help creghdfe:creghdfe}   - Replaces -reghdfe-"
    di as text "  {help cbinscatter:cbinscatter}   - Replaces -binscatter-"
    di as text "  {help cqreg:cqreg}   - Replaces -qreg-"
    di as text ""
    di as text "Utilities:"
    di as text "  {cmd:ctools, environment_check} - Check plugin status"
    di as text ""
    di as text "For help on a specific command, type:"
    di as text "  {cmd:help} {it:command}"
    di as text ""
end
