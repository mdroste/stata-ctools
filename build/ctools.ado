*! version 1.3.2 20261124 github.com/mdroste/stata-ctools

program define ctools, rclass
    version 14.1

    syntax [, Version ENVironment_check Verbose UPDATE FORCE THReads(string) PERManently]

    if "`permanently'" != "" & `"`threads'"' == "" {
        di as error "option permanently requires threads()"
        exit 198
    }
    if "`force'" != "" & "`update'" == "" {
        di as error "option force requires update"
        exit 198
    }

    local __ado_version "1.3.2"
    if "`update'" != "" {
        * GitHub's latest release excludes drafts and prereleases. Each
        * release tag holds the installable package in build/.
        local __api "https://api.github.com/repos/mdroste/stata-ctools/releases/latest"
        local __site "https://raw.githubusercontent.com/mdroste/stata-ctools"
        di as text ""
        di as text "Checking GitHub for the latest ctools release..."
        tempfile __reply
        capture copy "`__api'" "`__reply'"
        local __rc = _rc
        if `__rc' {
            * r(2) after a timeout and r(631) for an unknown host mean GitHub
            * could not be reached; r(601) includes GitHub refusing the request.
            di as error "Could not read the latest ctools release from GitHub (error code `__rc')."
            di as error "Check your internet connection and try again. GitHub limits how often one"
            di as error "network address may ask, so a refused check may work again within an hour."
            exit `__rc'
        }
        tempname __json
        scalar `__json' = fileread("`__reply'")
        local __tag
        if ustrregexm(`__json', `""tag_name"[ ]*:[ ]*"([^"]*)""') local __tag = ustrregexs(1)
        if !ustrregexm(`"`__tag'"', "^[vV]?[0-9]+\.[0-9]+\.[0-9]+$") {
            di as error "GitHub's reply did not name a ctools release version."
            exit 672
        }
        local __latest = ustrregexra("`__tag'", "^[vV]", "")

        * Compare major, minor and patch; the most significant difference wins.
        local __have = subinstr("`__ado_version'", ".", " ", .)
        local __want = subinstr("`__latest'", ".", " ", .)
        local __newer 0
        forvalues __i = 3(-1)1 {
            local __a : word `__i' of `__have'
            local __b : word `__i' of `__want'
            if `__b' > `__a' local __newer 1
            else if `__b' < `__a' local __newer -1
        }
        di as text "  Installed version" _col(22) as result "`__ado_version'"
        di as text "  Latest release" _col(22) as result "`__latest'"
        local __install = `__newer' > 0 | "`force'" != ""
        if !`__install' & `__newer' == 0 di as text "ctools is up to date."
        else if !`__install' {
            di as text "The installed ctools is newer than the latest release; nothing was installed."
            di as text "Type {cmd:ctools, update force} to install release `__latest' instead."
        }
        else {
            di as text ""
            di as text "Installing ctools `__latest'..."
            local __from "`__site'/`__tag'/build"
            capture noisily net install ctools, from("`__from'") replace
            local __rc = _rc
            if `__rc' {
                * r(631) is net's "host not found".
                di as error ""
                di as error "Update failed with error code `__rc'."
                di as error "Please check your internet connection and try again."
                exit `__rc'
            }

            * net keeps one package record per source address, so earlier
            * installs from GitHub (main or another tag) stay listed by ado
            * dir and make ado uninstall ctools ambiguous. Their records are
            * removed once the new record exists. ado uninstall keeps files
            * that another record lists, so only files a later release
            * dropped are deleted. Without the new record (net found every
            * file up to date and recorded nothing), nothing is removed.
            * Record lines are read through Mata so that quotes in package
            * descriptions are never parsed as Stata syntax. net records a
            * local folder with backslashes, so sources compare with slashes.
            local __stale
            local __found 0
            local __trk "`c(sysdir_plus)'stata.trk"
            capture confirm file "`__trk'"
            if !_rc {
                local __n 0
                tempname __fh
                file open `__fh' using "`__trk'", read text
                file read `__fh' __line
                while r(eof) == 0 {
                    mata: st_local("__key", substr(st_local("__line"), 1, 2))
                    if "`__key'" == "S " {
                        local ++__n
                        mata: st_local("__src", subinstr(substr(st_local("__line"), 3, .), char(92), "/"))
                    }
                    else if "`__key'" == "N " {
                        mata: st_local("__ours", strofreal(strtrim(st_local("__line")) == "N ctools.pkg" & strpos(st_local("__src"), "mdroste/stata-ctools") > 0))
                        mata: st_local("__this", strofreal(st_local("__src") == subinstr(st_local("__from"), char(92), "/")))
                        if `__ours' & `__this' local __found 1
                        else if `__ours' local __stale `__n' `__stale'
                    }
                    file read `__fh' __line
                }
                file close `__fh'
            }
            * Highest number first, so the numbers still to remove stay valid.
            if `__found' {
                foreach __i of local __stale {
                    capture ado uninstall [`__i']
                }
                if "`__stale'" != "" di as text "(records of earlier ctools installs from GitHub were removed from {cmd:ado dir})"
            }
            di as text ""
            di as text "{bf:ctools `__latest' has been installed.}"
            di as text "Restart Stata so the new plugin is loaded."
            di as text `"Release notes: {browse "https://github.com/mdroste/stata-ctools/releases/tag/`__tag'"}"'
        }
        * Returned last: commands run above may replace r().
        return local version "`__ado_version'"
        return local latest_version "`__latest'"
        return scalar updated = `__install'
        exit
    }

    if "`environment_check'" != "" | "`version'" != "" | `"`threads'"' != "" {
        * Load the plugin if not already loaded
        _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0

        if `"`threads'"' != "" {
            * threads(#) sets the session setting; permanently also saves it
            * in a per-user file that later Stata sessions read.
            local __request = strtrim(`"`threads'"')
            if lower(`"`__request'"') == "default" local __n 0
            else {
                capture confirm integer number `__request'
                if _rc | `: word count `__request'' != 1 {
                    di as error "threads() must be an integer from 1 through 256, or default"
                    exit 198
                }
                local __n = real("`__request'")
                if !inrange(`__n', 1, 256) {
                    di as error "threads() must be an integer from 1 through 256, or default"
                    exit 198
                }
            }
            if "`permanently'" != "" {
                capture noisily plugin call ctools_plugin, "threads_save `__n'"
                if _rc exit _rc
            }
            if `__n' > 0 global CTOOLS_THREADS `__n'
            else global CTOOLS_THREADS
        }

        if "`environment_check'" != "" | "`version'" != "" {
            capture noisily plugin call ctools_plugin, "version"
            local plugin_rc = _rc
            if `plugin_rc' == 0 {
                di as text "ctools ado `__ado_version'; plugin `ctools_plugin_version' (`ctools_plugin_revision'); protocol `ctools_plugin_protocol'"
                return local version "`__ado_version'"
                return local plugin_version "`ctools_plugin_version'"
                return local plugin_revision "`ctools_plugin_revision'"
                return local plugin_protocol "`ctools_plugin_protocol'"
            }
            else {
                di as error "ctools plugin identity could not be read."
                exit `plugin_rc'
            }
        }
        if "`environment_check'" == "" & `"`threads'"' == "" exit

        capture noisily plugin call ctools_plugin, "threads_info"
        if _rc exit _rc
        local __source "`ctools_threads_source'"
        if "`ctools_threads_source'" == "automatic" {
            local __source "automatic: usable CPUs, at most `ctools_threads_cap'"
        }
        if "`ctools_threads_session'" != "" & "`ctools_threads_session'" != "invalid" {
            local __source "session setting"
        }
        * An invalid session global makes every command fail, so no thread
        * count applies; r(threads_default) still reports the policy default.
        local __effective_show "`ctools_threads_effective'"
        if "`ctools_threads_session'" == "invalid" {
            local ctools_threads_effective .
            local __source "global CTOOLS_THREADS is invalid"
            local __effective_show "none (commands fail)"
        }
        foreach __item in session saved env omp omp_limit {
            local __show_`__item' "`ctools_threads_`__item''"
            if "`__show_`__item''" == "" local __show_`__item' "(not set)"
        }
        local __openmp = cond("`ctools_threads_openmp'" == "1", "enabled", "not available")
        di as text ""
        di as text "{ul:Threads}"
        di as text "  Commands without threads() use" _col(40) as result "`__effective_show'" as text "  (`__source')"
        di as text "  CPUs detected / usable" _col(40) as result "`ctools_threads_detected' / `ctools_threads_available'"
        di as text "  Session (global CTOOLS_THREADS)" _col(40) as result "`__show_session'"
        di as text "  Saved setting" _col(40) as result "`__show_saved'"
        if `"`ctools_threads_file'"' != "" {
            di as text "  Settings file" _col(40) as result `"`ctools_threads_file'"'
        }
        di as text "  CTOOLS_NUM_THREADS" _col(40) as result "`__show_env'"
        di as text "  OMP_NUM_THREADS" _col(40) as result "`__show_omp'"
        di as text "  OMP_THREAD_LIMIT" _col(40) as result "`__show_omp_limit'"
        di as text "  OpenMP" _col(40) as result "`__openmp'"
        if "`ctools_threads_session'" == "invalid" {
            di as error "  global CTOOLS_THREADS is invalid; ctools commands will fail until"
            di as error "  you run {cmd:ctools, threads(default)} or set it to an integer from 1 through 256"
        }
        if `"`threads'"' != "" {
            if `__n' > `ctools_threads_limit' {
                di as text "  Note: only `ctools_threads_limit' CPUs are usable here, so ctools uses at most `ctools_threads_limit' threads."
            }
            if "`permanently'" != "" & "`ctools_threads_env'" != "" {
                di as text "  Note: CTOOLS_NUM_THREADS takes precedence over the saved setting in new sessions."
            }
        }

        return scalar threads = `ctools_threads_effective'
        return scalar threads_default = `ctools_threads_default'
        return scalar threads_cap = `ctools_threads_cap'
        return scalar threads_limit = `ctools_threads_limit'
        return scalar cpus_detected = `ctools_threads_detected'
        return scalar cpus_available = `ctools_threads_available'
        return local threads_source "`__source'"
        return local threads_session "`ctools_threads_session'"
        return local threads_saved "`ctools_threads_saved'"
        return local threads_file `"`ctools_threads_file'"'
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
    di as text "  {it:C-accelerated tools for Stata}{col 56}v1.3.2"
    di as text "{hline 60}"
    di as text ""
    di as text "{ul:Data Management}"
    di as text "  {help cimport:cimport}      Import delimited/Excel data"
    di as text "  {help cexport:cexport}      Export delimited/Excel data"
    di as text "  {help csql:csql}         SQL on data, Parquet, CSV and .dta files"
    di as text "  {help cdistinct:cdistinct}    Report distinct values"
    di as text "  {help cegen:cegen}        Extensions to generate"
    di as text "  {help ccollapse:ccollapse}    Collapse to summary statistics"
    di as text "  {help chdfe:chdfe}        Remove fixed effects from variables"
    di as text "  {help csort:csort}        Sort dataset"
    di as text "  {help cmerge:cmerge}       Merge datasets"
    di as text "  {help cjoinby:cjoinby}      Join all pairs within groups"
    di as text "  {help creshape:creshape}     Reshape wide to long and back"
    di as text "  {help cfuzzyjoin:cfuzzyjoin}   Fuzzy joins and string groups"
    di as text "  {help csample:csample}      Sample without replacement"
    di as text "  {help cbsample:cbsample}     Sample with replacement"
    di as text "  {help cencode:cencode}      String to labeled numeric"
    di as text "  {help cdecode:cdecode}      Labeled numeric to string"
    di as text "  {help cdestring:cdestring}    String to numeric"
    di as text "  {help cwinsor:cwinsor}      Winsorize variables"
    di as text "  {help crangestat:crangestat}   Range statistics"
    di as text ""
    di as text "  {help cinterpolate:cinterpolate} Linear and cubic interpolation"
    di as text "  {help csplit:csplit}       Split strings into fields"
    di as text "  {help crangejoin:crangejoin}   Interval joins"
    di as text "  {help cgeonear:cgeonear}     Nearest neighbors and distances"
    di as text ""
    di as text "{ul:Estimation}"
    di as text "  {help creghdfe:creghdfe}     OLS with multi-way FE"
    di as text "  {help civreghdfe:civreghdfe}   2SLS/GMM with multi-way FE"
    di as text "  {help cqreg:cqreg}        Quantile regression"
    di as text "  {help cppmlhdfe:cppmlhdfe}    PPML with multi-way FE"
    di as text "  {help cpsmatch:cpsmatch}     Propensity score matching"
    di as text ""
    di as text "{ul:Visualization}"
    di as text "  {help cbinscatter:cbinscatter}  Binned scatter plots"
    di as text ""
    di as text "{ul:Utilities}"
    di as text "  {cmd:ctools, environment_check} - Check plugin status and threads"
    di as text "  {cmd:ctools, threads(#)}        - Set threads ({cmd:permanently} saves it)"
    di as text "  {cmd:ctools, update}            - Install the latest GitHub release"
    di as text ""
    di as text "For help: {cmd:help} {it:command}    Source: {browse github.com/mdroste/stata-ctools}"
    di as text "{hline 60}"
end

* ctools-helpers:begin (generated by scripts/sync_release.py from src/ado; edit those files)
* Shared plugin loader and identity contract
program define _ctools_load
    version 14.1
    * Not rclass, and findfile's r() is undone, so the caller's r() survives;
    * identity comes back through c_local.
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
        if `__is_arm' == 0 & !strpos(lower("`__machine'"), "intel") & !strpos(lower("`__machine'"), "x86") {
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

    * The ado/plugin argument protocol, stamped from config/release.json.
    local __ado_protocol "4"
    local __reinstall `"net install ctools, from("https://raw.githubusercontent.com/mdroste/stata-ctools/main/build") replace"'

    tempname __rhold
    _return hold `__rhold'
    capture findfile `__plugin'
    local __found = (_rc == 0)
    local __path `"`r(fn)'"'
    _return restore `__rhold'
    local __wanted "`__plugin'"
    if !`__found' & "`__plugin'" != "ctools.plugin" local __plugin "ctools.plugin"
    capture program ctools_plugin, plugin using("`__plugin'")
    local __load_rc = _rc
    if `__load_rc' != 0 & `__load_rc' != 110 {
        di as error "ctools: could not load the plugin `__plugin'"
        di as error `"  Stata reports c(os) "`__os'" and c(machine_type) "`__machine'"; ctools selected `__wanted'."'
        if `__found' di as error `"  The plugin file is `__path'."'
        else di as error "  `__wanted' was not found on the ado-path; reinstall ctools."
        di as error "  Stata's plugin loader reported:"
        capture noisily program ctools_plugin, plugin using("`__plugin'")
        if "`__os'" == "Unix" & !`__is_mac' {
            di as error "  The distribution Linux plugin needs x86-64 and glibc 2.17 or later; OpenMP is linked statically."
            di as error "  Run ldd -r on the plugin to inspect missing libraries or symbols."
        }
        else if `__is_mac' {
            di as error "  On macOS, check that the plugin is the notarized ctools release and is not quarantined."
        }
        di as error "  See https://github.com/mdroste/stata-ctools/blob/main/docs/INSTALLATION.md#troubleshooting"
        exit `__load_rc'
    }
    local ctools_plugin_hdfe_screen
    local ctools_plugin_hdfe_markout
    local ctools_plugin_hdfe_fvlevels
    capture plugin call ctools_plugin, "version"
    local __id_rc = _rc
    if `__id_rc' | "`ctools_plugin_protocol'" != "`__ado_protocol'" {
        local __has = cond("`ctools_plugin_protocol'" == "", "none", "`ctools_plugin_protocol'")
        di as error "ctools: the plugin `__plugin' does not match the installed ado files."
        local __ident = cond("`ctools_plugin_version'" == "", "", " (plugin version `ctools_plugin_version', revision `ctools_plugin_revision')")
        di as error "  The ado files need plugin protocol `__ado_protocol'; the loaded plugin reports `__has'`__ident'."
        di as error "  Reinstall ctools so its ado files and plugins come from the same release:"
        di as error `"    `__reinstall'"'
        di as error "  Then restart Stata, because a plugin stays loaded until Stata exits."
        exit 9
    }
    c_local __ctools_plugin "`__plugin'"
    c_local __ctools_plugin_version "`ctools_plugin_version'"
    c_local __ctools_plugin_revision "`ctools_plugin_revision'"
    c_local __ctools_plugin_protocol "`ctools_plugin_protocol'"
    c_local __ctools_plugin_hdfe_screen = ("`ctools_plugin_hdfe_screen'" == "1")
    c_local __ctools_plugin_hdfe_markout = ("`ctools_plugin_hdfe_markout'" == "1")
    c_local __ctools_plugin_hdfe_fvlevels = ("`ctools_plugin_hdfe_fvlevels'" == "1")
    capture confirm number 0
end
* ctools-helpers:end
