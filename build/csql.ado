*! version 1.4.0 20261127 github.com/mdroste/stata-ctools
*! SQL queries on Stata data and on Parquet, CSV, JSON and .dta files, through DuckDB (ctools)

program define csql, rclass
    version 14.1

    syntax [using/] [, Query(string) File(string) CLEAR SAVing(string asis) ///
        ROWs(integer 20) DATAbase(string) READonly MEMory(string) TEMPdir(string) ///
        noLABel noDATEtime STRL(integer 2045) LIBrary(string) THReads(integer 0) ///
        VERbose INSTall REPLACE FROM(string) VERsion]

    if `threads' < 0 {
        di as error "threads() must be nonnegative"
        exit 198
    }

    * =========================================================================
    * csql, install and csql, version
    * =========================================================================
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    _csql_platform `"`__ctools_plugin'"'
    local platform `"`r(platform)'"'

    if "`install'" != "" {
        if `"`using'`query'`file'`saving'`database'"' != "" | "`clear'" != "" {
            di as error "option {bf:install} may not be combined with a query"
            exit 198
        }
        _csql_install, platform(`platform') `replace' from(`"`from'"') library(`"`library'"')
        return add
        exit
    }
    if "`replace'" != "" | `"`from'"' != "" {
        di as error "options {bf:replace} and {bf:from()} belong to {bf:csql, install}"
        exit 198
    }

    _csql_library, platform(`platform') library(`"`library'"')
    local lib `"`r(path)'"'

    if "`version'" != "" {
        local __csql_lib `"`lib'"'
        plugin call ctools_plugin, "csql version"
        di as text "DuckDB " as result "`csql_version'" as text `" (`lib')"'
        return local version "`csql_version'"
        return local library `"`lib'"'
        exit
    }

    * =========================================================================
    * The query
    * =========================================================================
    if `"`query'"' != "" & `"`file'"' != "" {
        di as error "options {bf:query()} and {bf:file()} may not be combined"
        exit 198
    }
    if `"`query'"' == "" & `"`file'"' == "" {
        di as error "option {bf:query()} or {bf:file()} required"
        exit 198
    }
    if `rows' < 0 {
        di as error "rows() must be nonnegative"
        exit 198
    }
    if `strl' < 1 | `strl' > 2045 {
        di as error "strl() must be an integer from 1 through 2045"
        exit 198
    }
    if `"`file'"' != "" {
        confirm file `"`file'"'
    }
    if "`readonly'" != "" & `"`database'"' == "" {
        di as error "option {bf:readonly} requires {bf:database()}"
        exit 198
    }
    if `"`database'"' != "" & "`readonly'" != "" {
        confirm file `"`database'"'
    }

    * saving(filename [, replace]): .dta through the plugin, other formats
    * through DuckDB's COPY
    local savepath
    local copyfmt
    if `"`saving'"' != "" {
        _csql_saving `saving'
        local savepath `"`r(path)'"'
        local copyfmt "`r(format)'"
        if "`r(replace)'" == "" {
            confirm new file `"`savepath'"'
        }
        if "`copyfmt'" != "" & "`clear'" != "" {
            di as error "option {bf:clear} loads the result; with {bf:saving()} it requires a .dta file"
            exit 198
        }
    }

    local timing = ("`verbose'" != "")
    if `timing' quietly _ctools_clock __ct_start

    * The SQL goes to the plugin in a file: any length, any characters.
    tempfile sqlfile
    if `"`file'"' != "" {
        local sqlfile `"`file'"'
    }
    else {
        tempname fh
        quietly file open `fh' using `"`sqlfile'"', write text replace
        file write `fh' `"`macval(query)'"'
        file close `fh'
    }

    * The data in memory: the plugin reads it through the SPI; storage types,
    * formats and labels come from this description.
    local meta
    local varlist
    if c(k) > 0 & _N > 2147483647 {
        * Stata's plugin interface numbers observations with 32-bit integers
        di as text "{p 0 6 2}note: the data in memory have more than 2,147,483,647 observations;" ///
            " table {bf:stata} is not available{p_end}"
    }
    else if c(k) > 0 {
        tempfile meta
        mata: _csql_meta(st_local("meta"), "`label'" == "")
        local varlist "_all"
    }

    * Where the rows go
    local out
    local temporary
    if "`copyfmt'" != "" {
        local out `"`savepath'"'
    }
    else if `"`savepath'"' != "" {
        local out `"`savepath'"'
    }
    else if "`clear'" != "" {
        tempfile result
        local out `"`result'"'
        local temporary "temporary"
    }

    * DuckDB spills to disk past its memory limit; an in-memory database
    * would otherwise spill into the working directory.
    if `"`tempdir'"' == "" & `"`database'"' == "" {
        tempfile spill
        local tempdir `"`spill'_duckdb"'
        local own_tempdir 1
    }
    else local own_tempdir 0

    local __csql_lib `"`lib'"'
    local __csql_sqlfile `"`sqlfile'"'
    local __csql_meta `"`meta'"'
    local __csql_out `"`out'"'
    local __csql_using `"`using'"'
    local __csql_database `"`database'"'
    local __csql_memory `"`memory'"'
    local __csql_tempdir `"`tempdir'"'
    local labels_on = ("`label'" == "")
    local dates_on = ("`datetime'" == "")
    local opts "labels=`labels_on' datetime=`dates_on' rows=`rows' linesize=`c(linesize)' strl=`strl'"
    if "`readonly'" != "" local opts "`opts' readonly"
    if "`copyfmt'" != "" local opts "`opts' copy=`copyfmt'"
    if "`temporary'" != "" local opts "`opts' temporary"
    if `timing' local opts "`opts' verbose"
    local threadopt
    if `threads' > 0 local threadopt "threads(`threads') "

    if `timing' {
        quietly _ctools_clock __ct_now
        local __t_parse = `__ct_now' - `__ct_start'
        local __ct_start = `__ct_now'
    }
    local csql_renamed
    capture noisily plugin call ctools_plugin `varlist', "csql `threadopt'run `opts'"
    local rc = _rc
    if `own_tempdir' capture rmdir `"`tempdir'"'
    if `rc' {
        if "`temporary'" != "" capture erase `"`out'"'
        exit `rc'
    }
    local N = scalar(__csql_N)
    local k = scalar(__csql_k)
    local changed = scalar(__csql_changed)
    capture scalar drop __csql_N __csql_k __csql_changed
    if `timing' {
        quietly _ctools_clock __ct_now
        local __t_plugin = `__ct_now' - `__ct_start'
        local __ct_start = `__ct_now'
    }

    * =========================================================================
    * Results
    * =========================================================================
    if "`clear'" != "" {
        * use names the data after the file it read; a tempfile is no name
        mata: st_local("fn0", st_global("S_FN")); st_local("fndate0", st_global("S_FNDATE"))
        quietly use `"`out'"', clear
        if "`temporary'" != "" {
            capture erase `"`out'"'
            mata: st_global("S_FN", ""); st_global("S_FNDATE", "")
        }
        di as text "(" as result `k' as text " var" cond(`k' == 1, "", "s") ", " ///
            as result `N' as text " obs)"
    }
    else if `"`savepath'"' != "" {
        di as text `"(file `savepath' saved: "' as result `N' as text " row" cond(`N' == 1, "", "s") ")"
    }
    if `"`csql_renamed'"' != "" & `"`out'"' != "" & "`copyfmt'" == "" {
        di as text "{p 0 6 2}note: column names changed to valid Stata names: `csql_renamed'; " ///
            "the variable labels hold the original names{p_end}"
    }
    if `timing' {
        quietly _ctools_clock __ct_now
        local __t_load = `__ct_now' - `__ct_start'
        di as text "csql timing: parse " %6.3f `__t_parse' " s, query " %6.3f `__t_plugin' ///
            " s, load " %6.3f `__t_load' " s"
    }

    return scalar N = `N'
    return scalar k = `k'
    if `changed' < . return scalar changed = `changed'
    return local version "`csql_version'"
    return local library `"`lib'"'
    if `"`savepath'"' != "" return local saving `"`savepath'"'
    return local renamed `"`csql_renamed'"'
end

* =============================================================================
* Platform and DuckDB library
* =============================================================================

* The DuckDB release archive for the plugin this Stata loads.
program define _csql_platform, rclass
    args plugin
    local base = lower(substr(`"`plugin'"', strrpos(subinstr(`"`plugin'"', "\", "/", .), "/") + 1, .))
    if strpos("`base'", "windows") return local platform "windows"
    else if strpos("`base'", "mac_arm") return local platform "mac_arm"
    else if strpos("`base'", "mac_x86") return local platform "mac_x86"
    else if strpos("`base'", "linux") return local platform "linux"
    else return local platform ""
end

* The pinned DuckDB release: archive, its SHA-256 checksum, the library in it,
* and the name csql installs it under.
program define _csql_release, rclass
    syntax, platform(string)
    local v "2.0.1"
    return local version "`v'"
    local base "https://github.com/duckdb/duckdb/releases/download/v`v'"
    if "`platform'" == "windows" {
        return local archive "duckdb-shared-libs-windows-amd64.tar.gz"
        return local sha256 "d51f23831423b82418ee1d5bae324484f9fbe4a2c2424fa5b6be3392453c62ce"
        return local member "duckdb.dll"
        return local file "duckdb_windows.dll"
    }
    else if "`platform'" == "mac_arm" {
        return local archive "duckdb-shared-libs-osx-arm64.tar.gz"
        return local sha256 "d71e4800350f1135bc3f80a2e5cbcc01e199ca03759edf5a376db308b514b341"
        return local member "libduckdb.dylib"
        return local file "duckdb_mac_arm.dylib"
    }
    else if "`platform'" == "mac_x86" {
        return local archive "duckdb-shared-libs-osx-amd64.tar.gz"
        return local sha256 "e12f174ef01575590ec6e27a5d475308680a9e15292691623d9848c0fc5a9ba8"
        return local member "libduckdb.dylib"
        return local file "duckdb_mac_x86.dylib"
    }
    else if "`platform'" == "linux" {
        return local archive "duckdb-shared-libs-linux-amd64.tar.gz"
        return local sha256 "e25502612fed533a243b61f23a2a574daad030e7e88906222a5f0ab30d06769e"
        return local member "libduckdb.so"
        return local file "duckdb_linux.so"
    }
    else {
        di as error "csql: no DuckDB build is known for this platform"
        exit 198
    }
    return local url "`base'/`return(archive)'"
end

* The library to load: library(), then $CSQL_DUCKDB, then the installed file
* on the adopath.
program define _csql_library, rclass
    syntax, platform(string) [library(string)]
    if `"`library'"' == "" local library `"$CSQL_DUCKDB"'
    if `"`library'"' != "" {
        capture confirm file `"`library'"'
        if _rc {
            di as error `"DuckDB library `library' not found"'
            exit 601
        }
        return local path `"`library'"'
        exit
    }
    _csql_release, platform(`platform')
    local name "`r(file)'"
    capture findfile `name'
    if _rc {
        di as error "csql needs the DuckDB library, which is not installed"
        di as error "{p 4 4 2}Type {bf:csql, install} to download DuckDB `r(version)' from"
        di as error "github.com/duckdb/duckdb (about 20 MB), or see {help csql##install:help csql}.{p_end}"
        exit 601
    }
    return local path `"`r(fn)'"'
end

program define _csql_install, rclass
    syntax, platform(string) [replace from(string) library(string)]
    _csql_release, platform(`platform')
    local ver "`r(version)'"
    local url "`r(url)'"
    local sha "`r(sha256)'"
    local member "`r(member)'"
    local name "`r(file)'"
    local archive "`r(archive)'"

    if `"`library'"' != "" local dest `"`library'"'
    else {
        local plus = c(sysdir_plus)
        capture mkdir `"`plus'"'
        capture mkdir `"`plus'd"'
        local dest `"`plus'd/`name'"'
    }
    capture confirm file `"`dest'"'
    if !_rc & "`replace'" == "" {
        di as text `"DuckDB is already installed (`dest'); specify {bf:replace} to install it again"'
        return local library `"`dest'"'
        exit
    }

    local src `"`from'"'
    if `"`src'"' != "" {
        * a folder holding the archive, the archive itself, or the library
        mata: st_local("isdir", strofreal(direxists(st_local("src"))))
        if `isdir' local src `"`src'/`archive'"'
        confirm file `"`src'"'
        if !strpos(lower(`"`src'"'), ".tar.gz") {
            copy `"`src'"' `"`dest'"', replace
            di as text `"(DuckDB library copied to `dest')"'
            return local library `"`dest'"'
            exit
        }
    }
    else {
        tempfile download
        di as text "downloading DuckDB `ver' from github.com/duckdb/duckdb ..."
        capture noisily copy `"`url'"' `"`download'"'
        if _rc {
            di as error "{p 0 4 2}could not download `url'; download it with a web browser and"
            di as error "type {bf:csql, install from(}{it:folder}{bf:)}{p_end}"
            exit _rc
        }
        local src `"`download'"'
    }
    local __csql_archive `"`src'"'
    local __csql_sha256 "`sha'"
    local __csql_member "`member'"
    local __csql_dest `"`dest'"'
    plugin call ctools_plugin, "csql install"
    local __csql_lib `"`dest'"'
    plugin call ctools_plugin, "csql version"
    di as text "DuckDB " as result "`csql_version'" as text `" installed as `dest'"'
    return local version "`csql_version'"
    return local library `"`dest'"'
end

* saving(filename [, replace]): the format follows the extension.
program define _csql_saving, rclass
    syntax anything(name=path id="filename") [, REPLACE]
    local path `path'
    local low = lower(`"`path'"')
    local fmt
    if substr(`"`low'"', -4, .) == ".dta" local fmt ""
    else if substr(`"`low'"', -8, .) == ".parquet" | substr(`"`low'"', -5, .) == ".parq" local fmt "parquet"
    else if substr(`"`low'"', -4, .) == ".csv" | substr(`"`low'"', -7, .) == ".csv.gz" local fmt "csv"
    else if substr(`"`low'"', -5, .) == ".json" | substr(`"`low'"', -6, .) == ".jsonl" | ///
        substr(`"`low'"', -7, .) == ".ndjson" local fmt "json"
    else {
        * no extension: a Stata dataset, as save does
        if !strpos(substr(`"`path'"', strrpos(subinstr(`"`path'"', "\", "/", .), "/") + 1, .), ".") {
            local path `"`path'.dta"'
        }
        else {
            di as error "{p 0 4 2}saving(): the file must end in .dta, .parquet, .csv or .json{p_end}"
            exit 198
        }
    }
    return local path `"`path'"'
    return local format "`fmt'"
    return local replace "`replace'"
end

* =============================================================================
* The description of the data in memory, for the plugin
* =============================================================================
mata:
string scalar _csql_field(string scalar s)
{
    return(strofreal(strlen(s)) + ":" + s)
}

void _csql_meta(string scalar path, real scalar labels)
{
    real scalar fh, j, i, k
    string scalar vl, line
    string rowvector done
    real colvector vals
    string colvector txt

    k = st_nvar()
    fh = fopen(path, "w")
    for (j = 1; j <= k; j++) {
        line = "V " + _csql_field(st_varname(j)) + _csql_field(st_vartype(j)) +
            _csql_field(st_varformat(j)) +
            _csql_field(labels ? st_varvaluelabel(j) : "") +
            _csql_field(labels ? st_varlabel(j) : "")
        fwrite(fh, line + char(10))
    }
    if (labels) {
        done = J(1, 0, "")
        for (j = 1; j <= k; j++) {
            vl = st_varvaluelabel(j)
            if (vl == "") continue
            if (anyof(done, vl)) continue
            if (!st_vlexists(vl)) continue
            done = done, vl
            st_vlload(vl, vals, txt)
            for (i = 1; i <= rows(vals); i++) {
                if (vals[i] >= .) continue
                fwrite(fh, "L " + _csql_field(vl) + _csql_field(strofreal(vals[i], "%12.0f")) +
                    _csql_field(txt[i]) + char(10))
            }
        }
    }
    fclose(fh)
}
end

* ctools-helpers:begin (generated by scripts/sync_release.py from src/ado; edit those files)
* Private monotonic timestamp in seconds; no public Stata timer IDs.
program define _ctools_clock
    version 14.1
    args target
    tempname held
    _return hold `held'
    capture noisily {
        * Once registered, this plugin stays loaded for the Stata session.
        * Resolve its path and check its protocol only on first use (or
        * after discard), rather than searching the ado-path at every lap.
        capture plugin call ctools_clock_plugin, "clock"
        if _rc {
            _ctools_load
            capture program ctools_clock_plugin, plugin using("`__ctools_plugin'")
            if _rc != 0 & _rc != 110 exit _rc
            plugin call ctools_clock_plugin, "clock"
        }
    }
    local rc = _rc
    _return restore `held'
    if `rc' exit `rc'
    c_local `target' `ctools_clock_value'
end

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
