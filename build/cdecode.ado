*! version 1.3.2 20261124 github.com/mdroste/stata-ctools
*! original: Stata's decode (built in)
*! original: see help cdecode for details

* Hand one value label to the plugin as local __cd_map (format: see
* src/cdecode/cdecode_impl.c). Sets __cd_ok to 0 when decode itself must do the
* work: an undefined or empty label, a strL result, or a maxlength() cut inside
* text with multibyte characters, which decode never splits.
capture mata: mata drop _cdecode_pack()
mata:
void _cdecode_pack(string scalar label, real scalar maxlength)
{
    real colvector values, lengths
    string colvector text, entries
    string scalar packed
    real scalar i, width

    st_local("__cd_ok", "0")
    if (!st_vlexists(label)) return
    st_vlload(label, values, text)
    if (rows(values) == 0) return
    lengths = strlen(text)
    if (maxlength > 0) {
        for (i = 1; i <= rows(text); i++) {
            if (lengths[i] > maxlength) {
                if (max(ascii(text[i])) > 127) return
                text[i] = substr(text[i], 1, maxlength)
            }
        }
        lengths = strlen(text)
    }
    /* Like decode, the new variable is as wide as the longest label text. */
    width = max((1 \ lengths))
    if (width > 2045) return
    entries = strtrim(strofreal(values, "%21.0f")) :+ " " :+ strofreal(lengths) :+ " " :+ text
    packed = strofreal(rows(values)) + " " + invtokens(entries', " ") + " |"
    if (strlen(packed) != strlen(strofreal(rows(values))) + sum(strlen(entries)) + rows(values) + 2) return
    if (strlen(packed) >= c("macrolen")) return
    st_local("__cd_map", packed)
    st_local("__cd_len", strofreal(strlen(packed), "%21.0f"))
    st_local("__cd_width", strofreal(width))
    st_local("__cd_ok", "1")
}
end

* decode leaves the caller's r() in place, so cdecode keeps it too: it is
* held while cdecode works and returned along with r(N_vars) and
* r(N_plugin). It is in place while the if is evaluated, as it may name it.
program define cdecode, rclass
    version 14.1
    tempname caller_r
    _return hold `caller_r'
    _return restore `caller_r', hold
    capture noisily _cdecode_impl `0'
    local rc = _rc
    if !`rc' {
        local n_vars = r(N_vars)
        local n_plugin = r(N_plugin)
    }
    _return restore `caller_r'
    return add
    if `rc' exit `rc'
    return scalar N_vars = `n_vars'
    return scalar N_plugin = `n_plugin'
end

program define _cdecode_impl, rclass
    version 14.1
    syntax varlist [if] [in], [Generate(string) replace MAXLength(string) Verbose THReads(integer 0)]

    if ("`generate'" == "") == ("`replace'" == "") {
        di as error "cdecode: specify either generate() or replace"
        exit 198
    }
    * As decode: an explicit maxlength() is an integer from 1 to 32000;
    * 0 below means none
    if "`maxlength'" != "" {
        capture confirm integer number `maxlength'
        if _rc | !inrange(real("`maxlength'"), 1, 32000) {
            di as error "maxlength must be between 1 and 32000"
            exit 198
        }
    }
    else local maxlength 0
    if `threads' < 0 {
        di as error "cdecode: threads() must be nonnegative"
        exit 198
    }
    local nvars : word count `varlist'
    if "`replace'" != "" unab original_order : _all
    if "`generate'" != "" {
        local nout : word count `generate'
        local unique : list uniq generate
        local nunique : word count `unique'
        if `nout' != `nvars' | `nunique' != `nout' {
            di as error "cdecode: generate() requires one distinct name per input variable"
            exit 198
        }
        confirm new variable `generate'
    }
    * decode's messages and return codes
    foreach v of local varlist {
        capture confirm numeric variable `v'
        if _rc {
            di as error "not possible with string variable"
            exit 108
        }
        local label : value label `v'
        if "`label'" == "" {
            di as error "`v' not labeled"
            exit 182
        }
    }

    * Stata answers an if one observation at a time, so the plugin is passed
    * the sample as a marker variable, last, and reads it in parallel. An if
    * that names a numeric variable other than a source is that marker. Any
    * other if is marked here, once, before a helper can replace the r() or
    * e() results it may name.
    local sample `"`if' `in'"'
    local plugin_sample `"`sample'"'
    local ifflag = cond(`"`if'"' != "", "ifcond", "")
    local touse ""
    if `"`if'"' != "" {
        local touse = strtrim(substr(`"`if'"', 3, .))
        capture confirm numeric variable `touse', exact
        if _rc | wordcount(`"`touse'"') != 1 local touse ""
        else if `: list touse in varlist' local touse ""
        if "`touse'" != "" {
            local plugin_sample `"`in'"'
            local ifflag "ifcond touse"
        }
        else {
            marksample touse, novarlist
            local sample "if `touse'"
            local plugin_sample ""
            local ifflag "ifcond touse"
        }
    }
    * Without if/in the plugin writes every observation, so its destination
    * can be created without first being filled.
    local fill = cond(`"`if'`in'"' == "", "fillall", "")
    local maxopt ""
    if `maxlength' > 0 local maxopt "maxlength(`maxlength')"
    local threads_code ""
    if `threads' > 0 local threads_code "threads(`threads')"

    * The plugin is an accelerator only: without it every variable is decoded
    * by Stata, as are labels the plugin protocol does not cover.
    local have_plugin = 0
    capture _ctools_load
    * Above Stata's plugin limit of 2^31-1 observations Stata decodes everything.
    if !_rc & !(_N > 2147483647) {
        capture program ctools_plugin, plugin using("`__ctools_plugin'")
        if _rc == 0 | _rc == 110 local have_plugin = 1
    }
    capture confirm number 0

    * Every output is first built under a temporary name. The inputs are not
    * touched until all variables have been decoded, so a failure or Break
    * before then leaves the dataset as it was.
    local staged ""
    local n_plugin = 0
    foreach v of local varlist {
        tempvar decoded
        local done = 0
        local __cd_ok = 0
        if `have_plugin' & `maxlength' <= 32000 & _N > 0 {
            mata: _cdecode_pack("`: value label `v''", `maxlength')
        }
        if `__cd_ok' {
            mata: (void) st_addvar("str`__cd_width'", "`decoded'", "`fill'" != "")
            capture plugin call ctools_plugin `v' `decoded' `touse' `plugin_sample', ///
                "cdecode `threads_code' v2 maplen=`__cd_len' width=`__cd_width' `ifflag' `fill'"
            local plugin_rc = _rc
            if `plugin_rc' == 1 exit 1
            if `plugin_rc' {
                drop `decoded'
            }
            else {
                mata: st_varlabel("`decoded'", st_varlabel("`v'"))
                local done = 1
                local ++n_plugin
            }
        }
        local __cd_map ""
        if !`done' {
            quietly decode `v' `sample', generate(`decoded') `maxopt'
        }
        local staged `staged' `decoded'
    }

    * st_varrename() is built in; rename is an ado file whose first use in a
    * session costs about 10 ms. One order for all outputs moves the data once.
    nobreak {
        local i = 0
        foreach v of local varlist {
            local ++i
            local decoded : word `i' of `staged'
            * replace keeps the source's notes and characteristics
            if "`replace'" != "" {
                char rename `v' `decoded'
                drop `v'
            }
            else local v : word `i' of `generate'
            mata: st_varrename("`decoded'", "`v'")
        }
        if "`replace'" != "" order `original_order'
    }
    if "`verbose'" != "" {
        di as text "cdecode: decoded `nvars' variable(s): `n_plugin' in the ctools plugin, " ///
            "`=`nvars' - `n_plugin'' with Stata's decode"
        if `n_plugin' > 0 {
            capture local __threads_max = scalar(_cdecode_threads_max)
            if !_rc di as text "  threads available: `__threads_max'"
        }
    }
    return scalar N_vars = `nvars'
    return scalar N_plugin = `n_plugin'
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
