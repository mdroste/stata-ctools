*! version 1.2.3 20261030 github.com/mdroste/stata-ctools

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

program define cdecode, rclass
    version 14.1
    syntax varlist [if] [in], [Generate(string) replace MAXLength(integer 0) Verbose THReads(integer 0)]

    if ("`generate'" == "") == ("`replace'" == "") {
        di as error "cdecode: specify either generate() or replace"
        exit 198
    }
    if `maxlength' < 0 | `threads' < 0 {
        di as error "cdecode: maxlength() and threads() must be nonnegative"
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
    foreach v of local varlist {
        capture confirm numeric variable `v'
        if _rc exit 108
        local label : value label `v'
        if "`label'" == "" {
            di as error "cdecode: `v' has no value label attached"
            exit 182
        }
    }

    * Stata answers an if one observation at a time, so the plugin is passed
    * the sample as a marker variable, last, and reads it in parallel. An if
    * that names a numeric variable other than a source is that marker.
    * Several outputs otherwise share one marked sample, evaluated once. A
    * single output marks its if after its destination exists: the marker is
    * then the last variable, and dropping it after the call is free. With
    * in, which limits the range, a single output keeps any other if (marking
    * covers every observation).
    local sample `"`if' `in'"'
    local plugin_sample `"`sample'"'
    local ifflag = cond(`"`if'"' != "", "ifcond", "")
    local touse ""
    local mark_each = 0
    if `"`if'"' != "" {
        local touse = strtrim(substr(`"`if'"', 3, .))
        capture confirm numeric variable `touse', exact
        if _rc | wordcount(`"`touse'"') != 1 local touse ""
        else if `: list touse in varlist' local touse ""
        if "`touse'" != "" {
            local plugin_sample `"`in'"'
            local ifflag "ifcond touse"
        }
        else if `nvars' > 1 {
            marksample touse, novarlist
            local sample "if `touse'"
            local plugin_sample ""
            local ifflag "ifcond touse"
        }
        else if `"`in'"' == "" {
            local mark_each = 1
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
    if !_rc {
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
            local marker "`touse'"
            if `mark_each' {
                tempvar marker
                mark `marker' `if'
            }
            capture plugin call ctools_plugin `v' `decoded' `marker' `plugin_sample', ///
                "cdecode `threads_code' v2 maplen=`__cd_len' width=`__cd_width' `ifflag' `fill'"
            local plugin_rc = _rc
            if `mark_each' drop `marker'
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
            if "`replace'" != "" drop `v'
            else local v : word `i' of `generate'
            mata: st_varrename("`decoded'", "`v'")
        }
        if "`replace'" != "" order `original_order'
    }
    if "`verbose'" != "" {
        di as text "cdecode: decoded `nvars' variable(s): `n_plugin' in the ctools plugin, " ///
            "`=`nvars' - `n_plugin'' with Stata's decode"
        if `n_plugin' > 0 {
            capture local __threads_max = _cdecode_threads_max
            if !_rc di as text "  threads available: `__threads_max'"
        }
    }
    return scalar N_vars = `nvars'
    return scalar N_plugin = `n_plugin'
end
