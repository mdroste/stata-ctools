*! version 1.4.0 20261127 github.com/mdroste/stata-ctools
*! original: winsor2 1.1 (December 16, 2014) by Yujun Lian
*! original: gstats winsor from gtools 1.10.1 (December 5, 2022) by Mauricio Caceres Bravo
*! original: see help cwinsor for details

program define cwinsor, rclass
    version 14.1
    * A failed or interrupted call leaves the data as they were. New
    * variables (suffix(), prefix()) are only
    * added at the end and dropped again. In place, the plugin keeps the
    * original values of the observations it rewrites until cwinsor finishes
    * and writes them back after a failure; _ctools_undo then restores the
    * variable labels, the sort marker and c(changed). Option snapshot
    * preserves the whole dataset instead (on disk above c(max_preservemem));
    * option fast keeps no copy and undoes nothing. Replacing existing
    * generated targets uses a snapshot to retain their original contents.
    _ctools_mode `0'
    local __command `"`0'"'
    if "`mode'" == "" {
        syntax varlist(numeric) [if] [in] [aw fw iw pw], [REPLACE ///
            P(string) PREfix(string) Suffix(string) GENerate(namelist) *]
        local __target_prefix `"`prefix'"'
        if `"`p'"' != "" {
            capture confirm number `p'
            if _rc local __target_prefix `"`p'"'
        }
        if "`replace'" != "" & "`__target_prefix'`suffix'`generate'" != "" {
            local mode snapshot
            local __command `"`__command' snapshot"'
        }
        local 0 `"`__command'"'
    }
    * Break is recognized only inside _cwinsor_impl.
    nobreak {
        macro drop CTOOLS_cwinsor_plugin CTOOLS_cwinsor_targets CTOOLS_cwinsor_recast
        if "`mode'" == "snapshot" preserve
        else if "`mode'" == "" _ctools_undo record CTOOLS_u_cwinsor_
        capture noisily break _cwinsor_impl `0'
        local rc = _rc
        * gstats winsor's results, kept before anything else runs
        if !`rc' {
            foreach r in J minJ maxJ cutlow cuthigh {
                local r_`r' = r(`r')
            }
        }
        if `rc' {
            if "`mode'" == "" {
                * Values first, then labels, sort marker and c(changed)
                if "$CTOOLS_cwinsor_targets" != "" {
                    capture program ctools_plugin, plugin using("$CTOOLS_cwinsor_plugin")
                    capture noisily plugin call ctools_plugin $CTOOLS_cwinsor_targets, "cwinsor rollback"
                    if _rc di as error "cwinsor: the original values could not be restored"
                }
                * Storage promoted for fractional cutoffs goes back too (recast
                * refuses if a value would change)
                local recast "$CTOOLS_cwinsor_recast"
                while "`recast'" != "" {
                    gettoken target recast : recast
                    gettoken vtype recast : recast
                    capture recast `vtype' `target'
                }
                _ctools_undo rollback CTOOLS_u_cwinsor_
            }
            else if "`mode'" == "snapshot" restore
            if "$CTOOLS_cwinsor_plugin" != "" {
                capture program ctools_plugin, plugin using("$CTOOLS_cwinsor_plugin")
                capture plugin call ctools_plugin, "cwinsor commit"
            }
            macro drop CTOOLS_cwinsor_plugin CTOOLS_cwinsor_targets CTOOLS_cwinsor_recast
            exit `rc'
        }
        if "$CTOOLS_cwinsor_plugin" != "" {
            capture program ctools_plugin, plugin using("$CTOOLS_cwinsor_plugin")
            capture plugin call ctools_plugin, "cwinsor commit"
        }
        macro drop CTOOLS_cwinsor_plugin CTOOLS_cwinsor_targets CTOOLS_cwinsor_recast
        if "`mode'" == "snapshot" restore, not
        else if "`mode'" == "" _ctools_undo clear CTOOLS_u_cwinsor_
        foreach r in J minJ maxJ cutlow cuthigh {
            return scalar `r' = `r_`r''
        }
    }
end

program define _cwinsor_impl, rclass
    version 14.1

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    * Start total timing FIRST
    local __ct_elapsed90 = 0
    local __ct_elapsed91 = 0
    local __ct_elapsed92 = 0
    if "`__ct_elapsed90'" == "" local __ct_elapsed90 = 0
    quietly _ctools_clock __ct_start90

    * The abbreviations of winsor2 and gstats winsor: s() suffix, p() prefix,
    * t trim, c() cuts, l label, gen() generate. p() holding a number is
    * cwinsor's lower percentile (a prefix cannot be a number).
    syntax varlist(numeric) [if] [in] [aw fw iw pw], [Cuts(numlist min=2 max=2 >=0 <=100) P(string) PREfix(string) Q(real 99) Trim BY(varlist) Suffix(string) GENerate(namelist) REPLACE Label noMISSing Verbose THReads(integer 0) SNAPshot FAST ///
        compress forcestrl BENCHmark BENCHmarklevel(integer 0) HASHmethod(string) oncollision(string)]

    * Freeze incoming-r() qualifiers before parsing helpers replace r().
    if `"`if'"' != "" {
        tempvar __incoming_sample
        mark `__incoming_sample' `if' `in'
        local if "if `__incoming_sample'"
    }

    local __p_arg `"`p'"'
    local p = 1
    if `"`prefix'"' != "" {
        capture confirm number `prefix'
        if !_rc {
            local p = `prefix'
            local prefix
        }
    }
    if `"`__p_arg'"' != "" {
        capture confirm number `__p_arg'
        if !_rc local p = `__p_arg'
        else {
            if `"`prefix'"' != "" error 198
            local prefix `"`__p_arg'"'
        }
    }

    * gstats winsor's gtools options. compress and forcestrl concern reading
    * strL variables, which cwinsor reads without them; hashmethod() and
    * oncollision() choose gtools' hashing of the by() groups, which cwinsor
    * sorts instead. They are checked as gtools checks them and have no
    * effect; benchmark and benchmarklevel() show cwinsor's timing.
    if `"`hashmethod'"' != "" & !inlist(`"`hashmethod'"', "0", "1", "2", "default", "biject", "spooky") {
        di as error `"`hashmethod' not found"'
        exit 111
    }
    if `"`oncollision'"' != "" & !inlist(`"`oncollision'"', "fallback", "error") {
        di as error "option {opt oncollision()} must be 'fallback' or 'error'"
        exit 198
    }
    if "`benchmark'" != "" | `benchmarklevel' > 0 local verbose verbose

    local __do_timing = ("`verbose'" != "")

    * Handle percentile options
    if "`cuts'" != "" {
        local p : word 1 of `cuts'
        local q : word 2 of `cuts'
    }

    * Validate percentiles
    if `p' < 0 | `p' > 100 {
        di as error "cwinsor: lower percentile must be between 0 and 100"
        exit 198
    }
    if `q' < 0 | `q' > 100 {
        di as error "cwinsor: upper percentile must be between 0 and 100"
        exit 198
    }
    * Equal percentiles are allowed, as in winsor2
    if `p' > `q' {
        di as error "cwinsor: lower percentile (`p') must not exceed upper (`q')"
        exit 198
    }

    * Validate suffix/prefix/generate. As in gstats winsor, replace with them
    * overwrites targets that exist.
    * prefix() and suffix() together name prefix+source+suffix, as in gstats
    if "`generate'" != "" & "`suffix'`prefix'" != "" {
        di as error "cwinsor: specify only one of prefix()/suffix() or generate()"
        exit 198
    }
    if "`generate'" != "" {
        if `: word count `generate'' != `: word count `varlist'' {
            di as error "cwinsor: specify the same number of targets as sources with generate()"
            exit 198
        }
    }

    * Check numeric variables
    foreach v of local varlist {
        capture confirm numeric variable `v'
        if _rc != 0 {
            di as error "cwinsor: `v' is not numeric"
            exit 109
        }
    }

    local nvars : word count `varlist'

    * The plugin's stores into a frame alias variable have no effect; refuse
    * them as replace and winsor2 do
    local inplace = "`suffix'`prefix'`generate'" == ""
    if `inplace' & c(stata_version) >= 18 {
        foreach v of local varlist {
            mata: st_local("__cw_alias", strofreal(st_isalias(st_local("v"))))
            if `__cw_alias' {
                di as error "cannot change values of alias variables"
                di as error "    You tried to change values in alias variable `v'."
                exit 109
            }
        }
    }

    * Mark sample. With weights, as gstats winsor marks it: missing and zero
    * weights are excluded, negative ones refused (also iweights, which
    * gstats accepts; the weighted percentile needs positive weights), and
    * fweights must be integers.
    marksample touse, novarlist
    local weighted = ("`weight'" != "")
    if `weighted' {
        tempvar wvar
        quietly generate double `wvar' `exp' if `touse'
        quietly count if `touse' & `wvar' < 0
        if r(N) {
            di as error "negative weights encountered"
            exit 402
        }
        markout `touse' `wvar'
        quietly replace `touse' = 0 if `wvar' == 0
    }
    * nomissing: with by(), groups with a missing by-value are left out
    if "`missing'" == "nomissing" & "`by'" != "" markout `touse' `by', strok

    * Load plugin
    _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0

    * Count observations
    quietly count if `touse'
    local nobs_sample = r(N)
    if `nobs_sample' == 0 {
        di as error "cwinsor: no observations"
        exit 2000
    }

    * Count by-variables (sorting now done in C for speed)
    if "`__ct_elapsed91'" == "" local __ct_elapsed91 = 0
    quietly _ctools_clock __ct_start91
    local nby = 0
    if "`by'" != "" {
        local nby : word count `by'
    }
    quietly _ctools_clock __ct_now
    local __ct_elapsed91 = `__ct_elapsed91' + (`__ct_now' - `__ct_start91')

    * Floating-point outputs can be written directly from the plugin's double
    * working columns. Native recast otherwise copies and scans each column
    * again, then moves it back through the entire dataset.
    local direct_outputs = !`inplace' & "`trim'" == ""
    local overlap : list varlist & by
    if "`overlap'" != "" local direct_outputs = 0
    if "`replace'" != "" local direct_outputs = 0
    foreach v of local varlist {
        local output_type : type `v'
        if "`by'" != "" local output_type "`c(type)'"
        if !inlist("`output_type'", "float", "double") local direct_outputs = 0
    }

    * Handle generate new variables
    local target_varlist "`varlist'"
    if !`inplace' {
        local target_varlist ""
        local gen_index = 0
        * The names first, so nothing is created before they are checked
        foreach v of local varlist {
            local ++gen_index
            if "`suffix'`prefix'" != "" local newvar "`prefix'`v'`suffix'"
            else local newvar : word `gen_index' of `generate'
            confirm name `newvar'
            capture confirm variable `newvar', exact
            if _rc == 0 {
                if "`replace'" == "" {
                    di as error "cwinsor: variable `newvar' already exists"
                    exit 110
                }
                capture confirm numeric variable `newvar', exact
                if _rc {
                    di as error "cwinsor: variable `newvar' is not numeric"
                    exit 109
                }
                * A source as its own target is winsorized in place, which
                * an if/in or weighted sample would turn into missing values
                * outside the sample (gstats winsor refuses it with if/in)
                local source_pos : list posof "`newvar'" in varlist
                if `source_pos' & `source_pos' != `gen_index' {
                    di as error "cwinsor: target `newvar' is another source"
                    exit 198
                }
                if `source_pos' & (`"`if'`in'"' != "" | `weighted') {
                    di as error "cwinsor: source `newvar' as its own target is not allowed with if/in or weights"
                    exit 198
                }
            }
            local target_varlist "`target_varlist' `newvar'"
        }
        local dups : list dups target_varlist
        if "`dups'" != "" {
            di as error "cwinsor: target `: word 1 of `dups'' named twice"
            exit 198
        }
        local gen_index = 0
        foreach newvar of local target_varlist {
            local ++gen_index
            local v : word `gen_index' of `varlist'
            capture confirm variable `newvar', exact
            if _rc == 0 {
                * replace: the existing target is overwritten, keeping its
                * place in the data
                if "`newvar'" != "`v'" {
                    quietly recast double `newvar'
                    quietly replace `newvar' = cond(`touse' & !missing(`v'), `v', .)
                }
                continue
            }
            if `direct_outputs' {
                local output_type : type `v'
                if "`by'" != "" local output_type "`c(type)'"
                * Unfilled storage only when the plugin writes every
                * observation (no if/in, and no weights or nomissing that
                * leave observations out)
                local nofill = (`nobs_sample' == _N)
                mata: (void) st_addvar(st_local("output_type"), st_local("newvar"), `nofill')
            }
            else quietly gen double `newvar' = `v' if `touse' & !missing(`v')
        }
    }

    * In place, the default keeps an undo copy of the targets in the plugin,
    * and their variable labels in the record
    local undo_opt ""
    if `inplace' & "`snapshot'`fast'" == "" {
        local undo_opt "undo"
        _ctools_undo varlabels CTOOLS_u_cwinsor_ `varlist'
    }
    * Validation hooks: fail as Break would once new variables exist, or have
    * the plugin fail after writing half the observations
    if "$CTOOLS_cwinsor_fault" == "outputs" error 1
    if "$CTOOLS_cwinsor_fault" == "store" local undo_opt "`undo_opt' fault_store"

    * Checked SPI stores use positions in the plugin varlist, just like reads.
    * Only target variables need store indices; by-variables are read-only.
    * The targets are consecutive, so the plugin receives only the position
    * before the first one (a per-variable list overflows the plugin argument
    * on wide data).
    local store_indices = "storebase=" + string(cond(`direct_outputs', `nvars' + `nby' + `weighted', 0))

    * Build options
    local opts "p=`p' q=`q'"
    if `direct_outputs' local opts "`opts' newoutputs"
    * Every observation is in the sample without if/in (and without weights
    * or nomissing leaving some out): the plugin then reads no marker
    if `nobs_sample' == _N local opts "`opts' skipif"
    if "`trim'" != "" local opts "`opts' trim"
    if `weighted' local opts "`opts' weighted"
    * winsor2 expands non-grouped bounds through numeric local macros. This
    * rounding affects both clipping and trimming, including float midpoints.
    * Keep C's loaded columns and quickselect across this scalar-only step.
    if "`by'"=="" {
        tempname cwinsor_bounds
        matrix `cwinsor_bounds' = J(`nvars', 2, .)
        local opts "`opts' deferbounds"
    }
    else {
        * winsor2 stores by() cutoffs with egen in the source's type: the
        * plugin truncates them for byte, int and long (t) and rounds them to
        * float for float (f); double cutoffs are used as computed (d).
        local cwinsor_rounding
        foreach v of local varlist {
            local vtype : type `v'
            * With weights the cutoffs are gstats winsor's, used as computed
            if `weighted' local cwinsor_rounding = "`cwinsor_rounding'd"
            else local cwinsor_rounding = "`cwinsor_rounding'" + cond(inlist("`vtype'", "byte", "int", "long"), "t", cond("`vtype'" == "float", "f", "d"))
        }
        local opts "`opts' rounding"
    }

    if "`verbose'" != "" local opts "`opts' verbose"

    local threads_opt ""
    if `threads' > 0 {
        local threads_opt "threads(`threads')"
    }

    * Record pre-plugin time
    quietly _ctools_clock __ct_now
    local __ct_elapsed90 = `__ct_elapsed90' + (`__ct_now' - `__ct_start90')
    * Elapsed seconds are held in caller-local timing state.
    local __time_preplugin = `__ct_elapsed90'
    if "`__ct_elapsed90'" == "" local __ct_elapsed90 = 0
    quietly _ctools_clock __ct_start90

    * Load source + by columns; separate final-type outputs follow them when
    * requested. C uses plugin-local indices for both loading and writing.
    local nvars_target : word count `target_varlist'

    local plugin_varlist `target_varlist' `by' `wvar'
    if `direct_outputs' local plugin_varlist `varlist' `by' `wvar' `target_varlist'
    * Pass the sample marker directly, avoiding serial SF_ifobs evaluation.
    * Keep the same varlist for the deferred-bound application below.
    local plugin_varlist `plugin_varlist' `touse'
    local opts `opts' touse
    plugin call ctools_plugin `plugin_varlist', ///
        "cwinsor `threads_opt' `nvars_target' `nby' `store_indices' `opts' `undo_opt'"
    global CTOOLS_cwinsor_plugin "`__ctools_plugin'"
    if strpos("`undo_opt'", "undo") {
        * The plugin holds the original values until cwinsor commits
        global CTOOLS_cwinsor_targets "`target_varlist'"
    }

    * Include metadata, conversion and cleanup in the post-plugin measurement.
    quietly _ctools_clock __ct_start92

    if "`by'"=="" {
        if "$CTOOLS_cwinsor_fault" == "bounds" error 1
        local bound_index=0
        foreach target of local target_varlist {
            local ++bound_index
            local lower = `cwinsor_bounds'[`bound_index',1]
            local upper = `cwinsor_bounds'[`bound_index',2]
            * winsor2's macro expansion; gstats winsor (weights) uses the
            * cutoffs as computed
            if !`weighted' {
                matrix `cwinsor_bounds'[`bound_index',1] = `lower'
                matrix `cwinsor_bounds'[`bound_index',2] = `upper'
            }
            * winsor2's replace promotes an integer variable that receives a
            * fractional cutoff: byte and int to float, long to double.
            local vtype : type `target'
            if `inplace' & "`trim'" == "" & inlist("`vtype'", "byte", "int", "long") & ///
                (`lower' != floor(`lower') | `upper' != floor(`upper')) {
                local promoted = cond("`vtype'" == "long", "double", "float")
                global CTOOLS_cwinsor_recast "$CTOOLS_cwinsor_recast `target' `vtype'"
                recast `promoted' `target'
            }
        }
        capture noisily plugin call ctools_plugin `plugin_varlist', "cwinsor `threads_opt' apply"
        local apply_rc = _rc
        if `apply_rc' {
            * The plugin restores a failed store before releasing its cache.
            macro drop CTOOLS_cwinsor_targets
            exit `apply_rc'
        }
    }

    * winsor2's label text: -Winsor(p#,p#) or -Trim(p#,p#), a 0 before a
    * lower percentile below 1
    if "`label'" != "" {
        local label_low = `p'
        local label_high = `q'
        if `label_low' < 1 local label_low "0`label_low'"
        local label_note = cond("`trim'" != "", "-Trim", "-Winsor") + "(p`label_low',p`label_high')"
    }

    * gstats winsor's results: groups, the smallest and largest, and the cuts
    return scalar J = _cwinsor_ngroups
    return scalar minJ = _cwinsor_minJ
    return scalar maxJ = _cwinsor_maxJ
    return scalar cutlow = `p'
    return scalar cuthigh = `q'

    * Match winsor2's output metadata and final storage. Keep the working
    * copy double until after computing the bounds, so creation does not round
    * the input before percentile calculation.
    local target_index=0
    foreach v of local varlist {
        local ++target_index
        local target : word `target_index' of `target_varlist'
        local source_label : variable label `v'
        if `"`source_label'"'=="" local source_label "`v'"
        * label: winsor2's note of the cuts
        if "`label'" != "" local source_label `"`source_label'`label_note'"'
        label variable `target' `"`source_label'"'
        if "`target'"!="`v'" & !`direct_outputs' {
            local output_type : type `v'
            if "`by'"!="" & "`trim'"=="" local output_type "`c(type)'"
            recast `output_type' `target', force
            * recast retains the double working copy's format. gen/egen use
            * the default format for the final storage type instead.
            tempvar format_probe
            quietly gen `output_type' `format_probe' = .
            local output_format : format `format_probe'
            format `target' `output_format'
            drop `format_probe'
        }
        * prefix(), generate() and weights are gstats winsor's, whose
        * targets keep the source's display format (winsor2's gen gives the
        * default one)
        if "`target'" != "`v'" & ("`prefix'`generate'" != "" | `weighted') {
            format `target' `: format `v''
        }
    }

    * The plugin changed the variables in place: as with replace, drop Stata's
    * sort marker when it involves one of them (trimmed values sort last, and
    * clamping only the if/in or by() subsets can break the order)
    if "`suffix'" == "" & "`prefix'" == "" {
        local marked : sortedby
        local hit : list marked & target_varlist
        if "`hit'" != "" _cwinsor_clear_sortedby
    }
    * Validation hook: fail after the targets were written and relabeled
    if "$CTOOLS_cwinsor_fault" == "post" error 1
    * Plugin stores do not mark the data changed; replace and winsor2 do
    if `inplace' mata: st_updata(1)

    * Finish timing
    quietly _ctools_clock __ct_now
    local __ct_elapsed90 = `__ct_elapsed90' + (`__ct_now' - `__ct_start90')
    quietly _ctools_clock __ct_now
    local __ct_elapsed92 = `__ct_elapsed92' + (`__ct_now' - `__ct_start92')

    if `__do_timing' {
        * Elapsed seconds are held in caller-local timing state.
        local __time_total = `__ct_elapsed90'
        * Elapsed seconds are held in caller-local timing state.
        local __time_sort = `__ct_elapsed91'
        * Elapsed seconds are held in caller-local timing state.
        local __time_postplugin = `__ct_elapsed92'

        * Calculate pre-plugin time (excluding sort which is timed separately)
        local __time_preplugin_other = `__time_preplugin' - `__time_sort'

        di as text ""
        di as text "{hline 60}"
        di as text "cwinsor timing breakdown:"
        di as text "{hline 60}"
        di as text "  Stata pre-plugin:"
        di as text "    Syntax/validation:      " as result %8.4f `__time_preplugin_other' " sec"
        if `nby' > 0 {
            di as text "    Sort by-variables:      " as result %8.4f `__time_sort' " sec"
        }
        di as text "  {hline 58}"
        di as text "    Pre-plugin total:       " as result %8.4f `__time_preplugin' " sec"
        di as text "  {hline 58}"
        di as text "  C plugin phases:"
        di as text "    Data load:              " as result %8.4f _cwinsor_time_load " sec"
        if `nby' > 0 {
            di as text "    Sort by-variables:      " as result %8.4f _cwinsor_time_sort " sec"
        }
        di as text "    Group detection:        " as result %8.4f _cwinsor_time_groups " sec"
        di as text "    Winsorization:          " as result %8.4f _cwinsor_time_winsor " sec"
        di as text "    Data store:             " as result %8.4f _cwinsor_time_store " sec"
        di as text "  {hline 58}"
        di as text "    C plugin total:         " as result %8.4f _cwinsor_time_total " sec"
        di as text "  {hline 58}"
        di as text "  Stata post-plugin:        " as result %8.4f `__time_postplugin' " sec"
        di as text "{hline 60}"
        di as text "    Wall clock total:       " as result %8.4f `__time_total' " sec"
        di as text "{hline 60}"
        di as text ""
        di as text "  Dataset: " as result `nobs_sample' as text " obs x " as result `nvars_target' as text " vars"
        di as text "  Groups:  " as result %8.0f _cwinsor_ngroups
        di as text "  Mode:    " as result cond("`trim'"!="", "trim", "winsorize") as text " at [" as result `p' as text ", " as result `q' as text "] percentiles"
        di as text "  Threads: " as result %8.0f _cwinsor_threads
        di as text "{hline 60}"
    }

    * Summary
    if "`verbose'" == "" {
        local mode_str = cond("`trim'" != "", "trimmed", "winsorized")
        di as text "(`nvars_target' variables `mode_str' at [`p', `q'] percentiles)"
    }
end

program define _cwinsor_clear_sortedby
    local marked : sortedby
    if "`marked'" == "" | _N == 0 exit
    local v : word 1 of `marked'
    * Stata drops the sort marker when a sort variable changes: change one
    * value with replace, then put the original back through Mata (no copy
    * of the variable is needed)
    capture confirm string variable `v'
    if _rc {
        mata: __ctools_sb_hold = st_data(1, "`v'")
        quietly replace `v' = cond(`v' == 0, 1, 0) in 1
        mata: st_store(1, "`v'", __ctools_sb_hold)
    }
    else {
        mata: __ctools_sb_hold = st_sdata(1, "`v'")
        quietly replace `v' = cond(`v' == "", "x", "") in 1
        mata: st_sstore(1, "`v'", __ctools_sb_hold)
    }
    mata: mata drop __ctools_sb_hold
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

* Failure-handling mode in a ctools wrapper's options.
* Sets local mode in the caller to snapshot, fast or empty (the default: an
* exact undo without a snapshot of the dataset). Options are parsed after the
* first comma outside parentheses and quotes, so if-expressions such as
* inlist(g,1,2) are skipped. A command line syntax cannot parse leaves the
* default; the wrapper's own syntax then reports the error.
program define _ctools_mode
    version 14.1
    gettoken lhs 0 : 0, parse(",") bind
    capture syntax [, SNAPshot FAST *]
    if _rc {
        c_local mode
        exit
    }
    if "`snapshot'" != "" & "`fast'" != "" {
        di as error "options snapshot and fast may not be combined"
        exit 198
    }
    c_local mode `snapshot'`fast'
end

* Exact undo of a ctools wrapper's changes without a dataset snapshot.
* The wrappers that use it only add observations and variables at the end of
* the dataset, add value labels and _dta characteristics or notes, change
* variable labels they name, and lose the sort marker; the plugin restores any
* values it overwrote first. The record lives in globals named <prefix>*.
*
*   _ctools_undo record <prefix>
*       before the first change: _N, c(k), the sort marker, c(changed), the
*       value-label names and the _dta characteristics (notes included)
*   _ctools_undo varlabels <prefix> <varlist>
*       add the variable labels of varlist to the record
*   _ctools_undo rollback <prefix>
*       drop what was added and put the recorded state back (call under
*       nobreak, after the plugin has restored overwritten values)
*   _ctools_undo clear <prefix>
*       drop the record after success
program define _ctools_undo
    version 14.1
    gettoken sub 0 : 0
    gettoken p 0 : 0
    if "`sub'" == "record" _ctools_undo_record `p'
    else if "`sub'" == "varlabels" _ctools_undo_varlabels `p' `0'
    else if "`sub'" == "rollback" _ctools_undo_rollback `p'
    else if "`sub'" == "clear" macro drop `p'*
    else exit 198
end

program define _ctools_undo_record
    args p
    macro drop `p'*
    global `p'N = _N
    global `p'k = c(k)
    global `p'sortedby : sortedby
    global `p'changed = c(changed)
    tempname caller_r
    _return hold `caller_r'
    quietly label dir
    global `p'labels `r(names)'
    _return restore `caller_r'
    local names : char _dta[]
    global `p'dn `names'
    local i = 0
    foreach c of local names {
        local ++i
        mata: st_global("`p'd`i'", st_global("_dta[`c']"))
    }
    global `p'nl = 0
    global `p'active 1
end

program define _ctools_undo_varlabels
    gettoken p 0 : 0
    if "${`p'active}" == "" exit
    local i = ${`p'nl}
    foreach v of local 0 {
        local ++i
        global `p'lv`i' `v'
        mata: st_global("`p'll`i'", st_varlabel("`v'"))
    }
    global `p'nl = `i'
end

program define _ctools_undo_rollback
    args p
    if "${`p'active}" == "" {
        macro drop `p'*
        exit
    }
    local N0 = ${`p'N}
    local k0 = ${`p'k}
    if _N > `N0' quietly drop in `=`N0'+1'/l
    if c(k) > `k0' {
        mata: st_local("added", invtokens(st_varname((`k0'+1)..st_nvar())))
        drop `added'
    }
    quietly label dir
    local now `r(names)'
    local before ${`p'labels}
    local added : list now - before
    if "`added'" != "" label drop `added'
    local now : char _dta[]
    local old ${`p'dn}
    foreach c of local now {
        if !`: list c in old' char _dta[`c']
    }
    local i = 0
    foreach c of local old {
        local ++i
        mata: st_global("_dta[`c']", st_global("`p'd`i'"))
    }
    forvalues i = 1/${`p'nl} {
        mata: st_varlabel(st_global("`p'lv`i'"), st_global("`p'll`i'"))
    }
    local sortedby ${`p'sortedby}
    if "`: sortedby'" != "`sortedby'" {
        * The observations are back in their order: a stable sort only sets
        * the marker again
        if "`sortedby'" != "" sort `sortedby', stable
        else _ctools_undo_unsort
    }
    mata: st_updata(${`p'changed})
    macro drop `p'*
end

* Drop the sort marker without moving data: change one value of the first sort
* variable with replace, which Stata notices, then put it back through Mata.
program define _ctools_undo_unsort
    local marked : sortedby
    if "`marked'" == "" | _N == 0 exit
    local v : word 1 of `marked'
    capture confirm string variable `v'
    if _rc {
        mata: __ctools_undo_hold = st_data(1, "`v'")
        quietly replace `v' = cond(`v' == 0, 1, 0) in 1
        mata: st_store(1, "`v'", __ctools_undo_hold)
    }
    else {
        mata: __ctools_undo_hold = st_sdata(1, "`v'")
        quietly replace `v' = cond(`v' == "", "x", "") in 1
        mata: st_sstore(1, "`v'", __ctools_undo_hold)
    }
    mata: mata drop __ctools_undo_hold
end
* ctools-helpers:end
