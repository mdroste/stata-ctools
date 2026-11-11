* Generate the corpus with validation.support.parquet_native_fixtures.corpus.
local fx : environment CTOOLS_PARQUET_CORPUS
assert `"`fx'"' != ""
do "validation/support/validate_setup.do"
sysuse auto, clear
tempfile original
save `original'
local bad : dir `"`fx'"' files "bad_*.parquet"
assert `: word count `bad'' == 19
foreach name of local bad {
    local expected=691
    if inlist("`name'", "bad_magic.parquet", "bad_footer_size.parquet", "bad_truncated.parquet", "bad_row_total.parquet", "bad_row_overflow.parquet", "bad_negative_rows.parquet") local expected=3900
    if "`name'"=="bad_huge_rows.parquet" local expected=920
    foreach how in "" snapshot {
        capture noisily cimport parquet using `"`fx'/`name'"', clear `how'
        local rc=_rc
        di "malformed `name' `how': r(`rc'), expected r(`expected')"
        assert `rc'==`expected'
        cf _all using `original', all
        assert c(k)==12 & _N==74
        test_pass "malformed Parquet `name' preserves data `how'"
    }
}
print_summary "complete malformed Parquet rollback"
if $TESTS_FAILED > 0 exit 9
