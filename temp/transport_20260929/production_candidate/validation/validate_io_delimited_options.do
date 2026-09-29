do "validation/io_test_helpers.do"
args root
if `"`root'"' == "" local root "`c(pwd)'/temp/io_delimited_options"
foreach opt in `"delimiters("|;",collapse)"' `"delimiters("||",asstring)"' `"delimiters("|;")"' `"delimiters("|;") collapsedelimiters"' {
    cio_import_test using "`root'/multidelim.txt", kind(delimited) name("delimiter case") opts(varnames(nonames) `macval(opt)')
}
cio_import_test using "`root'/long.csv", kind(delimited) name("CSV long Unicode strL") opts(encoding(utf8))
cio_import_test using "`root'/long.csv", kind(delimited) name("CSV long Unicode favorstrfixed") opts(encoding(utf8) favorstrfixed)
foreach names in "one" "one two" {
    cio_import_test using "`root'/multiline.csv", kind(delimited) name("CSV explicit names `names'") select(`names')
}
foreach limit in 1 2 3 20 unlimited {
    capture noisily import delimited using "`root'/multiline.csv", bindquotes(strict) maxquotedrows(`limit') clear
    local nr=_rc
    capture noisily cimport delimited using "`root'/multiline.csv", bindquotes(strict) maxquotedrows(`limit') clear
    local cr=_rc
    local rc=cond(`nr'==`cr',0,9)
    cio_record "maxquotedrows(`limit') error" `rc'
}

foreach opt in "varnames(1) rowrange(3:4)" "varnames(1) rowrange(2:3)" "rowrange(3:4)" "varnames(1) rowrange(4)" "varnames(1) rowrange(1:1)" {
    cio_import_test using "`root'/range.csv", kind(delimited) name("CSV range inference `opt'") opts(`opt')
}
import delimited using "`root'/long.csv", encoding(utf8) clear
cio_export_test, kind(delimited) ext(csv) name("CSV long Unicode writer") opts(quote) iopts(encoding(utf8))
import delimited using "`root'/long.csv", encoding(utf8) clear
cio_export_test, kind(delimited) ext(csv) name("CSV long Unicode filtered writer") select(text if id==1) iopts(encoding(utf8))
foreach enc in utf-32le UTF-32BE utf-16be utf16 windows-1251 cp1251 SHIFT-JIS shift_jis sjis ISO-8859-2 latin2 {
    local f cp1251
    if strpos(lower("`enc'"),"32le") local f utf32le
    if strpos(lower("`enc'"),"32be") local f utf32be
    if strpos(lower("`enc'"),"16") local f utf16be
    if strpos(lower("`enc'"),"jis") | "`enc'"=="sjis" local f sjis
    if strpos(lower("`enc'"),"8859") | "`enc'"=="latin2" local f latin2
    cio_import_test using "`root'/`f'.csv", kind(delimited) name("CSV encoding `enc'") opts(encoding(`enc'))
}
cio_import_test using "`root'/utf32bom.csv", kind(delimited) name("CSV UTF-32 BOM detection")
foreach opt in "encoding(bogus)" "encoding(utf32le)" "encoding(utf16be)" "locale(bogus)" "charset(cp1252) encoding(utf8)" {
    capture noisily import delimited using "`root'/range.csv", `opt' clear
    local nr=_rc
    capture noisily cimport delimited using "`root'/range.csv", `opt' clear
    local cr=_rc
    cio_record "CSV option error `opt'" `=cond(`nr'==`cr',0,9)'
}
foreach width in 20 40 50 100 200 1000 {
    foreach pct in 1 28 29 34 35 50 63 64 65 66 67 70 90 93 94 95 {
        cio_import_test using "`root'/storagew`width'p`pct'.csv", kind(delimited) name("CSV storage width `width' fraction `pct'")
    }
    cio_import_test using "`root'/storagew`width'p1.csv", kind(delimited) name("CSV favorstrfixed width `width'") opts(favorstrfixed)
}
foreach pair in "us en_US" "german de" "german de_DE" "french fr_FR" "french_nbsp fr_FR" "russian ru_RU" "swiss de_CH" "arabic ar_EG" "arabic_minus ar_EG" "persian fa_IR" "persian_minus fa_IR" "hindi hi_IN" "bad_sign en_US" "scientific en_US" "infinity en_US" "bad_groups en_US" "groups en_US" "trailing_group en_US" {
    gettoken file locale : pair
    local locale=strtrim("`locale'")
    cio_import_test using "`root'/`file'.csv", kind(delimited) name("CSV locale `file' `locale'") opts(delimiters(";") parselocale(`locale') encoding(utf8))
}
foreach opt in "stripquotes(yes)" "stripquotes(no)" "stripquotes(default)" "bindquotes(nobind)" "rowrange(f:l)" "colrange(f:l)" {
    cio_import_test using "`root'/quotes.csv", kind(delimited) name("CSV quotes/ranges `opt'") opts(`opt')
}
foreach opt in "maxquotedrows(1.5)" "maxquotedrows(bogus)" "rowrange(0)" "rowrange(4:2)" "rowrange(:)" "colrange(0)" "colrange(4:2)" "parselocale(bogus)" "parselocale(en_XX)" {
    capture noisily import delimited using "`root'/range.csv", `opt' clear
    local nr=_rc
    capture noisily cimport delimited using "`root'/range.csv", `opt' clear
    local cr=_rc
    cio_record "CSV invalid range/locale `opt'" `=cond(`nr'==`cr',0,9)'
}
foreach opt in "stripquotes(yes)" "stripquotes(no)" "stripquotes(default)" "bindquotes(nobind)" "bindquotes(nobind) stripquotes(yes)" "numericcols(2) stripquotes(no)" {
    cio_import_test using "`root'/quote_numbers.csv", kind(delimited) name("CSV numeric quotes `opt'") opts(`opt')
}
foreach enc in utf8 UTF-8 windows-1252 cp1252 latin1 Latin1 LATIN1 latin2 latin9 ascii US-ASCII {
    cio_import_test using "`root'/range.csv", kind(delimited) name("CSV returned encoding `enc'") opts(encoding(`enc'))
}
foreach opt in "" "encoding(utf8)" "encoding(latin1)" {
    cio_import_test using "`root'/empty.csv", kind(delimited) name("CSV empty returned results `opt'") opts(`opt')
}
foreach opt in `"delimiters("tab")"' "delimiters(tab)" `"delimiters("comma")"' "delimiters(comma)" `"delimiters("space")"' "delimiters(space)" "delimiters(whitespace)" {
    cio_import_test using "`root'/words.txt", kind(delimited) name("CSV literal and shortcut delimiter") opts(`opt')
}
foreach kind in tab semi pipe colon space {
    cio_import_test using "`root'/autodelim_`kind'.csv", kind(delimited) name("CSV automatic delimiter `kind'")
}
foreach kind in western central cyrillic koi8 arabic greek hebrew turkish japanese eucjp korean chinese traditional iso2022 {
    cio_import_test using "`root'/detect_`kind'.csv", kind(delimited) name("CSV automatic regional encoding `kind'")
}
foreach kind in utf32le utf32be utf16be cp1251 sjis latin2 {
    cio_import_test using "`root'/`kind'.csv", kind(delimited) name("CSV automatic short encoding `kind'")
}
foreach opt in "encoding(bogus)" "delimiters(foo)" "delimiters(|)" {
    capture noisily import delimited using "`root'/empty.csv", `opt' clear
    local nr=_rc
    capture noisily cimport delimited using "`root'/empty.csv", `opt' clear
    local cr=_rc
    cio_record "CSV empty option error `opt'" `=cond(`nr'==`cr',0,9)'
}
di "CSV_OPTIONS_PASSED=$IO_FORMATS_PASSED FAILED=$IO_FORMATS_FAILED"
if $IO_FORMATS_FAILED exit 9
