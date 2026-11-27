/*******************************************************************************
 * cgeonear_helpers.do - Mata helpers of validation/stata/cgeonear.do
 ******************************************************************************/

* Metadata signature of the data in memory: data label, dataset and variable
* characteristics (in name order), types, formats, labels, value-label
* definitions, and optionally the sort marker
capture mata: mata drop _cgn_chars()
capture mata: mata drop _cgn_sig()
mata:
string scalar _cgn_chars(string scalar name)
{
    string rowvector c
    string scalar s
    real scalar k

    c = tokens(st_macroexpand("`" + ": char " + name + "[]" + "'"))
    if (cols(c)) c = sort(c', 1)'
    s = ""
    for (k = 1; k <= cols(c); k++) s = s + c[k] + "=" + st_global(name + "[" + c[k] + "]") + ";"
    return(s)
}

string scalar _cgn_sig(real scalar withsort)
{
    string scalar s
    string colvector labs, text
    real colvector values
    real scalar j

    s = "label=" + st_macroexpand("`" + ": data label" + "'") + "|dta " + _cgn_chars("_dta")
    for (j = 1; j <= st_nvar(); j++) {
        s = s + "|var " + st_varname(j) + " " + st_vartype(j) + " " + st_varformat(j) +
            " [" + st_varlabel(j) + "] (" + st_varvaluelabel(j) + ") " + _cgn_chars(st_varname(j))
    }
    stata("quietly label dir")
    labs = tokens(st_global("r(names)"))'
    if (rows(labs)) labs = sort(labs, 1)
    for (j = 1; j <= rows(labs); j++) {
        st_vlload(labs[j], values, text)
        s = s + "|vl " + labs[j] + ":" + invtokens(strofreal(values', "%21x")) + "=" + invtokens(text')
    }
    if (withsort) s = s + "|sorted " + st_macroexpand("`" + ": sortedby" + "'")
    return(s)
}
end
