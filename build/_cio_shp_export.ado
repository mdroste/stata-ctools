*! version 1.2.4 20261103 github.com/mdroste/stata-ctools
program define _cio_shp_export, rclass
    version 14.1
    syntax [anything] [using/] [, REPLACE SHX ID(varname numeric)]
    if `"`using'"'=="" gettoken using extra : anything
    if `"`using'"'=="" exit 198
    mata: st_local("base",pathbasename(st_local("using")))
    if !strpos(`"`base'"',".") local using `"`using'.shp"'
    if "`id'"=="" local id _ID
    confirm numeric variable `id' _X _Y
    confirm string variable rec_header
    * The characteristic may come from someone else's file: read it in Mata
    * and accept only a number, never re-expanding its text.
    mata: st_local("shape",strtrim(st_global("_dta[shp_type]")))
    mata: st_local("shape",st_local("shape")=="" ? "" : strofreal(strtoreal(st_local("shape")),"%12.0g"))
    if "`shape'"=="" | "`shape'"=="." exit 198
    local __cio_shp_type "`shape'"
    foreach field in xmin ymin xmax ymax zmin zmax mmin mmax {
        local __cio_shp_`field' : char _dta[shp_`field']
    }
    local __cio_filename `"`using'"'
    local __cio_replace = ("`replace'"!="")
    local __cio_shx = ("`shx'"!="")
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc!=0 & _rc!=110 exit 601
    local variables `id' _X _Y rec_header
    if inlist(`shape',11,13,15,18,31) local variables `variables' _Z
    if inlist(`shape',11,13,15,18,21,23,25,28,31) local variables `variables' _M
    capture noisily plugin call ctools_plugin `variables', "cio export shp"
    * A write-protected destination (plugin r(608)) gives native's code and text
    if _rc == 608 {
        di as error "could not copy shp file"
        exit 602
    }
    if _rc exit _rc
    return local filename `"`using'"'
end
