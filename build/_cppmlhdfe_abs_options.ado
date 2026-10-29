*! version 1.1.0 20261016 github.com/mdroste/stata-ctools
*! Private helper for cppmlhdfe.
program define _cppmlhdfe_abs_options
    syntax [, SAVEfe GENerate KEEPSINgletons TOLerance(string) ITERate(string) DOFadjustments(string) ACCELeration(string) TRANSform(string) POOLsize(string)]
    if "`acceleration'"!="" c_local acceleration `acceleration'
    if "`transform'"!="" c_local transform `transform'
    if "`poolsize'"!="" {
        confirm integer number `poolsize'
        if `poolsize'<0 | missing(`poolsize') error 198
        c_local poolsize `poolsize'
    }
    c_local save_all_fe = ("`savefe'`generate'"!="")
    if "`keepsingletons'"!="" c_local keepsingletons keepsingletons
    if "`tolerance'"!="" {
        confirm number `tolerance'
        if `tolerance'<=0 | missing(`tolerance') error 198
        c_local itolerance `tolerance'
    }
    if "`iterate'"!="" {
        confirm integer number `iterate'
        if `iterate'<1 | missing(`iterate') error 198
        c_local iterate `iterate'
    }
    if "`dofadjustments'"!="" c_local dofadjustments `dofadjustments'
end
