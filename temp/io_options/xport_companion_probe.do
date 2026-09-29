clear all
adopath ++ "build"
log using "temp/io_options/xport_companion_probe.log", text replace
foreach companion in fixture empty {
    if "`companion'"=="empty" {
        clear
        set obs 1
        gen byte x=1
        drop in 1
        export sasxport8 using "temp/io_allformats/native/empty.v8xpt", replace
    }
    capture noisily import sasxport8 using "temp/io_allformats/native/fixture.v8xpt", vlabfile("temp/io_allformats/native/`companion'.v8xpt") clear
    di "PROBE_`companion'_NATIVE_RC=" _rc
    capture noisily cimport sasxport8 using "temp/io_allformats/native/fixture.v8xpt", vlabfile("temp/io_allformats/native/`companion'.v8xpt") clear
    di "PROBE_`companion'_C_RC=" _rc
}
log close
exit, clear
