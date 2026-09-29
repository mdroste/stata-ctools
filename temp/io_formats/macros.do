clear all
adopath ++ "build"
log using "temp/io_formats/macros.log", text replace
program define check
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    local __cio_filename "temp/io_formats/points.shp"
    plugin call ctools_plugin, "cio scan shp"
    macro list __cio_nobs __cio_shp_type __cio_shp_xmin
    mata: st_global("__cio_nobs"),st_global("__cio_shp_type")
    plugin call ctools_plugin, "cio blob 4 1 0"
    macro list __cio_blobhex
    mata: st_global("__cio_blobhex")
    plugin call ctools_plugin, "cio clear"
end
check
log close
exit, clear
