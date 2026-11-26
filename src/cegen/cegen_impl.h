/*
    cegen_impl.h
    Extensions to generate: the functions of egen in C

    cegen runs in two plugin calls, so the new variable can be created with
    its final storage type (a string as wide as its longest value, a promoted
    or automatically chosen numeric type) before anything is written:

    plugin call ctools_plugin <sources...> <by...> [<touse>] [in],
        "cegen compute fn=<function> nsrc=# nby=# type=<storage type> ..."

        Loads the sample, computes the function and keeps the result, one
        value per loaded observation. Reports through locals of the caller
        (cegen_n, cegen_ngroups, cegen_kind, cegen_strlen, cegen_isint,
        cegen_min, cegen_max, ... and the phase times).

    plugin call ctools_plugin <newvar>, "cegen store type=<storage type> ..."

        Writes the kept result into the new variable, every observation of
        the dataset: the result where the sample has one, the function's
        fill value elsewhere. Reports cegen_nmiss, the number of missing
        values in the new variable.

    plugin call ctools_plugin, "cegen discard"

        Frees a kept result that will not be stored.

    License: MIT
*/

#ifndef CEGEN_IMPL_H
#define CEGEN_IMPL_H

#include "stplugin.h"

/*
    Main entry point for cegen.

    @param args  "compute ...", "store ..." or "discard"
    @return      0 on success, a Stata error code on failure
                 (198 invalid request, 920 no memory, 109 a value the
                 loaders cannot read)
*/
ST_retcode cegen_main(const char *args);

/* Free the result kept between compute and store. Idempotent. */
void cegen_cleanup_state(void);

#endif /* CEGEN_IMPL_H */
