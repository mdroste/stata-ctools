/*
    cjoinby_impl.h
    C-accelerated joinby: all pairwise combinations within groups

    Plugin calls, in order (see cjoinby_impl.c):
        using     stage the using data and index its keys
        plan      match the master keys, count and lay out the output
        write     write the joined observations
        clear     release everything (the join stands)
        rollback  write the master observations back, then release
*/

#ifndef CJOINBY_IMPL_H
#define CJOINBY_IMPL_H

#include "stplugin.h"

ST_retcode cjoinby_main(const char *args);

/* Idempotent: frees the staged data, the key index and the output plan. */
void cjoinby_cleanup_cache(void);

#endif /* CJOINBY_IMPL_H */
