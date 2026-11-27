/*
    cfuzzyjoin_impl.h
    C-accelerated fuzzy joins: record linkage (reclink), similarity joins of
    two datasets or two columns (matchit) and grouping of similar strings
    (strgroup)

    Plugin calls (see cfuzzyjoin_impl.c):
        match_master, match_weights, match_using, match_store, match_columns
        link_master, link_exclude, link_using, link_layout, link_store
        group
        clear
*/

#ifndef CFUZZYJOIN_IMPL_H
#define CFUZZYJOIN_IMPL_H

#include "stplugin.h"

ST_retcode cfuzzyjoin_main(const char *args);

/* Idempotent: frees the staged data, indexes and results. */
void cfuzzyjoin_cleanup_state(void);

#endif /* CFUZZYJOIN_IMPL_H */
