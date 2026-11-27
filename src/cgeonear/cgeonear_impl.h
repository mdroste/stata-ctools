#ifndef CGEONEAR_IMPL_H
#define CGEONEAR_IMPL_H
#include "stplugin.h"
/* cgeonear: base, using, compute, store and clear (see cgeonear_impl.c) */
ST_retcode cgeonear_main(const char *args);
/* Frees what cgeonear keeps between calls; idempotent */
void cgeonear_cleanup_state(void);
#endif
