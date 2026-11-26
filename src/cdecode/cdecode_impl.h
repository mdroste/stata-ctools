/* Retired label-save protocol; both entry points reject stale ado wrappers. */
#ifndef CDECODE_IMPL_H
#define CDECODE_IMPL_H
#include "stplugin.h"
ST_retcode cdecode_main(const char *args);
ST_retcode cdecode_scan_main(const char *args);

/* A value-label mapping in cdecode's wire format (see cdecode_impl.c), for
 * commands that write label text themselves. create takes ownership of
 * text, a malloc'd copy of the maplen-byte mapping with a terminator, whose
 * label texts are at most width bytes; it returns NULL and sets *rc to a
 * Stata return code when the mapping is malformed or memory runs out.
 * text returns the label of a stored value, or NULL when it has none. */
void *cdecode_label_map_create(char *text, size_t maplen, size_t width, ST_retcode *rc);
const char *cdecode_label_map_text(const void *map, double value);
void cdecode_label_map_free(void *map);
#endif
