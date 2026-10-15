#ifndef CTOOLS_TEXT_H
#define CTOOLS_TEXT_H
#include "ctools_data_io.h"
#include "stplugin.h"
/* Text only: binary strL is rejected. Zero limits mean unbounded; Excel uses
 * 32767 UTF-16 units, CSV uses neither limit. SPI reads run on the caller thread.
 * Each str# column packs its values into one block owned by the column's arena;
 * strL values are allocated one by one at their SF_sdatalen size.
 * On error fd owns any partially loaded data; the caller must free it. */
typedef struct
{
    size_t max_bytes, max_utf16_units;
} ctools_text_policy;
ST_retcode ctools_load_text_data(ctools_filtered_data *fd, ST_int nvars, ctools_text_policy policy);
#endif
