#ifndef CTOOLS_OUTPUT_H
#define CTOOLS_OUTPUT_H
#include <stdbool.h>
#include "stplugin.h"
typedef struct
{
    char *destination;
    char *temporary;
    bool replace;
    bool owns_temporary;
} ctools_output;

ST_retcode ctools_output_prepare(ctools_output *out, const char *filename, bool replace);
ST_retcode ctools_output_commit(ctools_output *out);
void ctools_output_cleanup(ctools_output *out);
#endif
