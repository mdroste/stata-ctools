#ifndef CTOOLS_OUTPUT_H
#define CTOOLS_OUTPUT_H
#include <stdbool.h>
#include <stddef.h>
#include "stplugin.h"
typedef struct
{
    char *destination;
    char *temporary;
    bool replace;
    bool owns_temporary;
    unsigned mode;
    void *creation_acl; /* macOS ACL retained while the temporary is private */
} ctools_output;

ST_retcode ctools_output_prepare(ctools_output *out, const char *filename, bool replace);
ST_retcode ctools_output_commit(ctools_output *out);
ST_retcode ctools_output_commit_group(ctools_output **outputs, size_t count);
void ctools_output_cleanup(ctools_output *out);
#endif
