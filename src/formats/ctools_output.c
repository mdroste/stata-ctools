#include "ctools_output.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <fcntl.h>
#ifdef _WIN32
#include <io.h>
#include <sys/stat.h>
#include <windows.h>
#else
#include <unistd.h>
#endif

/* Reserve a sibling file; publishing is a separate, atomic operation. */
ST_retcode ctools_output_prepare(ctools_output *out, const char *filename, bool replace)
{
    memset(out, 0, sizeof(*out));
    out->replace = replace;
    out->destination = strdup(filename);
    out->temporary = malloc(strlen(filename) + 20);
    if (!out->destination || !out->temporary)
    {
        ctools_output_cleanup(out);
        return 920;
    }
    snprintf(out->temporary, strlen(filename) + 20, "%s.ctools.XXXXXX", filename);
#ifdef _WIN32
    int fd = -1;
    if (_mktemp_s(out->temporary, strlen(out->temporary) + 1) == 0)
        fd =
            _open(out->temporary, _O_CREAT | _O_EXCL | _O_BINARY | _O_WRONLY, _S_IREAD | _S_IWRITE);
    if (fd >= 0)
        _close(fd);
#else
    int fd = mkstemp(out->temporary);
    if (fd >= 0)
        close(fd);
#endif
    if (fd < 0)
    {
        ctools_output_cleanup(out);
        return 603;
    }
    out->owns_temporary = true;
    return 0;
}

ST_retcode ctools_output_commit(ctools_output *out)
{
#ifdef _WIN32
    if (!MoveFileExA(out->temporary, out->destination,
                     MOVEFILE_WRITE_THROUGH | (out->replace ? MOVEFILE_REPLACE_EXISTING : 0)))
        return GetLastError() == ERROR_ALREADY_EXISTS || GetLastError() == ERROR_FILE_EXISTS ? 602
                                                                                             : 603;
#else
    if (out->replace)
    {
        if (rename(out->temporary, out->destination))
            return 603;
    }
    else
    {
        /* link() refuses an existing destination, including a racing writer. */
        if (link(out->temporary, out->destination))
            return errno == EEXIST ? 602 : 603;
        unlink(out->temporary);
    }
#endif
    return 0;
}

void ctools_output_cleanup(ctools_output *out)
{
    if (out->owns_temporary)
        remove(out->temporary);
    free(out->temporary);
    free(out->destination);
    memset(out, 0, sizeof(*out));
}
