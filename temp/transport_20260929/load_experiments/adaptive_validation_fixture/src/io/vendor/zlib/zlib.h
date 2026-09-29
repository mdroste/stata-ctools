/*
 * zlib.h - zlib-compatible API backed by the vendored miniz.
 *
 * ReadStat's SPSS .zsav code includes <zlib.h>. Rather than vendoring zlib,
 * this header maps the zlib names ReadStat uses (z_stream, deflateInit,
 * deflate, deflateEnd, deflateBound, uncompress, uLongf and the Z_*
 * constants) onto miniz, which ctools already compiles for XLSX I/O.
 *
 * The build puts this directory on the include path, so it is found before
 * any system zlib.h.
 */
#ifndef CTOOLS_ZLIB_SHIM_H
#define CTOOLS_ZLIB_SHIM_H

#include "../../../cimport/miniz/miniz.h"

#endif /* CTOOLS_ZLIB_SHIM_H */
