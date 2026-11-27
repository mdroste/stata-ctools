#ifndef CSQL_IMPL_H
#define CSQL_IMPL_H
#include "stplugin.h"

/*
 * csql: SQL queries through DuckDB (see csql_internal.h).
 *
 *   csql run [verbose] [labels=0|1] [datetime=0|1] [rows=#] [linesize=#]
 *            [strl=#] [readonly] [copy=parquet|csv|json]
 *       locals: __csql_lib, __csql_sqlfile, __csql_meta, __csql_out,
 *               __csql_using, __csql_database, __csql_memory, __csql_tempdir
 *   csql version      local __csql_lib; sets csql_version
 *   csql install      locals __csql_archive, __csql_sha256, __csql_member, __csql_dest
 */
ST_retcode csql_main(const char *args);

#endif /* CSQL_IMPL_H */
