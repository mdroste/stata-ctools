#ifndef CSQL_API_H
#define CSQL_API_H
/*
 * DuckDB through its C API, loaded at run time.
 *
 * ctools does not link DuckDB. csql opens the official DuckDB shared library
 * (duckdb.dll, libduckdb.so, libduckdb.dylib; `csql, install` downloads the
 * pinned release) by its full path and resolves the C API functions listed
 * below. The plugin stays a C build, every other command is unaffected when
 * the library is absent, and a newer DuckDB can be used without rebuilding
 * ctools: the C API is stable across releases.
 *
 * Types come from DuckDB's own header (src/vendor/duckdb/duckdb.h, MIT);
 * DUCKDB_STATIC_BUILD only stops it declaring the functions dllimport.
 */
#ifndef DUCKDB_STATIC_BUILD
#define DUCKDB_STATIC_BUILD
#endif
#include "vendor/duckdb/duckdb.h"

/* X(name, required): required functions exist in DuckDB 1.1 and later. */
#define CSQL_API_FUNCTIONS(X)                                                                      \
    X(duckdb_library_version, 1)                                                                   \
    X(duckdb_free, 1)                                                                              \
    X(duckdb_vector_size, 1)                                                                       \
    X(duckdb_create_config, 1)                                                                     \
    X(duckdb_set_config, 1)                                                                        \
    X(duckdb_destroy_config, 1)                                                                    \
    X(duckdb_open_ext, 1)                                                                          \
    X(duckdb_close, 1)                                                                             \
    X(duckdb_connect, 1)                                                                           \
    X(duckdb_disconnect, 1)                                                                        \
    X(duckdb_interrupt, 1)                                                                         \
    X(duckdb_extract_statements, 1)                                                                \
    X(duckdb_extract_statements_error, 1)                                                          \
    X(duckdb_prepare_extracted_statement, 1)                                                       \
    X(duckdb_destroy_extracted, 1)                                                                 \
    X(duckdb_prepare, 1)                                                                           \
    X(duckdb_prepare_error, 1)                                                                     \
    X(duckdb_destroy_prepare, 1)                                                                   \
    X(duckdb_prepared_statement_type, 1)                                                           \
    X(duckdb_prepared_statement_column_count, 0)                                                   \
    X(duckdb_prepared_statement_column_type, 0)                                                    \
    X(duckdb_prepared_statement_column_name, 0)                                                    \
    X(duckdb_pending_prepared, 1)                                                                  \
    X(duckdb_pending_prepared_streaming, 0)                                                        \
    X(duckdb_pending_execute_task, 1)                                                              \
    X(duckdb_pending_execute_check_state, 1)                                                       \
    X(duckdb_pending_execution_is_finished, 1)                                                     \
    X(duckdb_pending_error, 1)                                                                     \
    X(duckdb_execute_pending, 1)                                                                   \
    X(duckdb_destroy_pending, 1)                                                                   \
    X(duckdb_result_error, 1)                                                                      \
    X(duckdb_destroy_result, 1)                                                                    \
    X(duckdb_result_return_type, 1)                                                                \
    X(duckdb_rows_changed, 1)                                                                      \
    X(duckdb_column_count, 1)                                                                      \
    X(duckdb_column_name, 1)                                                                       \
    X(duckdb_column_logical_type, 1)                                                               \
    X(duckdb_fetch_chunk, 1)                                                                       \
    X(duckdb_result_chunk_count, 0)                                                                \
    X(duckdb_result_get_chunk, 0)                                                                  \
    X(duckdb_data_chunk_get_size, 1)                                                               \
    X(duckdb_data_chunk_get_vector, 1)                                                             \
    X(duckdb_data_chunk_set_size, 1)                                                               \
    X(duckdb_destroy_data_chunk, 1)                                                                \
    X(duckdb_vector_get_data, 1)                                                                   \
    X(duckdb_vector_get_validity, 1)                                                               \
    X(duckdb_vector_ensure_validity_writable, 1)                                                   \
    X(duckdb_vector_assign_string_element_len, 1)                                                  \
    X(duckdb_validity_set_row_invalid, 1)                                                          \
    X(duckdb_get_type_id, 1)                                                                       \
    X(duckdb_create_logical_type, 1)                                                               \
    X(duckdb_destroy_logical_type, 1)                                                              \
    X(duckdb_decimal_width, 1)                                                                     \
    X(duckdb_decimal_scale, 1)                                                                     \
    X(duckdb_decimal_internal_type, 1)                                                             \
    X(duckdb_enum_internal_type, 1)                                                                \
    X(duckdb_enum_dictionary_size, 1)                                                              \
    X(duckdb_enum_dictionary_value, 1)                                                             \
    X(duckdb_create_varchar, 1)                                                                    \
    X(duckdb_get_varchar, 1)                                                                       \
    X(duckdb_get_bool, 1)                                                                          \
    X(duckdb_get_int64, 1)                                                                         \
    X(duckdb_destroy_value, 1)                                                                     \
    X(duckdb_create_table_function, 1)                                                             \
    X(duckdb_destroy_table_function, 1)                                                            \
    X(duckdb_table_function_set_name, 1)                                                           \
    X(duckdb_table_function_add_parameter, 1)                                                      \
    X(duckdb_table_function_add_named_parameter, 1)                                                \
    X(duckdb_table_function_set_extra_info, 1)                                                     \
    X(duckdb_table_function_set_bind, 1)                                                           \
    X(duckdb_table_function_set_init, 1)                                                           \
    X(duckdb_table_function_set_local_init, 1)                                                     \
    X(duckdb_table_function_set_function, 1)                                                       \
    X(duckdb_table_function_supports_projection_pushdown, 1)                                       \
    X(duckdb_register_table_function, 1)                                                           \
    X(duckdb_bind_get_extra_info, 1)                                                               \
    X(duckdb_bind_add_result_column, 1)                                                            \
    X(duckdb_bind_get_parameter_count, 1)                                                          \
    X(duckdb_bind_get_parameter, 1)                                                                \
    X(duckdb_bind_get_named_parameter, 1)                                                          \
    X(duckdb_bind_set_bind_data, 1)                                                                \
    X(duckdb_bind_set_cardinality, 1)                                                              \
    X(duckdb_bind_set_error, 1)                                                                    \
    X(duckdb_init_get_extra_info, 1)                                                               \
    X(duckdb_init_get_bind_data, 1)                                                                \
    X(duckdb_init_set_init_data, 1)                                                                \
    X(duckdb_init_get_column_count, 1)                                                             \
    X(duckdb_init_get_column_index, 1)                                                             \
    X(duckdb_init_set_max_threads, 1)                                                              \
    X(duckdb_init_set_error, 1)                                                                    \
    X(duckdb_function_get_extra_info, 1)                                                           \
    X(duckdb_function_get_bind_data, 1)                                                            \
    X(duckdb_function_get_init_data, 1)                                                            \
    X(duckdb_function_get_local_init_data, 1)                                                      \
    X(duckdb_function_set_error, 1)                                                                \
    X(duckdb_add_replacement_scan, 1)                                                              \
    X(duckdb_replacement_scan_set_function_name, 1)                                                \
    X(duckdb_replacement_scan_add_parameter, 1)

typedef struct
{
#define CSQL_API_FIELD(name, required) __typeof__(name) *name;
    CSQL_API_FUNCTIONS(CSQL_API_FIELD)
#undef CSQL_API_FIELD
    char path[4096];   /* the library loaded */
    char version[64];  /* duckdb_library_version() */
    int loaded;
} csql_api;

/* The loaded library; valid after csql_api_load() returned 0. */
extern csql_api csql_db;

/*
 * Load the DuckDB library at path (UTF-8). The library stays loaded for the
 * rest of the Stata session (DuckDB keeps thread-local and static state);
 * a later call with another path loads that one instead. Returns 0, or 601
 * (cannot load the file) or 610 (not a usable DuckDB: a required function
 * is missing) with the reason in msg.
 */
int csql_api_load(const char *path, char *msg, size_t cap);

#endif /* CSQL_API_H */
