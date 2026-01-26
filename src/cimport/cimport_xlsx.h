/*
 * cimport_xlsx.h
 * High-performance XLSX import for Stata
 *
 * Part of the ctools suite - provides Excel import functionality
 * matching Stata's import excel command syntax and output.
 */

#ifndef CIMPORT_XLSX_H
#define CIMPORT_XLSX_H

#include "stplugin.h"
#include "cimport_context.h"
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

/* ============================================================================
 * Constants
 * ============================================================================ */

#define XLSX_MAX_SHEETS         256
#define XLSX_MAX_SHEET_NAME     64
#define XLSX_MAX_SHARED_STRINGS (1024 * 1024)  /* 1M strings max */
#define XLSX_MAX_COLUMNS        16384          /* Excel max columns */

/* ============================================================================
 * Cell Range Structure
 * ============================================================================ */

typedef struct {
    int start_col;    /* 0-based, -1 means from beginning */
    int start_row;    /* 1-based, -1 means from beginning */
    int end_col;      /* 0-based, -1 means to end */
    int end_row;      /* 1-based, -1 means to end */
} XLSXCellRange;

/* ============================================================================
 * Sheet Info Structure
 * ============================================================================ */

typedef struct {
    char name[XLSX_MAX_SHEET_NAME];
    int sheet_id;
    char rel_id[32];  /* Relationship ID (rId1, rId2, etc.) */
    int sheet_index;  /* 1-based index for worksheet filename */
} XLSXSheetInfo;

/* ============================================================================
 * Cell Data Structure
 * ============================================================================ */

typedef enum {
    XLSX_CELL_EMPTY = 0,
    XLSX_CELL_NUMBER,
    XLSX_CELL_STRING,
    XLSX_CELL_SHARED_STRING,
    XLSX_CELL_BOOLEAN,
    XLSX_CELL_ERROR,
    XLSX_CELL_DATE
} XLSXCellType;

typedef struct {
    int col;              /* 0-based column index */
    int row;              /* 1-based row index */
    XLSXCellType type;
    union {
        double number;
        int shared_string_idx;
        const char *inline_string;
        bool boolean;
    } value;
    int style_idx;        /* Style index for date detection */
} XLSXCell;

/* ============================================================================
 * Parsed Row Structure
 * ============================================================================ */

typedef struct {
    int row_num;          /* 1-based row number */
    XLSXCell *cells;
    int num_cells;
    int capacity;
} XLSXParsedRow;

/* ============================================================================
 * XLSX Context Structure
 * ============================================================================ */

typedef struct XLSXContext {
    /* File info */
    char *filename;
    void *zip_archive;        /* xlsx_zip_archive handle */

    /* Sheet metadata */
    XLSXSheetInfo sheets[XLSX_MAX_SHEETS];
    int num_sheets;
    int selected_sheet;       /* 0-based index of sheet to import */

    /* Shared strings table */
    char **shared_strings;
    uint32_t num_shared_strings;
    uint32_t shared_strings_capacity;
    char *shared_strings_pool;      /* Memory pool for string data */
    size_t shared_strings_pool_size;
    size_t shared_strings_pool_used;

    /* Date format detection (style indices that are dates) */
    bool *date_styles;
    int num_styles;

    /* Parsed worksheet data */
    XLSXParsedRow *rows;
    int num_rows;
    int rows_capacity;
    int max_col;              /* Maximum column index seen (0-based) */
    int min_row;              /* Minimum row number (1-based) */
    int max_row;              /* Maximum row number (1-based) */

    /* Column metadata (reuse CImportColumnInfo from CSV) */
    CImportColumnInfo *columns;
    int num_columns;

    /* Column caches for Stata output */
    CImportColumnCache *col_cache;
    bool cache_ready;

    /* Options */
    XLSXCellRange cell_range;
    bool firstrow;            /* First row contains variable names */
    bool allstring;           /* Import all as strings */
    int case_mode;            /* 0=preserve, 1=lower, 2=upper */
    bool verbose;

    /* Timing */
    double time_zip;
    double time_shared_strings;
    double time_worksheet;
    double time_type_infer;
    double time_cache;

    /* Error handling */
    char error_message[256];
} XLSXContext;

/* ============================================================================
 * Public API
 * ============================================================================ */

/*
 * Main entry point for xlsx import.
 * Called from cimport_main when .xlsx file detected.
 *
 * args format: "scan|load filename [options]"
 * Options: sheet=name, cellrange=A1:D100, firstrow, allstring, case=lower, verbose
 */
ST_retcode xlsx_import_main(const char *args);

/*
 * Initialize XLSX context with default values.
 */
void xlsx_context_init(XLSXContext *ctx);

/*
 * Free all resources in XLSX context.
 */
void xlsx_context_free(XLSXContext *ctx);

/*
 * Open and validate an XLSX file.
 * Returns 0 on success, Stata error code on failure.
 */
ST_retcode xlsx_open_file(XLSXContext *ctx, const char *filename);

/*
 * Parse workbook.xml to populate sheet list.
 */
ST_retcode xlsx_parse_workbook(XLSXContext *ctx);

/*
 * Parse sharedStrings.xml to build string table.
 */
ST_retcode xlsx_parse_shared_strings(XLSXContext *ctx);

/*
 * Parse styles.xml to detect date formats.
 */
ST_retcode xlsx_parse_styles(XLSXContext *ctx);

/*
 * Parse the selected worksheet.
 */
ST_retcode xlsx_parse_worksheet(XLSXContext *ctx);

/*
 * Infer column types from parsed data.
 */
ST_retcode xlsx_infer_types(XLSXContext *ctx);

/*
 * Build column caches for Stata output.
 */
ST_retcode xlsx_build_cache(XLSXContext *ctx);

/*
 * Select sheet by name. Returns 0-based index or -1 if not found.
 */
int xlsx_select_sheet_by_name(XLSXContext *ctx, const char *name);

/*
 * Parse a cell range string like "A1:D100" or "A1" or ":D100".
 */
bool xlsx_parse_cellrange(const char *range_str, XLSXCellRange *range);

#endif /* CIMPORT_XLSX_H */
