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
#include "../ctools_arena.h"
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
#define XLSX_MAX_ROWS           1048576        /* Excel max rows */
#define XLSX_MAX_PART_PATH      256            /* package part name */

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

/* Sheet part types, from the workbook relationship of each <sheet>. */
typedef enum {
    XLSX_SHEET_WORKSHEET = 0,
    XLSX_SHEET_CHART,
    XLSX_SHEET_DIALOG,
    XLSX_SHEET_MACRO,
    XLSX_SHEET_OTHER
} XLSXSheetKind;

typedef struct {
    char name[XLSX_MAX_SHEET_NAME];
    int sheet_id;
    char rel_id[32];  /* Relationship ID (rId1, rId2, etc.) */
    int sheet_index;  /* 1-based position, used only when relationships are absent */
    char path[XLSX_MAX_PART_PATH]; /* part resolved through workbook.xml.rels */
    int kind;         /* XLSXSheetKind */
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
    XLSX_CELL_DATE,
    XLSX_CELL_DATETIME
} XLSXCellType;

/* ============================================================================
 * XLSX Context Structure
 * ============================================================================ */

typedef struct XLSXContext {
    /* File info */
    char *filename;
    void *zip_archive;        /* xlsx_zip_archive handle */

    /* Package parts, resolved through the relationship parts */
    char workbook_path[XLSX_MAX_PART_PATH];
    char shared_strings_path[XLSX_MAX_PART_PATH];
    char styles_path[XLSX_MAX_PART_PATH];

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
    uint32_t *shared_string_lengths; /* Pre-computed lengths for O(1) lookup */

    /* Date format detection (style indices that are dates) */
    bool date1904;           /* workbookPr date system */
    uint8_t *date_styles;    /* 0 numeric, 1 daily date, 2 datetime */
    uint8_t *excel_formats;  /* native display-format mapping per style */
    int num_styles;

    /* Arena for parse-phase allocations (inline strings) */
    ctools_arena parse_arena;

    /* <dimension> of the worksheet (0 if absent): only a sizing hint */
    int dimension_rows;       /* Total rows from dimension ref */
    int dimension_cols;       /* Total columns from dimension ref */

    /* Used range of the worksheet */
    int max_col;              /* Maximum column index seen (0-based) */
    int min_row;              /* Minimum row number (1-based) */
    int max_row;              /* Maximum row number (1-based) */

    /* Column-major cell storage, grown to fit every cell with a value */
    double **cm_numeric;      /* cm_numeric[col][row] — values or (double)ss_index */
    uint8_t **cm_types;       /* cm_types[col][row] — XLSXCellType per cell */
    uint8_t **cm_formats;     /* cm_formats[col][row] — native format code */
    size_t cm_num_rows;       /* Number of data rows allocated */
    int cm_num_cols;          /* Number of columns allocated */
    bool cm_active;           /* true once storage has been allocated */

    /* Cells with a value inside the selected range, and cells stored */
    size_t cells_seen;
    size_t cells_stored;

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
 * Parse workbook.xml to populate the sheet list, and resolve each sheet,
 * sharedStrings and styles through workbook.xml.rels. Selects the first
 * worksheet by default.
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
 * Select sheet by name. Returns 0-based index or -1 if not found.
 */
int xlsx_select_sheet_by_name(XLSXContext *ctx, const char *name);

/*
 * Parse a cell range string like "A1:D100" or "A1" or ":D100".
 */
bool xlsx_parse_cellrange(const char *range_str, XLSXCellRange *range);

#endif /* CIMPORT_XLSX_H */
