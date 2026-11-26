#ifndef CRESHAPE_INTERNAL_H
#define CRESHAPE_INTERNAL_H
/*
 * creshape internals, shared by the files of src/creshape.
 *
 *   creshape_names.c   names and j values as text: no Stata calls
 *   creshape_plan.c    roles of the variables and the list of output variables
 *   creshape_keys.c    order, groups and checks on the staged observations
 *   creshape_schema.c  description of the result (labels, formats, characteristics)
 *   creshape_store.c   parallel writes of the result, and the restore
 *   creshape_impl.c    plugin calls, state lifetime, messages
 *
 * Index spaces: columns are zero-based positions in the dataset as the plan
 * call saw it (the plugin variable is column + 1, and the staged column has
 * the same number, because every variable is staged in order). Staged rows
 * are zero-based; Stata observations are one-based. Levels are positions in
 * levels[], the j values in the order of the result's variables (wide) or
 * observations (long).
 */
#include <stddef.h>
#include <stdint.h>
#include "ctools_types.h"
#include "ctools_rowstage.h"
#include "formats/ctools_dta_schema.h"

#define CRESHAPE_NAME_BYTES CTOOLS_DTA_NAME_BYTES
#define CRESHAPE_MAX_OBS ((size_t)2147483647)

typedef enum
{
    CRESHAPE_LONG = 1, /* wide data to long */
    CRESHAPE_WIDE = 2  /* long data to wide */
} creshape_direction;

/* One value of j. */
typedef struct
{
    double num; /* numeric j */
    char *text; /* owned: the string value, or the number as it appears in names */
    int found;  /* wide: an observation has this value; long: a variable carries it */
} creshape_level;

/* One xij entry: a stub such as inc, or inc@r. */
typedef struct
{
    char *spec;      /* owned: the entry as typed */
    char *prefix;    /* owned: text before @ (the whole stub when it has no @) */
    char *suffix;    /* owned: text after @ */
    char *long_name; /* owned: prefix, atwl() and suffix */
    int long_col;    /* wide: column of the long variable */
    int is_string;
    int *cols;    /* owned [nlevels]: long: column of the wide variable, -1 when absent */
    char **names; /* owned [nlevels]: names of the wide variables */
} creshape_stub;

typedef enum
{
    CRESHAPE_OUT_COPY = 0, /* an i or xi variable: column col of the source row */
    CRESHAPE_OUT_J,        /* long: the new j variable */
    CRESHAPE_OUT_STUB      /* long: stub col over the levels; wide: stub col at level */
} creshape_out_kind;

/* One variable of the result, in the result's order. */
typedef struct
{
    creshape_out_kind kind;
    int is_string;
    int col;   /* COPY: staged column; STUB: position in stubs[] */
    int level; /* wide STUB: position in levels[] */
    char name[CRESHAPE_NAME_BYTES];
} creshape_out;

typedef struct
{
    creshape_direction dir;
    int stage;    /* 0 nothing, 1 planned, 2 description written, 3 stored */
    int string_j; /* j holds strings */
    int verbose;
    int report;      /* print reshape's summary table after the store */
    int user_levels; /* the j values were listed in j() */
    int xi_listed;   /* the xi variables were listed (reshape xi); others are dropped */
    int maxvar;      /* c(maxvar) */

    /* the dataset as planned */
    size_t nvars;         /* variables = plugin variables of the plan call */
    size_t nobs;          /* observations */
    char *names_buf;      /* owned: the names, NUL-separated */
    char **names;         /* owned [nvars]: into names_buf */
    uint32_t *name_slots; /* owned: open-addressed name -> column + 1 */
    size_t name_cap;
    unsigned char *is_str;  /* owned [nvars]: 1 for a string variable */
    ctools_rowstage staged; /* every variable, in order; rows == NULL when nobs == 0 */

    /* roles */
    size_t ni;
    int *icols; /* owned [ni] */
    int jcol;   /* wide: column of j; long: -1 */
    char jname[CRESHAPE_NAME_BYTES];
    char *atwl; /* owned */
    size_t nstubs;
    creshape_stub *stubs; /* owned [nstubs] */
    size_t nxi;
    int *xicols; /* owned [nxi]: in dataset order */
    size_t nlevels;
    creshape_level *levels; /* owned [nlevels] */
    size_t nreport;         /* levels the summary and ReS_Xij_wide# list */
    size_t *report_levels;  /* owned [nreport]: positions in levels[], by value */

    /* order of the staged rows */
    perm_idx_t *order; /* owned (aligned) [nobs]: position -> row; NULL: in order */
    /* wide */
    unsigned char *first; /* owned [nobs]: 1 where a position starts a group of i */
    int scan_unsorted;    /* the rows are not in order of (i, j) */
    int scan_equal;       /* two neighbouring rows have the same (i, j) */
    size_t ngroups;
    perm_idx_t *group_start; /* owned [ngroups + 1]: positions in order */
    perm_idx_t *cell;        /* owned [ngroups * nlevels]: row, or CTOOLS_ROW_NONE */

    /* the result */
    size_t out_nobs, out_nvars;
    creshape_out *out; /* owned [out_nvars] */
    uint16_t j_type;   /* long: storage type of j as a .dta code (a width for str#) */
    size_t old_nobs, old_nvars;

    /* seconds: staging; order and checks (of which, for wide: the j values,
     * the sort, the groups and cells); description; store */
    double t_load, t_order, t_levels, t_sort, t_groups, t_schema, t_store;
} creshape_state;

/* ---- creshape_names.c ---------------------------------------------------- */

/* real() of text that may appear in a variable name: digits, with an
 * optional e, E, d or D exponent. 0 with *value set, or -1 (missing). */
int creshape_parse_real(const char *text, double *value);
/* A number as Stata puts it in a macro (integers below 1e16 in full). */
void creshape_format_number(double value, char *buf, size_t cap);
/* 1 when text can name a new variable. */
int creshape_name_ok(const char *text);
/* prefix + middle + suffix into a new string; NULL when out of memory. */
char *creshape_join(const char *prefix, const char *middle, const char *suffix);
char *creshape_strdup(const char *text);
/* Split a space-separated list in place; returns the number of words and
 * sets *words (owned by the caller, pointing into list), or -1. */
long creshape_split(char *list, char ***words);

/* Name table of the planned dataset. */
int creshape_names_index(creshape_state *st);
/* Column of name, or -1. */
int creshape_names_find(const creshape_state *st, const char *name);

/* ---- creshape_plan.c ----------------------------------------------------- */

/* A local macro of the calling ado program as a heap string ("" when it is
 * not set); NULL when out of memory. */
char *creshape_local(const char *name);
/* reshape's picture of the two forms, as SMCL for SF_error. */
extern const char creshape_picture[];
/* Messages in reshape's style: an error line, and a note in parentheses. */
void creshape_error(const char *fmt, ...);
void creshape_note(const char *fmt, ...);
/* Error text of several lines (SMCL), displayed a line at a time. */
void creshape_error_text(const char *text);
/* Read the request from the ado's locals and name the roles; 0 or a Stata
 * return code, with the message already displayed. */
int creshape_plan_read(creshape_state *st, const char *args);
/* Long: the j values from the variable names (when j() lists none), each
 * stub's variables and the xi variables. */
int creshape_plan_long(creshape_state *st);
/* Wide: the long variable of each stub and the xi variables. */
int creshape_plan_wide(creshape_state *st);
/* The result's variables in order (wide: after the j values are known), with
 * checks for invalid and repeated names and for c(maxvar). */
int creshape_plan_outputs(creshape_state *st);
/* Rebuild report_levels: the positions of levels[] in order of value. */
int creshape_levels_report(creshape_state *st);
/* Order of two j values: numbers by value, strings by their bytes. */
int creshape_level_cmp(const creshape_state *st, const creshape_level *a, const creshape_level *b);

/* ---- creshape_keys.c ----------------------------------------------------- */

/* Order the staged rows by the i variables (wide: then by j); *unique is 0
 * when two rows have the same key. */
int creshape_order_rows(creshape_state *st, int *unique);
/* Wide: the j values in the data, sorted. Fills levels unless j() listed
 * them, and marks which listed values occur. Displays its own messages. */
int creshape_wide_levels(creshape_state *st);
/* Wide: groups of equal i, the row of each (group, level), and the xi
 * variables that vary within a group (*varying owned, [nxi], or NULL). */
int creshape_wide_groups(creshape_state *st, unsigned char **varying);

/* ---- creshape_schema.c --------------------------------------------------- */

/* Read the template the ado saved and write the result's description. */
int creshape_schema_write(creshape_state *st, const char *template_path, const char *out_path,
                          char *msg, size_t cap);

/* ---- creshape_store.c ---------------------------------------------------- */

/* Fill the result: every variable, every observation. */
stata_retcode creshape_store(const creshape_state *st);

#endif /* CRESHAPE_INTERNAL_H */
