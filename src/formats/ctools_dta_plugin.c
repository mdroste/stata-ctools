/*
 * ctools_dta_plugin.c - internal plugin command "dta" for the .dta writer
 *
 *   dta save          the plugin varlist (all observations) as the .dta in
 *                     local __ctools_dta_path, with storage types, names and
 *                     formats from locals __ctools_dta_types, _names and
 *                     _formats, value label names from __ctools_dta_vallabs
 *                     ("." for none), variable label j from
 *                     __ctools_dta_label<j>, the data label from
 *                     __ctools_dta_datalabel, sorted-by variable numbers from
 *                     __ctools_dta_sort and a thread count from
 *                     __ctools_dta_threads. The round trip behind the tests.
 *   dta selftest R T  a fixed dataset of every storage type, missing value,
 *                     boundary value, strL kind, label and characteristic,
 *                     in format R (118 or 119) with T threads
 *   dta wide K        K byte variables of 3 observations (format 119 past
 *                     32,767 variables)
 *   dta error CASE    a spec the writer must refuse; returns its code
 *   dta producer N T KEEP   N rows from a synthetic producer (see below)
 *   dta bench N K KIND T   N x K synthetic values of KIND (double, byte,
 *                     str16 or strl40) with T threads; local
 *                     __ctools_dta_seconds holds the time the write took
 */
#include <inttypes.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "stplugin.h"
#include "ctools_runtime.h"
#include "formats/ctools_dta.h"

static int dta_report(int rc, const ctools_dta_result *res)
{
    if (rc && res->message[0]) {
        char line[300];
        snprintf(line, sizeof(line), "dta: %s\n", res->message);
        SF_error(line);
    }
    return rc;
}

/* A local macro into a malloc'd buffer; "" when undefined. */
static char *dta_local(const char *name)
{
    char macro[160];
    snprintf(macro, sizeof(macro), "_%s", name);
    for (size_t cap = (size_t)1 << 16; cap <= ((size_t)1 << 26); cap *= 4) {
        char *buf = calloc(cap, 1);
        if (!buf) return NULL;
        if (SF_macro_use(macro, buf, (int)cap) != 0) buf[0] = 0;
        if (strlen(buf) < cap - 1) return buf;
        free(buf);
    }
    return NULL;
}

/* Split s in place on spaces into at most n tokens; returns the count. */
static size_t dta_tokens(char *s, char **tok, size_t n)
{
    size_t k = 0;
    for (char *p = s; *p;) {
        while (*p == ' ') p++;
        if (!*p) break;
        char *e = p;
        while (*e && *e != ' ') e++;
        if (k < n) tok[k] = p;
        k++;
        if (*e) *e++ = 0;
        p = e;
    }
    return k;
}

/* ---------------------------------------------------------------- save */

static int dta_save(void)
{
    const size_t K = (size_t)SF_nvars();
    const uint64_t N = (uint64_t)SF_nobs();
    ctools_dta_result res;
    memset(&res, 0, sizeof(res));
    int rc = 0;
    char *path = dta_local("__ctools_dta_path"), *types = dta_local("__ctools_dta_types"),
         *names = dta_local("__ctools_dta_names"), *formats = dta_local("__ctools_dta_formats"),
         *vallabs = dta_local("__ctools_dta_vallabs"), *datalabel = dta_local("__ctools_dta_datalabel"),
         *sort = dta_local("__ctools_dta_sort"), *threads = dta_local("__ctools_dta_threads");
    char **ttok = calloc(K + 1, sizeof(char *)), **ntok = calloc(K + 1, sizeof(char *)),
         **ftok = calloc(K + 1, sizeof(char *)), **vtok = calloc(K + 1, sizeof(char *)),
         **stok = calloc(K + 1, sizeof(char *)), **labels = calloc(K + 1, sizeof(char *));
    ctools_dta_column *cols = calloc(K + 1, sizeof(*cols));
    void **data = calloc(K + 1, sizeof(void *));
    char **arena = calloc(K + 1, sizeof(char *));
    int *sortlist = calloc(K + 1, sizeof(int));
    if (!path || !types || !names || !formats || !vallabs || !datalabel || !sort || !threads || !ttok ||
        !ntok || !ftok || !vtok || !stok || !labels || !cols || !data || !arena || !sortlist) {
        rc = 920;
        goto done;
    }
    if (!*path || dta_tokens(types, ttok, K) != K || dta_tokens(names, ntok, K) != K ||
        dta_tokens(formats, ftok, K) != K || (*vallabs && dta_tokens(vallabs, vtok, K) != K)) {
        SF_error("dta save: the type, name and format lists must name every variable\n");
        rc = 198;
        goto done;
    }
    size_t nsort = dta_tokens(sort, stok, K);
    for (size_t i = 0; i < nsort && i < K; i++) sortlist[i] = atoi(stok[i]);

    for (size_t j = 0; j < K && !rc; j++) {
        ctools_dta_column *c = &cols[j];
        const int v = (int)j + 1;
        char lname[64];
        snprintf(lname, sizeof(lname), "__ctools_dta_label%zu", j + 1);
        if (!(labels[j] = dta_local(lname))) { rc = 920; break; }
        c->name = ntok[j];
        c->format = ftok[j];
        c->label = labels[j][0] ? labels[j] : NULL;
        c->value_label = vtok[j] && strcmp(vtok[j], ".") ? vtok[j] : NULL;
        const char *t = ttok[j];
        if (!strcmp(t, "strL")) {
            c->type = CTOOLS_DTA_STRL;
            c->source = CTOOLS_DTA_SRC_SPAN;
            ctools_dta_span *sp = calloc(N ? (size_t)N : 1, sizeof(*sp));
            uint64_t total = 0;
            if (!sp) { rc = 920; break; }
            data[j] = sp;
            for (uint64_t i = 0; i < N; i++) {
                int len = SF_sdatalen(v, (ST_int)(i + 1));
                if (len < 0) { rc = 459; break; }
                sp[i].len = (size_t)len;
                total += (uint64_t)len + 1;
            }
            if (rc) break;
            if (!(arena[j] = malloc(total ? (size_t)total : 1))) { rc = 920; break; }
            char *at = arena[j];
            for (uint64_t i = 0; i < N; i++) {
                if (sp[i].len && SF_strldata(v, (ST_int)(i + 1), at, (int)sp[i].len + 1) != (int)sp[i].len) {
                    rc = 459;
                    break;
                }
                sp[i].ptr = at;
                sp[i].binary = SF_var_is_binary(v, (ST_int)(i + 1)) ? 1 : 0;
                at += sp[i].len + 1;
            }
        } else if (!strncmp(t, "str", 3)) {
            int w = atoi(t + 3);
            if (w < 1 || w > 2045) { rc = 198; break; }
            c->type = CTOOLS_DTA_STR;
            c->width = w;
            c->source = CTOOLS_DTA_SRC_FLAT;
            c->stride = (size_t)w + 1;
            char *flat = calloc(N ? (size_t)N : 1, c->stride);
            if (!flat) { rc = 920; break; }
            data[j] = flat;
            for (uint64_t i = 0; i < N; i++)
                if (SF_sdata(v, (ST_int)(i + 1), flat + i * c->stride)) { rc = 459; break; }
        } else {
            c->type = !strcmp(t, "byte")    ? CTOOLS_DTA_BYTE
                      : !strcmp(t, "int")   ? CTOOLS_DTA_INT
                      : !strcmp(t, "long")  ? CTOOLS_DTA_LONG
                      : !strcmp(t, "float") ? CTOOLS_DTA_FLOAT
                                            : CTOOLS_DTA_DOUBLE;
            c->source = CTOOLS_DTA_SRC_F64;
            double *x = malloc((N ? (size_t)N : 1) * sizeof(double));
            if (!x) { rc = 920; break; }
            data[j] = x;
            for (uint64_t i = 0; i < N; i++)
                if (SF_vdata(v, (ST_int)(i + 1), &x[i])) { rc = 459; break; }
        }
        c->data = data[j];
    }
    if (!rc) {
        ctools_dta_spec spec;
        ctools_dta_spec_init(&spec);
        spec.nobs = N;
        spec.ncols = K;
        spec.cols = cols;
        spec.data_label = *datalabel ? datalabel : NULL;
        spec.sortlist = sortlist;
        spec.nsort = nsort;
        spec.threads = atoi(threads);
        spec.temporary = 1;
        rc = dta_report(ctools_dta_write(path, &spec, &res), &res);
    }
done:
    for (size_t j = 0; j < K; j++) {
        if (data) free(data[j]);
        if (arena) free(arena[j]);
        if (labels) free(labels[j]);
    }
    free(path), free(types), free(names), free(formats), free(vallabs), free(datalabel), free(sort);
    free(threads), free(ttok), free(ntok), free(ftok), free(vtok), free(stok), free(labels);
    free(cols), free(data), free(arena), free(sortlist);
    return rc;
}

/* ------------------------------------------------------------- selftest */

#define DTA_ST_COLS 9

static int dta_selftest(int release, int threads, uint64_t N)
{
    static const double bvals[5] = {-127, 100, 0, -1, 1}, ivals[5] = {-32767, 32740, 0, -1, 1},
                        lvals[5] = {-2147483647.0, 2147483620.0, 0, -1, 1},
                        fvals[5] = {-1.7014117331926443e38, 1.7014117331926443e38, 0.5, -0.25, 3},
                        dvals[5] = {-8.9884656743115785e307, 8.9884656743115785e307, 0.1, -1e-300, 1};
    const double *vals[5] = {bvals, ivals, lvals, fvals, dvals};
    ctools_dta_result res;
    memset(&res, 0, sizeof(res));
    int rc = 0;
    char *path = dta_local("__ctools_dta_path");
    double *num[5] = {0};
    char *s1 = NULL, *s2045 = NULL, *big = NULL;
    ctools_dta_span *L = NULL, *L2 = NULL;
    char *l2text = NULL;
    double *id = malloc((size_t)N * sizeof(double));
    for (int k = 0; k < 5; k++) num[k] = malloc((size_t)N * sizeof(double));
    s1 = calloc((size_t)N, 2);
    s2045 = calloc((size_t)N, 2046);
    big = malloc(100000);
    L = calloc((size_t)N, sizeof(*L));
    L2 = calloc((size_t)N, sizeof(*L2));
    l2text = calloc((size_t)N, 24);
    if (!path || !*path || !id || !num[0] || !num[1] || !num[2] || !num[3] || !num[4] || !s1 || !s2045 ||
        !big || !L || !L2 || !l2text) {
        rc = path && !*path ? 198 : 920;
        goto done;
    }
    static char b2045[2045], c2046[2046];
    memset(b2045, 'b', sizeof(b2045));
    memset(c2046, 'c', sizeof(c2046));
    memset(big, 'd', 100000);
    for (uint64_t r = 0; r < N; r++) {
        id[r] = (double)(r + 1);
        unsigned m = (unsigned)(r % 32);
        for (int k = 0; k < 5; k++) {
            if (m < 27) {
                uint64_t bits = UINT64_C(0x7FE0000000000000) + ((uint64_t)m << 40);
                memcpy(&num[k][r], &bits, 8);
            } else {
                num[k][r] = vals[k][m - 27];
            }
        }
        s1[r * 2] = r % 3 == 0 ? 0 : r % 3 == 1 ? 'a' : 'z';
        char *s = s2045 + r * 2046;
        if (r % 2) memset(s, 'A' + (int)(r % 26), 2045);
        else snprintf(s, 2046, "short%" PRIu64, r);
        switch (r % 7) {
        case 0: L[r] = (ctools_dta_span){"", 0, 0}; break;
        case 1: L[r] = (ctools_dta_span){"x", 1, 0}; break;
        case 2: L[r] = (ctools_dta_span){b2045, 2045, 0}; break;
        case 3: L[r] = (ctools_dta_span){c2046, 2046, 0}; break;
        case 4: L[r] = (ctools_dta_span){big, 100000, 0}; break;
        case 5: L[r] = (ctools_dta_span){"a\0b", 3, 0}; break; /* binary: holds a NUL */
        default: L[r] = (ctools_dta_span){"h\xc3\xa9llo w\xc3\xb6rld", 13, 0}; break;
        }
        if (r % 5 == 0) {
            int n = snprintf(l2text + r * 24, 24, "five%" PRIu64, r);
            L2[r] = (ctools_dta_span){l2text + r * 24, (size_t)n, r % 10 == 0};
        }
    }
    ctools_dta_column cols[DTA_ST_COLS];
    memset(cols, 0, sizeof(cols));
    cols[0] = (ctools_dta_column){.name = "id", .type = CTOOLS_DTA_LONG, .source = CTOOLS_DTA_SRC_F64, .data = id};
    cols[1] = (ctools_dta_column){.name = "b", .type = CTOOLS_DTA_BYTE, .source = CTOOLS_DTA_SRC_F64,
                                  .data = num[0], .value_label = "yn"};
    cols[2] = (ctools_dta_column){.name = "i", .type = CTOOLS_DTA_INT, .source = CTOOLS_DTA_SRC_F64, .data = num[1]};
    cols[3] = (ctools_dta_column){.name = "l", .type = CTOOLS_DTA_LONG, .source = CTOOLS_DTA_SRC_F64, .data = num[2]};
    cols[4] = (ctools_dta_column){.name = "f", .type = CTOOLS_DTA_FLOAT, .source = CTOOLS_DTA_SRC_F64, .data = num[3]};
    cols[5] = (ctools_dta_column){.name = "d", .type = CTOOLS_DTA_DOUBLE, .source = CTOOLS_DTA_SRC_F64,
                                  .data = num[4], .format = "%12.4f", .label = "a double"};
    cols[6] = (ctools_dta_column){.name = "s1", .type = CTOOLS_DTA_STR, .width = 0, .source = CTOOLS_DTA_SRC_FLAT,
                                  .data = s1, .stride = 2};
    cols[7] = (ctools_dta_column){.name = "s2045", .type = CTOOLS_DTA_STR, .width = 2045,
                                  .source = CTOOLS_DTA_SRC_FLAT, .data = s2045, .stride = 2046};
    cols[8] = (ctools_dta_column){.name = "L", .type = CTOOLS_DTA_STRL, .source = CTOOLS_DTA_SRC_SPAN, .data = L,
                                  .label = "a strL"};
    /* L2 follows as a tenth column, built from the same layout */
    ctools_dta_column all[DTA_ST_COLS + 1];
    memcpy(all, cols, sizeof(cols));
    all[DTA_ST_COLS] = (ctools_dta_column){.name = "L2", .type = CTOOLS_DTA_STRL, .source = CTOOLS_DTA_SRC_SPAN,
                                           .data = L2};
    static const int32_t yv[3] = {1, -1, 0};
    static const char *const yt[3] = {"one", "minus one", "zero"};
    const ctools_dta_value_label lab = {"yn", 3, yv, yt};
    const ctools_dta_char chars[5] = {{"_dta", "note0", "1"},
                                      {"_dta", "note1", "dataset note"},
                                      {"d", "note0", "1"},
                                      {"d", "note1", "double note"},
                                      {"s1", "mychar", "char text"}};
    static const int sortlist[1] = {1};
    ctools_dta_spec spec;
    ctools_dta_spec_init(&spec);
    spec.nobs = N;
    spec.ncols = DTA_ST_COLS + 1;
    spec.cols = all;
    spec.release = release;
    spec.data_label = "ctools dta selftest";
    spec.sortlist = sortlist;
    spec.nsort = 1;
    spec.chars = chars;
    spec.nchars = 5;
    spec.labels = &lab;
    spec.nlabels = 1;
    spec.threads = threads;
    spec.temporary = 1;
    rc = dta_report(ctools_dta_write(path, &spec, &res), &res);
    if (!rc) {
        char text[64];
        snprintf(text, sizeof(text), "%" PRIu64, res.strls);
        rc = SF_macro_save("___ctools_dta_strls", text);
    }
done:
    free(path), free(id), free(s1), free(s2045), free(big), free(L), free(L2), free(l2text);
    for (int k = 0; k < 5; k++) free(num[k]);
    return rc;
}

static int dta_wide(size_t K)
{
    ctools_dta_result res;
    memset(&res, 0, sizeof(res));
    char *path = dta_local("__ctools_dta_path");
    ctools_dta_column *cols = calloc(K ? K : 1, sizeof(*cols));
    char *names = calloc(K ? K : 1, 12);
    double *x = malloc((K ? K : 1) * 3 * sizeof(double));
    int rc = 0;
    if (!path || !cols || !names || !x) rc = 920;
    for (size_t j = 0; j < K && !rc; j++) {
        snprintf(names + 12 * j, 12, "v%zu", j + 1);
        for (int r = 0; r < 3; r++) x[3 * j + r] = (double)((j + (size_t)r) % 100);
        cols[j] = (ctools_dta_column){.name = names + 12 * j, .type = CTOOLS_DTA_BYTE,
                                      .source = CTOOLS_DTA_SRC_F64, .data = x + 3 * j};
    }
    if (!rc) {
        ctools_dta_spec spec;
        ctools_dta_spec_init(&spec);
        spec.nobs = 3;
        spec.ncols = K;
        spec.cols = cols;
        spec.temporary = 1;
        rc = dta_report(ctools_dta_write(path, &spec, &res), &res);
    }
    free(path), free(cols), free(names), free(x);
    return rc;
}

/* Specs the writer must refuse, without leaving a file behind. */
static int dta_error_case(const char *which)
{
    ctools_dta_result res;
    memset(&res, 0, sizeof(res));
    char *path = dta_local("__ctools_dta_path");
    if (!path) return 920;
    double x[2] = {1, 200};
    double frac[2] = {1, 2.5};
    const char *strs[2] = {"ab", "abcd"};
    ctools_dta_column c[2];
    memset(c, 0, sizeof(c));
    c[0] = (ctools_dta_column){.name = "x", .type = CTOOLS_DTA_DOUBLE, .source = CTOOLS_DTA_SRC_F64, .data = x};
    c[1] = (ctools_dta_column){.name = "y", .type = CTOOLS_DTA_DOUBLE, .source = CTOOLS_DTA_SRC_F64, .data = x};
    size_t ncols = 2;
    if (!strcmp(which, "dupname")) c[1].name = "x";
    else if (!strcmp(which, "badname")) c[1].name = "1y";
    else if (!strcmp(which, "reserved")) c[1].name = "_n";
    else if (!strcmp(which, "byte")) c[1].type = CTOOLS_DTA_BYTE;
    else if (!strcmp(which, "fraction")) c[1].type = CTOOLS_DTA_INT, c[1].data = frac;
    else if (!strcmp(which, "strwidth"))
        c[1] = (ctools_dta_column){.name = "y", .type = CTOOLS_DTA_STR, .width = 3, .source = CTOOLS_DTA_SRC_CSTR,
                                   .data = strs};
    else if (!strcmp(which, "strformat")) c[0].format = "%9s";
    else {
        free(path);
        return 198;
    }
    ctools_dta_spec spec;
    ctools_dta_spec_init(&spec);
    spec.nobs = 2;
    spec.ncols = ncols;
    spec.cols = c;
    int rc = ctools_dta_write(path, &spec, &res);
    char text[300];
    snprintf(text, sizeof(text), "%d", rc);
    SF_macro_save("___ctools_dta_rc", text);
    SF_macro_save("___ctools_dta_message", res.message);
    free(path);
    return 0;
}

/* ------------------------------------------------------------- producer */

/* A synthetic producer of groups of 1,000 rows next to a column held in
 * memory: pid (long), pnum (double, .b every 13th row), pstr (str8), pL
 * (strL: "", 3,000 bytes, short text, binary), m (double, in memory). keep
 * is the writer's keep budget (1 keeps nothing, so every group is filled
 * again for its strLs). Locals __ctools_dta_fills and _releases count calls. */
typedef struct
{
    uint64_t nobs;
    pthread_mutex_t lock;
    uint64_t fills, releases;
} dta_synth;

typedef struct
{
    double *pid, *pnum;
    ctools_dta_span *pstr, *pL;
    char *text;
} dta_synth_group;

#define DTA_SYNTH_ROWS 1000

static void dta_synth_release(void *arg, ctools_dta_group *g)
{
    dta_synth *s = arg;
    dta_synth_group *st = g->state;
    if (st) free(st->pid), free(st->pnum), free(st->pstr), free(st->pL), free(st->text), free(st);
    g->state = NULL;
    pthread_mutex_lock(&s->lock);
    s->releases++;
    pthread_mutex_unlock(&s->lock);
}

static int dta_synth_fill(void *arg, size_t group, int phase, ctools_dta_group *out)
{
    (void)phase;
    dta_synth *s = arg;
    uint64_t r0 = (uint64_t)group * DTA_SYNTH_ROWS, r1 = r0 + DTA_SYNTH_ROWS < s->nobs ? r0 + DTA_SYNTH_ROWS : s->nobs;
    size_t n = (size_t)(r1 - r0);
    dta_synth_group *st = calloc(1, sizeof(*st));
    if (!st) return 920;
    st->pid = malloc(n * sizeof(double));
    st->pnum = malloc(n * sizeof(double));
    st->pstr = calloc(n, sizeof(ctools_dta_span));
    st->pL = calloc(n, sizeof(ctools_dta_span));
    st->text = malloc(n * (3000 + 32));
    if (!st->pid || !st->pnum || !st->pstr || !st->pL || !st->text) {
        free(st->pid), free(st->pnum), free(st->pstr), free(st->pL), free(st->text), free(st);
        return 920;
    }
    for (size_t i = 0; i < n; i++) {
        uint64_t r = r0 + i;
        char *t = st->text + i * (3000 + 32);
        st->pid[i] = (double)(r + 1);
        if (r % 13 == 0) {
            uint64_t bits = UINT64_C(0x7FE0000000000000) + (UINT64_C(2) << 40);
            memcpy(&st->pnum[i], &bits, 8);
        } else {
            st->pnum[i] = (double)r * 0.25;
        }
        int len = snprintf(t, 32, "p%" PRIu64, r % 1000);
        st->pstr[i] = (ctools_dta_span){t, (size_t)len, 0};
        char *u = t + 32;
        switch (r % 4) {
        case 0: break;
        case 1:
            memset(u, 'a' + (int)(r % 26), 3000);
            st->pL[i] = (ctools_dta_span){u, 3000, 0};
            break;
        case 2:
            len = snprintf(u, 32, "short%" PRIu64, r);
            st->pL[i] = (ctools_dta_span){u, (size_t)len, 0};
            break;
        default:
            u[0] = 0;
            len = snprintf(u + 1, 31, "%" PRIu64, r);
            st->pL[i] = (ctools_dta_span){u, (size_t)len + 1, 0};
            break;
        }
    }
    out->columns[0] = st->pid;
    out->columns[1] = st->pnum;
    out->columns[2] = st->pstr;
    out->columns[3] = st->pL;
    out->state = st;
    pthread_mutex_lock(&s->lock);
    s->fills++;
    pthread_mutex_unlock(&s->lock);
    return 0;
}

static int dta_producer_test(uint64_t N, int threads, size_t keep)
{
    ctools_dta_result res;
    memset(&res, 0, sizeof(res));
    char *path = dta_local("__ctools_dta_path");
    size_t ngroups = (size_t)((N + DTA_SYNTH_ROWS - 1) / DTA_SYNTH_ROWS);
    uint64_t *starts = malloc((ngroups + 1) * sizeof(uint64_t)), *bytes = malloc((ngroups + 1) * sizeof(uint64_t));
    double *m = malloc((N ? (size_t)N : 1) * sizeof(double));
    int rc = path && starts && bytes && m ? 0 : 920;
    dta_synth s = {N, PTHREAD_MUTEX_INITIALIZER, 0, 0};
    if (!rc) {
        for (size_t g = 0; g <= ngroups; g++) starts[g] = (uint64_t)g * DTA_SYNTH_ROWS < N ? (uint64_t)g * DTA_SYNTH_ROWS : N;
        for (size_t g = 0; g < ngroups; g++) bytes[g] = (starts[g + 1] - starts[g]) * (3000 + 32 + 64);
        for (uint64_t r = 0; r < N; r++) m[r] = (double)r * 2;
        ctools_dta_column cols[5] = {
            {.name = "pid", .type = CTOOLS_DTA_LONG, .source = CTOOLS_DTA_SRC_PRODUCER},
            {.name = "pnum", .type = CTOOLS_DTA_DOUBLE, .source = CTOOLS_DTA_SRC_PRODUCER},
            {.name = "pstr", .type = CTOOLS_DTA_STR, .width = 8, .source = CTOOLS_DTA_SRC_PRODUCER},
            {.name = "pL", .type = CTOOLS_DTA_STRL, .source = CTOOLS_DTA_SRC_PRODUCER},
            {.name = "m", .type = CTOOLS_DTA_DOUBLE, .source = CTOOLS_DTA_SRC_F64, .data = m}};
        ctools_dta_producer pr = {ngroups, starts, bytes, &s, dta_synth_fill, dta_synth_release};
        ctools_dta_spec spec;
        ctools_dta_spec_init(&spec);
        spec.nobs = N;
        spec.ncols = 5;
        spec.cols = cols;
        spec.producer = &pr;
        spec.threads = threads;
        spec.keep_bytes = keep;
        spec.temporary = 1;
        rc = dta_report(ctools_dta_write(path, &spec, &res), &res);
        char text[64];
        snprintf(text, sizeof(text), "%" PRIu64, s.fills);
        if (!rc) rc = SF_macro_save("___ctools_dta_fills", text);
        snprintf(text, sizeof(text), "%" PRIu64, s.releases);
        if (!rc) rc = SF_macro_save("___ctools_dta_releases", text);
    }
    free(path), free(starts), free(bytes), free(m);
    return rc;
}

/* ---------------------------------------------------------------- bench */

static int dta_bench(uint64_t N, size_t K, const char *kind, int threads)
{
    ctools_dta_result res;
    memset(&res, 0, sizeof(res));
    char *path = dta_local("__ctools_dta_path");
    ctools_dta_column *cols = calloc(K ? K : 1, sizeof(*cols));
    char *names = calloc(K ? K : 1, 12);
    void **data = calloc(K ? K : 1, sizeof(void *));
    char *arena = NULL;
    int rc = path && cols && names && data ? 0 : 920;
    int strl = !strcmp(kind, "strl40"), str16 = !strcmp(kind, "str16"), byte = !strcmp(kind, "byte");
    if (!rc && strl && !(arena = malloc((size_t)N * 40 + 1))) rc = 920;
    for (uint64_t r = 0; !rc && strl && r < N; r++)
        snprintf(arena + r * 40, 41, "%020" PRIu64 "%020" PRIu64, r, r * 7);
    for (size_t j = 0; j < K && !rc; j++) {
        snprintf(names + 12 * j, 12, "v%zu", j + 1);
        cols[j].name = names + 12 * j;
        if (strl) {
            ctools_dta_span *sp = malloc((size_t)N * sizeof(*sp));
            if (!(data[j] = sp)) { rc = 920; break; }
            for (uint64_t r = 0; r < N; r++) sp[r] = (ctools_dta_span){arena + r * 40, 40, 0};
            cols[j] = (ctools_dta_column){.name = cols[j].name, .type = CTOOLS_DTA_STRL,
                                          .source = CTOOLS_DTA_SRC_SPAN, .data = sp, .flags = CTOOLS_DTA_COL_TEXT};
        } else if (str16) {
            char *flat = malloc((size_t)N * 17);
            if (!(data[j] = flat)) { rc = 920; break; }
            for (uint64_t r = 0; r < N; r++) snprintf(flat + r * 17, 17, "%016" PRIu64, r * (j + 1));
            cols[j] = (ctools_dta_column){.name = cols[j].name, .type = CTOOLS_DTA_STR, .width = 16,
                                          .source = CTOOLS_DTA_SRC_FLAT, .data = flat, .stride = 17};
        } else {
            double *x = malloc((size_t)N * sizeof(double));
            if (!(data[j] = x)) { rc = 920; break; }
            for (uint64_t r = 0; r < N; r++) x[r] = byte ? (double)((r + j) % 100) : (double)r * 0.5 + (double)j;
            cols[j] = (ctools_dta_column){.name = cols[j].name, .type = byte ? CTOOLS_DTA_BYTE : CTOOLS_DTA_DOUBLE,
                                          .source = CTOOLS_DTA_SRC_F64, .data = x};
        }
    }
    if (!rc) {
        ctools_dta_spec spec;
        ctools_dta_spec_init(&spec);
        spec.nobs = N;
        spec.ncols = K;
        spec.cols = cols;
        spec.threads = threads;
        spec.temporary = 1;
        double t0 = ctools_timer_seconds();
        rc = dta_report(ctools_dta_write(path, &spec, &res), &res);
        double t1 = ctools_timer_seconds();
        char text[64];
        snprintf(text, sizeof(text), "%.6f", t1 - t0);
        if (!rc) rc = SF_macro_save("___ctools_dta_seconds", text);
        snprintf(text, sizeof(text), "%.6f", res.write_seconds);
        if (!rc) rc = SF_macro_save("___ctools_dta_write_seconds", text);
    }
    for (size_t j = 0; data && j < K; j++) free(data[j]);
    free(path), free(cols), free(names), free(data), free(arena);
    return rc;
}

int ctools_dta_main(const char *args)
{
    char op[32] = "", a[64] = "", b[64] = "", c[64] = "", d[64] = "";
    if (!args || sscanf(args, "%31s %63s %63s %63s %63s", op, a, b, c, d) < 1) return 198;
    if (!strcmp(op, "save")) return dta_save();
    if (!strcmp(op, "selftest")) {
        int release = *a ? atoi(a) : 0, threads = *b ? atoi(b) : 0;
        uint64_t n = *c ? strtoull(c, NULL, 10) : 3000;
        return n ? dta_selftest(release, threads, n) : 198;
    }
    if (!strcmp(op, "wide")) return dta_wide(*a ? (size_t)strtoull(a, NULL, 10) : 0);
    if (!strcmp(op, "error")) return dta_error_case(a);
    if (!strcmp(op, "producer"))
        return dta_producer_test(strtoull(a, NULL, 10), *b ? atoi(b) : 0, *c ? (size_t)strtoull(c, NULL, 10) : 0);
    if (!strcmp(op, "bench"))
        return dta_bench(strtoull(a, NULL, 10), (size_t)strtoull(b, NULL, 10), c, *d ? atoi(d) : 0);
    return 198;
}
