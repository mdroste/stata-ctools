#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "stplugin.h"

static ST_plugin mock;
ST_plugin *_stata_ = &mock;
static size_t stream_threshold;
static uint64_t cache_budget;
#define CIMPORT_STREAM_THRESHOLD stream_threshold
#define CIMPORT_CACHE_BUDGET cache_budget
static char first[32], last[32], strip[16], locale[32];
static ST_retcode macro(char *key, char *out, ST_int cap)
{
    const char *value = "";
    if (!strcmp(key, "___cimport_startrow"))
        value = first;
    if (!strcmp(key, "___cimport_endrow"))
        value = last;
    if (!strcmp(key, "___cimport_stripquotes"))
        value = strip;
    if (!strcmp(key, "___cimport_parselocale"))
        value = locale;
    snprintf(out, (size_t)cap + 1, "%s", value);
    return 0;
}
static ST_retcode display(char *s)
{
    (void)s;
    return 0;
}
static ST_int poll(void)
{
    return 0;
}
double ctools_timer_ms(void)
{
    return 0;
}
#include "cimport/cimport_impl.c"

typedef struct
{
    double number;
    char *text;
} Cell;
static Cell *cells;
static size_t nobs;
static int ncols, fail_col;
static size_t released_sources, retained_sources;
static ST_retcode number_store(ST_int col, ST_int obs, double value)
{
    assert(col > 0 && col <= ncols && obs > 0 && (size_t)obs <= nobs);
    if (col == fail_col && obs == 2)
        return 459;
    cells[((size_t)obs - 1) * ncols + col - 1].number = value;
    return 0;
}
static ST_retcode string_store(ST_int col, ST_int obs, char *value)
{
    assert(col > 0 && col <= ncols && obs > 0 && (size_t)obs <= nobs);
    if (col == fail_col && obs == 2)
        return 459;
    Cell *cell = &cells[((size_t)obs - 1) * ncols + col - 1];
    assert(!cell->text);
    cell->text = strdup(value);
    assert(cell->text);
    return 0;
}

static Cell *load(CImportContext *ctx, int failure)
{
    nobs = ctx->total_rows;
    ncols = ctx->num_columns;
    cells = calloc((nobs ? nobs : 1) * ncols, sizeof(*cells));
    assert(cells);
    unsigned char *kind = malloc((size_t)ncols);
    assert(kind);
    for (int col = 0; col < ncols; col++)
        kind[col] = ctx->columns[col].type == CIMPORT_COL_STRING
                        ? (ctx->columns[col].use_strl ? CIMPORT_STORE_STRL : CIMPORT_STORE_STRING)
                        : CIMPORT_STORE_NUMERIC;
    fail_col = failure;
    size_t obs = 1;
    int errors = 0;
    for (int c = 0; c < ctx->num_chunks; c++)
    {
        size_t start = c == 0 && ctx->has_header ? 1 : 0;
        CImportStoreTile task = {ctx, kind, c, start, ctx->chunks[c].num_rows, (ST_int)obs, 0};
        if (ctx->compact_rows)
            cimport_store_compact_worker(&task);
        else if (ctx->stream_rows)
            cimport_store_stream_worker(&task);
        else
            cimport_store_tile_worker(&task);
        if (task.error)
        {
            assert(task.error == 459 && failure);
            errors++;
        }
        obs += ctx->chunks[c].num_rows - start;
    }
    assert(errors == (failure ? 1 : 0));
    free(kind);
    return cells;
}

static void free_cells(Cell *data, size_t count)
{
    for (size_t i = 0; i < count; i++)
        free(data[i].text);
    free(data);
}

static void compare(const char *path, int mode, int empty, bool header, int skip,
                    bool expected_stream, int threads)
{
    ctools_set_max_threads(threads);
    stream_threshold = SIZE_MAX;
    CImportContext *a =
        cimport_parse_csv(path, ',', header, skip, false, mode, CIMPORT_NUMTYPE_AUTO, '.', 0, empty,
                          0, CIMPORT_ENC_UNKNOWN, NULL, 0, NULL, 0);
    assert(a && !a->stream_rows);
    stream_threshold = 1;
    CImportContext *b =
        cimport_parse_csv(path, ',', header, skip, false, mode, CIMPORT_NUMTYPE_AUTO, '.', 0, empty,
                          0, CIMPORT_ENC_UNKNOWN, NULL, 0, NULL, 0);
    if (!b || b->stream_rows != expected_stream)
        fprintf(
            stderr,
            "stream mismatch mode=%d empty=%d header=%d skip=%d expected=%d actual=%d columns=%d\n",
            mode, empty, header, skip, expected_stream, b ? b->stream_rows : -1, a->num_columns);
    assert(b && b->stream_rows == expected_stream);
    assert(a->num_columns == b->num_columns && a->total_rows == b->total_rows);
    assert(a->num_unmatched_quote_warnings == b->num_unmatched_quote_warnings);
    assert(!memcmp(a->columns, b->columns, (size_t)a->num_columns * sizeof(*a->columns)));
    if (b->stream_rows)
    {
        cache_budget = UINT64_MAX;
        assert(!cimport_should_stream(b));
        cache_budget = 0;
        assert(cimport_should_stream(b));
        size_t allocated = 0;
        for (int c = 0; c < b->num_chunks; c++)
        {
            assert(b->chunks[c].capacity == 1);
            allocated += b->chunks[c].arena.total_allocated;
        }
        assert(allocated == (b->has_header ? cimport_row_size(b->num_columns) : 0));
    }
    Cell *left = load(a, 0), *right = load(b, 0);
    size_t count = a->total_rows * a->num_columns;
    for (size_t i = 0; i < count; i++)
    {
        assert(!memcmp(&left[i].number, &right[i].number, sizeof(double)));
        assert((left[i].text == NULL) == (right[i].text == NULL));
        if (left[i].text)
            assert(!strcmp(left[i].text, right[i].text));
    }
    free_cells(right, count);
    /* Additional compact/hybrid comparisons preserve the original retained vs
     * streaming assertions above. Compare exact final Stata representations:
     * the mock stores doubles, whereas real float columns round at the SPI. */
    const uint64_t budgets[] = {UINT64_MAX, 128 * 1024, 1};
    for (size_t trial = 0; trial < sizeof(budgets) / sizeof(*budgets); trial++)
    {
        cache_budget = budgets[trial];
        CImportContext *compact =
            cimport_parse_csv(path, ',', header, skip, false, mode, CIMPORT_NUMTYPE_AUTO, '.', 0,
                              empty, 0, CIMPORT_ENC_UNKNOWN, NULL, 0, NULL, 0);
        assert(compact && compact->compact_rows == expected_stream);
        assert(compact->total_rows == a->total_rows && compact->num_columns == a->num_columns);
        assert(!memcmp(a->columns, compact->columns, (size_t)a->num_columns * sizeof(*a->columns)));
        assert(atomic_load(&compact->cache_budget.used) <= cache_budget);
        cimport_release_cached_source(compact);
        if (compact->compact_rows)
        {
            if (compact->file_data)
                retained_sources++;
            else
                released_sources++;
        }
        Cell *actual = load(compact, 0);
        for (size_t i = 0; i < count; i++)
        {
            double want = left[i].number, got = actual[i].number;
            if (a->columns[i % a->num_columns].num_subtype == CIMPORT_NUM_FLOAT &&
                a->columns[i % a->num_columns].type == CIMPORT_COL_NUMERIC)
            {
                if (want < SV_missval)
                    want = (double)(float)want;
                if (got < SV_missval)
                    got = (double)(float)got;
            }
            assert(!memcmp(&want, &got, sizeof(double)));
            assert((left[i].text == NULL) == (actual[i].text == NULL));
            if (left[i].text)
                assert(!strcmp(left[i].text, actual[i].text));
        }
        free_cells(actual, count);
        if (compact->compact_rows && compact->total_rows >= 2)
        {
            for (int col = 1; col <= compact->num_columns; col++)
                if (!compact->columns[col - 1].use_strl)
                    free_cells(load(compact, col), count);
        }
        cimport_free_context(compact);
    }
    cache_budget = 0;
    free_cells(left, count);
    if (b->stream_rows && b->total_rows >= 2)
    {
        for (int col = 1; col <= b->num_columns; col++)
        {
            if (!b->columns[col - 1].use_strl)
                free_cells(load(b, col), count);
        }
    }
    cimport_free_context(a);
    cimport_free_context(b);
}

static void cache_boundaries(void)
{
    const double values[][4] = {{7, 7, 7, 7},          {-100, -99, 154, 0},
                                {-32768, 32766, 0, 0}, {0, 65535, 4294967294.0, 0},
                                {-0.0, .125, .25, 0},  {.12345678901234567, 1e100, -1e100, 0}};
    const unsigned char kinds[] = {CIMPORT_CACHE_CONSTANT, CIMPORT_CACHE_U8,
                                   CIMPORT_CACHE_U16,      CIMPORT_CACHE_U32,
                                   CIMPORT_CACHE_FLOAT,    CIMPORT_CACHE_DOUBLE};
    for (size_t trial = 0; trial < sizeof(kinds); trial++)
    {
        CImportCacheBuilder b;
        assert(!cimport_cache_init(&b, 1));
        CImportCacheBudget budget = {.limit = 1024 * 1024};
        atomic_init(&budget.used, 0);
        double expected[4];
        for (int r = 0; r < 4; r++)
        {
            expected[r] = r == 3 && trial != 0 ? SV_missval : values[trial][r];
            CImportFieldRef field = {UINT64_C(1) << 40, 4};
            cimport_cache_append(&b, &field, 1, &expected[r], field.offset, field.offset + 4);
        }
        b.numeric[0] = 1;
        CImportCachedBlock *block =
            cimport_cache_pack(&b, &budget, SV_missval, false, false, NULL, 0);
        assert(block && block->columns[0].kind == kinds[trial]);
        assert(atomic_load(&budget.used) == block->bytes);
        for (int r = 0; r < 4; r++)
        {
            double actual = cimport_cache_number(block, 0, r, SV_missval);
            assert(!memcmp(&actual, &expected[r], sizeof(double)));
        }
        cimport_cache_free(block, &budget);
        assert(!atomic_load(&budget.used));
        cimport_cache_destroy(&b);
    }
    for (int fixed = 0; fixed < 2; fixed++)
    {
        CImportCacheBuilder b;
        assert(!cimport_cache_init(&b, 1));
        CImportCacheBudget budget = {.limit = 1};
        atomic_init(&budget.used, 0);
        double unused = 0;
        CImportFieldRef expected[2] = {
            {UINT64_C(1) << 40, 4},
            {(UINT64_C(1) << 40) + UINT32_MAX, fixed ? 4 : (5 | CIMPORT_FIELD_QUOTED_FLAG)}};
        for (int r = 0; r < 2; r++)
            cimport_cache_append(&b, &expected[r], 1, &unused, expected[r].offset,
                                 expected[r].offset + 5);
        b.numeric[0] = 0;
        assert(!cimport_cache_pack(&b, &budget, SV_missval, false, false, NULL, 0));
        assert(!atomic_load(&budget.used) && b.rows == 2);
        budget.limit = 1024 * 1024;
        CImportCachedBlock *block =
            cimport_cache_pack(&b, &budget, SV_missval, false, false, NULL, 0);
        assert(block);
        for (int r = 0; r < 2; r++)
        {
            CImportFieldRef actual = cimport_cache_field(block, 0, r);
            assert(actual.offset == expected[r].offset && actual.length == expected[r].length);
        }
        cimport_cache_free(block, &budget);
        expected[1].offset++;
        for (int r = 0; r < 2; r++)
            cimport_cache_append(&b, &expected[r], 1, &unused, expected[r].offset,
                                 expected[r].offset + 5);
        assert(!cimport_cache_pack(&b, &budget, SV_missval, false, false, NULL, 0));
        assert(!atomic_load(&budget.used));
        cimport_cache_destroy(&b);
    }
}

int main(void)
{
    mock.macuse = macro;
    mock.spoutsml = display;
    mock.spouterr = display;
    mock.pollstd = poll;
    mock.safestore = number_store;
    mock.sstore = string_store;
    uint64_t bits = UINT64_C(0x7fe0000000000000);
    memcpy(&mock.missval, &bits, 8);
    cache_boundaries();
    char path[] = "/tmp/ctools-stream-XXXXXX";
    int fd = mkstemp(path);
    assert(fd >= 0);
    close(fd);
    const char *cases[] = {
        "id,value,text\n1,1.25,alpha\n2,-0,\"b,c\"\n3,1e-8,\"d\"\"e\"\n4,.,\n5,.a,x\n",
        "id,value,text\r\n1,\"12\",\"a\nb\"\r\n\r\n2,-12.5,z\r\n",
        "id,value,text\n1,2.5,a\n2\n\n3,,b",
        "id,value,text\n1,+2,a\n2,1e300,\"b\n3, 3 ,c\n",
        "id,value,text\n1,\"1,234.5\",a\n2,3,b\n",
        "id,value\n1,2,extra\n2,3,column\n",
        "id,text\n1,short\n2,"
        "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
        "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
        "xxxxxxxxxxxxxxxxxxxx\n",
        "id,value\n",
        "1,2\n3,4\n"};
    for (size_t i = 0; i < sizeof(cases) / sizeof(*cases); i++)
    {
        FILE *f = fopen(path, "wb");
        assert(f);
        fputs(cases[i], f);
        fclose(f);
        for (int mode = 0; mode < 3; mode++)
        {
            for (int empty = 0; empty < 2; empty++)
            {
                bool stream = i != 5 && i != 6 && !(mode == 2 && (i == 0 || i == 4));
                compare(path, mode, empty, i != 8, 0, stream, 1);
            }
        }
    }
    /* Multiple chunks, blanks at boundaries, CRLF, absent final newline. */
    FILE *f = fopen(path, "wb");
    assert(f);
    fputs("id,value,text\n", f);
    for (int i = 0; i < 100000; i++)
        fprintf(f, "%d,%.9f,\"quoted, text\"%s", i, i * .012345, i == 99999 ? "" : "\r\n\n");
    fclose(f);
    for (int mode = 0; mode < 3; mode++)
        compare(path, mode, 0, true, 0, mode != 2, 4);
    strcpy(first, "5");
    strcpy(last, "15");
    compare(path, 0, 0, true, 0, false, 4);
    first[0] = last[0] = 0;
    compare(path, 0, 0, true, 4, true, 4);
    strcpy(strip, "no");
    compare(path, 0, 0, true, 0, true, 4);
    strip[0] = 0;
    strcpy(locale, "en_US");
    compare(path, 0, 0, true, 0, true, 4);
    locale[0] = 0;
    /* Later blocks change the final storage type: preserve original numeric
     * precision on double promotion and spelling on numeric-to-text promotion. */
    for (int text = 0; text < 2; text++)
    {
        f = fopen(path, "wb");
        assert(f);
        fputs("id,value,text\n", f);
        for (int i = 0; i < 10000; i++)
            fprintf(f, "%d,0.12345678901234567,x\n", i);
        fprintf(f, "10000,%s,last\n", text ? "00012x" : "5000000000000.5");
        fclose(f);
        compare(path, 0, 0, true, 0, true, 1);
    }
    assert(released_sources && retained_sources);
    unlink(path);
    ctools_destroy_global_pool();
    return 0;
}
