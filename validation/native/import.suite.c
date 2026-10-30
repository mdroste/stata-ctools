/* Import tests. Select one CTOOLS_CASE_* per compilation.
 * Separate selectors preserve allocator/SPI interposition and independent main functions.
 */

#if defined(CTOOLS_CASE_cimport_options_native_source)
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "stplugin.h"
static ST_plugin mock;
ST_plugin *_stata_=&mock;
static char delimiters[40],collapse[8],literal[8],hex[16001],next[40],encoding[40];
static ST_retcode macro(char *key,char *out,ST_int cap) {
    const char *value="";
    if(!strcmp(key,"___cimport_delimiters"))value=delimiters;
    if(!strcmp(key,"___cimport_collapse"))value=collapse;
    if(!strcmp(key,"___cimport_asstring"))value=literal;
    if(!strcmp(key,"___cimport_encoding_name"))value=encoding;
    snprintf(out,(size_t)cap+1,"%s",value);return 0;
}
static ST_retcode save(char *key,char *value) {
    if(!strcmp(key,"_cimport_blobhex"))strcpy(hex,value);
    if(!strcmp(key,"_cimport_blobnext"))strcpy(next,value);
    return 0;
}
static ST_retcode display(char *s){(void)s;return 0;}
static ST_int poll(void){return 0;}
double ctools_timer_ms(void){return 0;}
#include "cimport/cimport_impl.c"
#include "cimport/cimport_cache.c"
static void normalize(const char *input,const char *expected,int mode) {
    CImportContext ctx={0};ctx.file_data=(char*)input;ctx.file_size=strlen(input);ctx.bindquotes=mode;
    assert(!cimport_normalize_delimiters(&ctx));
    for(size_t i=0;i<strlen(expected);i++)assert(ctx.file_data[i]==(expected[i]=='@'?ctx.delimiter:expected[i]));
    assert(ctx.file_size==strlen(expected));free(ctx.converted_data);
}
int main(void) {
    mock.macuse=macro;mock.macresave=save;mock.spoutsml=display;mock.spouterr=display;mock.pollstd=poll;
    uint64_t bits=UINT64_C(0x7fe0000000000000);memcpy(&mock.missval,&bits,8);
    strcpy(delimiters,"|;");strcpy(collapse,"1");
    normalize("||a;;\"b|c\";\n","@a@\"b|c\"@\n",CIMPORT_BINDQUOTES_LOOSE);
    strcpy(collapse,"0");normalize("a|b;c\n","a@b@c\n",CIMPORT_BINDQUOTES_LOOSE);
    strcpy(delimiters,"||");strcpy(literal,"1");
    normalize("a|b||c\n","a|b@c\n",CIMPORT_BINDQUOTES_LOOSE);
    strcpy(delimiters,"\xce\xa9|");strcpy(literal,"0");
    normalize("a\xce\xa9" "b|c\n","a@b@c\n",CIMPORT_BINDQUOTES_LOOSE);
    char *converted=NULL;size_t converted_size=0;
    const char utf32[]={0,0,0,'a',0,0,0,'b'};
    assert(!cimport_convert_named_to_utf8(utf32,sizeof(utf32),"UTF-32BE",&converted,&converted_size));
    assert(converted_size==2 && !strcmp(converted,"ab"));free(converted);
    assert(cimport_convert_named_to_utf8("x",1,"not-a-charset",&converted,&converted_size)==198);
    assert(!cimport_convert_named_to_utf8("\xff",1,"UTF-8",&converted,&converted_size));
    assert(!converted_size);free(converted);
    assert(!strcmp(cimport_detect_charset(NULL,0),"UTF-8"));
    assert(!strcmp(cimport_detect_charset("id,text\n1,first\n2,second\n",25),"ISO-8859-1"));
    const char utf32bom[]={ (char)0xff,(char)0xfe,0,0,'a',0,0,0 };
    assert(!strcmp(cimport_detect_charset(utf32bom,sizeof(utf32bom)),"UTF-32LE"));
    /* Prefixes and random byte buffers exercise every recognizer under ASan. */
    char sample[16000];uint32_t seed=17;
    for(size_t i=0;i<sizeof(sample);i++) {
        seed^=seed<<13;seed^=seed>>17;seed^=seed<<5;sample[i]=(char)seed;
    }
    for(size_t n=0;n<=8000;n+=31)assert(*cimport_detect_charset(sample,n));
    for(size_t n=0;n<64;n++)assert(*cimport_detect_charset(sample,n));
    const char *before=cimport_detect_charset(sample,8000);
    memset(sample+8000,'x',8000);
    assert(!strcmp(before,cimport_detect_charset(sample,sizeof(sample))));
    double value;
    assert(!cimport_configure_locale("ar_EG","",""));
    assert(cimport_parse_unquoted_number("١٬٢٣٤٫٥٦",(int)strlen("١٬٢٣٤٫٥٦"),&value,SV_missval,'.',0));
    assert(fabs(value-1234.56)<1e-10);
    assert(!cimport_configure_locale("fa_IR","",""));
    assert(cimport_parse_unquoted_number("‎−۱۲۳۴٫۵۶",(int)strlen("‎−۱۲۳۴٫۵۶"),&value,SV_missval,'.',0));
    assert(fabs(value+1234.56)<1e-10);
    assert(cimport_configure_locale("bogus","","")==198);
    assert(!cimport_configure_locale("en_US","",""));
    assert(!cimport_parse_unquoted_number("+123",4,&value,SV_missval,'.',0));
    assert(cimport_parse_unquoted_number("∞",3,&value,SV_missval,'.',0) && value==SV_missval);
    assert(!cimport_configure_locale("","",""));
    /* Empty files must still validate charsets; conversion receives NULL/0. */
    delimiters[0]=collapse[0]=literal[0]=0;
    char emptyfile[]="/tmp/ctools-cimport-empty-XXXXXX";
    int fd=mkstemp(emptyfile);assert(fd>=0);close(fd);
    CImportContext *empty=cimport_parse_csv(emptyfile,0,true,0,false,CIMPORT_BINDQUOTES_LOOSE,
        CIMPORT_NUMTYPE_AUTO,'.',0,CIMPORT_EMPTYLINES_SKIP,20,CIMPORT_ENC_UNKNOWN,NULL,0,NULL,0);
    assert(empty && empty->num_columns==0 && empty->delimiter==',');
    assert(!strcmp(empty->reported_encoding,"UTF-8"));cimport_free_context(empty);
    strcpy(encoding,"latin1");
    empty=cimport_parse_csv(emptyfile,0,true,0,false,CIMPORT_BINDQUOTES_LOOSE,
        CIMPORT_NUMTYPE_AUTO,'.',0,CIMPORT_EMPTYLINES_SKIP,20,CIMPORT_ENC_UNKNOWN,NULL,0,NULL,0);
    assert(empty && !strcmp(empty->reported_encoding,"ISO-8859-1"));cimport_free_context(empty);
    strcpy(encoding,"bogus");
    empty=cimport_parse_csv(emptyfile,0,true,0,false,CIMPORT_BINDQUOTES_LOOSE,
        CIMPORT_NUMTYPE_AUTO,'.',0,CIMPORT_EMPTYLINES_SKIP,20,CIMPORT_ENC_UNKNOWN,NULL,0,NULL,0);
    assert(!empty && g_cimport_parse_status==198);encoding[0]=0;unlink(emptyfile);
    /* Ignore decoded NULs without truncating fields or mutating the mapping.
     * Exercise both the mapped ASCII path and the owned conversion buffer. */
    const char zeros[]="i\0d,name\n1\0,te\0st\n2,\"n\0ormal\"\n";
    char nullfile[]="/tmp/ctools-cimport-nuls-XXXXXX";
    fd=mkstemp(nullfile);assert(fd>=0);
    assert(write(fd,zeros,sizeof(zeros)-1)==sizeof(zeros)-1);close(fd);
    for(int explicit=0;explicit<2;explicit++) {
        strcpy(encoding,explicit?"UTF-8":"");
        CImportContext *clean=cimport_parse_csv(nullfile,',',true,0,false,CIMPORT_BINDQUOTES_LOOSE,
            CIMPORT_NUMTYPE_AUTO,'.',0,CIMPORT_EMPTYLINES_SKIP,20,CIMPORT_ENC_UNKNOWN,NULL,0,NULL,0);
        assert(clean && clean->num_columns==2);
        const char expected[]="id,name\n1,test\n2,\"normal\"\n";
        assert(clean->file_size==sizeof(expected)-1);
        assert(!memcmp(clean->file_data,expected,sizeof(expected)-1));
        cimport_free_context(clean);
    }
    encoding[0]=0;unlink(nullfile);
    CImportContext numeric={0};numeric.file_data="1\"\"2";numeric.strip_quotes_all=true;
    CImportFieldRef number={0,4|CIMPORT_FIELD_QUOTED_FLAG};
    assert(cimport_context_number(&numeric,&number,&value) && value==12);
    numeric.file_data="NA";number=(CImportFieldRef){0,2};
    assert(!cimport_context_number(&numeric,&number,&value));
    numeric.file_data="1\"\"2";number=(CImportFieldRef){0,4|CIMPORT_FIELD_QUOTED_FLAG};
    numeric.strip_quotes_all=false;numeric.retain_quotes=true;
    assert(!cimport_context_number(&numeric,&number,&value));
    CImportContext limit={0};limit.file_data="1,\"a\nb\nc\"\n";limit.file_size=strlen(limit.file_data);
    limit.bindquotes=CIMPORT_BINDQUOTES_STRICT;limit.max_quoted_rows=2;
    assert(cimport_check_quoted_rows(&limit)==5101);limit.max_quoted_rows=3;assert(!cimport_check_quoted_rows(&limit));
    CImportContext *ctx=calloc(1,sizeof(*ctx));pthread_mutex_init(&ctx->warning_mutex,NULL);
    ctx->num_columns=1;ctx->total_rows=1;ctx->num_chunks=1;ctx->num_threads=1;ctx->quote_char='"';
    ctx->columns=calloc(1,sizeof(*ctx->columns));ctx->columns[0].type=CIMPORT_COL_STRING;ctx->columns[0].max_strlen=20000;
    ctx->columns[0].use_strl=true;
    ctx->chunks=calloc(1,sizeof(*ctx->chunks));ctx->chunks[0].num_rows=1;
    ctx->chunks[0].rows=calloc(1,sizeof(CImportParsedRow*));
    CImportParsedRow *row=ctools_arena_alloc(&ctx->chunks[0].arena,cimport_row_size(1));
    CImportFieldRef long_field={0,20000};
    assert(row && cimport_fill_row(row,&long_field,1));cimport_row_values(row)[0]=CIMPORT_NO_VALUE;
    ctx->chunks[0].rows[0]=row;
    ctx->file_data=ctx->converted_data=malloc(20001);ctx->file_size=20000;memset(ctx->file_data,'x',20000);ctx->file_data[20000]=0;
    atomic_init(&ctx->error_code,0);g_cimport_ctx=ctx;
    assert(!cimport_blob(1,1,0));assert(strlen(hex)==16000 && !strcmp(next,"8000"));
    assert(!cimport_blob(1,1,16000));assert(strlen(hex)==8000 && !strcmp(next,"0"));
    assert(strlen(ctx->col_cache[0].string_data[0])==20000);
    assert(cimport_blob(1,1,20001)==198);assert(cimport_blob(2,1,0)==198);
    cimport_free_context(ctx);g_cimport_ctx=NULL;
    ctools_destroy_global_pool();return 0;
}
#endif

#if defined(CTOOLS_CASE_cimport_stream_native)
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
#include "cimport/cimport_cache.c"

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
#endif

#if defined(CTOOLS_CASE_p2_native_2)
#include <assert.h>
#include "cimport/cimport_encoding.c"
int main(void) {
    const char le[] = "\xff\xfe\0\0a\0\0\0";
    const char be[] = "\0\0\xfe\xff\0\0\0a";
    CImportEncodingDetection a = cimport_detect_encoding(le, sizeof(le)-1);
    CImportEncodingDetection b = cimport_detect_encoding(be, sizeof(be)-1);
    assert(a.encoding == CIMPORT_ENC_UTF32LE && a.bom_length == 4);
    assert(b.encoding == CIMPORT_ENC_UTF32BE && b.bom_length == 4);
    return 0;
}
#endif

#if defined(CTOOLS_CASE_p2_native_4)
#include "cimport/cimport_impl.c"
#include "cimport/cimport_cache.c"
int main(void) {
    setup(); nrows = 3; bad_row = 2;
    /* One column of three rows ("1", "2", "3"): a store tile must stop at the
     * failing observation and report its code, for string and numeric columns. */
    CImportParsedRow *rows[3];
    for (int i = 0; i < 3; i++) {
        rows[i] = calloc(1, cimport_row_size(1));
        CImportFieldRef field = {.offset = (uint64_t)i, .length = 1};
        assert(cimport_fill_row(rows[i], &field, 1));
        double value = i + 1;
        memcpy(&cimport_row_values(rows[i])[0], &value, sizeof(value));
    }
    CImportParsedChunk chunk = {.rows = rows, .num_rows = 3};
    CImportContext ctx = {.chunks = &chunk, .num_chunks = 1, .total_rows = 3, .num_columns = 1,
                          .file_data = "123", .decimal_separator = '.', .quote_char = '"'};
    unsigned char kind[1] = {CIMPORT_STORE_STRING};
    string_var = 1; writes = 0;
    CImportStoreTile tile = {&ctx, kind, 0, 0, 3, 1, 0};
    assert(cimport_store_tile_worker(&tile) != NULL);
    assert(tile.error == 459 && writes == 2);
    string_var = 0; kind[0] = CIMPORT_STORE_NUMERIC; writes = 0;
    tile = (CImportStoreTile){&ctx, kind, 0, 0, 3, 1, 0};
    assert(cimport_store_tile_worker(&tile) != NULL);
    assert(tile.error == 459 && writes == 2);
    /* strL columns are transferred by the ado: the tile writes nothing */
    kind[0] = CIMPORT_STORE_STRL; writes = 0;
    tile = (CImportStoreTile){&ctx, kind, 0, 0, 3, 1, 0};
    assert(cimport_store_tile_worker(&tile) == NULL && tile.error == 0 && writes == 0);
    for (int i = 0; i < 3; i++) free(rows[i]);
    return 0;
}
#endif

#if defined(CTOOLS_CASE_p2_native_5)
#include "cimport/cimport_xlsx.c"
int main(void) {
    setup(); nrows = 3; string_var = 1; bad_row = 2;
    char *strings[] = {"a", "b", "c"}; double values[] = {1, 2, 3};
    CImportColumnInfo col = {.type = CIMPORT_COL_STRING};
    CImportColumnCache cache = {.count = 3, .string_data = strings, .numeric_data = values};
    xlsx_store_task task = {.col = &col, .cache = &cache, .var = 1, .stata_nobs = 3};
    assert(xlsx_store_worker(&task) != NULL && task.error == 459);
    string_var = 0; col.type = CIMPORT_COL_NUMERIC;
    assert(xlsx_store_worker(&task) != NULL && task.error == 459);
    task.stata_nobs = 2; writes = 0;
    assert(xlsx_store_worker(&task) != NULL && task.error == 198 && writes == 0);
    return 0;
}
#endif

#if defined(CTOOLS_CASE_io_perf_csv)
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "stplugin.h"
static ST_plugin mock;
ST_plugin *_stata_=&mock;
static char explicit_encoding[40];
static ST_retcode macro(char *key,char *out,ST_int cap) {
    snprintf(out,(size_t)cap+1,"%s",!strcmp(key,"___cimport_encoding_name")?explicit_encoding:"");
    return 0;
}
static ST_retcode display(char *s){(void)s;return 0;}
static ST_int poll(void){return 0;}
double ctools_timer_ms(void){return 0;}
#include "cimport/cimport_impl.c"
#include "cimport/cimport_cache.c"
int main(void) {
    mock.macuse=macro;mock.spoutsml=display;mock.spouterr=display;mock.pollstd=poll;
    uint64_t bits=UINT64_C(0x7fe0000000000000);memcpy(&mock.missval,&bits,8);
    ctools_set_max_threads(1);
    char data[32000]="v1,v2\n";
    for(int i=0;i<1000;i++)strcat(data,"123.456789,-678.901234\n");
    size_t n=strlen(data);
    const char *charset=cimport_detect_charset(data,n);
    assert(!strcmp(charset,"ISO-8859-2"));
    for(int nonascii=0;nonascii<2;nonascii++) {
        if(nonascii)strcpy(data+n,"1,\xa3\n");
        size_t size=strlen(data);
        char path[]="/tmp/ctools-ascii-import-XXXXXX";
        int fd=mkstemp(path);assert(fd>=0);
        assert(write(fd,data,size)==(ssize_t)size);close(fd);
        CImportContext *ctx=cimport_parse_csv(path,',',true,0,false,CIMPORT_BINDQUOTES_LOOSE,
            CIMPORT_NUMTYPE_AUTO,'.',0,CIMPORT_EMPTYLINES_SKIP,20,CIMPORT_ENC_UNKNOWN,NULL,0,NULL,0);
        assert(ctx && !strcmp(ctx->reported_encoding,charset));
        char *expected=NULL;size_t expected_size;
        assert(!cimport_convert_named_to_utf8(data,size,charset,&expected,&expected_size));
        assert(ctx->file_size==expected_size && !memcmp(ctx->file_data,expected,expected_size));
        assert(ctx->total_rows==(size_t)(1000+nonascii));
        if(!nonascii)assert(ctx->converted_data==NULL);
        else assert(ctx->converted_data!=NULL);
        free(expected);cimport_free_context(ctx);
        strcpy(explicit_encoding,"bogus");
        ctx=cimport_parse_csv(path,',',true,0,false,CIMPORT_BINDQUOTES_LOOSE,
            CIMPORT_NUMTYPE_AUTO,'.',0,CIMPORT_EMPTYLINES_SKIP,20,CIMPORT_ENC_UNKNOWN,NULL,0,NULL,0);
        assert(!ctx && g_cimport_parse_status==198);explicit_encoding[0]=0;
        unlink(path);
    }
    ctools_destroy_global_pool();
}
#endif
