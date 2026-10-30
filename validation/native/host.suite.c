/* Host tests. Select one CTOOLS_CASE_* per compilation.
 * Separate selectors preserve allocator/SPI interposition and independent main functions.
 */

#if defined(CTOOLS_CASE_contracts_native_source)
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "stplugin.h"
#include "ctools_config.h"
#include "ctools_parse.h"
#include "ctools_runtime.h"
#include "ctools_types.h"
#include "ctools_ols.h"
static ST_plugin mock;
ST_plugin *_stata_ = &mock;
static int freed[32];
void cmerge_cleanup_cache(void) { freed[CTOOLS_CMD_CMERGE]++; }
void cimport_cleanup_cache(void) { freed[CTOOLS_CMD_CIMPORT]++; }
void cio_cleanup(void) { freed[CTOOLS_CMD_CIO]++; }
void csplit_cleanup_cache(void) { freed[CTOOLS_CMD_CSPLIT]++; }
void crangejoin_cleanup_cache(void) { freed[CTOOLS_CMD_CRANGEJOIN]++; }
void creghdfe_cleanup_state(void) { freed[CTOOLS_CMD_CREGHDFE]++; }
void civreghdfe_cleanup_state(void) { freed[CTOOLS_CMD_CIVREGHDFE]++; }
void cppmlhdfe_cleanup_state(void) { freed[CTOOLS_CMD_CPPMLHDFE]++; }
void cexport_cleanup_state(void) { freed[CTOOLS_CMD_CEXPORT]++; }
void cwinsor_cleanup_state(void) { freed[CTOOLS_CMD_CWINSOR]++; }
static const ctools_command_descriptor merge={"cmerge",CTOOLS_CMD_CMERGE,NULL,cmerge_cleanup_cache};
static const ctools_command_descriptor import={"cimport",CTOOLS_CMD_CIMPORT,NULL,cimport_cleanup_cache};
static const ctools_command_descriptor io={"cio",CTOOLS_CMD_CIO,NULL,cio_cleanup};
static const ctools_command_descriptor join={"crangejoin",CTOOLS_CMD_CRANGEJOIN,NULL,crangejoin_cleanup_cache};
static const ctools_command_descriptor winsor={"cwinsor",CTOOLS_CMD_CWINSOR,NULL,cwinsor_cleanup_state};
static void invoke(const ctools_command_descriptor *d,const char *phase,int rc) {
 assert(ctools_command_begin(d,phase)==0); ctools_command_finish(rc);
}
static void lifecycle(void) {
 invoke(&merge,"load_using",0); invoke(&merge,"plan",0); invoke(&merge,"execute",0);
 assert(freed[CTOOLS_CMD_CMERGE]==0); // execute keeps its undo copy for commit or rollback
 invoke(&merge,"commit",0);
 assert(freed[CTOOLS_CMD_CMERGE]==1);
 invoke(&merge,"load_using",0); invoke(&merge,"plan",0); invoke(&merge,"execute",0);
 invoke(&merge,"rollback",0);
 assert(freed[CTOOLS_CMD_CMERGE]==2);
 invoke(&merge,"load_using",0); assert(ctools_command_begin(&merge,"execute")==198);
 assert(freed[CTOOLS_CMD_CMERGE]==3); // execute requires plan
 invoke(&merge,"load_using",0); invoke(&import,"scan",0); invoke(&merge,"load_using",0);
 assert(freed[CTOOLS_CMD_CIMPORT]==1); // Regression: A -> B -> A.
 invoke(&import,"scan",0); invoke(&io,"scan",0);
 assert(freed[CTOOLS_CMD_CIMPORT]==2); // A -> B -> C, distinct CSV and generic IO.
 invoke(&io,"column 1",0); invoke(&io,"load",0); invoke(&io,"blob 1 1 0",459);
 assert(freed[CTOOLS_CMD_CIO]==1);
 assert(ctools_command_begin(&io,"column 1")==198);
 invoke(&join,"using 1",0); invoke(&join,"prepare 1",920);
 assert(ctools_command_begin(&join,"write")==198);
 invoke(&join,"using 1",0); invoke(&join,"prepare 1",0); invoke(&join,"write",0);
 invoke(&import,"scan",459); assert(ctools_command_begin(&import,"blob 1")==198);
 invoke(&import,"scan",0); invoke(&import,"scan",0); // Interrupted/restarted scan frees old state.
 int count=freed[CTOOLS_CMD_CIMPORT];
 ctools_command_reset(); ctools_command_reset(); assert(freed[CTOOLS_CMD_CIMPORT]==count+1);
 invoke(&import,"scan",0); assert(ctools_command_begin(NULL,"")==198);
 assert(freed[CTOOLS_CMD_CIMPORT]==count+2);
 // crangejoin write keeps the master copy for rollback or clear; cwinsor keeps
 // an in-place undo copy until commit or rollback, or until another command.
 int joins=freed[CTOOLS_CMD_CRANGEJOIN], winsors=freed[CTOOLS_CMD_CWINSOR];
 invoke(&join,"using 1",0); invoke(&join,"prepare 1",0); invoke(&join,"write",0);
 assert(freed[CTOOLS_CMD_CRANGEJOIN]==joins);
 invoke(&join,"rollback",0); assert(freed[CTOOLS_CMD_CRANGEJOIN]==joins+1);
 invoke(&winsor,"1 0 1 undo",0); assert(freed[CTOOLS_CMD_CWINSOR]==winsors);
 invoke(&winsor,"commit",0); assert(freed[CTOOLS_CMD_CWINSOR]==winsors+1);
 invoke(&winsor,"1 0 1 undo",0); invoke(&winsor,"rollback",0);
 assert(freed[CTOOLS_CMD_CWINSOR]==winsors+2);
 invoke(&winsor,"1 0 1 undo",0); invoke(&join,"using 1",0);
 assert(freed[CTOOLS_CMD_CWINSOR]==winsors+3);
 ctools_command_reset();
 assert(ctools_command_begin(&winsor,"apply")==198);
 winsors=freed[CTOOLS_CMD_CWINSOR];
 invoke(&winsor,"1 0 1 deferbounds undo",0);
 invoke(&winsor,"apply",0);
 assert(freed[CTOOLS_CMD_CWINSOR]==winsors);
 invoke(&winsor,"rollback",0);
 assert(freed[CTOOLS_CMD_CWINSOR]==winsors+1);
 invoke(&winsor,"1 0 1 deferbounds",0); invoke(&winsor,"apply",0);
 assert(ctools_command_begin(&winsor,"apply")==198);
 invoke(&winsor,"1 0 1 deferbounds",0); invoke(&winsor,"apply",459);
 assert(ctools_command_begin(&winsor,"apply")==198);
}
static void parsing(void) {
 char b[32]="untouched";int n=47;uint64_t u=47;double d=47;
 assert(ctools_parse_bool_option("notverbose verbose","verbose"));
 assert(!ctools_parse_bool_option("notverbose","verbose"));
 assert(ctools_parse_string_option("nolabel=wrong label=right","label",b,sizeof b)==1&&!strcmp(b,"right"));
 assert(ctools_parse_string_option("label=\"with spaces\"","label",b,sizeof b)==1&&!strcmp(b,"with spaces"));
 assert(ctools_parse_string_option("label=\"unterminated","label",b,sizeof b)==-1);
 assert(ctools_parse_string_option("label=x label=y","label",b,sizeof b)==-1);
 assert(ctools_parse_string_option("label=1234","label",b,4)==-1);
 assert(ctools_parse_int_option("n=2147483647","n",&n)==1&&n==INT_MAX);
 assert(ctools_parse_int_option("n=-2147483648","n",&n)==1&&n==INT_MIN);
 for(size_t i=0;i<6;i++){
  const char *bad[]={"n=2147483648","n=-2147483649","n=1x","n=","n","n=1 n=2"};
  n=47;assert(ctools_parse_int_option(bad[i],"n",&n)==-1&&n==47);
 }
 assert(ctools_parse_int_option("other=3","n",&n)==0&&n==47);
 assert(ctools_parse_u64_option("n=18446744073709551615","n",&u)==1&&u==UINT64_MAX);
 assert(ctools_parse_u64_option("n=18446744073709551616","n",&u)==-1);
 assert(ctools_parse_u64_option("n=-1","n",&u)==-1);
 assert(ctools_parse_double_checked("d=nan","d",&d)==-1&&d==47);
 assert(ctools_parse_double_checked("d=1e999","d",&d)==-1);
 assert(ctools_parse_double_checked("d=1.5junk","d",&d)==-1);
 assert(ctools_parse_seed_option("seedhi=2147483647 seedlo=2147483647",&u)==1&&u==((UINT64_C(1)<<62)-1));
 assert(ctools_parse_seed_option("seedhi=1",&u)==-1);
 assert(ctools_parse_seed_option("seedhi=2147483648 seedlo=0",&u)==-1);
 const char *p=" \t2147483648",*saved=p;assert(ctools_parse_next_int(&p,&n)==-1&&p==saved);
 const char bounded[]={'a','=','1',' ','b'};ctools_arg_cursor c;ctools_arg_token t;
 ctools_args_init(&c,bounded,3);assert(ctools_args_next(&c,&t)==1&&t.length==3);assert(ctools_args_next(&c,&t)==0);
}
/* Model Windows allocator separation on every host: a wrong deallocator aborts.
 * Real Windows builds use _aligned_free through the same explicit storage mode. */
static void *aligned[32];static int live;
static void *alloc_aligned(size_t n) {void *p=ctools_cacheline_alloc(n);assert(p);aligned[live++]=p;return p;}
static void release_aligned(void *p) {if(!p)return;int found=0;for(int i=0;i<live;i++)if(aligned[i]==p){aligned[i]=aligned[--live];found=1;break;}assert(found);ctools_aligned_free(p);}
static void release_heap(void *p) {for(int i=0;i<live;i++)assert(aligned[i]!=p);free(p);}
#define ctools_aligned_free release_aligned
#define free release_heap
#include "ctools_types.c"
#include "ctools_platform.c"
#undef free
#undef ctools_aligned_free
static void ownership(void) {
 for(int stage=0;stage<5;stage++){
  stata_data d={.nobs=2,.nvars=2};
  if(stage>=1)d.sort_order=alloc_aligned(2*sizeof(perm_idx_t));
  if(stage>=2){d.vars=alloc_aligned(2*sizeof(stata_variable));memset(d.vars,0,2*sizeof(stata_variable));}
  if(stage>=3){stata_variable *v=d.vars;v->nobs=2;v->data.dbl=alloc_aligned(2*sizeof(double));}
  if(stage>=4){stata_variable *v=d.vars+1;v->type=STATA_TYPE_STRING;v->nobs=2;v->buffer_storage=CTOOLS_BUFFER_HEAP;v->string_storage=CTOOLS_STRINGS_BORROWED;v->data.str=calloc(2,sizeof(char*));v->data.str[0]="borrowed";v->data.str[1]="literal";}
  stata_data_free(&d);stata_data_free(&d);assert(live==0&&!d.vars&&!d.sort_order);
 }
 char *owner=strdup("survives"); char *view[]={owner};
 stata_variable v={.type=STATA_TYPE_STRING,.nobs=1,.data.str=view,.buffer_storage=CTOOLS_BUFFER_BORROWED,.string_storage=CTOOLS_STRINGS_BORROWED};
 perm_idx_t order[]={0};stata_data view_data={.nobs=1,.nvars=1,.vars=&v,.sort_order=order};
 assert(ctools_apply_permutation(&view_data)==STATA_ERR_UNSUPPORTED_TYPE);
 stata_variable_free(&v);assert(!strcmp(owner,"survives"));free(owner);
}
int main(void) {
 lifecycle();parsing();ownership();
 int rank=41;assert(ctools_cholesky(NULL,0)==STATA_OK);
 assert(ctools_invsym(NULL,0,NULL,&rank)==STATA_OK&&rank==0);
 double bad[]={1,2,2,1}, b[]={1,1}, x[]={41,42};
 assert(ctools_solve_cholesky(NULL,b,2,x)==STATA_ERR_INVALID_INPUT);
 assert(ctools_solve_cholesky(bad,b,2,x)==STATA_ERR_SINGULAR);
 assert(x[0]==41&&x[1]==42);
 assert(ctools_spi_status(1,STATA_ERR_STATA_READ)==STATA_ERR_CANCELLED);
 assert(ctools_spi_status(459,STATA_ERR_STATA_WRITE)==STATA_ERR_STATA_WRITE);
 assert(ctools_stata_rc(STATA_ERR_MEMORY)==920);
 assert(ctools_stata_rc(STATA_ERR_INVALID_INPUT)==198);
 assert(ctools_stata_rc(STATA_ERR_STATA_READ)==459);
 assert(ctools_stata_rc(STATA_ERR_STATA_WRITE)==459);
 assert(ctools_stata_rc(STATA_ERR_CANCELLED)==1);
 assert(ctools_stata_rc(STATA_ERR_UNSUPPORTED_TYPE)==109);
 puts("PASS parser, lifecycle, partial ownership, allocator pairing, and status contracts");
}
#endif

#if defined(CTOOLS_CASE_refactor_boundaries)
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdint.h>
#include <string.h>
#include "stplugin.h"
#include "ctools_rng.h"
#include "ctools_group.h"
#include "ctools_platform.h"
#include "ctools_runtime.h"
#include "formats/excel_address.h"

static ST_plugin mock;
ST_plugin *_stata_ = &mock;
static char message[1024];
static ST_retcode capture(char *text)
{
    assert(strlen(text) < sizeof(message));
    strcpy(message, text);
    return 0;
}

static atomic_int ready, go;
static ctools_prefetch_distances snapshots[16];
static void *worker(void *argument)
{
    size_t i = (size_t)argument;
    atomic_fetch_add(&ready, 1);
    while (!atomic_load(&go))
    {
    }
    snapshots[i] = ctools_get_prefetch_distances();
    double first = ctools_timer_seconds();
    for (int j = 0; j < 1000; j++)
        assert(ctools_timer_seconds() >= first);
    return NULL;
}

int main(void)
{
    /* Simultaneous first use must publish fully initialized process-wide state. */
    pthread_t threads[16];
    for (size_t i = 0; i < 16; i++)
        assert(!pthread_create(&threads[i], NULL, worker, (void *)i));
    while (atomic_load(&ready) != 16)
    {
    }
    atomic_store(&go, 1);
    for (int i = 0; i < 16; i++)
        assert(!pthread_join(threads[i], NULL));
    for (int i = 0; i < 16; i++)
    {
        assert(snapshots[i].near_dist > 0);
        assert(snapshots[i].near_dist == snapshots[0].near_dist);
        assert(snapshots[i].far_dist == snapshots[0].far_dist);
        assert(snapshots[i].stream == snapshots[0].stream);
        assert(snapshots[i].general == snapshots[0].general);
    }

    /* Published xoshiro256** transition for state {1,2,3,4}; no seeding ambiguity. */
    ctools_rng rng = {{1, 2, 3, 4}};
    assert(ctools_rng_uniform(&rng, 0) == 0);
    assert(ctools_rng_next(&rng) == UINT64_C(11520));
    assert(ctools_rng_next(&rng) == UINT64_C(0));
    assert(ctools_rng_next(&rng) == UINT64_C(1509978240));

    int col = -1, row = -1;
    assert(ctools_excel_address("XFD1048576", 16384, 1048576, false, &col, &row, NULL));
    assert(col == 16384 && row == 1048576);
    const char *bad[] = {"XFE1",
                         "A1048577",
                         "A0",
                         "A",
                         "1",
                         "$A$1",
                         "A1suffix",
                         "ZZZZZZZZZZZZZZ1",
                         "A9999999999999999999"};
    for (size_t i = 0; i < sizeof(bad) / sizeof(*bad); i++)
    {
        col = row = -1;
        assert(!ctools_excel_address(bad[i], 16384, 1048576, false, &col, &row, NULL));
        assert(col == -1 && row == -1);
    }
    const char *end;
    assert(ctools_excel_address("$az$12:B20", 16384, 1048576, true, &col, &row, &end));
    assert(col == 52 && row == 12 && !strcmp(end, ":B20"));

    double values[] = {0.0, -0.0, 8.98846567431158e307, 8.98846567431158e307,
                       nextafter(8.98846567431158e307, INFINITY)};
    char *strings[] = {NULL, "", "a", "a", "a"};
    stata_variable vars[2] = {{.type = STATA_TYPE_DOUBLE, .data.dbl = values},
                              {.type = STATA_TYPE_STRING, .data.str = strings}};
    stata_data data = {.nobs = 5, .nvars = 2, .vars = vars};
    ctools_group_range groups[5];
    size_t count;
    ctools_group_ranges(&data, NULL, NULL, 2, groups, &count);
    assert(count == 3 && groups[0].count == 2 && groups[1].count == 2);
    perm_idx_t order[] = {4, 3, 2, 1, 0};
    ctools_group_ranges(&data, order, NULL, 2, groups, &count);
    assert(count == 3 && groups[0].count == 1 && groups[1].start == 1 && groups[2].count == 2);

    mock.spouterr = capture;
    mock.spoutsml = capture;
    char oversized[4096];
    memset(oversized, 'x', sizeof(oversized) - 1);
    oversized[sizeof(oversized) - 1] = 0;
    ctools_error(oversized, "%s", oversized);
    assert(strlen(message) == 1023 && message[1022] == '\n');
    ctools_msg("unit", "%s", oversized);
    assert(strlen(message) == 1023 && !strncmp(message, "unit: ", 6) && message[1022] == '\n');
    ctools_verbose("unit", 0, "hidden");
    assert(strlen(message) == 1023);
    ctools_msg(NULL, "complete\n");
    assert(!strcmp(message, "complete\n"));
    return 0;
}
#endif

#if defined(CTOOLS_CASE_spi_matrix_native_source)
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "ctools_spi.h"
static ST_plugin host;
ST_plugin *_stata_ = &host;
static double *matrix;
static int checked_calls, unchecked_calls, fail_write;
static ST_retcode message(char *s) { (void)s; return 0; }
static ST_retcode save(char *s, double v) { (void)s; (void)v; return 0; }
static ST_retcode unchecked(char *s, ST_int row, ST_int col, double value) {
    (void)s;
    unchecked_calls++;
    matrix[(row - 1) * 2 + col - 1] = value;
    return 0;
}
static ST_retcode checked(char *s, ST_int row, ST_int col, double value) {
    checked_calls++;
    if (strcmp(s, "result")) return 111;
    if (row < 1 || row > 2 || col < 1 || col > 2) return 503;
    if (fail_write) return 459;
    matrix[(row - 1) * 2 + col - 1] = value;
    return 0;
}
int main(void) {
    matrix = calloc(4, sizeof(double));
    assert(matrix);
    host.matstore = unchecked;
    host.safematstore = checked;
    host.spouterr = message;
    host.scalsave = save;
    /* Under SD_FASTMODE the old wrapper overran this buffer (ASan). */
    assert(ctools_mat_store("result", 3, 1, 42) == 503);
    assert(ctools_mat_store("result", 1, 3, 42) == 503);
    assert(ctools_mat_store("result", 0, 1, 42) == 503);
    assert(ctools_mat_store("result", 1, 0, 42) == 503);
    assert(ctools_mat_store("absent", 1, 1, 42) == 111);
    fail_write = 1;
    assert(ctools_mat_store("result", 1, 1, 42) == 459);
    for (int i = 0; i < 4; i++) assert(matrix[i] == 0);
    fail_write = 0;
    assert(ctools_mat_store("result", 2, 2, 42) == 0);
    assert(matrix[3] == 42 && checked_calls == 7 && unchecked_calls == 0);
    free(matrix);
    return 0;
}
#endif

#if defined(CTOOLS_CASE_spi_matrix_native_publication)
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "ctools_spi.h"
static ST_plugin host;
ST_plugin *_stata_ = &host;
static int failed_alloc, write_calls, fail_at, host_K = 3;
static double got_b[3], got_V[9];
static void *map_malloc(size_t size) { return failed_alloc ? NULL : malloc(size); }
#define malloc map_malloc
#include "civreghdfe/civreghdfe_impl.c"
#include "civreghdfe/civreghdfe_result.c"
#include "civreghdfe/civreghdfe_input.c"
#undef malloc
#include "cqreg/cqreg_regress.c"
static ST_retcode message(char *s) { (void)s; return 0; }
static ST_retcode scalar_save(char *s, double v) { (void)s; (void)v; return 0; }
static ST_retcode unsafe_matrix(char *s, ST_int row, ST_int col, double v) {
    (void)s; (void)row; (void)col; (void)v;
    assert(!"unchecked matrix callback"); return 999;
}
static ST_retcode checked_matrix(char *s, ST_int row, ST_int col, double v) {
    write_calls++;
    if (write_calls == fail_at) return 459;
    int b = !strcmp(s, "__civreghdfe_b");
    if (row < 1 || row > (b ? 1 : host_K) || col < 1 || col > host_K) return 503;
    if (b) got_b[col - 1] = v;
    else got_V[(row - 1)*host_K + col - 1] = v;
    return 0;
}
int main(void) {
    host.matstore = unsafe_matrix;
    host.safematstore = checked_matrix;
    host.spouterr = message;
    host.scalsave = scalar_save;
    /* Internal [exog, omitted exog, endog] -> posted [endog, exog, omitted]. */
    double b[] = {2, 3}, V[] = {4, 1, 1, 9};
    int omitted[] = {0, 1, 0};
    assert(store_iv_matrices(b, V, 2, 1, 2, omitted) == 0);
    const double expected_b[] = {3, 2, 0}, expected_V[] = {9,1,0, 1,4,0, 0,0,0};
    assert(!memcmp(got_b, expected_b, sizeof got_b));
    assert(!memcmp(got_V, expected_V, sizeof got_V));
    for (int f = 1; f <= 12; f++) {
        write_calls = 0; fail_at = f;
        assert(store_iv_matrices(b, V, 2, 1, 2, omitted) == 459);
        assert(write_calls == f);
    }
    fail_at = 0; failed_alloc = 1; write_calls = 0;
    assert(store_iv_matrices(b, V, 2, 1, 2, omitted) == 920);
    assert(write_calls == 0);
    failed_alloc = 0; host_K = 2;
    assert(store_iv_matrices(b, V, 2, 1, 2, omitted) == 503);
    host_K = 3;
    assert(store_iv_matrices(b, V, 2, 1, 1, omitted) == 503);
    assert(store_iv_matrices(b, V, 2, 1, 2, omitted) == 0);
    cqreg_state state = {0};
    state.K = 2; state.N = 12; state.beta = b; state.V = V;
    for (int f = 1; f <= 4; f++) {
        write_calls = 0; fail_at = f;
        assert(store_results(&state) == 459);
        assert(write_calls == f);
    }
    fail_at = 0; host_K = 1;
    assert(store_results(&state) == 503);
    host_K = 2;
    assert(store_results(&state) == 0);
    assert(!memcmp(got_V, V, sizeof V));
    return 0;
}
#endif

#if defined(CTOOLS_CASE_thread_policy)
/* Thread policy: parsing, precedence, CPU limits, cgroup quotas and the
 * per-user settings file. Modes:
 *   unit <tmpdir> <expected settings path>
 *   report    prints the policy resolved from this process's environment */
#undef NDEBUG  /* some compilers define it at -O1; these asserts are the test */
#include "ctools_threads.c"
#include <assert.h>

static ctools_thread_policy resolve(int available, int omp_limit, int env, int saved, int omp)
{
    ctools_thread_policy p;
    memset(&p, 0, sizeof(p));
    p.detected = available;
    p.available = available;
    p.omp_limit = omp_limit;
    p.env_ctools = env;
    p.saved = saved;
    p.env_omp = omp;
    ctools_thread_policy_resolve(&p);
    return p;
}

static void put(const char *path, const char *text)
{
    FILE *fp = ctools_fopen(path, "wb");
    assert(fp != NULL);
    fputs(text, fp);
    assert(fclose(fp) == 0);
}

/* Text mode: the file uses the platform's line endings. */
static void get(const char *path, char *text, size_t size)
{
    FILE *fp = ctools_fopen(path, "r");
    size_t n;
    assert(fp != NULL);
    n = fread(text, 1, size - 1, fp);
    text[n] = '\0';
    fclose(fp);
}

/* Whether the directory holding a file still exists. */
static int parent_exists(const char *path)
{
    char dir[4096];
    snprintf(dir, sizeof(dir), "%s", path);
    *strrchr(dir, CTOOLS_PATH_SEP) = '\0';
#if defined(_WIN32)
    wchar_t *wide = ctools_wpath(dir);
    DWORD attributes = wide ? GetFileAttributesW(wide) : INVALID_FILE_ATTRIBUTES;
    free(wide);
    return attributes != INVALID_FILE_ATTRIBUTES;
#else
    struct stat st;
    return stat(dir, &st) == 0;
#endif
}

static void parsing(void)
{
    assert(ctools_parse_thread_value(NULL, 0) == 0);
    assert(ctools_parse_thread_value("", 0) == 0);
    assert(ctools_parse_thread_value(" \t", 0) == 0);
    assert(ctools_parse_thread_value("8", 0) == 8);
    assert(ctools_parse_thread_value(" 8 \r\n", 0) == 8);
    assert(ctools_parse_thread_value("0", 0) == -1);
    assert(ctools_parse_thread_value("-2", 0) == -1);
    assert(ctools_parse_thread_value("+2", 0) == -1);
    assert(ctools_parse_thread_value("8x", 0) == -1);
    assert(ctools_parse_thread_value("8 9", 0) == -1);
    assert(ctools_parse_thread_value("8,4", 0) == -1);
    assert(ctools_parse_thread_value("8,4", 1) == 8);
    assert(ctools_parse_thread_value("99999999999999999999", 0) == INT_MAX);
}

static void precedence(void)
{
    ctools_thread_policy p;
    p = resolve(200, 0, 0, 0, 0);
    assert(p.threads == CTOOLS_DEFAULT_THREAD_CAP && p.limit == 200);
    assert(strcmp(p.source, "automatic") == 0);
    p = resolve(8, 0, 0, 0, 0);
    assert(p.threads == 8 && strcmp(p.source, "automatic") == 0);
    p = resolve(200, 0, 0, 0, 64);
    assert(p.threads == 64 && strcmp(p.source, "OMP_NUM_THREADS") == 0);
    p = resolve(200, 0, 0, 32, 64);
    assert(p.threads == 32 && strcmp(p.source, "saved setting") == 0);
    p = resolve(200, 0, 48, 32, 64);
    assert(p.threads == 48 && strcmp(p.source, "CTOOLS_NUM_THREADS") == 0);
    /* Every source is clamped to usable CPUs, OMP_THREAD_LIMIT and 256. */
    p = resolve(4, 0, 48, 0, 0);
    assert(p.threads == 4 && strcmp(p.source, "CTOOLS_NUM_THREADS") == 0);
    p = resolve(200, 2, 0, 32, 0);
    assert(p.limit == 2 && p.threads == 2);
    p = resolve(1000, 0, 999, 0, 0);
    assert(p.limit == CTOOLS_MAX_THREADS && p.threads == CTOOLS_MAX_THREADS);
    /* Invalid inputs are ignored. */
    p = resolve(200, -1, -1, -1, -1);
    assert(p.limit == 200 && p.threads == CTOOLS_DEFAULT_THREAD_CAP);
    assert(strcmp(p.source, "automatic") == 0);
    p = resolve(0, 0, 0, 0, 0);
    assert(p.limit == 1 && p.threads == 1);
}

#if defined(__linux__)
static void cgroups(const char *tmp)
{
    char root[4096], path[4200];
    snprintf(root, sizeof(root), "%s/cg", tmp);
    assert(mkdir(root, 0700) == 0);
    snprintf(path, sizeof(path), "%s/a", root);
    assert(mkdir(path, 0700) == 0);
    snprintf(path, sizeof(path), "%s/a/b", root);
    assert(mkdir(path, 0700) == 0);

    /* v2: the tightest quota on the path wins, rounded up to whole CPUs. */
    assert(cgroup_walk(root, "/a/b", 1) == 0);
    snprintf(path, sizeof(path), "%s/a/b/cpu.max", root);
    put(path, "max 100000\n");
    snprintf(path, sizeof(path), "%s/a/cpu.max", root);
    put(path, "400000 100000\n");
    snprintf(path, sizeof(path), "%s/cpu.max", root);
    put(path, "150000 100000\n");
    assert(cgroup_walk(root, "/a/b", 1) == 2);
    assert(cgroup_walk(root, "/a/b/", 1) == 2);
    put(path, "max 100000\n");
    assert(cgroup_walk(root, "/a/b", 1) == 4);
    assert(cgroup_walk(root, "/", 1) == 0);
    assert(cgroup_walk(root, "/../a", 1) == 0);  /* never leaves the root */
    assert(cgroup_walk(root, "/missing/x", 1) == 0);
    snprintf(path, sizeof(path), "%s/a/cpu.max", root);
    put(path, "garbage\n");
    assert(cgroup_walk(root, "/a/b", 1) == 0);

    /* v1: CFS quota over period; -1 means unlimited. */
    snprintf(path, sizeof(path), "%s/cpu.cfs_quota_us", root);
    put(path, "-1\n");
    snprintf(path, sizeof(path), "%s/cpu.cfs_period_us", root);
    put(path, "100000\n");
    snprintf(path, sizeof(path), "%s/a/cpu.cfs_quota_us", root);
    put(path, "250000\n");
    snprintf(path, sizeof(path), "%s/a/cpu.cfs_period_us", root);
    put(path, "100000\n");
    assert(cgroup_walk(root, "/a/b", 0) == 3);
    assert(cgroup_walk(root, "/", 0) == 0);
}
#endif

static void settings(const char *expected)
{
    const ctools_thread_policy *p = ctools_thread_policy_get();
    char text[1024];
    const char *path = p->settings_path;

    assert(strcmp(path, expected) == 0);
    assert(p->saved == 0);
    ctools_init_thread_count();

    /* Save creates missing directories and becomes this process's default. */
    assert(ctools_save_thread_setting(3) == 0);
    assert(p->saved == 3 && read_thread_setting(path) == 3);
    assert(strcmp(p->source, "saved setting") == 0 || p->env_ctools > 0);
    assert(ctools_get_max_threads() == p->threads);
    get(path, text, sizeof(text));
    assert(strcmp(text, CTOOLS_SETTINGS_HEADER "threads=3\n") == 0);

    /* Unknown keys and comments survive a rewrite; the last threads= wins. */
    put(path, "future = 1\n# note\nthreads = 6\nthreads=7\n");
    assert(read_thread_setting(path) == 7);
    assert(ctools_save_thread_setting(12) == 0);
    get(path, text, sizeof(text));
    assert(strcmp(text, CTOOLS_SETTINGS_HEADER "future = 1\n# note\nthreads=12\n") == 0);
    assert(ctools_save_thread_setting(0) == 0);
    get(path, text, sizeof(text));
    assert(strcmp(text, CTOOLS_SETTINGS_HEADER "future = 1\n# note\n") == 0);
    assert(read_thread_setting(path) == 0 && p->saved == 0);

    /* Hand-edited files: BOM and CRLF are accepted, bad values are invalid. */
    put(path, "\xEF\xBB\xBFthreads=5\r\n");
    assert(read_thread_setting(path) == 5);
    put(path, "threads=abc\n");
    assert(read_thread_setting(path) == -1);
    put(path, "threadsx=4\nthreads 4\n");
    assert(read_thread_setting(path) == 0);

    /* Removing the only setting removes the file and its now-empty
     * directory; removing twice is fine, and saving recreates both. */
    put(path, CTOOLS_SETTINGS_HEADER "threads=5\n");
    assert(ctools_save_thread_setting(0) == 0);
    assert(ctools_fopen(path, "r") == NULL && !parent_exists(path));
    assert(ctools_save_thread_setting(0) == 0);
    assert(ctools_save_thread_setting(-4) == 0);
    assert(ctools_save_thread_setting(100000) == 0 && read_thread_setting(path) == CTOOLS_MAX_THREADS);
    assert(ctools_save_thread_setting(0) == 0);
}

int main(int argc, char **argv)
{
    if (argc == 4 && strcmp(argv[1], "unit") == 0) {
        parsing();
        precedence();
#if defined(__linux__)
        cgroups(argv[2]);
#endif
        settings(argv[3]);
        puts("unit ok");
        return 0;
    }
    if (argc == 2 && strcmp(argv[1], "report") == 0) {
        const ctools_thread_policy *p = ctools_thread_policy_get();
        ctools_init_thread_count();
        ctools_reset_max_threads();
        printf("threads=%d source=%s limit=%d detected=%d available=%d saved=%d "
               "env=%d omp=%d omp_limit=%d max=%d clamp=%d\n",
               p->threads, p->source, p->limit, p->detected, p->available, p->saved,
               p->env_ctools, p->env_omp, p->omp_limit, ctools_get_max_threads(),
               ctools_clamp_thread_request(100000));
        return 0;
    }
    return 2;
}
#endif

#if defined(CTOOLS_CASE_p2_native_3)
#include "ctools_config.h"
#undef MIN_OBS_PER_THREAD
#define MIN_OBS_PER_THREAD 256
#include "ctools_data_io.c"
#include "ctools_transport.c"
#include "ctools_data_plan.c"
#include "ctools_data_store.c"
int main(void) {
    setup();
    int indices[] = {1, 2}, widths[] = {8, 8};
    ctools_filtered_data fd;
    for (int pool = 0; pool <= 1; pool++) {
        ctools_set_max_threads(pool ? 4 : 1);
        for (int columns = 1; columns <= 2; columns++)
        for (int strings = 0; strings <= 1; strings++)
        for (int hints = 0; hints <= 1; hints++)
        for (filter = 0; filter <= 1; filter++) {
            string_var = strings ? 1 : 0;
            bad_row = 3;
            assert(ctools_data_load_ex(&fd, indices, columns, 1, nrows, 0, hints ? widths : NULL) == STATA_ERR_STATA_READ);
            assert(fd.data.vars == NULL && fd.obs_map == NULL);
            bad_row = 0;
            assert(ctools_data_load_ex(&fd, indices, columns, 1, nrows, 0, hints ? widths : NULL) == 0);
            bad_row = 3;
            assert(ctools_data_store_ex(&fd.data, indices, columns, 1) == STATA_ERR_STATA_WRITE);
            if (!strings) {
                assert(ctools_store_filtered(fd.data.vars[0].data.dbl, fd.data.nobs, 1, fd.obs_map) == STATA_ERR_STATA_WRITE);
                assert(ctools_store_filtered_rowpar(fd.data.vars[0].data.dbl, fd.data.nobs, 1, fd.obs_map) == STATA_ERR_STATA_WRITE);
            } else {
                assert(ctools_store_filtered_str(fd.data.vars[0].data.str, fd.data.nobs, 1, fd.obs_map) == STATA_ERR_STATA_WRITE);
            }
            ctools_filtered_data_free(&fd);
        }
        ctools_destroy_global_pool();
    }
    bad_row = filter = 0; string_var = 1; strl = 1; reads = 0;
    assert(ctools_data_load(&fd, indices, 1, 1, 10, 0) == 0 && reads == 10);
    ctools_filtered_data_free(&fd);
    reads = 0; strl_length = 2046;
    assert(ctools_data_load(&fd, indices, 1, 1, 10, 0) == STATA_ERR_UNSUPPORTED_TYPE && reads == 0);
    strl_length = 4; binary_value = 1;
    assert(ctools_data_load(&fd, indices, 1, 1, 10, 0) == STATA_ERR_UNSUPPORTED_TYPE && reads == 0);
    binary_value = 0; strl = 0; widths[0] = 1;
    assert(ctools_data_load_ex(&fd, indices, 1, 1, 10, 0, widths) == STATA_ERR_STATA_READ);
    indices[0] = 3;
    assert(ctools_data_load(&fd, indices, 1, 1, 10, 0) == STATA_ERR_INVALID_INPUT);
    indices[0] = 1;
    assert(ctools_data_load(&fd, indices, 1, 10, 1, 0) == STATA_ERR_INVALID_INPUT);
    assert(ctools_data_load(&fd, indices, 1, 1, nrows+1, 0) == STATA_ERR_INVALID_INPUT);
    double values[] = {1, 2}; perm_idx_t map[] = {1, 6001}; writes = 0; string_var = 0;
    assert(ctools_store_filtered_rowpar(values, 2, 1, map) == STATA_ERR_INVALID_INPUT && writes == 0);
    int64_t rows[] = {0, 2, 1}; bad_row = 3; writes = 0;
    assert(ctools_stream_var_permuted(1, rows, 3, 1) == STATA_ERR_STATA_READ && writes == 0);
    string_var = 1; writes = 0;
    assert(ctools_stream_var_permuted(1, rows, 3, 1) == STATA_ERR_STATA_READ && writes == 0);
    strl = 1; reads = 0;
    assert(ctools_stream_var_permuted(1, rows, 3, 1) == STATA_ERR_UNSUPPORTED_TYPE && reads == 0 && writes == 0);
    ctools_destroy_global_pool();
    return 0;
}
#endif

#if defined(CTOOLS_CASE_p2_native_6)
#include "cbinscatter/cbinscatter_impl.c"
static ST_retcode scalars(char *s, ST_double d) { (void)s; (void)d; return 0; }
static ST_retcode matrix(char *s, ST_int row, ST_int col, ST_double d) {
    (void)s; (void)row; (void)d; return col == 4 ? 503 : 0;
}
int main(void) {
    setup(); mock.scalsave = scalars; mock.safematstore = matrix;
    BinStats bin = {.n_obs = 3};
    ByGroupResult group = {.num_bins = 1, .bins = &bin};
    BinscatterResults results = {.num_by_groups = 1, .groups = &group};
    BinscatterConfig config = {0};
    assert(store_results(&results, &config) == 503);
    return 0;
}
#endif

#if defined(CTOOLS_CASE_sep24_native_pool)
#include <assert.h>
#include <errno.h>
#include <pthread.h>
#include "ctools_threads.h"
static int calls, fail_at;
static int create(pthread_t *t,const pthread_attr_t *a,void *(*f)(void*),void *p) {
 if(++calls==fail_at)return EAGAIN;return pthread_create(t,a,f,p);
}
#define pthread_create create
#include "ctools_threads.c"
int main(void) {
 for(int run=0;run<100;run++)for(int f=1;f<=4;f++) {
  ctools_persistent_pool pool;calls=0;fail_at=f;
  assert(ctools_persistent_pool_init(&pool,4)!=0);
  assert(!pool.initialized && !pool.workers && pool.num_workers==0);
  fail_at=0;assert(ctools_persistent_pool_init(&pool,2)==0);
  ctools_persistent_pool_destroy(&pool);
 }
 return 0;
}
#endif

#if defined(CTOOLS_CASE_memory_threads)
#include "ctools_threads.c"
#include <assert.h>
#include <limits.h>
int main(void) {
    ctools_set_max_threads(INT_MAX);
    assert(ctools_get_max_threads() == CTOOLS_MAX_THREADS);
    assert(ctools_workspace_threads(64, 1, 80000000) == 1);
    assert(ctools_workspace_threads(64, 100, 80000000) == 3);
    assert(ctools_workspace_threads(64, 100, 0) == 64);
    ctools_set_max_threads(2);
    assert(ctools_get_max_threads() == 2);
    ctools_reset_max_threads();
    assert(ctools_get_max_threads() > 0);
    return 0;
}
#endif
