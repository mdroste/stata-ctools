/* Commands tests. Select one CTOOLS_CASE_* per compilation.
 * Separate selectors preserve allocator/SPI interposition and independent main functions.
 */

#if defined(CTOOLS_CASE_cpsmatch_paths)
/* Independent scan/enumeration oracles for the optimized matching paths.
 * The tracked allocator is prepended by the Python runner. */
#include "cpsmatch/cpsmatch_internal.h"
#include "cpsmatch/cpsmatch_index.h"
#include "cpsmatch/cpsmatch_kernels.c"
#include "cpsmatch/cpsmatch_radius.c"

static ST_plugin host;
ST_plugin *_stata_ = &host;
int ctools_get_max_threads(void)
{
    return 4;
}

static void availability(void)
{
    enum
    {
        N = 1024
    };
    perm_idx_t previous[N + 1], next[N + 1];
    int live[N];
    for (size_t i = 0; i <= N; i++)
        previous[i] = next[i] = (perm_idx_t)i;
    for (size_t i = 0; i < N; i++)
        live[i] = 1;
    for (size_t step = 0; step <= N; step++)
    {
        for (size_t q = 0; q <= N; q++)
        {
            size_t r = q;
            while (r < N && !live[r])
                r++;
            size_t l = q;
            while (l > 0 && !live[l - 1])
                l--;
            assert(cpsmatch_available_find(next, (perm_idx_t)q) == r);
            assert(cpsmatch_available_find(previous, (perm_idx_t)q) == l);
        }
        if (step == N)
            break;
        size_t deleted = step * 733 % N;
        live[deleted] = 0;
        previous[deleted + 1] = cpsmatch_available_find(previous, (perm_idx_t)deleted);
        next[deleted] = cpsmatch_available_find(next, (perm_idx_t)(deleted + 1));
    }
}

static void radius(int shape, int inject)
{
    enum
    {
        NC = 513,
        NT = 777,
        N = NC + NT
    };
    sorted_control_t controls[NC];
    double ps[NT], weights[N] = {0}, ids[N], nn[N] = {0}, support[N];
    double expected[N] = {0};
    size_t ci[NC], ti[NT];
    perm_idx_t obs[N];
    for (size_t i = 0; i < N; i++)
    {
        obs[i] = (perm_idx_t)(3 * i + 1);
        ids[i] = support[i] = host.missval;
    }
    for (size_t j = 0; j < NC; j++)
    {
        ci[j] = j;
        controls[j] = (sorted_control_t){0};
        controls[j].orig_idx = NC - 1 - j;
        controls[j].pscore = shape == 0 ? .5 : (double)(j / 3) / 171;
    }
    for (size_t t = 0; t < NT; t++)
    {
        ti[t] = NC + t;
        ps[t] = shape == 0 ? .5 : shape == 4 ? 2.0 : (double)((t * 173) % 1000) / 1000;
    }
    double width = shape == 2 ? 1e-12 : shape == 3 ? .03 : .25;
    int common = shape == 3;
    fail_at = inject;
    calls = 0;
    int matched = 0, off = 0;
    int rc = cpsmatch_match_radius(controls, NC, ps, ti, NT, ci, obs, width, common, .2, .8,
                                   weights, ids, nn, support, &matched, &off);
    if (inject)
    {
        assert(rc == 920 && matched == 0 && off == 0);
        for (size_t i = 0; i < N; i++)
            assert(weights[i] == 0 && ids[i] == host.missval);
        assert(alive == 0);
        return;
    }
    assert(rc == 0);
    int expected_matched = 0;
    for (size_t t = 0; t < NT; t++)
    {
        size_t start = NC, count = 0;
        if (!common || (ps[t] >= .2 && ps[t] <= .8))
        {
            for (size_t j = 0; j < NC; j++)
                if (controls[j].pscore >= ps[t] - width && controls[j].pscore <= ps[t] + width)
                {
                    if (!count)
                        start = j;
                    count++;
                }
        }
        assert(nn[ti[t]] == (double)count);
        assert(support[ti[t]] == (count ? 1 : 0));
        if (!count)
        {
            assert(ids[ti[t]] == host.missval);
            continue;
        }
        expected_matched++;
        assert(ids[ti[t]] == obs[ci[controls[start].orig_idx]]);
        expected[ti[t]] = 1;
        for (size_t j = start; j < start + count; j++)
            expected[ci[controls[j].orig_idx]] += 1.0 / count;
    }
    assert(matched == expected_matched && off == NT - expected_matched);
    for (size_t i = 0; i < N; i++)
    {
        if (expected[i] == 0)
            assert(weights[i] == 0);
        else
            assert(weights[i] > 0 && fabs(weights[i] - expected[i]) / expected[i] < 1e-7);
    }
    assert(alive == 0);
}

int main(void)
{
    host.missval = 8.98846567431158e307;
    int counts[] = {1, 2, 3, 17, 1023, 2147483647};
    for (size_t i = 0; i < sizeof(counts) / sizeof(counts[0]); i++)
    {
        double weight = 1.0 / counts[i];
        cpsmatch_radius_sum small = cpsmatch_radius_encode(weight);
        assert(cpsmatch_radius_value(small) == weight);
        cpsmatch_radius_sum big = cpsmatch_radius_encode(1.0);
        for (int bit = 0; bit < 31; bit++)
            big = cpsmatch_radius_add(big, big);
        cpsmatch_radius_sum recovered = cpsmatch_radius_sub(cpsmatch_radius_add(big, small), big);
        assert(recovered.hi == small.hi && recovered.lo == small.lo);
        recovered = cpsmatch_radius_add(cpsmatch_radius_sub(small, big), big);
        assert(recovered.hi == small.hi && recovered.lo == small.lo);
    }
    availability();
    for (int shape = 0; shape < 5; shape++)
        for (int injection = 0; injection <= 2; injection++)
            radius(shape, injection);
    puts("CPSMATCH_PATHS_NATIVE_COMPLETE");
    return 0;
}
#endif

#if defined(CTOOLS_CASE_cpsmatch_paths_oom)
#include "cpsmatch/cpsmatch_internal.h"
#include "cpsmatch/cpsmatch_kernels.c"
#include "cpsmatch/cpsmatch_options.c"
#include "cpsmatch/cpsmatch_stats.c"
#include "cpsmatch/cpsmatch_radius.c"
#include "cpsmatch/cpsmatch_impl.c"

static ST_plugin host;
ST_plugin *_stata_ = &host;
enum
{
    N = 192
};
static double treatment[N], score[N];
static perm_idx_t observations[N];
static stata_variable columns[2];
static int stores;
static ST_int first(void)
{
    return 1;
}
static ST_int last(void)
{
    return N;
}
static ST_boolean missing(double x)
{
    return x >= host.missval;
}
static ST_retcode store(ST_int v, ST_int row, double value)
{
    (void)v;
    (void)row;
    (void)value;
    stores++;
    return 0;
}
static ST_retcode scalar(char *name, double value)
{
    (void)name;
    (void)value;
    return 0;
}
static ST_retcode message(char *text)
{
    (void)text;
    return 0;
}
double ctools_timer_seconds(void)
{
    return 0;
}
int ctools_get_max_threads(void)
{
    return 4;
}
void ctools_filtered_data_init(ctools_filtered_data *d)
{
    memset(d, 0, sizeof(*d));
}
void ctools_filtered_data_free(ctools_filtered_data *d)
{
    (void)d;
}
stata_retcode ctools_data_load(ctools_filtered_data *d, int *v, size_t nv, size_t start, size_t end,
                               int flags)
{
    (void)v;
    (void)nv;
    (void)start;
    (void)end;
    (void)flags;
    columns[0].data.dbl = treatment;
    columns[1].data.dbl = score;
    d->data.nobs = N;
    d->data.nvars = 2;
    d->data.vars = columns;
    d->obs_map = observations;
    return STATA_OK;
}
int main(void)
{
    host.missval = 8.98846567431158e307;
    host.nobs1 = first;
    host.nobs2 = last;
    host.ismissing = missing;
    host.safestore = store;
    host.store = store;
    host.scalsave = scalar;
    host.spouterr = message;
    for (int i = 0; i < N; i++)
    {
        observations[i] = (perm_idx_t)(i + 1);
        treatment[i] = i % 2;
        score[i] = (double)(i % 17) / 17;
    }
    const char *args[] = {"1 2 0 3 4 5 192 0 1 0 0 0 .06 0 0 0",
                          "1 2 0 3 4 5 192 1 1 .1 0 0 .06 1 0 0"};
    for (int method = 0; method < 2; method++)
    {
        calls = fail_at = stores = 0;
        assert(cpsmatch_main(args[method]) == 0 && stores > 0 && alive == 0);
        int total = calls;
        for (int fail = 1; fail <= total; fail++)
        {
            calls = stores = 0;
            fail_at = fail;
            assert(cpsmatch_main(args[method]) == 920);
            assert(stores == 0 && alive == 0);
        }
        printf("path=%d allocation failures checked=%d\n", method, total);
    }
    puts("CPSMATCH_PATHS_OOM_COMPLETE");
}
#endif

#if defined(CTOOLS_CASE_encode_label_file)
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "ctools_hash.c"

int main(void)
{
    char path[] = "/tmp/ctools-labels-XXXXXX";
    int fd = mkstemp(path);
    assert(fd >= 0);
    close(fd);
    char unicode_path[256];
    snprintf(unicode_path,sizeof(unicode_path),"%s-é東京",path);
    assert(ctools_label_write_file(NULL,NULL,0,unicode_path)==0);
    FILE *unicode_file=ctools_fopen(unicode_path,"rb");
    assert(unicode_file);
    fclose(unicode_file);
    ctools_remove(unicode_path);
    char longest[2046];
    memset(longest, 'x', 2045);
    longest[2045] = '\0';
    const char *texts[] = {"", "\r\n\t", "`macro' $global \"quoted\"", "é東京", longest};
    const int codes[] = {-10, 0, 1, 7, 65536};
    assert(ctools_label_write_file(texts, codes, 5, path) == 0);
    FILE *fp = fopen(path, "rb");
    assert(fp);
    char magic[32], actual[2046];
    assert(fgets(magic, sizeof(magic), fp));
    assert(strcmp(magic, "CTOOLS_LABELS_1\n") == 0);
    for (size_t i = 0; i < 5; i++)
    {
        int code;
        size_t length;
        assert(fscanf(fp, "%d %zu", &code, &length) == 2);
        assert(code == codes[i] && length == strlen(texts[i]));
        assert(fgetc(fp) == '\n');
        assert(fread(actual, 1, length, fp) == length);
        assert(memcmp(actual, texts[i], length) == 0);
    }
    assert(fgetc(fp) == EOF && !ferror(fp));
    assert(fclose(fp) == 0);
    unlink(path);
    assert(ctools_label_write_file(texts, codes, 5, "/nonexistent/ctools-labels") == -1);
#ifdef __linux__
    assert(ctools_label_write_file(texts, codes, 5, "/dev/full") == -1);
#endif
    puts("PASS label file framing and failure handling");
    return 0;
}
#endif

#if defined(CTOOLS_CASE_group_moments)
#include "crangestat/crangestat_internal.h"
#include <assert.h>

static ST_plugin host;
ST_plugin *_stata_ = &host;

int main(void)
{
    host.missval = 8.98846567431158e307;
    double x[] = {1.0, 2.0, 3.0, host.missval};
    double constant[] = {2.0, 2.0};
    for (int stat = STAT_SKEWNESS; stat <= STAT_KURTOSIS; stat++)
    {
        assert(crangestat_compute_window_stat_simple(stat, x, 0, 0, 0, 0) == host.missval);
        assert(crangestat_compute_window_stat_simple(stat, x, 0, 1, 0, 0) == host.missval);
        assert(crangestat_compute_window_stat_simple(stat, constant, 0, 2, 0, 0) == host.missval);
        double pair = stat == STAT_SKEWNESS ? 0.0 : 1.0;
        assert(crangestat_compute_window_stat_simple(stat, x, 0, 2, 0, 0) == pair);
        assert(crangestat_compute_window_stat_simple(stat, x, 0, 4, 1, 1) == pair);
        double triple = stat == STAT_SKEWNESS ? 0.0 : 1.5;
        assert(crangestat_compute_window_stat_simple(stat, x, 0, 4, 0, 0) == triple);
    }
    return 0;
}
#endif

#if defined(CTOOLS_CASE_match_winsor_optimizations)
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "ctools_types.h"
#include "cpsmatch/cpsmatch_stats.c"

static int fail_boundary;
static void *boundary_malloc(size_t n)
{
    return fail_boundary ? NULL : malloc(n);
}
#define malloc boundary_malloc
#include "ctools_order.c"
#undef malloc

static ST_plugin host;
ST_plugin *_stata_ = &host;
int ctools_get_max_threads(void)
{
    return 4;
}
static int scalar_failure, scalar_calls;
static ST_retcode scalar_save(char *name, double value)
{
    (void)name;
    (void)value;
    return ++scalar_calls == scalar_failure ? 459 : 0;
}

static void statistics(void)
{
    double treatment[] = {0, 0, 0, 1, 1, 1};
    double y[] = {1, 2, 4, 3, 6, 9};
    double weight[] = {2, 1, .5, 1, 1, 1};
    double support[] = {1, 1, 1, 1, 1, 0};
    cpsmatch_statistics s;
    cpsmatch_compute_statistics(6, treatment, y, weight, support, &s);
    assert(s.y1 == 4.5 && s.n1 == 2 && s.var1 == 4.5);
    assert(fabs(s.y0 - 12.0 / 7) < 1e-15 && s.n0 == 3);
    assert(fabs(s.var0 - 7.0 / 3) < 1e-15 && s.weight_squared == 5.25);
    assert(s.y1_u == 6 && s.var1_u == 9);
    assert(s.n0_on == 3 && s.n1_on == 2 && s.n1_off == 1 && s.n0_off == 0);
    for (scalar_failure = 1; scalar_failure <= 17; scalar_failure++)
    {
        scalar_calls = 0;
        assert(cpsmatch_save_statistics(&s) == 459);
        assert(scalar_calls == scalar_failure);
    }
    scalar_failure = 0;
    scalar_calls = 0;
    assert(cpsmatch_save_statistics(&s) == 0 && scalar_calls == 17);
    cpsmatch_compute_statistics(6, treatment, NULL, weight, support, &s);
    assert(s.y1 == host.missval && s.var1 == host.missval && s.n1 == 0);
    assert(s.n1_on == 2 && s.n1_off == 1);
    cpsmatch_compute_statistics(0, NULL, NULL, NULL, NULL, &s);
    assert(s.y0 == host.missval && s.y1 == host.missval && s.n0_on == 0);
    cpsmatch_compute_statistics(1, treatment, y, weight, support, &s);
    assert(s.var0 == host.missval && s.var0_u == host.missval);
}

static void grouping(void)
{
    const size_t n = 70001;
    stata_variable columns[2] = {0};
    stata_data d = {0};
    d.nobs = n;
    d.nvars = 2;
    d.vars = columns;
    columns[0].type = STATA_TYPE_STRING;
    columns[1].type = STATA_TYPE_DOUBLE;
    char **strings = malloc(n * sizeof(*strings));
    double *numbers = malloc(n * sizeof(*numbers));
    perm_idx_t *order = malloc(n * sizeof(*order));
    size_t *expected = malloc((n + 1) * sizeof(*expected)),
           *actual = malloc((n + 1) * sizeof(*actual));
    assert(strings && numbers && order && expected && actual);
    columns[0].data.str = strings;
    columns[1].data.dbl = numbers;
    int keys[] = {0, 1};
    for (int shape = 0; shape < 4; shape++)
    {
        for (size_t i = 0; i < n; i++)
        {
            strings[i] = shape == 0 ? "same" : i < 32768 ? "prefix_a" : "prefix_b";
            numbers[i] = shape < 2 ? 0 : shape == 2 ? (double)(i / 31) : (double)i;
            if (shape == 2 && i >= 65536)
                numbers[i] = host.missval;
            order[i] = (perm_idx_t)(n - i - 1);
        }
        for (int indirect = 0; indirect < 2; indirect++)
        {
            const perm_idx_t *p = indirect ? order : NULL;
            size_t groups = ctools_group_starts(&d, p, keys, 2, expected);
            for (fail_boundary = 0; fail_boundary < 2; fail_boundary++)
            {
                size_t got = ctools_group_starts_parallel(&d, p, keys, 2, actual);
                assert(got == groups && !memcmp(expected, actual, (groups + 1) * sizeof(*actual)));
            }
        }
    }
    fail_boundary = 0;
    free(strings);
    free(numbers);
    free(order);
    free(expected);
    free(actual);
}

int main(void)
{
    host.missval = 8.98846567431158e307;
    host.scalsave = scalar_save;
    statistics();
    grouping();
    puts("MATCH_WINSOR_OPTIMIZATIONS_NATIVE_COMPLETE");
    return 0;
}
#endif

#if defined(CTOOLS_CASE_newcommands_native_source)
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "stplugin.h"
#include "ctools_types.h"
static ST_plugin mock;
ST_plugin *_stata_ = &mock;
static int nr, nv, types[16], calls, fail_at, alive, write_error, load_error;
static double numbers[16][16];
static char strings[16][16][2046], widths[16000], nout[64], nmissing[64];
static void *owned[1000];
static void track(void *p) { if(p) {for(int i=0;i<1000;i++)if(!owned[i]){owned[i]=p;alive++;return;}abort();} }
static void *test_malloc(size_t n) {if(++calls==fail_at)return NULL;void *p=malloc(n);track(p);return p;}
static void *test_calloc(size_t n,size_t s) {if(++calls==fail_at)return NULL;void *p=calloc(n,s);track(p);return p;}
static void test_free(void *p) {if(!p)return;for(int i=0;i<1000;i++)if(owned[i]==p){owned[i]=NULL;alive--;free(p);return;}assert(!"unowned free");}
static ST_int nobs(void){return nr;}
static ST_int nvars(void){return nv;}
static ST_boolean isstr(ST_int v){assert(v>=1&&v<=nv);return types[v-1];}
static ST_boolean isstrl(ST_int v){(void)v;return 0;}
static ST_retcode putnum(ST_int v,ST_int r,ST_double x){assert(v>=1&&v<=nv&&r>=1&&r<=nr);if(write_error)return 459;numbers[v-1][r-1]=x;return 0;}
static ST_retcode putstr(ST_int v,ST_int r,char *x){assert(v>=1&&v<=nv&&r>=1&&r<=nr);assert(strlen(x)<=2045);if(write_error)return 459;strcpy(strings[v-1][r-1],x);return 0;}
static ST_retcode macro_save(char *key,char *value){if(!strcmp(key,"_csplit_widths"))strcpy(widths,value);else if(!strcmp(key,"_crangejoin_nout"))strcpy(nout,value);else if(!strcmp(key,"_cinterpolate_nmissing"))strcpy(nmissing,value);else abort();return 0;}
static ST_retcode macro_get(char *key,char *value,ST_int size){const char *s=!strcmp(key,"_csplit_delim1")?"::":":";assert(size>=(int)strlen(s)+1);strcpy(value,s);return 0;}
int ctools_get_max_threads(void){return 1;}
double ctools_timer_seconds(void){return 0;}
void ctools_verbose(const char *m,int v,const char *f,...){(void)m;(void)v;(void)f;}
void ctools_error(const char *m,const char *f,...){(void)m;(void)f;}
#define malloc test_malloc
#define calloc test_calloc
#define free test_free
void ctools_filtered_data_init(ctools_filtered_data *fd){memset(fd,0,sizeof(*fd));}
void ctools_filtered_data_free(ctools_filtered_data *fd){
 for(size_t k=0;k<fd->data.nvars;k++){
  stata_variable *v=fd->data.vars+k;
  if(v->type==STATA_TYPE_STRING){if(v->data.str)for(size_t i=0;i<v->nobs;i++)free(v->data.str[i]);free(v->data.str);}
  else free(v->data.dbl);
 }
 free(fd->data.vars);free(fd->data.sort_order);free(fd->obs_map);memset(fd,0,sizeof(*fd));
}
stata_retcode ctools_data_load(ctools_filtered_data *fd,int *vars,size_t n,size_t start,size_t end,int flags){
 (void)start;(void)end;(void)flags;if(load_error)return load_error;if(!vars)n=nv;
 fd->data.vars=calloc(n,sizeof(stata_variable));if(!fd->data.vars)return STATA_ERR_MEMORY;
 fd->data.nvars=n;fd->data.nobs=nr;
 fd->data.sort_order=malloc(nr*sizeof(perm_idx_t));fd->obs_map=malloc(nr*sizeof(perm_idx_t));
 if(!fd->data.sort_order||!fd->obs_map)return STATA_ERR_MEMORY;
 for(int i=0;i<nr;i++){fd->data.sort_order[i]=i;fd->obs_map[i]=i+1;}
 for(size_t k=0;k<n;k++){
  int source=vars?vars[k]-1:(int)k;stata_variable *v=fd->data.vars+k;v->nobs=nr;v->type=types[source]?STATA_TYPE_STRING:STATA_TYPE_DOUBLE;
  if(types[source]){
   v->data.str=calloc(nr,sizeof(char*));if(!v->data.str)return STATA_ERR_MEMORY;
   for(int i=0;i<nr;i++){size_t len=strlen(strings[source][i])+1;v->data.str[i]=malloc(len);if(!v->data.str[i])return STATA_ERR_MEMORY;memcpy(v->data.str[i],strings[source][i],len);}
  }else{v->data.dbl=malloc(nr*sizeof(double));if(!v->data.dbl)return STATA_ERR_MEMORY;memcpy(v->data.dbl,numbers[source],nr*sizeof(double));}
 }
 return STATA_OK;
}
stata_retcode ctools_store_filtered(double *v,size_t n,int col,perm_idx_t *map){for(size_t i=0;i<n;i++)if(putnum(col,map[i],v[i]))return STATA_ERR_STATA_WRITE;return STATA_OK;}
#include "ctools_order.c"
#include "cinterpolate/cinterpolate_impl.c"
#define state split_state
#include "csplit/csplit_impl.c"
#undef state
#define state join_state
#include "crangejoin/crangejoin_impl.c"
#undef state

static void reset(void){assert(alive==0);memset(numbers,0,sizeof(numbers));memset(strings,0,sizeof(strings));memset(types,0,sizeof(types));nr=nv=0;calls=0;write_error=0;}
static int interpolation(void){
 reset();nr=4;nv=3;numbers[0][0]=numbers[0][2]=mock.missval;numbers[0][1]=10;numbers[0][3]=30;
 for(int i=0;i<nr;i++)numbers[1][i]=i;
 int rc=cinterpolate_main("0 1 0");if(!rc){for(int i=0;i<nr;i++)assert(numbers[2][i]==10*i);assert(!strcmp(nmissing,"0"));}
 assert(alive==0);return rc;
}
static int cubic_interpolation(void){
 reset();nr=9;nv=3;
 for(int i=0;i<nr;i++){numbers[1][i]=i;numbers[0][i]=(i%2)?mock.missval:i*i*i;}
 int rc=cinterpolate_main("0 0 0 cubic");
 if(!rc){
  for(int i=0;i<nr;i++)assert(numbers[2][i]==((i==1||i==7)?mock.missval:i*i*i));
  assert(!strcmp(nmissing,"2"));
 }
 assert(alive==0);return rc;
}
static int splitting(void){
 reset();nr=2;nv=1;types[0]=1;strcpy(strings[0][0],"a::b:c");strcpy(strings[0][1],"é:東京");
 int rc=csplit_main("scan 2 2046 1 0");
 if(!rc){assert(!strcmp(widths,"2 6 1 1 "));nv=4;for(int k=0;k<nv;k++)types[k]=1;
  rc=csplit_main("write");if(!rc){assert(!strcmp(strings[0][0],"a"));assert(!strcmp(strings[1][0],""));assert(!strcmp(strings[2][0],"b"));assert(!strcmp(strings[3][0],"c"));}}
 csplit_cleanup_cache();assert(alive==0);return rc;
}
static int joining(void){
 reset();nr=3;nv=2;numbers[0][0]=3;numbers[0][1]=1;numbers[0][2]=1;types[1]=1;
 strcpy(strings[1][0],"three");strcpy(strings[1][1],"first");strcpy(strings[1][2],"second");
 int rc=crangejoin_main("using 0 1 0");
 if(!rc){nr=2;nv=5;memset(types,0,sizeof(types));types[4]=1;numbers[0][0]=100;numbers[0][1]=200;numbers[1][0]=1;numbers[1][1]=9;numbers[2][0]=3;numbers[2][1]=10;
  rc=crangejoin_main("prepare 1");
  if(!rc){assert(!strcmp(nout,"4"));nr=4;rc=crangejoin_main("write");
   if(!rc){assert(numbers[0][0]==100&&numbers[0][2]==100&&numbers[0][3]==200);assert(numbers[3][0]==1&&numbers[3][2]==3&&numbers[3][3]==mock.missval);assert(!strcmp(strings[4][0],"first")&&!strcmp(strings[4][1],"second")&&!strcmp(strings[4][3],""));}}}
 crangejoin_cleanup_cache();assert(alive==0);return rc;
}
#ifndef CTOOLS_HARNESS_MAIN
#define CTOOLS_HARNESS_MAIN main
#endif
int CTOOLS_HARNESS_MAIN(void){
 mock.nobs=nobs;mock.nvars=nvars;mock.nvar=nvars;mock.isstr=isstr;mock.isstrl=isstrl;
 mock.store=putnum;mock.safestore=putnum;mock.sstore=putstr;mock.macresave=macro_save;mock.macuse=macro_get;mock.missval=8.98846567431158e307;
 int (*tests[])(void)={interpolation,cubic_interpolation,splitting,joining};
 for(int k=0;k<4;k++){fail_at=0;assert(tests[k]()==0);int allocations=calls;
  for(int f=1;f<=allocations;f++){fail_at=f;assert(tests[k]()!=0);}
  fail_at=0;assert(tests[k]()==0);
 }
 assert(csplit_main("write")==198);assert(crangejoin_main("write")==198);assert(alive==0);
 for(int e=STATA_ERR_MEMORY;e<=STATA_ERR_CANCELLED;e++) {
  reset();nr=4;nv=3;load_error=e;fail_at=0;
  assert(cinterpolate_main("0 1 0")==ctools_stata_rc(e));assert(alive==0);
  for(int i=0;i<nr;i++)assert(numbers[2][i]==0);
 }
 load_error=0;
 /* No-fill destinations must be completely overwritten, including empty
  * source strings and rows with fewer tokens than the maximum. */
 reset();nr=3;nv=1;types[0]=1;
 strcpy(strings[0][0],"a::b");strcpy(strings[0][2],"c");
 assert(csplit_main("scan 1 2046 1 0 skipif fillall")==0);
 nv=2;types[0]=types[1]=1;
 for(int k=0;k<2;k++)for(int i=0;i<3;i++)strcpy(strings[k][i],"dirty");
 assert(csplit_main("write")==0);
 assert(!strcmp(strings[0][0],"a")&&!strcmp(strings[1][0],"b"));
 assert(!strcmp(strings[0][1],"")&&!strcmp(strings[1][1],""));
 assert(!strcmp(strings[0][2],"c")&&!strcmp(strings[1][2],""));
 reset();nr=2;nv=1;types[0]=1;
 assert(csplit_main("scan 1 2046 1 0 skipif fillall")==2000);assert(alive==0);
 reset();nr=3;nv=3;
 numbers[0][0]=10;numbers[0][1]=mock.missval;numbers[0][2]=30;
 numbers[1][0]=0;numbers[1][1]=mock.missval;numbers[1][2]=2;
 assert(cinterpolate_main("0 1 0 skipif")==0);
 assert(numbers[2][0]==10&&numbers[2][1]==mock.missval&&numbers[2][2]==30);
 assert(!strcmp(nmissing,"1"));assert(alive==0);
 /* Exact missing-key grouping and stable ties, independent of Stata's sorter. */
 reset();nr=4;nv=1;numbers[0][0]=mock.missval;numbers[0][1]=nextafter(mock.missval,INFINITY);numbers[0][2]=mock.missval;numbers[0][3]=2;
 ctools_filtered_data fd={0};assert(ctools_data_load(&fd,NULL,0,0,0,0)==0);int keys[]={0};assert(ctools_order_stable(&fd.data,keys,1)==0);
 assert(fd.data.sort_order[0]==3&&fd.data.sort_order[1]==0&&fd.data.sort_order[2]==2&&fd.data.sort_order[3]==1);ctools_filtered_data_free(&fd);assert(alive==0);
 puts("PASS new command kernels, allocation failures, stale-cache rejection, exact key ordering");return 0;
}
#endif

#if defined(CTOOLS_CASE_split_tokens_native_extra)
static const char *single_delimiter;
static ST_retcode single_macro(char *key,char *value,ST_int size) {
    (void)key;assert(size>(int)strlen(single_delimiter));strcpy(value,single_delimiter);return 0;
}
static void single_case(const char *input,const char *delimiter,int trim,int limit,int count,const char **expected) {
    reset();nr=1;nv=1;types[0]=1;strcpy(strings[0][0],input);
    single_delimiter=delimiter;mock.macuse=single_macro;
    char args[64];snprintf(args,sizeof(args),"scan 1 %d %d 0",limit,trim);
    assert(csplit_main(args)==0);nv=count;
    for(int j=0;j<count;j++)types[j]=1;
    assert(csplit_main("write")==0);
    for(int j=0;j<count;j++)assert(!strcmp(strings[j][0],expected[j]));
    assert(alive==0);
}
int main(void) {
    assert(original_main()==0);
    const char *a[]={"a","b","","c"},*b[]={"a","b"},*empty[]={""},*x[]={"x"},*utf[]={"é","東京"};
    single_case("  a::b::::c:: ","::",1,2046,4,a);
    single_case("  a   b  "," ",2,2046,2,b);
    single_case("  a   b  "," ",2,1,1,b);
    single_case("     "," ",2,2046,1,empty);
    single_case("     x","::",1,2046,1,x);
    single_case("é::東京","::",0,2046,2,utf);
    char longest[2046];memset(longest,'x',2045);longest[2045]=0;
    const char *long_expected[]={longest};single_case(longest,"::",0,2046,1,long_expected);
    puts("PASS in-place split: delimiter widths, trimming, empty fields, limits, UTF-8, maximum-width input");
    return 0;
}
#endif

#if defined(CTOOLS_CASE_oct22_match_oom)
#include <assert.h>
#include "cpsmatch/cpsmatch_internal.h"
static int fail_growth;
static void *probe_realloc(void *p,size_t n) {if(fail_growth) return NULL; return realloc(p,n);}
#define realloc probe_realloc
#include "cpsmatch/cpsmatch_kernels.c"
#undef realloc
#include "cpsmatch/cpsmatch_options.c"
#include "cpsmatch/cpsmatch_stats.c"
#include "cpsmatch/cpsmatch_radius.c"
#include "cpsmatch/cpsmatch_impl.c"
static ST_plugin mock;
ST_plugin *_stata_=&mock;
static int stores;
static double tr[301], ps[301], outputs[7][302], matched;
static perm_idx_t map[301];
static stata_variable columns[2];
static ST_int one(void){return 1;} static ST_int count(void){return 301;}
static ST_retcode store(ST_int v,ST_int i,ST_double x){stores++;outputs[v][i]=x;return 0;}
static ST_retcode scal(char *name,double x){if(!strcmp(name,"_cpsmatch_n_matched")) matched=x;return 0;}
static ST_boolean missing(double x){return x>=8.98846567431158e307;}
static ST_retcode message(char *s){fputs(s,stderr);return 0;}
double ctools_timer_seconds(void){return 0;}
int ctools_get_max_threads(void){return 1;}
void ctools_filtered_data_init(ctools_filtered_data *fd){memset(fd,0,sizeof(*fd));}
void ctools_filtered_data_free(ctools_filtered_data *fd){(void)fd;}
stata_retcode ctools_data_load(ctools_filtered_data *fd,int *v,size_t nv,size_t s,size_t e,int flags){
 (void)v;(void)nv;(void)s;(void)e;(void)flags;
 columns[0].data.dbl=tr;columns[1].data.dbl=ps;
 fd->data.nobs=301;fd->data.nvars=2;fd->data.vars=columns;fd->obs_map=map;return STATA_OK;
}
int main(void){
 mock.ismissing=missing;mock.nobs1=one;mock.nobs2=count;mock.safestore=store;mock.store=store;mock.scalsave=scal;mock.spouterr=message;mock.missval=8.98846567431158e307;
 tr[0]=1;for(int i=0;i<301;i++){ps[i]=0.5;map[i]=i+1;}
 const char *args="1 2 0 3 4 5 301 2 1 0 0 0 1 1 0 0";
 for(fail_growth=0;fail_growth<2;fail_growth++){
 stores=0; int rc=cpsmatch_main(args);
 assert(rc==(fail_growth ? 920 : 0));
 assert(!fail_growth || stores==0);
 printf("fail_growth=%d rc=%d matched=%.0f treated_support=%.0f\n",fail_growth,rc,matched,outputs[5][1]);
 }
}
#endif

#if defined(CTOOLS_CASE_oct22_range)
#include <stdio.h>
#include <assert.h>
#include "crangestat/crangestat_kernels.c"
ST_plugin mock;
ST_plugin *_stata_ = &mock;
int main(void) {
 mock.missval = 8.98846567431158e307;
 double data[] = {1e16,1,2,3};
 double sum[5],sq[5]; size_t count[5];
 prefix_arrays pa = {.prefix_sum=sum,.prefix_sum2=sq,.prefix_count=count};
 crangestat_build_prefix_arrays(data,0,4,&pa);
 for(size_t i=0;i<4;i++) {
  assert(crangestat_compute_stat_with_prefix(STAT_MEAN,&pa,data,0,i,i+1,i,0)==data[i]);
  assert(crangestat_compute_stat_with_prefix(STAT_SUM,&pa,data,0,i,i+1,i,0)==data[i]);
 }
 assert(crangestat_compute_stat_with_prefix(STAT_SUM,&pa,data,0,1,4,2,0)==6);
 assert(crangestat_compute_stat_with_prefix(STAT_MEAN,&pa,data,0,1,4,2,1)==2);
 puts("oct22 prefix cancellation passed");
}
#endif

#if defined(CTOOLS_CASE_p2_native_1)
#include <assert.h>
#include "ctools_hash.c"
int main(void) {
    ctools_str_hash_table ht;
    assert(ctools_str_hash_init(&ht, 16) == 0);
    const char *names[] = {"negative", "zero", "positive"};
    int codes[] = {-10, 0, 7};
    for (int i = 0; i < 3; i++) {
        assert(ctools_str_hash_insert_value(&ht, names[i], codes[i]) == 0);
        int value = 999;
        assert(ctools_str_hash_lookup(&ht, names[i], &value) == 1 && value == codes[i]);
    }
    int value = 999;
    assert(ctools_str_hash_lookup(&ht, "absent", &value) == 0 && value == 999);
    ctools_str_hash_free(&ht);
    return 0;
}
#endif

#if defined(CTOOLS_CASE_winsor_transport)
/* Production grouping, quantiles, gather and stores; mocked host and loader.
 * The tracked allocator from memory_safety/alloc.h is prepended by the runner. */
#include "ctools_types.h"
#include "ctools_parse.c"
#include "ctools_order.c"
#include "cwinsor/cwinsor_impl.c"

static ST_plugin host;
ST_plugin *_stata_ = &host;
static int targets, byvars, store_at, store_calls, return_bounds, load_calls, read_at, read_calls;
static int host_rows = 256;

static ST_int host_nobs(void) { return host_rows; }

static ST_int host_nvars(void)
{
    return targets + byvars;
}
static ST_retcode save_scalar(char *name, double value)
{
    return 0;
}
static double input_value(int v, int obs);
static ST_retcode read_macro(char *name, char *out, ST_int length)
{
    assert(strcmp(name, "_cwinsor_bounds") == 0);
    snprintf(out, length, "bounds");
    return 0;
}
static ST_retcode store_bound(char *name, ST_int row, ST_int col, double value)
{
    assert(strcmp(name, "bounds") == 0);
    assert(value == input_value(row, col == 1 ? 1 : 255));
    return ++store_calls == store_at ? 459 : 0;
}
static ST_retcode read_bound(char *name, ST_int row, ST_int col, double *value)
{
    assert(strcmp(name, "bounds") == 0);
    *value = input_value(row, col == 1 ? 1 : 255);
    return ++read_calls == read_at ? 459 : 0;
}
int ctools_get_max_threads(void)
{
    return 4;
}
double ctools_timer_seconds(void)
{
    return 0;
}
void ctools_error(const char *a, const char *b, ...) {}
void ctools_error_alloc(const char *a) {}

void ctools_filtered_data_init(ctools_filtered_data *d)
{
    memset(d, 0, sizeof(*d));
}
void ctools_filtered_data_free(ctools_filtered_data *d)
{
    for (size_t j = 0; j < d->data.nvars; j++)
        free(d->data.vars[j].data.dbl);
    free(d->data.vars);
    free(d->data.sort_order);
    free(d->obs_map);
    memset(d, 0, sizeof(*d));
}

static double input_value(int v, int obs)
{
    return obs % 7 == 0 ? host.missval : obs + 1000.0 * v;
}

stata_retcode ctools_data_load(ctools_filtered_data *d, int *vars, size_t n, size_t start,
                               size_t end, int flags)
{
    load_calls++;
    const int count = 128;
    d->data.vars = calloc(n, sizeof(stata_variable));
    if (!d->data.vars)
        return STATA_ERR_MEMORY;
    d->data.nvars = n;
    d->data.nobs = count;
    d->obs_map = malloc(count * sizeof(perm_idx_t));
    if (!d->obs_map)
        goto failure;
    if (!(flags & CTOOLS_LOAD_NO_SORT_ORDER))
    {
        d->data.sort_order = malloc(count * sizeof(perm_idx_t));
        if (!d->data.sort_order)
            goto failure;
    }
    for (int i = 0; i < count; i++)
    {
        d->obs_map[i] = 2 * i + 1;
        if (d->data.sort_order)
            d->data.sort_order[i] = i;
    }
    for (size_t j = 0; j < n; j++)
    {
        d->data.vars[j].type = STATA_TYPE_DOUBLE;
        d->data.vars[j].nobs = count;
        double *x = d->data.vars[j].data.dbl = malloc(count * sizeof(double));
        if (!x)
            goto failure;
        for (int i = 0; i < count; i++)
            x[i] = vars[j] <= targets ? input_value(vars[j], 2 * i + 1) : i % 4;
    }
    return STATA_OK;
failure:
    /* The real loader frees partial allocations before returning an error. */
    ctools_filtered_data_free(d);
    return STATA_ERR_MEMORY;
}

stata_retcode ctools_sort_dispatch(stata_data *d, int *keys, size_t n, sort_algorithm_t alg)
{
    /* Replace the permutation allocation, as a sorting backend is allowed to. */
    perm_idx_t *order = malloc(d->nobs * sizeof(perm_idx_t));
    if (!order)
        return STATA_ERR_MEMORY;
    size_t next = 0;
    for (size_t group = 0; group < 4; group++)
        for (size_t i = group; i < d->nobs; i += 4)
            order[next++] = i;
    free(d->sort_order);
    d->sort_order = order;
    return STATA_OK;
}

stata_retcode ctools_store_filtered(double *values, size_t n, int var, perm_idx_t *map)
{
    for (size_t i = 0; i < n; i++)
    {
        assert(values[i] == input_value(var, map[i]));
        if (__atomic_add_fetch(&store_calls, 1, __ATOMIC_RELAXED) == store_at)
            return 459;
    }
    return STATA_OK;
}
stata_retcode ctools_store_filtered_rowpar(double *v, size_t n, int var, perm_idx_t *map)
{
    return ctools_store_filtered(v, n, var, map);
}

static int run_command(void)
{
    char args[128];
    snprintf(args, sizeof(args), "%d %d %s p=0 q=100 %s", targets, byvars,
             targets == 1 ? "1" : "1 2 3",
             return_bounds == 2 ? "deferbounds" : return_bounds ? "trim returnbounds" : "");
    store_calls = read_calls = load_calls = 0;
    int rc = cwinsor_main(args);
    if (!rc && return_bounds == 2) {
        assert(pending.nvars == (size_t)targets);
        assert(store_calls == 2*targets); /* Only the cutoff matrix was written. */
        rc = cwinsor_main("apply");
        assert(!pending.nvars);
        assert(load_calls == 1);
    }
    cwinsor_cleanup_state();
    leftovers();
    assert(alive == 0);
    return rc;
}

int main(void)
{
    host.nvars = host_nvars;
    host.nobs = host_nobs;
    host.scalsave = save_scalar;
    host.macuse = read_macro;
    host.matstore = store_bound;
    host.safematstore = store_bound;
    host.matel = read_bound;
    host.safematel = read_bound;
    host.missval = 8.98846567431158e307;
#ifdef _OPENMP
    omp_set_num_threads(4);
#endif
    int checked = 0;
    for (return_bounds = 0; return_bounds <= 2; return_bounds++)
        for (targets = 1; targets <= 3; targets += 2)
            for (byvars = 0; byvars <= 2; byvars += 2)
            {
                if (return_bounds && byvars)
                    continue;
                calls = fail_at = store_at = 0;
                assert(run_command() == 0);
                int allocations = calls;
                for (int failure = 1; failure <= allocations; failure++)
                {
                    calls = 0;
                    fail_at = failure;
                    int rc = run_command();
                    assert(rc == 0 || rc == 920);
                    fail_at = 0;
                    assert(run_command() == 0);
                    checked++;
                }
                for (store_at = 1; store_at <= (return_bounds == 1 ? 2 * targets : 128); store_at *= 2)
                {
                    assert(run_command() == 459);
                    checked++;
                }
                store_at = 0;
                if (return_bounds == 2) {
                    for (read_at = 1; read_at <= 2*targets; read_at++) {
                        assert(run_command() == 459);
                        assert(store_calls == 2*targets); /* No data stores on a cutoff read failure. */
                        checked++;
                    }
                    read_at = 0;
                }
            }
    targets=1; byvars=0; fail_at=store_at=read_at=0;
    assert(cwinsor_main("apply") == 198);
    assert(cwinsor_main("1 0 1 p=0 q=100 deferbounds") == 0);
    host_rows++;
    assert(cwinsor_main("apply") == 198);
    host_rows--;
    assert(!pending.nvars);
    assert(cwinsor_main("1 0 1 p=0 q=100 deferbounds") == 0);
    assert(cwinsor_main("commit") == 0);
    assert(cwinsor_main("apply") == 198);
    leftovers(); assert(alive == 0);
    printf("PASS winsor transport: %d allocation/store injections, no leaks\n", checked);
    return 0;
}
#endif

#if defined(CTOOLS_CASE_sep22_native_source)
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "stplugin.h"
static ST_plugin mock;
ST_plugin *_stata_ = &mock;
static int calls, fail_at, live;
static ST_boolean ismissing(ST_double x) { return x >= 8.98846567431158e307; }
static void *alloc(size_t n) {
    if (++calls == fail_at) return NULL;
    void *p = malloc(n); if (p) live++; return p;
}
static void *zero(size_t n, size_t k) {
    if (++calls == fail_at) return NULL;
    void *p = calloc(n,k); if (p) live++; return p;
}
static void release(void *p) { if (p) live--; free(p); }
#define malloc alloc
#define calloc zero
#define free release
#include "ctools_hdfe_utils.c"
#include "cppmlhdfe/cppmlhdfe_separation.c"
int main(void) {
    mock.missval = 8.98846567431158e307;
    mock.ismissing = ismissing;
    double codes[] = {-1, -.1, .1, .2, 1e12, 1e12+1, -.1, 8.98846567431158e307};
    int ids[8], groups;
    assert(ctools_numeric_to_cluster_ids(codes,8,ids,&groups) == 0);
    assert(groups == 6 && ids[7] == -1 && ids[1] == ids[6]);
    for (int i=0; i<6; i++) for (int j=0; j<i; j++) assert(ids[i] != ids[j]);
    assert(live == 0);
    calls=0; fail_at=1;
    assert(ctools_numeric_to_cluster_ids(codes,8,ids,&groups) < 0 && live == 0);
    fail_at=0;
    /* Removed blocks at start/middle/end, with an initial excluded row and
       a singleton. Allocations are exact-sized, so ASan sees boundary errors. */
    for (int block=0; block<3; block++) {
        int *level = malloc(14*sizeof(int));
        int *mask = malloc(14*sizeof(int));
        double *y = malloc(14*sizeof(double));
        for (int i=0; i<14; i++) { level[i]=i/4+1; mask[i]=1; y[i]=1; }
        mask[13]=0;
        for (int i=block*4; i<block*4+4; i++) y[i]=0;
        int *levels[]={level}, nlevels[]={4}, singletons, separated;
        int base=live;
        calls=0;
        assert(ppml_select_fe_sample(y,levels,nlevels,1,14,mask,&singletons,&separated)==0);
        int allocations=calls;
        assert(singletons==1 && separated==4 && live==base);
        for (int i=0; i<14; i++) assert(mask[i] == (i<12 && i/4!=block));
        for (int failure=1; failure<=allocations; failure++) {
            for(int i=0;i<14;i++) mask[i]=i!=13;
            calls=0; fail_at=failure;
            assert(ppml_select_fe_sample(y,levels,nlevels,1,14,mask,&singletons,&separated)<0);
            assert(live==base);
        }
        fail_at=0;
        free(level); free(mask); free(y);
        assert(live==0);
    }
    /* Removing a zero group creates a singleton in a crossing FE. */
    int a[]={1,1,2,2,2}, b[]={1,1,1,2,2}, mask[]={1,1,1,1,1};
    double y[]={0,0,1,1,1};
    int *levels[]={a,b}, counts[]={2,2}, ns, nz;
    assert(ppml_select_fe_sample(y,levels,counts,2,5,mask,&ns,&nz)==0);
    assert(nz==2 && ns==1 && !mask[2] && mask[3] && mask[4]);
    assert(live==0);
    return 0;
}
#endif

#if defined(CTOOLS_CASE_memory_rangestat)
/* Full command with mocked Stata loading/output; numerical code is production code. */
#include "alloc.h"
#include "ctools_types.h"
#include "ctools_parse.c"
#include "stplugin.h"
#include <math.h>
static int store_at, store_calls;
static ST_plugin mock;
ST_plugin *_stata_ = &mock;
static ST_int one(void) { return 1; }
static ST_int nvars(void) { return 6; }
static ST_int obs(void) { return 12; }
static ST_boolean selected(int i) { return 1; }
static ST_retcode save(char *k, double v) { return 0; }
static ST_retcode store(ST_int v, ST_int i, double x) {
    assert(x == (i - 1) % 10);
    return __atomic_add_fetch(&store_calls, 1, __ATOMIC_RELAXED) == store_at ? 459 : 0;
}
int ctools_get_max_threads(void) { return 8; }
double ctools_timer_seconds(void) { return 0; }
void ctools_error(const char *a, const char *b, ...) {}
void ctools_error_alloc(const char *a) {}
void ctools_msg(const char *a, const char *b, ...) {}
void ctools_filtered_data_init(ctools_filtered_data *d) { memset(d, 0, sizeof(*d)); }
void ctools_filtered_data_free(ctools_filtered_data *d) {
    for (size_t j = 0; j < d->data.nvars; j++)
        free(d->data.vars[j].data.dbl);
    free(d->data.vars);
    free(d->data.sort_order);
    free(d->obs_map);
    memset(d, 0, sizeof(*d));
}
stata_retcode ctools_data_load(ctools_filtered_data *d, int *v, size_t n, size_t s, size_t e,
                               int f) {
    int N = 512;
    d->data.vars = calloc(n, sizeof(stata_variable));
    if (!d->data.vars)
        return STATA_ERR_MEMORY;
    d->data.nvars = n;
    d->data.nobs = N;
    d->data.sort_order = malloc(N * sizeof(perm_idx_t));
    d->obs_map = malloc(N * sizeof(perm_idx_t));
    if (!d->obs_map || !d->data.sort_order)
        return STATA_ERR_MEMORY;
    for (int i = 0; i < N; i++) {
        d->obs_map[i] = i + 1;
        d->data.sort_order[i] = i;
    }
    for (size_t j = 0; j < n; j++) {
        stata_variable *q = d->data.vars + j;
        q->type = STATA_TYPE_DOUBLE;
        q->nobs = N;
        q->data.dbl = malloc(N * sizeof(double));
        if (!q->data.dbl)
            return STATA_ERR_MEMORY;
        for (int i = 0; i < N; i++)
            q->data.dbl[i] = j == 0 ? i / 64 : j == 1 ? i % 64 : i % 10;
    }
    return STATA_OK;
}
stata_retcode ctools_sort_dispatch(stata_data *d, int *v, size_t n, sort_algorithm_t a) {
    return STATA_OK;
}
stata_retcode ctools_apply_permutation(stata_data *d) { return STATA_OK; }

#include "crangestat/crangestat_impl.c"
#include "crangestat/crangestat_kernels.c"
#include "ctools_order.c"
int main(int argc, char **argv) {
    store_at = argc > 2 ? atoi(argv[2]) : 0;
    fail_at = argc > 1 ? atoi(argv[1]) : 0;
    mock.nobs1 = one;
    mock.nobs2 = obs;
    mock.nobs = obs;
    mock.nvars = nvars;
    mock.selobs = selected;
    mock.scalsave = save;
    mock.safestore = store;
    mock.store = store;
    mock.missval = 8.98846567431158e307;
    int rc = crangestat_main("2 2 1 2 3 4 5 6 1 3 4 0 1 low=0 high=0");
    int first_rc = rc, first_calls = calls;
    assert(alive == 0);
    fail_at = 0;
    store_at = 0;
    calls = 0;
    assert(crangestat_main("2 2 1 2 3 4 5 6 1 3 4 0 1 low=0 high=0") == 0);
    assert(alive == 0);
    rc = first_rc;
    calls = first_calls;
    printf("rc=%d calls=%d alive=%d\n", rc, calls, alive);
    leftovers();
    assert(alive == 0);
    return 0;
}
#endif

#if defined(CTOOLS_CASE_memory_sampling)
/* Full command with mocked Stata loading/output; numerical code is production code. */
#include "alloc.h"
#include "ctools_types.h"
#include "ctools_parse.c"
#include "ctools_group.c"
#include "stplugin.h"
#include <math.h>
static int store_at, store_calls;
static ST_plugin mock;
ST_plugin *_stata_ = &mock;
static ST_int one(void) { return 1; }
static ST_int obs(void) { return 12; }
static ST_boolean selected(int i) { return 1; }
static ST_retcode save(char *k, double v) { return 0; }
static ST_retcode store(ST_int v, ST_int i, double x) {
    return ++store_calls == store_at ? 459 : 0;
}
int ctools_get_max_threads(void) { return 8; }
double ctools_timer_seconds(void) { return 0; }
void ctools_error(const char *a, const char *b, ...) {}
void ctools_error_alloc(const char *a) {}
void ctools_msg(const char *a, const char *b, ...) {}
void ctools_filtered_data_init(ctools_filtered_data *d) { memset(d, 0, sizeof(*d)); }
void ctools_filtered_data_free(ctools_filtered_data *d) {}
stata_retcode ctools_data_load(ctools_filtered_data *d, int *v, size_t n, size_t s, size_t e,
                               int f) {
    abort();
}
stata_retcode ctools_sort_dispatch(stata_data *d, int *v, size_t n, sort_algorithm_t a) { abort(); }
stata_retcode ctools_apply_permutation(stata_data *d) { abort(); }
#ifdef BOOTSTRAP
#include "cbsample/cbsample_impl.c"
#else
#include "csample/csample_impl.c"
#endif
int main(int argc, char **argv) {
    store_at = argc > 2 ? atoi(argv[2]) : 0;
    fail_at = argc > 1 ? atoi(argv[1]) : 0;
    mock.nobs1 = one;
    mock.nobs2 = obs;
    mock.nobs = obs;
    mock.nvars = one;
    mock.selobs = selected;
    mock.scalsave = save;
    mock.safestore = store;
    mock.store = store;
#ifdef BOOTSTRAP
    int rc = cbsample_main("1 0 0 n=6 seedhi=1 seedlo=2");
#else
    int rc = csample_main("1 0 count=6 seedhi=1 seedlo=2");
#endif
    assert(alive == 0);
    /* Even with eight configured workers, one group gets one RNG/workspace. */
#ifdef BOOTSTRAP
    assert(peak_bytes <= 700);
#else
    assert(peak_bytes <= 512);
#endif
    if (store_at)
        assert(rc == 459);
    printf("rc=%d allocations=%d alive=%d\n", rc, calls, alive);
    leftovers();
    assert(alive == 0);
    return 0;
}
#endif

#if defined(CTOOLS_CASE_memory_remap)
#include "alloc.h"
#include "stplugin.h"
#include <math.h>
static ST_plugin mock;
ST_plugin *_stata_ = &mock;
#include "hdfe/hdfe_remap.c"
int main(void) {
    double cases[][6] = {{1, 1, 2, 2, -1, -1},
                         {1.5, 1.5, 2.5, 2.5, -1, -1},
                         {1e30, 1e30, -1e30, -1e30, 0, 0},
                         {-0x1p63, 0, 0x1p62, 1, -0x1p63, 1},
                         {0x1p63, 0x1p63, -0x1p64, 2, 2, -0x1p64}};
    double weights[] = {1, 2, 3, 4, 5, 6};
    for (size_t c = 0; c < sizeof(cases) / sizeof(*cases); c++) {
        for (int weighted = 0; weighted < 2; weighted++) {
            int allocations = 0;
            for (int f = 0; f <= allocations; f++) {
                ST_int levels[6], nlevels = 123, *counts = (void *)1;
                double *wcounts = (void *)1;
                calls = 0;
                fail_at = f;
                int rc = remap_and_count(cases[c], 6, levels, &nlevels, &counts,
                                         weighted ? weights : NULL, &wcounts);
                if (!f)
                    allocations = calls;
                if (f) {
                    assert(rc && nlevels == 0 && !counts && !wcounts);
                } else {
                    assert(!rc && nlevels > 0 && counts);
                    int total = 0;
                    double sum = 0;
                    for (int j = 0; j < nlevels; j++) {
                        total += counts[j];
                        if (weighted)
                            sum += wcounts[j];
                    }
                    assert(total == 6 && (!weighted || sum == 21));
                    for (int i = 0; i < 6; i++)
                        for (int j = 0; j < 6; j++)
                            assert((levels[i] == levels[j]) == (cases[c][i] == cases[c][j]));
                    free(counts);
                    free(wcounts);
                }
                assert(alive == 0);
            }
        }
    }
    return 0;
}
#endif

#if defined(CTOOLS_CASE_sep24_native_6)
#include "ctools_arena.c"
int main(void) {
 ctools_arena a;ctools_arena_init(&a,32);
 size_t capacity=0;
 for(int run=0;run<50;run++) {
  assert(ctools_arena_alloc(&a,24));assert(ctools_arena_alloc(&a,64));
  void *aligned=ctools_arena_alloc_aligned(&a,53,64);assert(aligned && (uintptr_t)aligned%64==0);
  assert(!ctools_arena_alloc(&a,SIZE_MAX));assert(!ctools_arena_alloc_aligned(&a,SIZE_MAX,64));
  if(run==0)capacity=live_count;else assert(capacity==(size_t)live_count);
  ctools_arena_reset(&a);
 }
 ctools_arena_free(&a);assert(live_count==0);return 0;
}
#endif

#if defined(CTOOLS_CASE_sep24_native_8)
#include "ctools_types.h"
#include "cmerge/cmerge_keys.c"
#include "cmerge/cmerge_group_search.c"
#include "cmerge/cmerge_join.c"
int main(void) {
 double number[]={1};char *text[]={"1"};
 stata_variable a[2]={{.type=STATA_TYPE_DOUBLE,.data.dbl=number},{.type=STATA_TYPE_STRING,.data.str=text}};
 stata_variable b[2]={{.type=STATA_TYPE_DOUBLE,.data.dbl=number},{.type=STATA_TYPE_DOUBLE,.data.dbl=number}};
 stata_data da={.nobs=1,.nvars=2,.vars=a},db={.nobs=1,.nvars=2,.vars=b};
 cmerge_output_spec_t *output=(void*)1;
 assert(cmerge_sorted_join(&da,&db,2,MERGE_1_1,&output)==-4 && !output);
 assert(cmerge_sorted_join(&db,&da,2,MERGE_M_M,&output)==-4 && !output);
 assert(live_count==0);return 0;
}
#endif

#if defined(CTOOLS_CASE_sep24_native_atomic)
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stplugin.h"
static ST_plugin mock;ST_plugin *_stata_=&mock;
#include "cexport/cexport_parse.c"
#include "formats/ctools_output.c"
static void put(const char *p,const char *s){FILE *f=fopen(p,"w");assert(f);assert(fputs(s,f)>=0);assert(!fclose(f));}
static void check(const char *p,const char *s){char buf[30]={0};FILE *f=fopen(p,"r");assert(f);fread(buf,1,29,f);fclose(f);assert(!strcmp(buf,s));}
int main(int argc,char **argv){
 assert(argc==2);ctools_output out;const char *p=argv[1];put(p,"original");
 assert(!ctools_output_prepare(&out,p,1));put(out.temporary,"partial");ctools_output_cleanup(&out);check(p,"original");
 assert(!ctools_output_prepare(&out,p,0));put(out.temporary,"new");assert(ctools_output_commit(&out)==602);ctools_output_cleanup(&out);check(p,"original");
 assert(!ctools_output_prepare(&out,p,1));put(out.temporary,"complete");assert(!ctools_output_commit(&out));ctools_output_cleanup(&out);check(p,"complete");
 remove(p);return 0;
}
#endif
