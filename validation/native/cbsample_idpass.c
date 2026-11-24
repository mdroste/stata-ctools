/* Expanded bootstrap rows: allocation/read/store failures must propagate. */
#include "newcommands_mock.h"
#include "ctools_runtime.h"
void ctools_error_alloc(const char *m) { (void)m; }
#include "cbsample/cbsample_impl.c"

static int reads,fail_read,read_rc;
static ST_int one(void) { return 1; }
static ST_retcode getnum(ST_int v,ST_int r,ST_double *x)
{
    assert(v>=1&&v<=3&&r>=1&&r<=nr);
    if(++reads==fail_read)return read_rc;
    *x=numbers[v-1][r-1];return 0;
}
static void init(void)
{
    assert(!live);clear_fixture();resize_fixture(6);nv=5;
    reads=fail_read=0;read_rc=459;
    /* Two copies of a two-row cluster, then one two-row cluster. */
    double positions[]={1,2,1,2,3,4};
    for(int r=0;r<6;r++) {
        numbers[0][r]=positions[r];numbers[1][r]=r<4?1:2;
        numbers[2][r]=r<4?1:2;
        numbers[3][r]=numbers[4][r]=-1;
    }
}
int main(void)
{
    mock.nobs1=one;mock.nobs2=nobs;mock.nobs=nobs;mock.nvars=nvars;
    mock.vdata=getnum;mock.safevdata=getnum;mock.store=putnum;mock.safestore=putnum;
    mock.missval=8.98846567431158e307;
    for(int random=0;random<2;random++) {
        const char *args=random?"idpass iduniq idrandom seedhi=1 seedlo=2":"idpass iduniq seedhi=1 seedlo=2";
        init();assert(!cbsample_idcluster_pass(args));assert(!live);
        int count=allocations;
        for(int fault=1;fault<=count;fault++) {
            init();fail_alloc=fault;
            assert(cbsample_idcluster_pass(args)==920&&live==0&&stores==0);
        }
        for(int fault=1;fault<=18;fault++) for(int code=0;code<2;code++) {
            init();fail_read=fault;read_rc=code?1:459;
            assert(cbsample_idcluster_pass(args)==read_rc&&live==0&&stores==0);
        }
        for(int fault=1;fault<=12;fault++) {
            init();fail_store=fault;
            assert(cbsample_idcluster_pass(args)==459&&live==0);
        }
        init();assert(!cbsample_idcluster_pass(args)&&live==0);
    }
    clear_fixture();puts("PASS bootstrap id pass: every allocation, read, store and Break fault");
    return 0;
}
