/* Complete file/column pipeline, allocation faults, and bounded mutation corpus. */
#include "parquet_alloc.h"
#include "io/parquet_thrift.c"
#include "io/parquet_decode.c"
#include "io/parquet_codec.c"
#include "io/parquet_zstd.c"
#include "io/parquet_brotli.c"
#include PARQUET_CASES
static unsigned read_u32(FILE *f)
{ unsigned char b[4];assert(fread(b,1,4,f)==4);return b[0]|(unsigned)b[1]<<8|(unsigned)b[2]<<16|(unsigned)b[3]<<24; }
static void expected(FILE *f,const pq_values *v,int64_t first,int64_t last,int rows)
{
    int kind=fgetc(f);assert(kind>=0&&kind<=2);
    for(int i=0;i<rows;i++) {
        int valid=fgetc(f);assert(valid==0||valid==1);int selected=i>=first&&i<=last;
        size_t k=(size_t)(i-first);
        if(selected)assert(v->valid[k]==valid);
        if(!valid)continue;
        if(kind==2) {
            unsigned n=read_u32(f);unsigned char s[4096];assert(n<=sizeof(s)&&fread(s,1,n,f)==n);
            if(selected)assert(v->off&&v->off[k+1]-v->off[k]==n&&(!n||!memcmp(v->bytes+v->off[k],s,n)));
        } else if(kind==1) {
            int64_t x;assert(fread(&x,8,1,f)==1);if(selected)assert(v->i64&&v->i64[k]==x);
        } else {
            double x;assert(fread(&x,8,1,f)==1);if(selected)assert(v->num&&v->num[k]==x);
        }
    }
}
static int run_file(const char *name,int valid,int failure,int subset)
{
    assert(!pq_live&&!pq_live_bytes);pq_allocations=0;pq_fail_at=failure;pq_denied=0;
    char path[2048],err[512];snprintf(path,sizeof(path),"%s/%s.parquet",PARQUET_FIXTURES,name);
    pq_file file;int rc=pq_open(&file,path,err,sizeof(err));FILE *ref=NULL;int rows=0,cols=0;
    if(!rc&&valid) {
        snprintf(path,sizeof(path),"%s/%s.ref",PARQUET_FIXTURES,name);ref=fopen(path,"rb");assert(ref);
        rows=(int)read_u32(ref);cols=(int)read_u32(ref);assert(file.num_rows==rows&&file.nfields>=cols);
    }
    for(int j=0;!rc&&j<file.nfields;j++) {
        if(file.fields[j].leaf<0) { rc=pq_check_field(&file,file.fields+j,err,sizeof(err));continue; }
        int64_t first=subset?2:0,last=valid?rows-1:file.num_rows-1;
        if(subset&&last>7)last=7;
        if(!valid&&last>31)last=31;
        pq_values values;rc=pq_read_column(&file,file.fields+j,first,last,2,&values,err,sizeof(err));
        if(!rc&&valid&&j<cols)expected(ref,&values,first,last,rows);
        pq_values_free(&values);
    }
    if(ref)fclose(ref);pq_close(&file);assert(!pq_live&&!pq_live_bytes);
    if(failure&&rc==0)fprintf(stderr,"allocation %d unexpectedly succeeded in %s\n",failure,name);
    return rc;
}
static void mutations(void)
{
    char path[2048];snprintf(path,sizeof(path),"%s/good_0_1_plain.parquet",PARQUET_FIXTURES);
    FILE *f=fopen(path,"rb");assert(f);unsigned char data[4096],original[4096];
    size_t n=fread(original,1,sizeof(original),f);assert(n<sizeof(original));fclose(f);
    snprintf(path,sizeof(path),"%s/mutant.parquet",PARQUET_FIXTURES);
    for(size_t i=0;i<n;i++)for(int mode=0;mode<2;mode++) {
        memcpy(data,original,n);data[i]=mode?(unsigned char)(data[i]^0x80):0xff;
        f=fopen(path,"wb");assert(f&&fwrite(data,1,n,f)==n);fclose(f);
        (void)run_file("mutant",0,0,0);
    }
    for(size_t cut=0;cut<n;cut++) {
        f=fopen(path,"wb");assert(f&&fwrite(original,1,cut,f)==cut);fclose(f);
        assert(run_file("mutant",0,0,0)!=0);
    }
    printf("PASS %zu bounded full-file byte mutations/truncations\n",n*3);
}
int main(void)
{
    libdeflate_set_memory_allocator(pq_test_malloc,pq_test_free);
    for(size_t i=0;i<sizeof(cases)/sizeof(*cases);i++) {
        int rc=run_file(cases[i].name,cases[i].valid,0,0);
        if(cases[i].valid) {
            if(rc)fprintf(stderr,"valid %s rc=%d\n",cases[i].name,rc);
            assert(!rc);int count=pq_allocations;
            assert(!run_file(cases[i].name,1,0,1));
            if(cases[i].sweep)for(int f=1;f<=count;f++)assert(run_file(cases[i].name,1,f,0)!=0);
        } else {
            if(!rc)fprintf(stderr,"malformed %s accepted\n",cases[i].name);
            assert(rc);
            if(!strcmp(cases[i].name,"bad_page_allocation"))assert(!pq_denied);
        }
    }
    mutations();puts("PASS complete Parquet files, independent cells, row ranges, codecs and allocation faults");return 0;
}
