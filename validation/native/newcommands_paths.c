#include "newcommands_mock.h"
#include "ctools_order.c"
#include "cinterpolate/cinterpolate_impl.c"
#define state split_state
#include "csplit/csplit_impl.c"
#undef state
#define state join_state
#include "crangejoin/crangejoin_impl.c"
#undef state

static int last_allocations,last_stores;
static void group_key(int col,int row,int group,int string)
{
    if(string) { char s[24];snprintf(s,sizeof(s),"group-%d",group);text_cell(col,row,s);types[col]=1; }
    else numbers[col][row]=group;
}

static void interpolation_cases(int groups,int by,int string,int cubic)
{
    clear_fixture();resize_fixture(groups*12);nv=by+3;
    const int xs[]={0,1,2,3,4,4,5,6,7,8,9,-1};
    size_t missing=0,loaded=0;
    for(int i=0;i<nr;i++) {
        int g=i%groups,t=i/groups,x=xs[t];
        numbers[0][i]=mock.missval;numbers[1][i]=x<0?mock.missval:x;
        group_key(2,i,g,string);
        if(by==2)group_key(3,i,g%7,!string);
        if(g%5!=0) {
            if(t==4)numbers[0][i]=g+(cubic?62:6);
            if(t==5)numbers[0][i]=g+(cubic?66:10);
            if(g%5!=1 && (cubic ? (x>=0&&x%2==0&&x!=4) : (t==0||t==9)))
                numbers[0][i]=g+(cubic?x*x*x:2*x);
        }
        selected[i]=(g%13!=12);
        loaded+=selected[i];
        if(!selected[i]||x<0||g%5==0||(g%5==1&&x!=4)||(cubic&&(x==1||x==7||x==9)))missing++;
    }
    char args[24];snprintf(args,sizeof(args),"%d %d 0 %s",by,!cubic,cubic?"cubic":"");
    assert(!cinterpolate_main(args));assert(!live);
    int allocs=allocations;last_stores=stores;
    assert((size_t)last_stores==loaded && strtoull(nmissing,NULL,10)==missing);
    for(int i=0;i<nr;i++) {
        int g=i%groups,x=xs[i/groups];
        double expected=(!selected[i]||x<0||g%5==0||(g%5==1&&x!=4)||
            (cubic&&(x==1||x==7||x==9)))?mock.missval:g+(cubic?x*x*x:2*x);
        assert(numbers[nv-1][i]==expected);
    }
    int points[]={1,(int)loaded/2,(int)loaded};
    for(int k=0;k<3;k++) {
        allocations=stores=0;fail_store=points[k];
        assert(cinterpolate_main(args)==459 && !live);
        fail_store=0;assert(!cinterpolate_main(args) && !live);
    }
    fail_macro=1;assert(cinterpolate_main(args)==459 && !live);
    fail_macro=0;assert(!cinterpolate_main(args) && !live);
    if(groups<10)for(int f=1;f<=allocs;f++) {
        allocations=0;fail_alloc=f;assert(cinterpolate_main(args)==920 && !live);
    }
    fail_alloc=0;assert(!cinterpolate_main(args) && !live);
    clear_fixture();
}

static char longest[2046];
static const char *split_input(int pattern,int mode)
{
    static const char *input[]={"  alpha::beta  ","","é::東京","a::b:c"," lonely ","   ",NULL,"a::::b"};
    if(mode==2)return pattern%2?"  a   b  ":"   ";
    return pattern==6?longest:input[pattern];
}
static const char *split_expected(int pattern,int mode,int trim,int col)
{
    static const char *single[][4]={{"alpha","beta","",""},{"","","",""},
        {"é","東京","",""},{"a","b:c","",""},{"lonely","","",""},
        {"","","",""},{NULL,"","",""},{"a","","b",""}};
    static const char *multi[][4]={{"alpha","","beta",""},{"","","",""},
        {"é","","東京",""},{"a","","b","c"},{"lonely","","",""},
        {"","","",""},{NULL,"","",""},{"a","","",""}};
    if(mode==2)return pattern%2?(col==0?"a":col==1?"b":""):"";
    if(pattern==6&&col==0)return longest;
    if(!trim) {
        if(pattern==0 && col==0)return "  alpha";
        if(pattern==0 && col==(mode?2:1))return "beta  ";
        if(pattern==4&&col==0)return " lonely ";
        if(pattern==5&&col==0)return "   ";
    }
    return mode?multi[pattern][col]:single[pattern][col];
}
static int split_case(int rows,int mode,int trim,int limit,int filter,int allocation,int publication,int macro)
{
    clear_fixture();resize_fixture(rows);nv=1;types[0]=1;
    delimiters[0]=mode==2?" ":"::";delimiters[1]=":";
    for(int i=0;i<rows;i++) { text_cell(0,i,split_input(i%8,mode));selected[i]=!filter||i%11; }
    fail_alloc=allocation;fail_store=publication;fail_macro=macro;
    char args[96];snprintf(args,sizeof(args),"scan %d %d %d 0 %s",mode==1?2:1,limit,trim,filter?"":"skipif fillall");
    int rc=csplit_main(args);
    if(!rc) {
        int fields=mode==0?3:mode==1?4:2;if(fields>limit)fields=limit;
        assert(split_state.fields==fields);nv=fields;
        const char *width=widths;
        for(int c=0;c<fields;c++) {
            types[c]=1;size_t maximum=1;
            for(int i=0;i<rows;i++) {
                if(selected[i]) { size_t len=strlen(split_expected(i%8,mode,trim,c));if(len>maximum)maximum=len; }
                text_cell(c,i,filter?"":"dirty");
            }
            char *end;assert(strtoul(width,&end,10)==maximum);width=end;
        }
        rc=csplit_main("write");
        if(!rc)for(int c=0;c<fields;c++)for(int i=0;i<rows;i++)
            assert(!strcmp(string_cell(c,i),selected[i]?split_expected(i%8,mode,trim,c):""));
#ifdef _OPENMP
        if(!rc&&threads>1&&rows>10000)assert(worker_store);
#endif
    }
    last_allocations=allocations;last_stores=stores;
    csplit_cleanup_cache();assert(!live);clear_fixture();return rc;
}

static double join_group(int group)
{ return group==4?mock.missval:group==5?nextafter(mock.missval,INFINITY):group; }
static int master_group(int row)
{ return row<4000?(row/503)%8:row%8; }
static void join_key(int col,int row,int group,int string)
{
    if(string) { char s[16];if(group)snprintf(s,sizeof(s),"g%d",group);else s[0]=0;text_cell(col,row,s);types[col]=1; }
    else numbers[col][row]=join_group(group);
}
static void interval(int row,double *low,double *high)
{
    switch(row%6) {
    case 0:*low=1;*high=3;break;
    case 1:*low=3;*high=3;break;
    case 2:*low=9;*high=10;break;
    case 3:*low=5;*high=1;break;
    case 4:*low=-mock.missval;*high=mock.missval;break;
    default:*low=mock.missval;*high=mock.missval;
    }
}
static void check_master(int row,int output,int by,int string)
{
    assert(numbers[0][output]==row+1);
    int g=master_group(row);
    for(int c=0;c<by;c++) {
        int str=c? !string:string;
        if(str) { char s[16];if(g)snprintf(s,sizeof(s),"g%d",g);else s[0]=0;assert(!strcmp(string_cell(c+1,output),s)); }
        else assert(numbers[c+1][output]==join_group(g));
    }
}
static int join_case(int rows,int by,int string,int allocation,int publication,int macro,int using_rows)
{
    clear_fixture();resize_fixture(using_rows);nv=by+3;
    const int key[]={1,1,3,5};
    for(int i=0;i<using_rows;i++) {
        for(int c=0;c<by;c++)join_key(c,i,i%6,c?!string:string);
        numbers[by][i]=key[(i/6)%4];numbers[by+1][i]=100+i;
        char tag[24];snprintf(tag,sizeof(tag),"u%d",i);text_cell(by+2,i,tag);
    }
    types[by+2]=1;fail_alloc=allocation;fail_store=publication;
    char args[64];snprintf(args,sizeof(args),"using %d %d 0",by,by+1);
    int rc=crangejoin_main(args);if(rc)goto done;
    clear_data();resize_fixture(rows);int nm=by+1;nv=nm+5;
    types[nm+4]=1;
    for(int i=0;i<rows;i++) {
        numbers[0][i]=i+1;
        for(int c=0;c<by;c++)join_key(c+1,i,master_group(i),c?!string:string);
        interval(i,&numbers[nm][i],&numbers[nm+1][i]);
    }
    snprintf(args,sizeof(args),"prepare %d%s",nm,by==2?" 2 3":by==1?" 2":"");
    fail_macro=macro;rc=crangejoin_main(args);if(rc)goto done;
    size_t expected=0;
    for(int i=0;i<rows;i++) {
        double lo,hi;interval(i,&lo,&hi);int count=0;
        for(int u=0;u<using_rows;u++)if((!by||master_group(i)==u%6)&&lo<=key[(u/6)%4]&&key[(u/6)%4]<=hi)count++;
        expected+=count?count:1;
    }
    assert(strtoull(nout,NULL,10)==expected);
    resize_fixture((int)expected);rc=crangejoin_main("write");
    if(rc) { for(int i=0;i<rows;i++)check_master(i,i,by,string);goto done; }
    size_t out=0;const int ordered[]={1,3,5};
    for(int i=0;i<rows;i++) {
        double lo,hi;interval(i,&lo,&hi);int count=0;
        for(int k=0;k<3;k++)for(int u=0;u<using_rows;u++)
            if((!by||master_group(i)==u%6)&&key[(u/6)%4]==ordered[k]&&lo<=key[(u/6)%4]&&key[(u/6)%4]<=hi) {
                check_master(i,(int)out,by,string);
                assert(numbers[nm+2][out]==key[(u/6)%4]&&numbers[nm+3][out]==u+100);
                char tag[24];snprintf(tag,sizeof(tag),"u%d",u);assert(!strcmp(string_cell(nm+4,(int)out),tag));
                out++;count++;
            }
        if(!count) {
            check_master(i,(int)out,by,string);
            assert(numbers[nm+2][out]==mock.missval&&numbers[nm+3][out]==mock.missval);
            assert(!string_cell(nm+4,(int)out)[0]);out++;
        }
    }
    assert(out==expected);
#ifdef _OPENMP
    if(threads>1&&expected>10000)assert(worker_store);
#endif
done:
    last_allocations=allocations;last_stores=stores;
    crangejoin_cleanup_cache();assert(!live);clear_fixture();return rc;
}

int main(void)
{
    mock.nobs=nobs;mock.nvars=nvars;mock.nvar=nvars;mock.isstr=isstr;mock.isstrl=isstrl;
    mock.store=putnum;mock.safestore=putnum;mock.sstore=putstr;mock.macresave=macro_save;mock.macuse=macro_get;
    mock.missval=8.98846567431158e307;
    memset(longest,'x',2045);longest[2045]=0;
    for(threads=1;threads<=4;threads*=4) {
        for(int by=1;by<=2;by++)for(int string=0;string<=1;string++) {
            for(int cubic=0;cubic<=1;cubic++) {
                interpolation_cases(3,by,string,cubic);interpolation_cases(1001,by,string,cubic);
            }
        }
        for(int mode=0;mode<3;mode++)for(int filter=0;filter<=1;filter++) {
            int trim=mode==2?2:1;
            assert(!split_case(12001,mode,trim,4,filter,0,0,0));
            int total=last_stores;
            int points[]={1,total/2,total};
            for(int k=0;k<3;k++)assert(split_case(12001,mode,trim,4,filter,0,points[k],0)==459);
            assert(split_case(12001,mode,trim,4,filter,0,0,1)==459);
            assert(!split_case(12001,mode,trim,2,filter,0,0,0));
            if(mode<2)assert(!split_case(12001,mode,0,4,filter,0,0,0));
        }
        for(int by=0;by<=2;by++)for(int string=0;string<=(by!=0);string++) {
            assert(!join_case(12001,by,string,0,0,0,24));int total=last_stores;
            int points[]={1,total/2,total};
            for(int k=0;k<3;k++)assert(join_case(12001,by,string,0,points[k],0,24)==459);
            assert(join_case(12001,by,string,0,0,1,24)==459);
            assert(!join_case(12001,by,string,0,0,0,24));
        }
    }
    for(threads=1;threads<=4;threads*=4) {
        assert(!join_case(31,2,0,0,0,0,36024));
        assert(!join_case(31,2,1,0,0,0,36024));
    }
    for(int mode=0;mode<3;mode++) {
        assert(!split_case(8,mode,mode==2?2:1,4,0,0,0,0));int total=last_allocations;
        for(int f=1;f<=total;f++)assert(split_case(8,mode,mode==2?2:1,4,0,f,0,0)==920);
    }
    assert(!join_case(31,2,1,0,0,0,24));int total=last_allocations;
    for(int f=1;f<=total;f++)assert(join_case(31,2,1,f,0,0,24)==920);
    size_t lo[]={0,0,0},hi[]={INT_MAX-2,0,1},offset[4];
    assert(!rangejoin_offsets(lo,hi,3,offset)&&offset[3]==INT_MAX);
    hi[2]=2;assert(rangejoin_offsets(lo,hi,3,offset)==920);
    hi[0]=SIZE_MAX;assert(rangejoin_offsets(lo,hi,1,offset)==920);
    assert(csplit_main("write")==198&&crangejoin_main("write")==198&&!live);
    puts("PASS grouped/parallel interpolation, split, join, publication failures and output limits");
    return 0;
}
