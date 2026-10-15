
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
