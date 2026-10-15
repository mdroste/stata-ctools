
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
