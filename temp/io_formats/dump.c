#include <stdio.h>
#include "readstat.h"
int var(int i,readstat_variable_t *v,const char *l,void *c) {printf("VAR %d %s %s %zu %s\n",i,v->name,v->format,v->storage_width,l?l:"");return 0;}
int val(int i,readstat_variable_t *v,readstat_value_t x,void *c) {if (i<5) {if (readstat_value_type_class(x)==READSTAT_TYPE_CLASS_STRING) printf("VAL %d %s [%s]\n",i,v->name,readstat_string_value(x));else printf("VAL %d %s %.17g %.17g %.17g\n",i,v->name,readstat_double_value(x),(readstat_double_value(x)-11903760000.0)*1000,readstat_double_value(x)*1000-11903760000000.0);}return 0;}
int main(int n,char **a) {readstat_parser_t *p=readstat_parser_init();readstat_set_variable_handler(p,var);readstat_set_value_handler(p,val);readstat_error_t e;if(*a[1]=='s')e=readstat_parse_sav(p,a[2],NULL);else if(*a[1]=='x')e=readstat_parse_xport(p,a[2],NULL);else e=readstat_parse_sas7bdat(p,a[2],NULL);printf("RC %d\n",e);readstat_parser_free(p);}
