
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stplugin.h"
static ST_plugin mock;ST_plugin *_stata_=&mock;
#include "cexport/cexport_parse.c"
static void put(const char *p,const char *s){FILE *f=fopen(p,"w");assert(f);assert(fputs(s,f)>=0);assert(!fclose(f));}
static void check(const char *p,const char *s){char buf[30]={0};FILE *f=fopen(p,"r");assert(f);fread(buf,1,29,f);fclose(f);assert(!strcmp(buf,s));}
int main(int argc,char **argv){
 assert(argc==2);ctools_output out;const char *p=argv[1];put(p,"original");
 assert(!ctools_output_prepare(&out,p,1));put(out.temporary,"partial");ctools_output_cleanup(&out);check(p,"original");
 assert(!ctools_output_prepare(&out,p,0));put(out.temporary,"new");assert(ctools_output_commit(&out)==602);ctools_output_cleanup(&out);check(p,"original");
 assert(!ctools_output_prepare(&out,p,1));put(out.temporary,"complete");assert(!ctools_output_commit(&out));ctools_output_cleanup(&out);check(p,"complete");
 remove(p);return 0;
}
