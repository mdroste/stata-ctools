#include <stdio.h>
#include "xls.h"
int main(int argc,char **argv){xls_error_t e;xlsWorkBook*w=xls_open_file(argv[1],"UTF-8",&e);xlsWorkSheet*s=xls_getWorkSheet(w,0);xls_parseWorkSheet(s);for(int r=1;r<=5;r++)printf("%.17g\n",xls_cell(s,r,4)->d);xls_close_WS(s);xls_close_WB(w);}
