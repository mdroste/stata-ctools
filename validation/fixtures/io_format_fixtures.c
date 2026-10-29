/* Independent format fixtures made with the pinned codec's public writer API.
   The native differential test compares its reader with the C integration. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "readstat.h"
static ssize_t output(const void *data,size_t n,void *context) {
    size_t wrote=fwrite(data,1,n,context);return wrote==n?(ssize_t)n:-1;
}
static void checked(readstat_error_t error) {
    if(error!=READSTAT_OK){fprintf(stderr,"%s\n",readstat_error_message(error));exit(1);}
}
static void dataset(const char *path,int spss) {
    FILE *file=fopen(path,"wb");if(!file)exit(1);
    readstat_writer_t *w=readstat_writer_init();readstat_set_data_writer(w,output);
    readstat_writer_set_file_label(w,"External codec fixture");
    if(spss)checked(readstat_writer_set_compression(w,READSTAT_COMPRESS_BINARY));
    readstat_variable_t *id=readstat_add_variable(w,"id",READSTAT_TYPE_DOUBLE,0),*value=readstat_add_variable(w,"value",READSTAT_TYPE_DOUBLE,0),*text=readstat_add_variable(w,"text",READSTAT_TYPE_STRING,2100);
    readstat_variable_set_format(value,spss?"F12.2":"BEST12");
    readstat_variable_set_label(value,"Exact numeric fixture");
    readstat_label_set_t *labels=readstat_add_label_set(w,READSTAT_TYPE_DOUBLE,"status");
    readstat_label_double_value(labels,1,"Domestic");readstat_label_double_value(labels,2,"Foreign");
    readstat_variable_set_label_set(id,labels);
    checked(spss?readstat_begin_writing_sav(w,file,3):readstat_begin_writing_sas7bdat(w,file,3));
    char longtext[2101];memset(longtext,'a',2100);longtext[2100]=0;
    const char *strings[]={longtext," leading ",""};
    for(int i=0;i<3;i++) {
        checked(readstat_begin_row(w));checked(readstat_insert_double_value(w,id,i+1));
        checked(i==2?readstat_insert_missing_value(w,value):readstat_insert_double_value(w,value,1.1234567890123456+i));
        checked(readstat_insert_string_value(w,text,strings[i]));checked(readstat_end_row(w));
    }
    checked(readstat_end_writing(w));readstat_writer_free(w);if(fclose(file))exit(1);
}
static void catalog(const char *path) {
    FILE *file=fopen(path,"wb");if(!file)exit(1);
    readstat_writer_t *w=readstat_writer_init();readstat_set_data_writer(w,output);
    readstat_label_set_t *labels=readstat_add_label_set(w,READSTAT_TYPE_DOUBLE,"status");
    readstat_label_double_value(labels,1,"Domestic");readstat_label_double_value(labels,2,"Foreign");
    checked(readstat_begin_writing_sas7bcat(w,file));checked(readstat_end_writing(w));
    readstat_writer_free(w);if(fclose(file))exit(1);
}
static void formats(const char *path,int spss) {
    const char *sav[]={"F8.2","COMMA10.2","DOT10.2","DOLLAR10.2","E12.2","PCT10.2","N8.2","Z8.2","DATE11","ADATE10","EDATE10","SDATE10","JDATE7","DATETIME23.2","TIME12.2","DTIME12.2","WKDAY8","MONTH8","MOYR8","QYR8","WKYR10"};
    const char *sas[]={"BEST12","DATE","DDMMYY","YYMMDD","JULIAN","MONYY","YEAR","WEEKDATE","TIME","DATETIME","E8601DA","E8601DT","DATEAMPM","B8601DT","TOD","HHMM"};
    const char **fmts=spss?sav:sas;int k=spss?21:16;
    FILE *file=fopen(path,"wb");if(!file)exit(1);
    readstat_writer_t *w=readstat_writer_init();readstat_set_data_writer(w,output);
    readstat_variable_t *v[24];char name[20];
    for(int j=0;j<k;j++) {
        snprintf(name,sizeof(name),"v%d",j+1);v[j]=readstat_add_variable(w,name,READSTAT_TYPE_DOUBLE,0);
        readstat_variable_set_format(v[j],fmts[j]);
    }
    v[k]=readstat_add_variable(w,"missing",READSTAT_TYPE_DOUBLE,0);
    v[k+1]=readstat_add_variable(w,"strmiss",READSTAT_TYPE_STRING,4);
    if(spss) {
        checked(readstat_variable_add_missing_double_value(v[k],99));
        checked(readstat_variable_add_missing_double_range(v[k],999,1001));
        checked(readstat_variable_add_missing_string_value(v[k+1],"NA"));
    }
    checked(spss?readstat_begin_writing_sav(w,file,5):readstat_begin_writing_sas7bdat(w,file,5));
    double values[]={12.5,-0.5,1234567890,0,999};double missing[]={99,999,1000,1001,1};
    for(int i=0;i<5;i++) {
        checked(readstat_begin_row(w));
        for(int j=0;j<k;j++) checked(readstat_insert_double_value(w,v[j],values[i]));
        checked(readstat_insert_double_value(w,v[k],missing[i]));
        checked(readstat_insert_string_value(w,v[k+1],i==0?"NA":"abc"));
        checked(readstat_end_row(w));
    }
    checked(readstat_end_writing(w));readstat_writer_free(w);if(fclose(file))exit(1);
}
static void companion(const char *path) {
    FILE *file=fopen(path,"wb");if(!file)exit(1);
    readstat_writer_t *w=readstat_writer_init();readstat_set_data_writer(w,output);
    checked(readstat_writer_set_file_format_version(w,8));
    readstat_variable_t *id=readstat_add_variable(w,"id",READSTAT_TYPE_DOUBLE,0);
    readstat_variable_t *text=readstat_add_variable(w,"text",READSTAT_TYPE_STRING,12);
    checked(readstat_begin_writing_xport(w,file,3));
    for(int i=0;i<3;i++) {
        checked(readstat_begin_row(w));checked(readstat_insert_double_value(w,id,i+1));
        checked(readstat_insert_string_value(w,text,i==0?"abc":""));checked(readstat_end_row(w));
    }
    checked(readstat_end_writing(w));readstat_writer_free(w);if(fclose(file))exit(1);
}
int main(int argc,char **argv) {
    if(argc!=2)return 1;char path[4096];
    snprintf(path,sizeof(path),"%s/fixture.zsav",argv[1]);dataset(path,1);
    snprintf(path,sizeof(path),"%s/long.sas7bdat",argv[1]);dataset(path,0);
    snprintf(path,sizeof(path),"%s/labels.sas7bcat",argv[1]);catalog(path);
    snprintf(path,sizeof(path),"%s/formats.sav",argv[1]);formats(path,1);
    snprintf(path,sizeof(path),"%s/formats.sas7bdat",argv[1]);formats(path,0);
    snprintf(path,sizeof(path),"%s/raw.v8xpt",argv[1]);companion(path);
    return 0;
}
