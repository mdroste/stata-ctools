"""ASan/UBSan checks for lossless worksheet/style merge operations."""
import os
import platform
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = r'''
#include <assert.h>
#include <string.h>
#include "stplugin.h"
ST_plugin *_stata_;
#include "cexport/cexport_xlsx.c"
int main(void) {
    const char *old = "<x:worksheet xmlns:x=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\">"
        "<!-- > <row> --> <x:dimension ref='A1:F9'/><x:sheetData>"
        "<x:row r='2' ht='30' customHeight='1'><x:c r='A2' s='8'><x:f>SUM(A1)&gt;0</x:f><x:v>9</x:v></x:c>"
        "<x:c r='B2' s='7'><x:v>3</x:v></x:c><x:c r='F2'><x:is><x:t>old</x:t></x:is></x:c></x:row>"
        "<x:row r='9'><x:c r='F9'><x:v>12</x:v></x:c></x:row></x:sheetData>"
        "<x:mergeCells count='1'><x:mergeCell ref='G10:H10'/></x:mergeCells></x:worksheet>";
    const char *fresh="<worksheet><sheetData><row r='2'><c r='B2' s='1'><v>20</v></c><c r='C2'/></row>"
        "<row r='3'><c r='B3'><is><t>new &amp; text</t></is></c></row></sheetData></worksheet>";
    char *got=xedit_worksheet(old,fresh,1,20);
    assert(got && strstr(got,"s=\"7\"") && strstr(got,"s=\"20\""));
    assert(strstr(got,"<x:c r='A2' s='8'><x:f>SUM(A1)&gt;0</x:f><x:v>9</x:v></x:c>"));
    assert(strstr(got,"<x:c r='F9'><x:v>12</x:v></x:c>"));
    assert(strstr(got,"ht='30' customHeight='1'"));
    assert(strstr(got,"ref=\"A2:F9\"") && strstr(got,"ref='G10:H10'"));
    assert(!strstr(got,"<x:v>3</x:v>")); free(got);
    got=xedit_worksheet(old,fresh,0,20);assert(got && strstr(got,"s=\"21\""));free(got);
    got=xedit_worksheet("<worksheet><sheetData/></worksheet>",fresh,0,0);
    assert(got && strstr(got,"new &amp; text"));free(got);
    xedit_node n;
    const char *entities="<sheet name='A &amp; B &#x3a9;'/>";
    assert(xedit_find(entities,entities+strlen(entities),"sheet",&n));
    got=xedit_attr(&n,"name");assert(got && !strcmp(got,"A & B \xce\xa9"));free(got);
    const char *styles="<styleSheet><numFmts count='1'><numFmt numFmtId='181' formatCode='custom'/></numFmts>"
        "<cellXfs count='1'><xf numFmtId='181' fontId='2' fillId='3'/></cellXfs></styleSheet>";
    unsigned base=99;got=xedit_styles((char*)styles,&base);
    assert(got && base==1 && strstr(got,"numFmtId=\"182\"") && strstr(got,"count=\"4\""));
    assert(strstr(got,"<xf numFmtId='181' fontId='2' fillId='3'/>"));free(got);
    styles="<styleSheet><cellXfs count='0'/></styleSheet>";
    got=xedit_styles((char*)styles,&base);assert(got && base==0 && strstr(got,"count=\"3\""));free(got);
    got=xedit_translate_formula("SUM(A2,$A2,A$2,$A$2,'A2 sheet'!B3,LOG10(A2),\"A2\",Table1[A2],A2:B4)",1,2);
    assert(got && !strcmp(got,"SUM(B4,$A4,B$2,$A$2,'A2 sheet'!C5,LOG10(B4),\"A2\",Table1[A2],B4:C6)"));free(got);
    got=xedit_translate_formula("A1+$A1+A$1+$A$1",-1,-1);
    assert(got && !strcmp(got,"#REF!+#REF!+#REF!+$A$1"));free(got);
    got=xedit_translate_formula("SUM($1:2,A:$B,'A1:A3'!C2,A1:A3!B2)",1,2);
    assert(got && !strcmp(got,"SUM($1:2,A:$B,'A1:A3'!D4,A1:A3!C4)"));free(got);
    const char *shared="<worksheet><sheetData><row r='2'><c r='C2'><f t='shared' ref='C2:C4' si='0'>A2+B2</f><v>5</v></c></row>"
        "<row r='3'><c r='C3'><f t='shared' si='0'/><v>7</v></c></row><row r='4'><c r='C4'><f t='shared' si='0'/><v>9</v></c></row></sheetData></worksheet>";
    const char *overwritten="<worksheet><sheetData><row r='2'><c r='C2'><v>99</v></c></row>"
        "<row r='3'><c r='C3'><f t='shared' si='0'/><v>7</v></c></row><row r='4'><c r='C4'><f t='shared' si='0'/><v>9</v></c></row></sheetData></worksheet>";
    got=xedit_shared_anchors(shared,overwritten);
    assert(got && strstr(got,"<f t='shared' ref='C2:C4' si='0'>A3+B3</f>"));
    assert(strstr(got,"<v>99</v>") && strstr(got,"<c r='C4'><f t='shared' si='0'/>"));free(got);
    char chain[100];
    got=xedit_remove_chain("<Relationships><Relationship Id='one' Type='url/worksheet' Target='s.xml'/><Relationship Id='two' Type='url/calcChain' Target='calcChain.xml'/></Relationships>","Relationship","Type","/calcChain",chain,sizeof(chain));
    assert(got && !strcmp(chain,"xl/calcChain.xml") && strstr(got,"Id='one'") && !strstr(got,"Id='two'"));free(got);
    got=xedit_recalculate("<x:workbook xmlns:x='urn:x'><x:calcPr fullCalcOnLoad='0' calcId='42'/></x:workbook>");
    assert(got && strstr(got,"fullCalcOnLoad=\"1\"") && strstr(got,"calcId='42'"));free(got);
    got=xedit_recalculate("<workbook/>");
    assert(got && strstr(got,"fullCalcOnLoad=\"1\"") && strstr(got,"</workbook>"));free(got);
    const char *lex="SUM($1:2,A:$B,'A1:A3'!C2,A1:A3!B2,Table1[[A2]:[B3]],\"A2\"\"B3\",a2)";
    for(size_t i=0;i<=strlen(lex);i++) {
        char *truncated=malloc(i+1);memcpy(truncated,lex,i);truncated[i]=0;
        got=xedit_translate_formula(truncated,1,1);assert(got);free(got);free(truncated);
    }
    /* Truncated tags must be bounded: no reads past supplied lengths. */
    for(size_t i=0;i<strlen(old);i++) {
        char *truncated=malloc(i+1);memcpy(truncated,old,i);truncated[i]=0;
        xedit_find(truncated,truncated+i,"sheetData",&n);
        got=xedit_worksheet(truncated,fresh,1,0);free(got);free(truncated);
    }
    return 0;
}
'''

def main():
    with tempfile.TemporaryDirectory(prefix='ctools-xlsx-edit-') as tmp:
        directory = Path(tmp)
        source = directory / 'test.c'
        source.write_text(SOURCE)
        binary = directory / 'test'
        flags = ['-std=c11', '-O1', '-g', '-DSD_FASTMODE', '-fno-fast-math',
                 '-ffp-contract=off', '-fsanitize=address,undefined',
                 '-fno-omit-frame-pointer', '-pthread', '-I', str(ROOT / 'src')]
        if platform.system() == 'Darwin':
            flags += ['-DSYSTEM=APPLEMAC', '-D_DARWIN_C_SOURCE', '-Wl,-dead_strip',
                      '-Wl,-undefined,dynamic_lookup']
        else:
            flags += ['-DSYSTEM=STUNIX', '-D_GNU_SOURCE', '-ffunction-sections',
                      '-fdata-sections', '-Wl,--gc-sections', '-lm']
        subprocess.run([os.environ.get('CC', 'cc'), *flags, str(source),
                        str(ROOT / 'src/cimport/cimport_xlsx_xml.c'), '-o', str(binary)], check=True)
        subprocess.run([str(binary)], check=True, timeout=60,
                       env={**os.environ, 'UBSAN_OPTIONS': 'halt_on_error=1'})
        print('PASS XLSX edits preserve styles/formulas/XML, namespaces, truncated tags (ASan/UBSan)')

if __name__ == '__main__':
    main()
