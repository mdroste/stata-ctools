"""Independent Excel format fixtures for native/C differential readers."""
from pathlib import Path
from xml.sax.saxutils import escape
from zipfile import ZipFile, ZIP_DEFLATED
import struct


def fixtures(directory: Path):
    builtin = list(range(14, 23)) + list(range(27, 37)) + list(range(45, 48)) + list(range(50, 59))
    custom = ['yyyy-mm-dd', 'yyyy-mm-dd hh:mm:ss', '[h]:mm:ss', 'hh:mm:ss', '0.00 "days"']
    formats = builtin + list(range(165, 165+len(custom)))
    def column(index):
        out=''
        while index:
            index, digit = divmod(index-1,26);out=chr(65+digit)+out
        return out
    ns='http://schemas.openxmlformats.org/spreadsheetml/2006/main'
    types='http://schemas.openxmlformats.org/package/2006/content-types'
    parts={
        '[Content_Types].xml':f'<Types xmlns="{types}"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/><Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/></Types>',
        '_rels/.rels':'<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>',
        'xl/workbook.xml':f'<workbook xmlns="{ns}" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="Formats" sheetId="1" r:id="rId1"/></sheets></workbook>',
        'xl/_rels/workbook.xml.rels':'<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/><Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/></Relationships>',
    }
    numfmts=''.join(f'<numFmt numFmtId="{165+i}" formatCode="{escape(code, {chr(34):"&quot;"})}"/>' for i,code in enumerate(custom))
    xfs='<xf numFmtId="0" fontId="0" fillId="0" borderId="0"/>'+''.join(f'<xf numFmtId="{fmt}" fontId="0" fillId="0" borderId="0"/>' for fmt in formats)
    parts['xl/styles.xml']=f'<styleSheet xmlns="{ns}"><numFmts count="{len(custom)}">{numfmts}</numFmts><fonts count="1"><font><sz val="11"/><name val="Calibri"/></font></fonts><fills count="2"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill></fills><borders count="1"><border/></borders><cellStyleXfs count="1"><xf numFmtId="0"/></cellStyleXfs><cellXfs count="{len(formats)+1}">{xfs}</cellXfs></styleSheet>'
    rows=['<row r="1">'+''.join(f'<c r="{column(j+1)}1" t="inlineStr"><is><t>fmt{fmt}</t></is></c>' for j,fmt in enumerate(formats))+'</row>']
    for row,value in [(2,'44300.5123456'),(3,'0.5123456'),(4,'44300.125'),(5,'44197.001')]:
        rows.append(f'<row r="{row}">'+''.join(f'<c r="{column(j+1)}{row}" s="{j+1}"><v>{value}</v></c>' for j in range(len(formats)))+'</row>')
    parts['xl/worksheets/sheet1.xml']=f'<worksheet xmlns="{ns}"><dimension ref="A1:{column(len(formats))}5"/><sheetData>'+''.join(rows)+'</sheetData></worksheet>'
    with ZipFile(directory/'display_formats.xlsx','w',compression=ZIP_DEFLATED) as book:
        for name,text in parts.items():book.writestr(name,text)

    # Independent BIFF8/OLE container with the same values and format IDs.
    def record(kind, body=b''):
        return struct.pack('<HH',kind,len(body))+body
    def bof(kind):
        return record(0x809,struct.pack('<4H2I',0x600,kind,0xdbb,1997,0x41,6))
    sheet=bof(0x10)+record(0x200,struct.pack('<IIHHH',0,5,0,len(formats),0))
    for j,fmt in enumerate(formats):
        sheet+=record(0xfd,struct.pack('<HHHI',0,j,0,j))
    for row,value in [(1,44300.5123456),(2,0.5123456),(3,44300.125),(4,44197.001)]:
        for j in range(len(formats)):
            sheet+=record(0x203,struct.pack('<HHHd',row,j,j+1,value))
    sheet+=record(0x0a)
    globals_=bof(5)+record(0x42,struct.pack('<H',1200))+record(0x22,struct.pack('<H',0))
    for j,code in enumerate(custom):
        globals_+=record(0x41e,struct.pack('<HHB',165+j,len(code),1)+code.encode('utf-16le'))
    for fmt in [0]+formats:
        globals_+=record(0xe0,struct.pack('<HH',0,fmt)+bytes(16))
    strings=b''.join(struct.pack('<HB',len('fmt'+str(fmt)),0)+('fmt'+str(fmt)).encode('ascii') for fmt in formats)
    globals_+=record(0xfc,struct.pack('<II',len(formats),len(formats))+strings)
    bound=record(0x85,struct.pack('<IBBBB',0,0,0,len('Formats'),0)+b'Formats')
    offset=len(globals_)+len(bound)+4
    bound=record(0x85,struct.pack('<IBBBB',offset,0,0,len('Formats'),0)+b'Formats')
    stream=globals_+bound+record(0x0a)+sheet
    size=max(4096,(len(stream)+511)//512*512);n=size//512
    header=bytearray(512);header[:8]=bytes.fromhex('d0cf11e0a1b11ae1')
    struct.pack_into('<HHHHH',header,24,0x3e,3,0xfffe,9,6)
    struct.pack_into('<IIIII',header,40,0,1,n,0,4096)
    struct.pack_into('<4I',header,60,0xfffffffe,0,0xfffffffe,0)
    struct.pack_into('<109I',header,76,n+1,*([0xffffffff]*108))
    def entry(name,kind,start,length,child=0xffffffff):
        out=bytearray(128);encoded=(name+'\0').encode('utf-16le');out[:len(encoded)]=encoded
        struct.pack_into('<HBBIII',out,64,len(encoded),kind,1,0xffffffff,0xffffffff,child)
        struct.pack_into('<IQ',out,116,start,length);return out
    ole_directory=entry('Root Entry',5,0xfffffffe,0,1)+entry('Workbook',2,0,size)+bytes(256)
    fat=[i+1 for i in range(n)];fat[-1]=0xfffffffe
    fat += [0xfffffffe,0xfffffffd]+[0xffffffff]*(128-n-2)
    (directory/'display_formats.xls').write_bytes(bytes(header)+stream.ljust(size,b'\0')+ole_directory+struct.pack('<128I',*fat))
