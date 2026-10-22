"""Small deterministic XLSX/SHP regression fixtures; no external data."""
from pathlib import Path
import struct
import zipfile


def create(folder):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    ns = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
    rns = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
    pns = 'http://schemas.openxmlformats.org/package/2006/relationships'
    book = f'<workbook xmlns="{ns}" xmlns:r="{rns}"><sheets><sheet name="Sheet1" sheetId="1" r:id="rId1"/></sheets></workbook>'
    rels = f'<Relationships xmlns="{pns}"><Relationship Id="rId1" Type="{rns}/worksheet" Target="worksheets/sheet1.xml"/><Relationship Id="rId2" Type="{rns}/sharedStrings" Target="sharedStrings.xml"/></Relationships>'
    base = {'xl/workbook.xml': book, 'xl/_rels/workbook.xml.rels': rels,
            '_rels/.rels': f'<Relationships xmlns="{pns}"><Relationship Id="rId1" Type="{rns}/officeDocument" Target="xl/workbook.xml"/></Relationships>',
            '[Content_Types].xml': '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/><Override PartName="/xl/sharedStrings.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sharedStrings+xml"/></Types>',
            'xl/sharedStrings.xml': f'<sst xmlns="{ns}"><si><t>first</t></si></sst>',
            'xl/styles.xml': f'<styleSheet xmlns="{ns}"><cellXfs count="1"><xf numFmtId="0"/></cellXfs></styleSheet>'}
    cell = '<c r="A1" t="s"><v>0</v></c>'
    def write(name, text=cell, **parts):
        content = dict(base)
        content['xl/worksheets/sheet1.xml'] = f'<worksheet xmlns="{ns}"><sheetData><row r="1">{text}</row></sheetData></worksheet>'
        content.update(parts)
        with zipfile.ZipFile(folder / (name+'.xlsx'), 'w', zipfile.ZIP_DEFLATED) as archive:
            for key, value in content.items(): archive.writestr(key,value)
    write('good')
    for name, part in [('sst','xl/sharedStrings.xml'),('styles','xl/styles.xml'),('book','xl/workbook.xml'),('rels','xl/_rels/workbook.xml.rels')]:
        write('truncated_'+name, **{part:base[part].rsplit('</',1)[0]})
    write('invalid_index', '<c r="A1" t="s"><v>1</v></c>')
    for name, text in [('cdata','<![CDATA[customer & name]]>'),('comment','before<!-- with > and < -->after'),('escaped','a_x000D_b')]:
        write(name, f'<c r="A1" t="inlineStr"><is><t>{text}</t></is></c>')
        write('shared_'+name, **{'xl/sharedStrings.xml':f'<sst xmlns="{ns}"><si><t>{text}</t></si></sst>'})
    records=b''.join(struct.pack('>II',i,10)+struct.pack('<Idd',1,float(i),float(i)*2) for i in [1,2])
    header=struct.pack('>7I',9994,0,0,0,0,0,(100+len(records))//2)+struct.pack('<2I8d',1000,1,1,2,2,4,0,0,0,0)
    (folder/'good.shp').write_bytes(header+records)
    (folder/'truncated.shp').write_bytes(header+records[:28])
