from pathlib import Path
import zipfile,struct
D=Path('/tmp/ctools-io-audit')
ns='http://schemas.openxmlformats.org/spreadsheetml/2006/main'
rns='http://schemas.openxmlformats.org/officeDocument/2006/relationships'
pns='http://schemas.openxmlformats.org/package/2006/relationships'
book=f'<workbook xmlns="{ns}" xmlns:r="{rns}"><sheets><sheet name="Sheet1" sheetId="1" r:id="rId1"/></sheets></workbook>'
rels=f'<Relationships xmlns="{pns}"><Relationship Id="rId1" Type="{rns}/worksheet" Target="worksheets/sheet1.xml"/><Relationship Id="rId2" Type="{rns}/sharedStrings" Target="sharedStrings.xml"/></Relationships>'
types='<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/><Override PartName="/xl/sharedStrings.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sharedStrings+xml"/></Types>'
def write(name, sheet, sst=None, workbook=book):
    with zipfile.ZipFile(D/name,'w',zipfile.ZIP_DEFLATED) as z:
        z.writestr('[Content_Types].xml',types)
        z.writestr('_rels/.rels',f'<Relationships xmlns="{pns}"><Relationship Id="rId1" Type="{rns}/officeDocument" Target="xl/workbook.xml"/></Relationships>')
        z.writestr('xl/workbook.xml',workbook)
        z.writestr('xl/_rels/workbook.xml.rels',rels)
        z.writestr('xl/worksheets/sheet1.xml',f'<worksheet xmlns="{ns}"><sheetData>{sheet}</sheetData></worksheet>')
        if sst is not None:z.writestr('xl/sharedStrings.xml',sst)
sheet='<row r="1"><c r="A1" t="s"><v>0</v></c><c r="B1"><v>1</v></c></row><row r="2"><c r="A2" t="s"><v>1</v></c><c r="B2"><v>2</v></c></row>'
write('truncated_sst.xlsx',sheet,f'<sst xmlns="{ns}"><si><t>first</t></si><si><t>second')
write('good_sst.xlsx',sheet,f'<sst xmlns="{ns}"><si><t>first</t></si><si><t>second</t></si></sst>')
for name,text in [('cdata','<![CDATA[customer & name]]>'),('comment','before<!-- comment -->after'),('escaped','a_x000D_b')]:
    write(name+'.xlsx',f'<row r="1"><c r="A1" t="inlineStr"><is><t>{text}</t></is></c><c r="B1"><v>1</v></c></row>')
records=b''.join(struct.pack('>II',i,10)+struct.pack('<Idd',1,float(i),float(i)*2) for i in [1,2])
head=struct.pack('>7I',9994,0,0,0,0,0,(100+len(records))//2)+struct.pack('<2I8d',1000,1,1,2,2,4,0,0,0,0)
(D/'good.shp').write_bytes(head+records)
(D/'truncated.shp').write_bytes(head+records[:28])
