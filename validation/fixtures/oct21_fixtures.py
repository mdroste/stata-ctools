"""Small corrupt-container and geometry fixtures; no external data needed."""
from pathlib import Path
import struct
import zipfile

def fixtures(p):
    p.mkdir(exist_ok=True)
    ns = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
    rel = 'http://schemas.openxmlformats.org/package/2006/relationships'
    orel = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
    files = {
        '[Content_Types].xml': (
            '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
            '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
            '<Default Extension="xml" ContentType="application/xml"/>'
            '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
            '<Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
            '</Types>'
        ),
        '_rels/.rels': (
            f'<Relationships xmlns="{rel}"><Relationship Id="rId1" '
            f'Type="{orel}/officeDocument" Target="xl/workbook.xml"/></Relationships>'
        ),
        'xl/workbook.xml': (
            f'<workbook xmlns="{ns}" xmlns:r="{orel}"><sheets>'
            '<sheet name="Sheet1" sheetId="1" r:id="rId1"/></sheets></workbook>'
        ),
        'xl/_rels/workbook.xml.rels': (
            f'<Relationships xmlns="{rel}"><Relationship Id="rId1" '
            f'Type="{orel}/worksheet" Target="worksheets/sheet1.xml"/></Relationships>'
        ),
        'xl/worksheets/sheet1.xml': (
            f'<worksheet xmlns="{ns}"><dimension ref="A1:A2"/><sheetData>'
            '<row r="1"><c r="A1"><v>123</v></c></row>'
            '<row r="2"><c r="A2"><v>456</v></c></row>'
            '</sheetData></worksheet>'
        ),
    }

    def write_xlsx(name, parts, level=6):
        with zipfile.ZipFile(p / name, 'w', compression=zipfile.ZIP_DEFLATED,
                             compresslevel=level) as z:
            for key, value in parts.items():
                z.writestr(key, value)

    write_xlsx('valid.xlsx', files)
    write_xlsx('bitflip.xlsx', files, level=0)
    b = bytearray((p / 'bitflip.xlsx').read_bytes())
    b[b.index(b'<v>123</v>') + 3] ^= 8
    (p / 'bitflip.xlsx').write_bytes(b)
    try:
        with zipfile.ZipFile(p / 'bitflip.xlsx') as z:
            z.read('xl/worksheets/sheet1.xml')
    except zipfile.BadZipFile as error:
        print(error)  # Expected CRC failure

    truncated = dict(files)
    truncated['xl/worksheets/sheet1.xml'] = (
        files['xl/worksheets/sheet1.xml'].split('<row r="2">')[0]
        + '<row r="2"><c r="A2"><v>45'
    )
    write_xlsx('truncated-xml.xlsx', truncated)
    with zipfile.ZipFile(p / 'truncated-xml.xlsx') as z:
        assert z.testzip() is None  # ZIP is intact; contained XML is incomplete

    h = bytearray(100)
    struct.pack_into('>I', h, 0, 9994)
    struct.pack_into('>I', h, 24, 64)
    struct.pack_into('<II', h, 28, 1000, 1)
    struct.pack_into('<dddd', h, 36, 1, 2, 1, 2)
    (p / 'point.shp').write_bytes(
        h + struct.pack('>II', 1, 10) + struct.pack('<Idd', 1, 1, 2)
    )
    (p / 'partial.shx').write_text('old index sentinel')
    (p / 'corrupt.dbf').write_bytes(b'not a dbf')
    (p / 'transport.sas').write_text('old SAS sentinel')
