"""Independent input fixtures and structural output checks for Stata and native tests."""

import os
import platform
import struct
import subprocess
import tempfile
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET
from xml.sax.saxutils import escape
from zipfile import ZIP_DEFLATED, ZipFile

from validation import ROOT


def multichunk():
    """Over 2 MB, so cimport parses it in parallel chunks: quoted delimiters,
    padded and leading-zero numbers, exponents, ragged rows, and a column
    that becomes text only near the end of the file."""
    rows = ["id,num,padded,lead0,sci,text,late,ragged"]
    for i in range(1, 60001):
        text = '"has, comma"' if i % 7 == 0 else '"q ""x"" y"' if i % 11 == 0 else f"plain{i}"
        # id stays in int range: native's %12.0g-for-long heuristic is a separate gap.
        row = [
            str(i % 30000 + 1),
            f"{(i * 7919 % 100003) / 97:.6f}",
            f" {i % 500} ",
            "000" + str(i % 1000),
            f"{i % 97 + 1}.5e-{i % 30 + 1}",
            text,
            "word" if i == 59995 else str(i % 1000),
            str(i % 3),
        ]
        if i % 1000 == 0:
            row = row[:5]
        elif i % 1500 == 0:
            row.append("extra")
        rows.append(",".join(row))
    return "\n".join(rows) + "\n"


def parity_fixtures(p):
    p.mkdir(parents=True, exist_ok=True)
    cases = {
        "headers": "First Name,1,,a-b,a b,a_b,UPPER,if,_all,über\n1,2,3,4,5,6,7,8,9,10\n11,12,13,14,15,16,17,18,19,20\n",
        "noheaders": "1,2,abc\n3,4,def\n",
        "numbers": "small,decimal,precise,missing,na,empty\n1,1.1,1.12345678901234,.a,NA,\n2,1.2,2.12345678901234,.z,1,\n",
        "quoted": 'value,text\n"1.25","hello"\n"2.5","a""b"\n',
        "exponents": "value\n1e-2\n2e-3\n",
        "outlier": "id,value\n"
        + "".join(f"{i},{1000000000 if i == 10001 else 1}\n" for i in range(1, 20002)),
        # Whitespace-only fields are text; padded numbers and missings are numeric
        # (padding turns .a into system missing).
        "whitespace": "sp1,lead,trail,both,sp2,tab,qsp,qempty,spa,atr,xa,qa,trailtab\n"
        "1,1,1,1,1,1,1,1,1,1,1,1,1\n"
        ' , 2,2 , 2 ,  ,\t," ",""," .a",.a ,.a,".a",2\t\n'
        "3,3,3,3,3,3,3,3,3,3,3,3,3\n",
        "digits": "lead,sci\n0000000000000000007.9770,2.59462e-18\n00000000000000.0851289740626,7.7651552348e-14\n"
        "1167.411680745998,-4.356865436322338e-21\n0.00000000000000000012,3.474157e-30\n",
        "multichunk": multichunk(),
    }
    for name, data in cases.items():
        (p / (name + ".csv")).write_text(data)
    return p


def shapefiles(directory):
    """Small ESRI records, built independently of either reader/writer."""

    def write(name, kind, bodies):
        records = b"".join(
            struct.pack(">II", i + 1, len(body) // 2) + body for i, body in enumerate(bodies)
        )
        header = struct.pack(">7I", 9994, 0, 0, 0, 0, 0, (100 + len(records)) // 2)
        header += struct.pack("<2I8d", 1000, kind, 1, 2, 3, 4, 7, 8, 9, 10)
        (directory / (name + ".shp")).write_bytes(header + records)

    xy = struct.pack("<4d", 1, 2, 3, 4)
    for kind in (0, 1, 11, 21):
        body = struct.pack("<I", kind)
        if kind:
            body += struct.pack("<2d", 1, 2)
            if kind == 11:
                body += struct.pack("<d", 7)
            if kind in (11, 21):
                body += struct.pack("<d", 9)
        write("shape" + str(kind), kind, [body])
    for kind in (3, 5, 8, 13, 15, 18, 23, 25, 28, 31):
        multi = kind in (8, 18, 28)
        body = struct.pack("<I4d", kind, 1, 2, 3, 4)
        body += struct.pack("<I", 2) if multi else struct.pack("<3I", 1, 2, 0)
        if kind == 31:
            body += struct.pack("<I", 2)
        body += xy
        if kind in (13, 15, 18, 31):
            body += struct.pack("<4d", 7, 8, 7, 8)
        if kind in (13, 15, 18, 23, 25, 28, 31):
            body += struct.pack("<4d", 9, 10, 9, 10)
        write("shape" + str(kind), kind, [body])
    body = struct.pack("<I4d4I", 3, 1, 2, 3, 4, 2, 3, 0, 1)
    body += struct.pack("<6d", 1, 2, 3, 4, 2, 3)
    write("multipart", 3, [body])
    write("nullrecord", 1, [struct.pack("<I", 0), struct.pack("<I2d", 1, 1, 2)])
    (directory / "corrupt.sav").write_bytes(b"$FL2" + b"\0" * 30)
    (directory / "corrupt.shp").write_bytes(b"\0" * 100)
    (directory / "corrupt.dbf").write_bytes(b"\x03\0")
    (directory / "corrupt.xls").write_bytes(b"\xd0\xcf\x11\xe0")
    (directory / "memory.csv").write_text("id,text\n1,first\n2,second\n")


def codec_fixtures(directory):
    vendor = ROOT / "src/vendor"
    sources = [
        *sorted((vendor / "readstat").rglob("*.c")),
        *sorted((vendor / "zlib").glob("*.c")),
        # vendor/zlib/zlib.h maps ReadStat's zlib calls onto miniz.
        *sorted((ROOT / "src/vendor/miniz").glob("*.c")),
    ]
    with tempfile.TemporaryDirectory(prefix="ctools-format-fixtures-") as tmp:
        binary = Path(tmp) / "fixtures"
        command = [
            os.environ.get("CC", "cc"),
            "-std=c11",
            "-O2",
            "-DHAVE_ZLIB=1",
            "-D_GNU_SOURCE",
            "-D_DARWIN_C_SOURCE",
            "-I",
            str(vendor / "readstat"),
            "-I",
            str(vendor / "zlib"),
            str(ROOT / "validation/fixtures/io_format_fixtures.c"),
            *map(str, sources),
            "-lm",
        ]
        if platform.system() == "Darwin":
            command += ["-liconv"]
        subprocess.run([*command, "-o", str(binary)], check=True)
        subprocess.run([str(binary), str(directory)], check=True)


def failure_fixtures(directory):
    headers = {
        "sas.sas7bdat": (directory / "long.sas7bdat").read_bytes()[:16],
        "spss.sav": b"$FL2" + bytes(12),
        "sasxport5.xpt": b"HEADER RECORD***",
        "sasxport8.v8xpt": b"HEADER RECORD***",
        "dbase.dbf": b"\x03" + bytes(15),
        "xls.xls": bytes.fromhex("d0cf11e0a1b11ae1") + bytes(8),
        "xlsx.xlsx": b"PK\x03\x04" + bytes(12),
        "shp.shp": struct.pack(">I", 9994) + bytes(12),
    }
    for name, header in headers.items():
        kind, ext = name.split(".")
        (directory / f"{kind}_invalid.{ext}").write_bytes(b"invalid file header\n")
        (directory / f"{kind}_truncated.{ext}").write_bytes(header)
        (directory / f"{kind}_longinvalid.{ext}").write_bytes(b"not valid\0" * 200)


def delimited_fixtures(directory):
    text = "ID,Text\n1,café Ω 日本語\n2,é\n"
    for name, encoding in [
        ("utf32le", "utf-32le"),
        ("utf32be", "utf-32be"),
        ("utf16be", "utf-16be"),
    ]:
        (directory / (name + ".csv")).write_bytes(text.encode(encoding))
    (directory / "utf32bom.csv").write_bytes(text.encode("utf-32"))
    for name, encoding, value in [
        ("cp1251", "cp1251", "Привет"),
        ("sjis", "shift_jis", "日本語"),
        ("latin2", "iso8859_2", "Łódź"),
    ]:
        (directory / (name + ".csv")).write_bytes(
            ("ID,Text\n1," + value + "\n2,abc\n").encode(encoding)
        )
    (directory / "multidelim.txt").write_text('1||2;3\n4;5||6\n7;"8|9";10\n')
    (directory / "multiline.csv").write_text('id,text\n1,"a\nb\nc"\n2,d\n')
    (directory / "range.csv").write_text("id,value\n1,textual\n2,1.5\n3,2.25\n4,99\n")
    (directory / "empty.csv").write_bytes(b"")
    (directory / "words.txt").write_text("idtabvalue\n1tab2\n3tab4\n")
    for name, separator in [
        ("tab", "\t"),
        ("semi", ";"),
        ("pipe", "|"),
        ("colon", ":"),
        ("space", " "),
    ]:
        (directory / ("autodelim_" + name + ".csv")).write_text(
            separator.join(("id", "value"))
            + "\n"
            + separator.join(("1", "2"))
            + "\n"
            + separator.join(("3", "4"))
            + "\n"
        )
    for width in (20, 40, 50, 100, 200, 1000):
        for pct in (1, 28, 29, 34, 35, 50, 63, 64, 65, 66, 67, 70, 90, 93, 94, 95):
            (directory / f"storagew{width}p{pct}.csv").write_text(
                "ID,Text\n"
                + "".join(
                    f"{i}," + ("x" * (width - 3) + f"{i:03d}" if i <= pct else f"s{i}") + "\n"
                    for i in range(1, 101)
                )
            )
    numbers = {
        "us": "1,234.56",
        "german": "1.234,56",
        "french": "1\u202f234,56",
        "french_nbsp": "1\u00a0234,56",
        "russian": "1\u00a0234,56",
        "swiss": "1’234.56",
        "arabic": "١٬٢٣٤٫٥٦",
        "arabic_minus": "\u061c-١٢٣٤٫٥٦",
        "persian": "۱۲۳۴٫۵۶",
        "persian_minus": "\u200e−۱۲۳۴٫۵۶",
        "hindi": "१२३४.५६",
        "bad_sign": "+1,234.56",
        "scientific": "1.23e3",
        "infinity": "∞",
        "bad_groups": "1.2,3",
        "groups": "1,,23",
        "trailing_group": "123,",
    }
    for name, value in numbers.items():
        (directory / (name + ".csv")).write_text("id;value\n1;" + value + "\n2;1\n")
    (directory / "quote_numbers.csv").write_text('id,value\n1,"123"\n2,1""2\n3,12"\n4,"34\n')
    (directory / "quotes.csv").write_text('id,text\n1,a"b"c\n2,"a""b"\n3,c"\n4,"d\n')
    (directory / "long.csv").write_text(
        'id,text\n1,"' + ('éΩ""' * 3000) + '"\n2,"' + ("x" * 5000) + '"\n'
    )
    regional = [
        ("western", "cp1252", "La France et les États-Unis ont une économie européenne. "),
        ("central", "iso8859_2", "Polska gospodarka i społeczeństwo rozwijają się. "),
        ("cyrillic", "cp1251", "Российская экономика и международная торговля. "),
        ("koi8", "koi8_r", "Российская экономика и международная торговля. "),
        ("arabic", "cp1256", "الاقتصاد الدولي والتجارة العربية في العالم. "),
        ("greek", "iso8859_7", "Η οικονομία και το διεθνές εμπόριο στην Ελλάδα. "),
        ("hebrew", "iso8859_8", "כלכלה בינלאומית ומסחר בישראל. "),
        ("turkish", "iso8859_9", "Türkiye ekonomisi ve uluslararası ticaret. "),
        ("japanese", "shift_jis", "日本の経済と世界の貿易についての研究。"),
        ("eucjp", "euc_jp", "日本の経済と世界の貿易についての研究。"),
        ("korean", "euc_kr", "한국의 경제와 국제 무역에 대한 연구. "),
        ("chinese", "gb18030", "中国经济与国际贸易的研究。"),
        ("traditional", "big5", "國際經濟與貿易的研究。"),
        ("iso2022", "iso2022_jp", "日本の経済と世界の貿易についての研究。"),
    ]
    for name, encoding, value in regional:
        (directory / ("detect_" + name + ".csv")).write_bytes(
            ('id,text\n1,"' + value * 100 + '"\n2,abc\n').encode(encoding)
        )


def excel_fixtures(directory: Path):
    builtin = list(range(14, 23)) + list(range(27, 37)) + list(range(45, 48)) + list(range(50, 59))
    custom = ["yyyy-mm-dd", "yyyy-mm-dd hh:mm:ss", "[h]:mm:ss", "hh:mm:ss", '0.00 "days"']
    formats = builtin + list(range(165, 165 + len(custom)))

    def column(index):
        out = ""
        while index:
            index, digit = divmod(index - 1, 26)
            out = chr(65 + digit) + out
        return out

    ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    types = "http://schemas.openxmlformats.org/package/2006/content-types"
    parts = {
        "[Content_Types].xml": f'<Types xmlns="{types}"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/><Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/></Types>',
        "_rels/.rels": '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>',
        "xl/workbook.xml": f'<workbook xmlns="{ns}" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="Formats" sheetId="1" r:id="rId1"/></sheets></workbook>',
        "xl/_rels/workbook.xml.rels": '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/><Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/></Relationships>',
    }
    numfmts = "".join(
        f'<numFmt numFmtId="{165+i}" formatCode="{escape(code, {chr(34):"&quot;"})}"/>'
        for i, code in enumerate(custom)
    )
    xfs = '<xf numFmtId="0" fontId="0" fillId="0" borderId="0"/>' + "".join(
        f'<xf numFmtId="{fmt}" fontId="0" fillId="0" borderId="0"/>' for fmt in formats
    )
    parts["xl/styles.xml"] = (
        f'<styleSheet xmlns="{ns}"><numFmts count="{len(custom)}">{numfmts}</numFmts><fonts count="1"><font><sz val="11"/><name val="Calibri"/></font></fonts><fills count="2"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill></fills><borders count="1"><border/></borders><cellStyleXfs count="1"><xf numFmtId="0"/></cellStyleXfs><cellXfs count="{len(formats)+1}">{xfs}</cellXfs></styleSheet>'
    )
    rows = [
        '<row r="1">'
        + "".join(
            f'<c r="{column(j+1)}1" t="inlineStr"><is><t>fmt{fmt}</t></is></c>'
            for j, fmt in enumerate(formats)
        )
        + "</row>"
    ]
    for row, value in [(2, "44300.5123456"), (3, "0.5123456"), (4, "44300.125"), (5, "44197.001")]:
        rows.append(
            f'<row r="{row}">'
            + "".join(
                f'<c r="{column(j+1)}{row}" s="{j+1}"><v>{value}</v></c>'
                for j in range(len(formats))
            )
            + "</row>"
        )
    parts["xl/worksheets/sheet1.xml"] = (
        f'<worksheet xmlns="{ns}"><dimension ref="A1:{column(len(formats))}5"/><sheetData>'
        + "".join(rows)
        + "</sheetData></worksheet>"
    )
    with ZipFile(directory / "display_formats.xlsx", "w", compression=ZIP_DEFLATED) as book:
        for name, text in parts.items():
            book.writestr(name, text)

    # Independent BIFF8/OLE container with the same values and format IDs.
    def record(kind, body=b""):
        return struct.pack("<HH", kind, len(body)) + body

    def bof(kind):
        return record(0x809, struct.pack("<4H2I", 0x600, kind, 0xDBB, 1997, 0x41, 6))

    sheet = bof(0x10) + record(0x200, struct.pack("<IIHHH", 0, 5, 0, len(formats), 0))
    for j, fmt in enumerate(formats):
        sheet += record(0xFD, struct.pack("<HHHI", 0, j, 0, j))
    for row, value in [(1, 44300.5123456), (2, 0.5123456), (3, 44300.125), (4, 44197.001)]:
        for j in range(len(formats)):
            sheet += record(0x203, struct.pack("<HHHd", row, j, j + 1, value))
    sheet += record(0x0A)
    globals_ = bof(5) + record(0x42, struct.pack("<H", 1200)) + record(0x22, struct.pack("<H", 0))
    for j, code in enumerate(custom):
        globals_ += record(
            0x41E, struct.pack("<HHB", 165 + j, len(code), 1) + code.encode("utf-16le")
        )
    for fmt in [0] + formats:
        globals_ += record(0xE0, struct.pack("<HH", 0, fmt) + bytes(16))
    strings = b"".join(
        struct.pack("<HB", len("fmt" + str(fmt)), 0) + ("fmt" + str(fmt)).encode("ascii")
        for fmt in formats
    )
    globals_ += record(0xFC, struct.pack("<II", len(formats), len(formats)) + strings)
    bound = record(0x85, struct.pack("<IBBBB", 0, 0, 0, len("Formats"), 0) + b"Formats")
    offset = len(globals_) + len(bound) + 4
    bound = record(0x85, struct.pack("<IBBBB", offset, 0, 0, len("Formats"), 0) + b"Formats")
    stream = globals_ + bound + record(0x0A) + sheet
    size = max(4096, (len(stream) + 511) // 512 * 512)
    n = size // 512
    header = bytearray(512)
    header[:8] = bytes.fromhex("d0cf11e0a1b11ae1")
    struct.pack_into("<HHHHH", header, 24, 0x3E, 3, 0xFFFE, 9, 6)
    struct.pack_into("<IIIII", header, 40, 0, 1, n, 0, 4096)
    struct.pack_into("<4I", header, 60, 0xFFFFFFFE, 0, 0xFFFFFFFE, 0)
    struct.pack_into("<109I", header, 76, n + 1, *([0xFFFFFFFF] * 108))

    def entry(name, kind, start, length, child=0xFFFFFFFF):
        out = bytearray(128)
        encoded = (name + "\0").encode("utf-16le")
        out[: len(encoded)] = encoded
        struct.pack_into("<HBBIII", out, 64, len(encoded), kind, 1, 0xFFFFFFFF, 0xFFFFFFFF, child)
        struct.pack_into("<IQ", out, 116, start, length)
        return out

    ole_directory = (
        entry("Root Entry", 5, 0xFFFFFFFE, 0, 1) + entry("Workbook", 2, 0, size) + bytes(256)
    )
    fat = [i + 1 for i in range(n)]
    fat[-1] = 0xFFFFFFFE
    fat += [0xFFFFFFFE, 0xFFFFFFFD] + [0xFFFFFFFF] * (128 - n - 2)
    (directory / "display_formats.xls").write_bytes(
        bytes(header) + stream.ljust(size, b"\0") + ole_directory + struct.pack("<128I", *fat)
    )


NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"


CASES = [
    ("simple", location, 1, 1, "modify") for location in ("C2", "C3", "C4", "D2", "D3", "D4", "A2")
] + [
    ("grid", "C2", 1, 1, "modify"),
    ("grid", "D3", 1, 1, "modify"),
    ("grid", "C2", 2, 2, "modify"),
    ("grid", "B2", 2, 3, "modify"),
    ("grid", "C2", 3, 3, "modify"),
    ("grid", "C2", 1, 1, "replace"),
    ("grid", "C2", 1, 1, "new"),
    ("grid", "C2", 1, 3, "modify"),
]


def formula_fixtures(directory: Path):
    with ZipFile(directory / "display_formats.xlsx") as book:
        base = {name: book.read(name).decode() for name in book.namelist()}
    for kind in ("simple", "grid"):
        parts = base.copy()
        parts["xl/workbook.xml"] = (
            f'<workbook xmlns="{NS}" xmlns:r="http://schemas.openxmlformats.org/'
            'officeDocument/2006/relationships"><sheets><sheet name="Data" '
            'sheetId="1" r:id="rId1"/></sheets><calcPr calcId="124519" '
            f'fullCalcOnLoad="{1 if kind == "simple" else 0}"/></workbook>'
        )
        parts["xl/_rels/workbook.xml.rels"] = parts["xl/_rels/workbook.xml.rels"].replace(
            "</Relationships>",
            '<Relationship Id="rId3" Type="http://schemas.openxmlformats.org/'
            'officeDocument/2006/relationships/calcChain" Target="calcChain.xml"/>'
            "</Relationships>",
        )
        parts["[Content_Types].xml"] = parts["[Content_Types].xml"].replace(
            "</Types>",
            '<Override PartName="/xl/calcChain.xml" ContentType="application/'
            'vnd.openxmlformats-officedocument.spreadsheetml.calcChain+xml"/></Types>',
        )
        rows = ['<row r="1"><c r="A1" t="inlineStr"><is><t>Input</t></is></c></row>']
        chain = []
        for row in range(2, 5):
            cells = [f'<c r="A{row}"><v>{row}</v></c>', f'<c r="B{row}"><v>{row+1}</v></c>']
            for col in ("C",) if kind == "simple" else ("C", "D", "E"):
                ref = f"{col}{row}"
                if ref == "C2":
                    formula = (
                        "A2+B2"
                        if kind == "simple"
                        else 'IF(A2>1,A2+$A2+A$2+$A$2+SUM(A2:B4)+SUM(A:A)+SUM(2:2),LEN("A2"))'
                    )
                    f = (
                        f'<f t="shared" ref="C2:{"C" if kind == "simple" else "E"}4" '
                        f'si="0">{escape(formula)}</f>'
                    )
                else:
                    f = '<f t="shared" si="0"/>'
                cells.append(f'<c r="{ref}">{f}<v>{2*row+1}</v></c>')
                chain.append(f'<c r="{ref}" i="1"/>')
            arraycol = "D" if kind == "simple" else "G"
            f = f'<f t="array" ref="{arraycol}2:{arraycol}4">A2:A4*B2:B4</f>' if row == 2 else ""
            cells.append(f'<c r="{arraycol}{row}">{f}<v>{row*(row+1)}</v></c>')
            rows.append(f'<row r="{row}">' + "".join(cells) + "</row>")
        chain.append(f'<c r="{arraycol}2" i="1"/>')
        parts["xl/worksheets/sheet1.xml"] = (
            f'<worksheet xmlns="{NS}"><dimension ref="A1:{arraycol}4"/>'
            "<sheetData>" + "".join(rows) + "</sheetData></worksheet>"
        )
        parts["xl/calcChain.xml"] = f'<calcChain xmlns="{NS}">' + "".join(chain) + "</calcChain>"
        with ZipFile(directory / f"formulas_{kind}.xlsx", "w", compression=ZIP_DEFLATED) as book:
            for name, xml in parts.items():
                book.writestr(name, xml)


def formula_state(path: Path):
    """Compare formula nodes and their cached values, independent of XML layout."""
    formulas = []
    with ZipFile(path) as book:
        for name in book.namelist():
            if name.endswith(".xml") or name.endswith(".rels"):
                root = ET.fromstring(book.read(name))
                if name == "xl/workbook.xml":
                    calc = root.find(f"{{{NS}}}calcPr")
                    assert calc is not None and calc.get("fullCalcOnLoad") in ("true", "1"), path
                for element in root.iter():
                    assert not element.tag.endswith("}calcChain"), path
                    assert not element.get("Type", "").endswith("/calcChain"), path
                    assert not element.get("ContentType", "").endswith("calcChain+xml"), path
                if name.startswith("xl/worksheets/") and name.endswith(".xml"):
                    for cell in root.findall(f".//{{{NS}}}c"):
                        f = cell.find(f"{{{NS}}}f")
                        if f is not None:
                            formulas.append(
                                (
                                    name,
                                    cell.get("r"),
                                    sorted(f.attrib.items()),
                                    f.text or "",
                                    cell.findtext(f"{{{NS}}}v", ""),
                                )
                            )
    return sorted(formulas)


def compare_formula_outputs(directory: Path):
    for index, _ in enumerate(CASES, 1):
        native = formula_state(directory / f"formula_native_{index}.xlsx")
        result = formula_state(directory / f"formula_result_{index}.xlsx")
        if native != result:
            raise AssertionError(f"Formula structure case {index}:\nNative: {native}\nC: {result}")
    return len(CASES)


def integrity_fixtures(p):
    p.mkdir(exist_ok=True)
    ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    rel = "http://schemas.openxmlformats.org/package/2006/relationships"
    orel = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
    files = {
        "[Content_Types].xml": (
            '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
            '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
            '<Default Extension="xml" ContentType="application/xml"/>'
            '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
            '<Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
            "</Types>"
        ),
        "_rels/.rels": (
            f'<Relationships xmlns="{rel}"><Relationship Id="rId1" '
            f'Type="{orel}/officeDocument" Target="xl/workbook.xml"/></Relationships>'
        ),
        "xl/workbook.xml": (
            f'<workbook xmlns="{ns}" xmlns:r="{orel}"><sheets>'
            '<sheet name="Sheet1" sheetId="1" r:id="rId1"/></sheets></workbook>'
        ),
        "xl/_rels/workbook.xml.rels": (
            f'<Relationships xmlns="{rel}"><Relationship Id="rId1" '
            f'Type="{orel}/worksheet" Target="worksheets/sheet1.xml"/></Relationships>'
        ),
        "xl/worksheets/sheet1.xml": (
            f'<worksheet xmlns="{ns}"><dimension ref="A1:A2"/><sheetData>'
            '<row r="1"><c r="A1"><v>123</v></c></row>'
            '<row r="2"><c r="A2"><v>456</v></c></row>'
            "</sheetData></worksheet>"
        ),
    }

    def write_xlsx(name, parts, level=6):
        with zipfile.ZipFile(
            p / name, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=level
        ) as z:
            for key, value in parts.items():
                z.writestr(key, value)

    write_xlsx("valid.xlsx", files)
    write_xlsx("bitflip.xlsx", files, level=0)
    b = bytearray((p / "bitflip.xlsx").read_bytes())
    b[b.index(b"<v>123</v>") + 3] ^= 8
    (p / "bitflip.xlsx").write_bytes(b)
    try:
        with zipfile.ZipFile(p / "bitflip.xlsx") as z:
            z.read("xl/worksheets/sheet1.xml")
    except zipfile.BadZipFile as error:
        print(error)  # Expected CRC failure

    truncated = dict(files)
    truncated["xl/worksheets/sheet1.xml"] = (
        files["xl/worksheets/sheet1.xml"].split('<row r="2">')[0] + '<row r="2"><c r="A2"><v>45'
    )
    write_xlsx("truncated-xml.xlsx", truncated)
    with zipfile.ZipFile(p / "truncated-xml.xlsx") as z:
        assert z.testzip() is None  # ZIP is intact; contained XML is incomplete

    h = bytearray(100)
    struct.pack_into(">I", h, 0, 9994)
    struct.pack_into(">I", h, 24, 64)
    struct.pack_into("<II", h, 28, 1000, 1)
    struct.pack_into("<dddd", h, 36, 1, 2, 1, 2)
    (p / "point.shp").write_bytes(h + struct.pack(">II", 1, 10) + struct.pack("<Idd", 1, 1, 2))
    (p / "partial.shx").write_text("old index sentinel")
    (p / "corrupt.dbf").write_bytes(b"not a dbf")
    (p / "transport.sas").write_text("old SAS sentinel")


def regression_fixtures(folder):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    rns = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
    pns = "http://schemas.openxmlformats.org/package/2006/relationships"
    book = f'<workbook xmlns="{ns}" xmlns:r="{rns}"><sheets><sheet name="Sheet1" sheetId="1" r:id="rId1"/></sheets></workbook>'
    rels = f'<Relationships xmlns="{pns}"><Relationship Id="rId1" Type="{rns}/worksheet" Target="worksheets/sheet1.xml"/><Relationship Id="rId2" Type="{rns}/sharedStrings" Target="sharedStrings.xml"/></Relationships>'
    base = {
        "xl/workbook.xml": book,
        "xl/_rels/workbook.xml.rels": rels,
        "_rels/.rels": f'<Relationships xmlns="{pns}"><Relationship Id="rId1" Type="{rns}/officeDocument" Target="xl/workbook.xml"/></Relationships>',
        "[Content_Types].xml": '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/><Override PartName="/xl/sharedStrings.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sharedStrings+xml"/></Types>',
        "xl/sharedStrings.xml": f'<sst xmlns="{ns}"><si><t>first</t></si></sst>',
        "xl/styles.xml": f'<styleSheet xmlns="{ns}"><cellXfs count="1"><xf numFmtId="0"/></cellXfs></styleSheet>',
    }
    cell = '<c r="A1" t="s"><v>0</v></c>'

    def write(name, text=cell, **parts):
        content = dict(base)
        content["xl/worksheets/sheet1.xml"] = (
            f'<worksheet xmlns="{ns}"><sheetData><row r="1">{text}</row></sheetData></worksheet>'
        )
        content.update(parts)
        with zipfile.ZipFile(folder / (name + ".xlsx"), "w", zipfile.ZIP_DEFLATED) as archive:
            for key, value in content.items():
                archive.writestr(key, value)

    write("good")
    for name, part in [
        ("sst", "xl/sharedStrings.xml"),
        ("styles", "xl/styles.xml"),
        ("book", "xl/workbook.xml"),
        ("rels", "xl/_rels/workbook.xml.rels"),
    ]:
        write("truncated_" + name, **{part: base[part].rsplit("</", 1)[0]})
    write("invalid_index", '<c r="A1" t="s"><v>1</v></c>')
    for name, text in [
        ("cdata", "<![CDATA[customer & name]]>"),
        ("comment", "before<!-- with > and < -->after"),
        ("escaped", "a_x000D_b"),
    ]:
        write(name, f'<c r="A1" t="inlineStr"><is><t>{text}</t></is></c>')
        write(
            "shared_" + name,
            **{"xl/sharedStrings.xml": f'<sst xmlns="{ns}"><si><t>{text}</t></si></sst>'},
        )
    records = b"".join(
        struct.pack(">II", i, 10) + struct.pack("<Idd", 1, float(i), float(i) * 2) for i in [1, 2]
    )
    header = struct.pack(">7I", 9994, 0, 0, 0, 0, 0, (100 + len(records)) // 2) + struct.pack(
        "<2I8d", 1000, 1, 1, 2, 2, 4, 0, 0, 0, 0
    )
    (folder / "good.shp").write_bytes(header + records)
    (folder / "truncated.shp").write_bytes(header + records[:28])
