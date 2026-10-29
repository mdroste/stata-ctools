"""Independent shared/array formula fixtures and structural comparisons."""
from pathlib import Path
from xml.etree import ElementTree as ET
from xml.sax.saxutils import escape
from zipfile import ZipFile, ZIP_DEFLATED

NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
CASES = [
    ("simple", location, 1, 1, "modify")
    for location in ("C2", "C3", "C4", "D2", "D3", "D4", "A2")
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
        parts["xl/_rels/workbook.xml.rels"] = parts[
            "xl/_rels/workbook.xml.rels"
        ].replace(
            "</Relationships>",
            '<Relationship Id="rId3" Type="http://schemas.openxmlformats.org/'
            'officeDocument/2006/relationships/calcChain" Target="calcChain.xml"/>'
            '</Relationships>',
        )
        parts["[Content_Types].xml"] = parts["[Content_Types].xml"].replace(
            "</Types>",
            '<Override PartName="/xl/calcChain.xml" ContentType="application/'
            'vnd.openxmlformats-officedocument.spreadsheetml.calcChain+xml"/></Types>',
        )
        rows = ['<row r="1"><c r="A1" t="inlineStr"><is><t>Input</t></is></c></row>']
        chain = []
        for row in range(2, 5):
            cells = [f'<c r="A{row}"><v>{row}</v></c>',
                     f'<c r="B{row}"><v>{row+1}</v></c>']
            for col in ("C",) if kind == "simple" else ("C", "D", "E"):
                ref = f"{col}{row}"
                if ref == "C2":
                    formula = ("A2+B2" if kind == "simple" else
                               'IF(A2>1,A2+$A2+A$2+$A$2+SUM(A2:B4)+SUM(A:A)+SUM(2:2),LEN("A2"))')
                    f = (f'<f t="shared" ref="C2:{"C" if kind == "simple" else "E"}4" '
                         f'si="0">{escape(formula)}</f>')
                else:
                    f = '<f t="shared" si="0"/>'
                cells.append(f'<c r="{ref}">{f}<v>{2*row+1}</v></c>')
                chain.append(f'<c r="{ref}" i="1"/>')
            arraycol = "D" if kind == "simple" else "G"
            f = (f'<f t="array" ref="{arraycol}2:{arraycol}4">A2:A4*B2:B4</f>'
                 if row == 2 else "")
            cells.append(f'<c r="{arraycol}{row}">{f}<v>{row*(row+1)}</v></c>')
            rows.append(f'<row r="{row}">'+"".join(cells)+"</row>")
        chain.append(f'<c r="{arraycol}2" i="1"/>')
        parts["xl/worksheets/sheet1.xml"] = (
            f'<worksheet xmlns="{NS}"><dimension ref="A1:{arraycol}4"/>'
            '<sheetData>'+"".join(rows)+"</sheetData></worksheet>"
        )
        parts["xl/calcChain.xml"] = f'<calcChain xmlns="{NS}">'+"".join(chain)+"</calcChain>"
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
                            formulas.append((name, cell.get("r"), sorted(f.attrib.items()),
                                             f.text or "", cell.findtext(f"{{{NS}}}v", "")))
    return sorted(formulas)


def compare_formula_outputs(directory: Path):
    for index, _ in enumerate(CASES, 1):
        native = formula_state(directory / f"formula_native_{index}.xlsx")
        result = formula_state(directory / f"formula_result_{index}.xlsx")
        if native != result:
            raise AssertionError(f"Formula structure case {index}:\nNative: {native}\nC: {result}")
    return len(CASES)
