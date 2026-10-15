"""Generate numeric locale facts; JDK runs only here, never in the C reader.

Use the same OpenJDK/CLDR version as the native comparison installation:
python3 tools/generation/generate_cimport_locales.py --java /path/to/jdk17/bin/java
"""

import argparse
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
JAVA = r"""
import java.util.*;
import java.text.*;
public class LocaleFacts {
 public static void main(String[] args) {
  var seen=new TreeSet<String>();
  for (Locale l:Locale.getAvailableLocales()) seen.add(l.toLanguageTag().toLowerCase(Locale.ROOT));
  for(String tag:seen) {
   Locale l=Locale.forLanguageTag(tag);DecimalFormat f=(DecimalFormat)NumberFormat.getNumberInstance(l);
   var s=f.getDecimalFormatSymbols();
   System.out.printf("%s\t%d\t%d\t%s\t%s%n",tag,(int)s.getDecimalSeparator(),(int)s.getGroupingSeparator(),f.getNegativePrefix(),s.getExponentSeparator());
  }
  for(int cp=0;cp<65536;cp++)
   if(Character.getType(cp)==Character.DECIMAL_DIGIT_NUMBER && Character.digit(cp,10)==0)
    System.out.println("DIGIT\t"+cp);
 }
}
"""


def c_string(value):
    return (
        '"'
        + "".join(
            chr(b) if 32 <= b < 127 and b not in (34, 92) else "\\%03o" % b for b in value.encode()
        )
        + '"'
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--java", help="verify pinned facts against the exact recorded JDK")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    inputs = Path(__file__).with_name("inputs")
    metadata = json.loads((inputs / "locales.json").read_text())
    facts = (inputs / "locale_facts.tsv").read_text()
    version = metadata["runtime"]
    if hashlib.sha256(facts.encode()).hexdigest() != metadata["facts_sha256"]:
        raise SystemExit("Locale input checksum mismatch")
    if hashlib.sha256(JAVA.encode()).hexdigest() != metadata["generator_java_sha256"]:
        raise SystemExit("Java extraction source changed; review pinned inputs")
    if args.java:
        with tempfile.TemporaryDirectory(prefix="ctools-locale-facts-") as tmp:
            source = Path(tmp) / "LocaleFacts.java"
            source.write_text(JAVA)
            actual = subprocess.check_output([args.java, str(source)], text=True)
            runtime = subprocess.check_output(
                [args.java, "-version"], stderr=subprocess.STDOUT, text=True
            ).splitlines()[0]
        if actual != facts or runtime != version:
            raise SystemExit("Runtime or locale facts differ from pinned inputs")
    rows = []
    digits = []
    for line in facts.splitlines():
        if line.startswith("DIGIT\t"):
            digits.append(hex(int(line.split("\t")[1])))
            continue
        tag, decimal, group, negative, exponent = line.split("\t")
        rows.append(
            "    {"
            + ",".join((c_string(tag), decimal, group, c_string(negative), c_string(exponent)))
            + "},"
        )
    text = "/* Numeric locale facts generated from " + version.replace("*/", "") + ".\n"
    text += "   No Java code executes at runtime; see tools/generation/generate_cimport_locales.py. */\n"
    text += (
        "static const cimport_locale_profile cimport_locale_profiles[] = {\n"
        + "\n".join(rows)
        + "\n};\n"
    )
    text += "static const uint32_t cimport_digit_starts[] = {\n    " + ",".join(digits) + "\n};\n"
    output = ROOT / "src/cimport/generated/cimport_locale_profiles.inc"
    if args.check:
        if output.read_text() != text:
            raise SystemExit("Stale locale table")
    else:
        output.write_text(text)
    print(f"Generated {len(rows)} locale profiles and {len(digits)} BMP digit blocks")


if __name__ == "__main__":
    main()
