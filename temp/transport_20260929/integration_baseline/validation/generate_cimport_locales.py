"""Generate numeric locale facts; JDK runs only here, never in the C reader.

Use the same OpenJDK/CLDR version as the native comparison installation:
python3 validation/generate_cimport_locales.py --java /path/to/jdk17/bin/java
"""
import argparse
import subprocess
import tempfile
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
JAVA = r'''
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
'''
def c_string(value):
    return '"'+''.join(chr(b) if 32<=b<127 and b not in (34,92) else '\\%03o'%b
                       for b in value.encode())+'"'
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--java',default='java');args=parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='ctools-locale-facts-') as tmp:
        source=Path(tmp)/'LocaleFacts.java';source.write_text(JAVA)
        facts=subprocess.check_output([args.java,str(source)],text=True)
        version=subprocess.check_output([args.java,'-version'],stderr=subprocess.STDOUT,text=True).splitlines()[0]
    rows=[];digits=[]
    for line in facts.splitlines():
        if line.startswith('DIGIT\t'): digits.append(hex(int(line.split('\t')[1])));continue
        tag,decimal,group,negative,exponent=line.split('\t')
        rows.append('    {'+','.join((c_string(tag),decimal,group,c_string(negative),c_string(exponent)))+'},')
    text='/* Numeric locale facts generated from '+version.replace('*/','')+'.\n'
    text+='   No Java code executes at runtime; see validation/generate_cimport_locales.py. */\n'
    text+='static const cimport_locale_profile cimport_locale_profiles[] = {\n'+'\n'.join(rows)+'\n};\n'
    text+='static const uint32_t cimport_digit_starts[] = {\n    '+','.join(digits)+'\n};\n'
    (ROOT/'src/cimport/cimport_locale_profiles.inc').write_text(text)
    print(f'Generated {len(rows)} locale profiles and {len(digits)} BMP digit blocks')
if __name__=='__main__':main()
