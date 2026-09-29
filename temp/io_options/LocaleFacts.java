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
 }
}
