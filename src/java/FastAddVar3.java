import com.stata.sfi.*;

public class FastAddVar3 {
    public static int run(String[] args) {
        String typesStr = Macro.getGlobal("ADDVAR_TYPES");
        String namesStr = Macro.getGlobal("ADDVAR_NAMES");

        SFIToolkit.displayln("Types: " + typesStr);
        SFIToolkit.displayln("Names: " + namesStr);

        if (typesStr == null || typesStr.isEmpty()) {
            SFIToolkit.errorln("ADDVAR_TYPES not set");
            return 198;
        }

        String[] types = typesStr.trim().split(",");
        String[] names = namesStr.trim().split(",");

        for (int i = 0; i < types.length && i < names.length; i++) {
            String type = types[i].trim().toLowerCase();
            String name = names[i].trim();

            if (type.equals("double")) {
                Data.addVarDouble(name);
            } else if (type.equals("long")) {
                Data.addVarLong(name);
            } else if (type.equals("byte")) {
                Data.addVarByte(name);
            } else if (type.startsWith("str")) {
                int len = 244;
                try { len = Integer.parseInt(type.substring(3)); } catch (Exception e) {}
                Data.addVarStr(name, len);
            }
        }
        return 0;
    }
}
