import com.stata.sfi.*;

public class FastAddVar2 {
    public static int addvars(String[] args) {
        // Read from globals: ADDVAR_TYPES and ADDVAR_NAMES (comma-separated)
        String typesStr = Macro.getGlobal("ADDVAR_TYPES");
        String namesStr = Macro.getGlobal("ADDVAR_NAMES");

        SFIToolkit.displayln("Types: " + typesStr);
        SFIToolkit.displayln("Names: " + namesStr);

        if (typesStr == null || typesStr.isEmpty()) {
            SFIToolkit.errorln("ADDVAR_TYPES not set");
            return 198;
        }
        if (namesStr == null || namesStr.isEmpty()) {
            SFIToolkit.errorln("ADDVAR_NAMES not set");
            return 198;
        }

        String[] types = typesStr.trim().split(",");
        String[] names = namesStr.trim().split(",");

        // Simple sequential add
        for (int i = 0; i < types.length && i < names.length; i++) {
            String type = types[i].trim().toLowerCase();
            String name = names[i].trim();

            SFIToolkit.displayln("Adding: " + type + " " + name);

            int rc = -1;
            if (type.equals("byte")) {
                rc = Data.addVarByte(name);
            } else if (type.equals("int")) {
                rc = Data.addVarInt(name);
            } else if (type.equals("long")) {
                rc = Data.addVarLong(name);
            } else if (type.equals("float")) {
                rc = Data.addVarFloat(name);
            } else if (type.equals("double")) {
                rc = Data.addVarDouble(name);
            } else if (type.startsWith("str")) {
                int len = 1;
                try {
                    len = Integer.parseInt(type.substring(3));
                } catch (Exception e) {
                    len = 244;
                }
                rc = Data.addVarStr(name, len);
            }

            SFIToolkit.displayln("  Result: " + rc);
        }

        return 0;
    }
}
