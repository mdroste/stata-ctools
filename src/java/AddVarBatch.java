import com.stata.sfi.*;

public class AddVarBatch {
    // Takes comma-separated types and names as two args
    public static int create(String[] args) {
        if (args.length < 2) {
            SFIToolkit.errorln("Need types and names (comma-separated)");
            return 198;
        }

        String[] types = args[0].split("\\+");
        String[] names = args[1].split("\\+");

        if (types.length != names.length) {
            SFIToolkit.errorln("Number of types must match number of names");
            return 198;
        }

        for (int i = 0; i < types.length; i++) {
            String type = types[i].trim().toLowerCase();
            String name = names[i].trim();

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
                int len = 244;
                if (type.length() > 3) {
                    try { len = Integer.parseInt(type.substring(3)); } catch (Exception e) {}
                }
                rc = Data.addVarStr(name, len);
            }

            if (rc < 0) {
                SFIToolkit.errorln("Failed to create variable: " + name);
                return 920;
            }
        }

        return 0;
    }
}
