import com.stata.sfi.*;

public class AddVar {
    public static int create(String[] args) {
        SFIToolkit.displayln("AddVar.create called with " + args.length + " args");
        for (int i = 0; i < args.length; i++) {
            SFIToolkit.displayln("  args[" + i + "] = " + args[i]);
        }

        if (args.length < 2) {
            SFIToolkit.errorln("Usage: javacall AddVar create, args(type name)");
            return 198;
        }

        String type = args[0].trim().toLowerCase();
        String name = args[1].trim();
        SFIToolkit.displayln("Creating " + type + " variable: " + name);

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

        return 0;
    }
}
