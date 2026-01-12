import com.stata.sfi.*;

public class AddVar3 {
    public static int create(String[] args) {
        if (args.length < 2) {
            SFIToolkit.errorln("Need type and name");
            return 198;
        }

        String type = args[0].trim().toLowerCase();
        String name = args[1].trim();

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

        return (rc >= 0) ? 0 : 920;
    }
}
