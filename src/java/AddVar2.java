import com.stata.sfi.*;

public class AddVar2 {
    public static int create(String[] args) {
        SFIToolkit.displayln("AddVar2.create called with " + args.length + " args");
        for (int i = 0; i < args.length; i++) {
            SFIToolkit.displayln("  args[" + i + "] = " + args[i]);
        }
        return 0;
    }
}
