import com.stata.sfi.*;

public class FastAddVar4 {
    public static int run(String[] args) {
        SFIToolkit.displayln("Starting...");
        int rc = Data.addVarDouble("testvar");
        SFIToolkit.displayln("Added variable, rc=" + rc);
        return 0;
    }
}
