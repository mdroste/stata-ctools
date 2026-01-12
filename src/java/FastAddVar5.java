import com.stata.sfi.*;

public class FastAddVar5 {
    public static int run(String[] args) {
        try {
            SFIToolkit.displayln("Starting with Frame API...");

            // Get the current working frame
            Frame frame = Frame.connect("default");
            if (frame == null) {
                SFIToolkit.errorln("Could not connect to default frame");
                return 1;
            }

            int rc = frame.addVarDouble("testvar");
            SFIToolkit.displayln("Added variable via Frame, rc=" + rc);

            return 0;
        } catch (Exception e) {
            SFIToolkit.errorln("Error: " + e.getMessage());
            e.printStackTrace();
            return 1;
        }
    }
}
