import com.stata.sfi.*;
import java.util.concurrent.*;

public class FastAddVar {

    public static int addvars(String[] args) {
        // Read from globals: ADDVAR_TYPES and ADDVAR_NAMES (comma-separated)
        String typesStr = Macro.getGlobal("ADDVAR_TYPES");
        String namesStr = Macro.getGlobal("ADDVAR_NAMES");

        if (typesStr == null || typesStr.isEmpty() || namesStr == null || namesStr.isEmpty()) {
            SFIToolkit.errorln("Usage: set globals ADDVAR_TYPES and ADDVAR_NAMES before calling");
            return 198;
        }

        String[] types = typesStr.trim().split(",");
        String[] names = namesStr.trim().split(",");

        if (types.length != names.length) {
            SFIToolkit.errorln("Number of types must match number of names");
            return 198;
        }

        int n = types.length;
        if (n == 0) return 0;

        // Use parallel streams for variable creation
        int nThreads = Math.min(n, Runtime.getRuntime().availableProcessors());
        ExecutorService executor = Executors.newFixedThreadPool(nThreads);
        CountDownLatch latch = new CountDownLatch(n);

        int[] results = new int[n];  // Store variable indices or -1 for error

        for (int i = 0; i < n; i++) {
            final int idx = i;
            final String type = types[i].toLowerCase();
            final String name = names[i];

            executor.submit(() -> {
                try {
                    int varIdx = addVariable(type, name);
                    results[idx] = varIdx;
                } catch (Exception e) {
                    results[idx] = -1;
                } finally {
                    latch.countDown();
                }
            });
        }

        try {
            latch.await();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            return 1;
        }

        executor.shutdown();

        // Check for errors
        for (int i = 0; i < n; i++) {
            if (results[i] < 0) {
                SFIToolkit.errorln("Failed to create variable: " + names[i]);
                return 920;
            }
        }

        return 0;
    }

    private static int addVariable(String type, String name) {
        if (type.equals("byte")) {
            return Data.addVarByte(name);
        } else if (type.equals("int")) {
            return Data.addVarInt(name);
        } else if (type.equals("long")) {
            return Data.addVarLong(name);
        } else if (type.equals("float")) {
            return Data.addVarFloat(name);
        } else if (type.equals("double")) {
            return Data.addVarDouble(name);
        } else if (type.startsWith("str")) {
            // Parse string length: str20 -> 20
            int len = 1;
            if (type.length() > 3) {
                try {
                    len = Integer.parseInt(type.substring(3));
                } catch (NumberFormatException e) {
                    len = 244;  // Default max str length
                }
            }
            return Data.addVarStr(name, len);
        } else {
            // Default to double
            return Data.addVarDouble(name);
        }
    }
}
