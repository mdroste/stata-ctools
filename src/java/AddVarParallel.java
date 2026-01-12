import com.stata.sfi.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

public class AddVarParallel {
    public static int create(String[] args) {
        if (args.length < 2) {
            SFIToolkit.errorln("Need types and names (+-separated)");
            return 198;
        }

        String[] types = args[0].split("\\+");
        String[] names = args[1].split("\\+");

        if (types.length != names.length) {
            SFIToolkit.errorln("Number of types must match number of names");
            return 198;
        }

        int n = types.length;
        int nThreads = Math.min(n, Runtime.getRuntime().availableProcessors());

        SFIToolkit.displayln("AddVarParallel: " + n + " vars, " + nThreads + " threads");

        ExecutorService executor = Executors.newFixedThreadPool(nThreads);
        CountDownLatch latch = new CountDownLatch(n);
        AtomicInteger errors = new AtomicInteger(0);

        for (int i = 0; i < n; i++) {
            final String type = types[i].trim().toLowerCase();
            final String name = names[i].trim();

            executor.submit(() -> {
                try {
                    int rc = addVariable(type, name);
                    if (rc < 0) {
                        errors.incrementAndGet();
                    }
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

        return errors.get() > 0 ? 920 : 0;
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
            int len = 244;
            if (type.length() > 3) {
                try { len = Integer.parseInt(type.substring(3)); } catch (Exception e) {}
            }
            return Data.addVarStr(name, len);
        }
        return -1;
    }
}
