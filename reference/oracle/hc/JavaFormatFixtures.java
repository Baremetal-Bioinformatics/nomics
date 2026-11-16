// Fixture generator: Java's formatting of doubles as GATK/htsjdk use it.
// Run in the reference container's JVM (OpenJDK 17.0.17):
//   java JavaFormatFixtures.java N SEED > fixtures.tsv
// Each line: hex bits of the double, then String.format(Locale.US, f, d) for
// f in %.2f %.3f %.4f %.3e, then Double.toString(d), tab-separated.
import java.util.Locale;
import java.util.Random;

public class JavaFormatFixtures {
    public static void main(String[] a) {
        int n = Integer.parseInt(a[0]);
        Random r = new Random(Long.parseLong(a[1]));
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < n; ++i) {
            double d;
            switch (i % 8) {
                case 0: d = r.nextDouble() * 100; break;                       // QD, MQ, AF-like
                case 1: d = (r.nextDouble() - 0.5) * 20; break;                // rank sums, incl. negatives
                case 2: d = r.nextDouble() * 1e5; break;                       // QUAL
                case 3: { // decimal ties at 2-4 places, +/- 0-2 ulps
                    int places = 2 + r.nextInt(3);
                    long m = r.nextInt(2000000) - 1000000;
                    d = (m + 0.5) / Math.pow(10, places);
                    int u = r.nextInt(5) - 2;
                    for (int k = 0; k < Math.abs(u); ++k) d = u > 0 ? Math.nextUp(d) : Math.nextDown(d);
                    break;
                }
                case 4: d = Math.pow(10, r.nextDouble() * 12 - 8) * (r.nextBoolean() ? 1 : -1); break; // wide magnitudes
                case 5: d = Double.longBitsToDouble(r.nextLong()); if (Double.isNaN(d) || Double.isInfinite(d) || Math.abs(d) > 1e15 || Math.abs(d) < 1e-300) d = r.nextDouble(); break;
                case 6: d = r.nextInt(1000) / 1000.0 + r.nextInt(3) * Math.ulp(1.0); break;
                default: d = -10 * Math.log10(r.nextDouble()) + 0.0; break;    // phred-like
            }
            sb.setLength(0);
            sb.append(Long.toHexString(Double.doubleToRawLongBits(d)));
            for (String f : new String[] {"%.2f", "%.3f", "%.4f", "%.3e"}) sb.append('\t').append(String.format(Locale.US, f, d));
            sb.append('\t').append(Double.toString(d));
            System.out.println(sb);
        }
        // Edge cases.
        double[] edge = {0.0, -0.0, 0.005, 0.015, 0.125, 0.375, 1.005, 2.675, 1e-20, 9.995, 99.995, 0.0005, 0.00049999999999999999,
            Double.MIN_VALUE, Double.MAX_VALUE, 1e7, 1e-3, 1e-4, 123456789.125, 0.1, 0.2, 0.3, 1.0 / 3, 2.0 / 3};
        for (double d : edge) {
            sb.setLength(0);
            sb.append(Long.toHexString(Double.doubleToRawLongBits(d)));
            for (String f : new String[] {"%.2f", "%.3f", "%.4f", "%.3e"}) sb.append('\t').append(String.format(Locale.US, f, d));
            sb.append('\t').append(Double.toString(d));
            System.out.println(sb);
        }
    }
}
