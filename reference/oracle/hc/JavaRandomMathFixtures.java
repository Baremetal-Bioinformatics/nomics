// Fixture generator: java.util.Random call streams and StrictMath functions,
// from the reference container's JVM (OpenJDK 17.0.17).
//   java JavaRandomMathFixtures.java random N > random.tsv
//   java JavaRandomMathFixtures.java math N > math.tsv
// random.tsv: one line per call on a single Random(47382911):
//   "i <bound> <value>" (nextInt), "g <hexbits>" (nextGaussian), "d <hexbits>" (nextDouble)
// math.tsv: "<fn> <hex x> [<hex y>] <hex result>" for StrictMath log, log10, exp, pow, log1p, expm1.
import java.util.Random;

public class JavaRandomMathFixtures {
    static String h(double d) { return Long.toHexString(Double.doubleToRawLongBits(d)); }

    public static void main(String[] a) {
        int n = Integer.parseInt(a[1]);
        StringBuilder sb = new StringBuilder();
        if (a[0].equals("random")) {
            Random r = new Random(47382911L);
            Random pick = new Random(7);
            for (int i = 0; i < n; ++i) {
                int op = pick.nextInt(10);
                if (op < 6) {
                    int bound = 51 + pick.nextInt(2000);
                    if (pick.nextInt(20) == 0) bound = 1 << pick.nextInt(31); // powers of two
                    if (pick.nextInt(50) == 0) bound = Integer.MAX_VALUE - pick.nextInt(1000); // rejection loop
                    sb.append("i ").append(bound).append(' ').append(r.nextInt(bound)).append('\n');
                } else if (op < 9) {
                    sb.append("g ").append(h(r.nextGaussian())).append('\n');
                } else {
                    sb.append("d ").append(h(r.nextDouble())).append('\n');
                }
            }
        } else {
            Random r = new Random(11);
            for (int i = 0; i < n; ++i) {
                double x;
                switch (i % 4) {
                    case 0: x = r.nextDouble(); break;
                    case 1: x = Math.pow(10, r.nextDouble() * 40 - 20); break;
                    case 2: x = Double.longBitsToDouble(r.nextLong() & 0x7fffffffffffffffL); if (Double.isNaN(x) || Double.isInfinite(x)) x = 1.5; break;
                    default: x = 1 + (r.nextDouble() - 0.5) * 1e-6; break;
                }
                sb.append("log ").append(h(x)).append(' ').append(h(StrictMath.log(x))).append('\n');
                sb.append("log10 ").append(h(x)).append(' ').append(h(StrictMath.log10(x))).append('\n');
                sb.append("log1p ").append(h(x)).append(' ').append(h(StrictMath.log1p(x))).append('\n');
                double e = (r.nextDouble() - 0.5) * 1400 * (i % 3 == 0 ? 1 : 0.01);
                sb.append("exp ").append(h(e)).append(' ').append(h(StrictMath.exp(e))).append('\n');
                sb.append("expm1 ").append(h(e)).append(' ').append(h(StrictMath.expm1(e))).append('\n');
                double y = (r.nextDouble() - 0.5) * 40;
                sb.append("pow ").append(h(10.0)).append(' ').append(h(y)).append(' ').append(h(StrictMath.pow(10.0, y))).append('\n');
                double b = r.nextDouble() * 4;
                sb.append("pow ").append(h(b)).append(' ').append(h(y)).append(' ').append(h(StrictMath.pow(b, y))).append('\n');
            }
        }
        System.out.print(sb);
    }
}
