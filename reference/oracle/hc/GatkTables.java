// Dump GATK 4.6.2.0 numeric tables/constants from the reference jar, for
// bit-exact comparison with nomics-hc. Run in the reference container:
//   java -cp /opt/conda/share/gatk4-4.6.2.0-0/gatk-package-4.6.2.0-local.jar GatkTables.java > gatk_tables.tsv
// Lines: "<name> <args...> <hex result>" (hex = raw double bits).
import java.util.Random;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.CombinatoricsUtils;

public class GatkTables {
    static String h(double d) { return Long.toHexString(Double.doubleToRawLongBits(d)); }

    public static void main(String[] a) {
        StringBuilder sb = new StringBuilder();
        for (int q = 0; q < 255; ++q) { // the qualToProbLog10 table has 255 entries
            sb.append("qualToProbLog10 ").append(q).append(' ').append(h(QualityUtils.qualToProbLog10((byte) q))).append('\n');
            sb.append("qualToErrorProbLog10 ").append(q).append(' ').append(h(QualityUtils.qualToErrorProbLog10((byte) q))).append('\n');
            sb.append("qualToErrorProb ").append(q).append(' ').append(h(QualityUtils.qualToErrorProb((byte) q))).append('\n');
            sb.append("qualToProb ").append(q).append(' ').append(h(QualityUtils.qualToProb((byte) q))).append('\n');
        }
        for (double x : new double[] {0.01, 1.01, 2.01, 10.0, 11.0, 12.0, 0.00125, 1.00125, 2.00125, 1.0, 2.0, 3.0})
            sb.append("logGamma ").append(h(x)).append(' ').append(h(Gamma.logGamma(x))).append('\n');
        for (int n = 0; n <= 2; ++n)
            sb.append("binomialCoefficientLog 2 ").append(n).append(' ').append(h(CombinatoricsUtils.binomialCoefficientLog(2, n))).append('\n');
        // approximateLog10SumLog10 on the table grid, near grid midpoints, and random pairs.
        Random r = new Random(5);
        for (int k = 0; k <= 80000; ++k) {
            double b = -r.nextDouble() * 50, d = k * 0.0001;
            sb.append("approxSum ").append(h(b - d)).append(' ').append(h(b)).append(' ').append(h(MathUtils.approximateLog10SumLog10(b - d, b))).append('\n');
        }
        for (int i = 0; i < 200000; ++i) {
            double x = -r.nextDouble() * 60, y = x - r.nextDouble() * (i % 2 == 0 ? 9 : 0.01);
            if (i % 7 == 0) y = x - (Math.floor(r.nextDouble() * 80000) + 0.5) * 0.0001; // ties of fastRound
            sb.append("approxSum ").append(h(x)).append(' ').append(h(y)).append(' ').append(h(MathUtils.approximateLog10SumLog10(x, y))).append('\n');
        }
        for (int i = 0; i <= 100; ++i)
            sb.append("normalDistribution 50 17 ").append(i).append(' ').append(h(MathUtils.normalDistribution(50, 17, i))).append('\n');
        // AlleleFrequencyCalculator.makeCalculator(default GenotypeCalculationArgumentCollection):
        // the exact pseudocounts and, for the biallelic single-sample posterior (ploidy 2), the
        // per-genotype prior term logToLog10(binomialCoefficientLog(2, n) + logGamma(n + snp) + logGamma(2 - n + ref)).
        {
            double snpHet = 0.001, hetSd = 0.01, indelHet = 1.25e-4;
            double refPc = snpHet / Math.pow(hetSd, 2);
            double snpPc = snpHet * refPc;
            double indelPc = indelHet * refPc;
            sb.append("refPseudocount ").append(h(refPc)).append('\n');
            sb.append("snpPseudocount ").append(h(snpPc)).append('\n');
            sb.append("indelPseudocount ").append(h(indelPc)).append('\n');
            for (int n = 0; n <= 2; ++n) {
                double term = MathUtils.logToLog10(CombinatoricsUtils.binomialCoefficientLog(2, n) + Gamma.logGamma(n + snpPc) + Gamma.logGamma(2 - n + refPc));
                sb.append("biallelicPriorTerm ").append(n).append(' ').append(h(term)).append('\n');
                sb.append("logGammaSnp ").append(n).append(' ').append(h(Gamma.logGamma(n + snpPc))).append('\n');
                sb.append("logGammaRef ").append(n).append(' ').append(h(Gamma.logGamma(2 - n + refPc))).append('\n');
            }
        }
        sb.append("LOG10_ONE_THIRD ").append(h(MathUtils.LOG10_ONE_THIRD)).append('\n');
        sb.append("LOG10_E ").append(h(MathUtils.LOG10_E)).append('\n');
        sb.append("log10_2 ").append(h(Math.log10(2))).append('\n');
        System.out.print(sb);
    }
}
