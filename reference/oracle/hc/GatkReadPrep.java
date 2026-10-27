// Oracle for nomics-hc readprep: GATK 4.6.2.0's ReadLikelihoodCalculationEngine.findTandemRepeatUnits
// and PairHMMLikelihoodCalculationEngine.getErrorModelAdjustedQual (CONSERVATIVE, rate 3.0).
//   javac -cp gatk.jar -d cls GatkReadPrep.java && java -cp gatk.jar:cls org.broadinstitute.hellbender.tools.walkers.haplotypecaller.GatkReadPrep N SEED
// Lines: "C i qual" for the PCR cache; "T bases rl0,rl1,..." per read (repeat length at every offset).
package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import java.util.*;

public class GatkReadPrep {
    public static void main(String[] a) {
        StringBuilder out = new StringBuilder();
        for (int i = 0; i <= 20; ++i) out.append("C ").append(i).append(' ').append(PairHMMLikelihoodCalculationEngine.getErrorModelAdjustedQual(i, 3.0)).append('\n');
        int n = Integer.parseInt(a[0]);
        Random r = new Random(Long.parseLong(a[1]));
        byte[] acgt = {'A', 'C', 'G', 'T'};
        for (int k = 0; k < n; ++k) {
            int len = 1 + r.nextInt(160);
            byte[] b = new byte[len];
            int i = 0;
            while (i < len) {
                if (r.nextInt(3) == 0) { // tandem repeat block
                    int ul = 1 + r.nextInt(r.nextInt(4) == 0 ? 9 : 4);
                    byte[] unit = new byte[ul];
                    for (int u = 0; u < ul; ++u) unit[u] = acgt[r.nextInt(4)];
                    int reps = 1 + r.nextInt(12);
                    for (int t = 0; t < reps * ul && i < len; ++t) b[i++] = unit[t % ul];
                } else {
                    b[i++] = r.nextInt(100) == 0 ? (byte) 'N' : acgt[r.nextInt(4)];
                }
            }
            out.append("T ").append(new String(b)).append(' ');
            for (int o = 0; o < len; ++o) {
                if (o > 0) out.append(',');
                out.append(ReadLikelihoodCalculationEngine.findTandemRepeatUnits(b, o).getRight());
            }
            out.append('\n');
        }
        System.out.print(out);
    }
}
