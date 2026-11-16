// Test-only independent inputs and raw-bit outputs from the unmodified jar.
import java.util.*;
import java.lang.reflect.*;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.broadinstitute.hellbender.tools.walkers.mutect.*;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.*;

public class MutectSomaticNumerics {
    static void bits(double x) { System.out.printf(" %016x", Double.doubleToRawLongBits(x)); }
    public static void main(String[] args) throws Exception {
        Random random = new Random(174728);
        for (int k = 0; k < 80; k++) {
            int a = 1 + k % 5, n = 1 + k * 3 % 71;
            double[][] data = new double[a][n];
            for (int i=0;i<a;i++) for (int j=0;j<n;j++) {
                data[i][j] = k < 5 ? (i == j % a ? -0.01 : -10.0) : -random.nextDouble() * (k % 4 == 0 ? 1000 : 30);
                if (k % 7 == 0 && i > 0 && j % 3 == 0) data[i][j] = Double.NEGATIVE_INFINITY;
            }
            double[] prior = new double[a];
            for (int i=0;i<a;i++) prior[i] = k % 3 == 0 ? 1.0 : 0.25 + random.nextDouble() * 10;
            var matrix = new Array2DRowRealMatrix(data);
            System.out.print(a + " " + n);
            for (double[] row : data) for (double x : row) bits(x);
            for (double x : prior) bits(x);
            for (double x : SomaticLikelihoodsEngine.alleleFractionsPosterior(matrix, prior)) bits(x);
            bits(SomaticLikelihoodsEngine.logEvidence(matrix, prior));
            bits(SomaticLikelihoodsEngine.logDirichletNormalization(prior));

            // Exercise the actual protected/private genotyper methods, rather
            // than reimplementing their formulas in the fixture generator.
            var engine = new SomaticGenotypingEngine(new M2ArgumentCollection(), Set.of(), null, null, null);
            List<Allele> alleles = new ArrayList<>();
            for (int i=0;i<a;i++) alleles.add(Allele.create("A" + "C".repeat(i), i == 0));
            List<GATKRead> reads = new ArrayList<>();
            for (int j=0;j<n;j++) {
                var read = new SAMRecordToGATKReadAdapter(new htsjdk.samtools.SAMRecord(null));
                read.setName("r" + j); read.setPosition("chr1", 100+j);
                read.setBases("A".getBytes()); read.setBaseQualities(new byte[]{30}); read.setCigar("1M");
                reads.add(read);
            }
            var lk = new AlleleLikelihoods<GATKRead,Allele>(new IndexedSampleList("S"),new IndexedAlleleList<>(alleles),Map.of("S",reads));
            for (int i=0;i<a;i++) for(int j=0;j<n;j++) lk.sampleMatrix(0).set(i,j,data[i][j]);
            for (String method : List.of("somaticLogOdds", "diploidAltLogOdds")) {
                Method m = SomaticGenotypingEngine.class.getDeclaredMethod(method, LikelihoodMatrix.class); m.setAccessible(true);
                @SuppressWarnings("unchecked") var result = (PerAlleleCollection<Double>)m.invoke(engine,lk.sampleMatrix(0));
                for(int i=1;i<a;i++) bits(result.getAlt(alleles.get(i)));
            }
            System.out.println();
        }
    }
}
