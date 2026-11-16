// Oracle for nomics-hc genotyping: GATK 4.6.2.0's HaplotypeCallerGenotypingEngine on synthetic
// per-event read x allele likelihood matrices (single sample, ploidy 2, HC defaults).
//   javac -cp gatk.jar -d cls GatkGenotype.java
//   java -cp gatk.jar:cls org.broadinstitute.hellbender.tools.walkers.haplotypecaller.GatkGenotype N SEED
// Per case: "CASE alleles(ref,alt1,...) nReads", "M a r hex" rows (likelihood matrix),
// "GL hex..." (calculateGLsForThisEvent's genotype likelihoods), "PL ...",
// then "NULL" or "CALL qualLog10Hex keptAlleleIdx GT PL GQ MLEAC MLEAF".
package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.*;
import java.util.*;

import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.*;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

public class GatkGenotype {
    static String h(double d) { return Long.toHexString(Double.doubleToRawLongBits(d)); }
    static final byte[] ACGT = {'A', 'C', 'G', 'T'};

    public static void main(String[] a) {
        int n = Integer.parseInt(a[0]);
        Random r = new Random(Long.parseLong(a[1]));
        SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        IndexedSampleList samples = new IndexedSampleList("S");
        HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
        GenotypePriorCalculator gpc = GenotypePriorCalculator.assumingHW(Math.log10(0.001), Math.log10(1.25e-4));
        StringBuilder out = new StringBuilder();
        for (int c = 0; c < n; ++c) {
            // Alleles: a ref of 1-6 bases; 1-3 alts: SNP-like (same length), indels, rarely '*'.
            int refLen = 1 + r.nextInt(r.nextInt(3) == 0 ? 6 : 1);
            byte[] rb = new byte[refLen];
            for (int i = 0; i < refLen; ++i) rb[i] = ACGT[r.nextInt(4)];
            List<Allele> alleles = new ArrayList<>();
            alleles.add(Allele.create(rb, true));
            int nAlt = 1 + (r.nextInt(4) == 0 ? 1 + r.nextInt(2) : 0);
            Set<String> seen = new HashSet<>();
            seen.add(new String(rb));
            while (alleles.size() < 1 + nAlt) {
                Allele alt;
                if (r.nextInt(30) == 0 && !seen.contains("*")) { alt = Allele.SPAN_DEL; seen.add("*"); alleles.add(alt); continue; }
                int len = r.nextInt(2) == 0 ? refLen : 1 + r.nextInt(8);
                byte[] ab = new byte[len];
                ab[0] = rb[0];
                for (int i = 1; i < len; ++i) ab[i] = ACGT[r.nextInt(4)];
                if (len == refLen) ab[len - 1] = ACGT[r.nextInt(4)];
                if (seen.add(new String(ab))) alleles.add(Allele.create(ab, false));
            }
            int nA = alleles.size();
            int nReads = r.nextInt(5) == 0 ? r.nextInt(4) : 1 + r.nextInt(80);
            List<GATKRead> reads = new ArrayList<>();
            for (int i = 0; i < nReads; ++i)
                reads.add(ArtificialReadUtils.createArtificialRead(header, "r" + i, 0, 1000, new byte[] {'A'}, new byte[] {30}, "1M"));
            Map<String, List<GATKRead>> ev = new HashMap<>();
            ev.put("S", reads);
            AlleleLikelihoods<GATKRead, Allele> lk = new AlleleLikelihoods<>(samples, new IndexedAlleleList<>(alleles), ev);
            LikelihoodMatrix<GATKRead, Allele> m = lk.sampleMatrix(0);
            boolean[] unsupported = new boolean[nA];
            for (int al = 1; al < nA; ++al) unsupported[al] = r.nextInt(12) == 0;
            // Truth genotype: hom-ref, het or hom-alt mixtures.
            int g1 = r.nextInt(nA), g2 = r.nextInt(3) == 0 ? g1 : r.nextInt(nA);
            out.append("CASE ");
            for (int al = 0; al < nA; ++al) out.append(al > 0 ? "," : "").append(alleles.get(al).getDisplayString());
            out.append(' ').append(nReads).append('\n');
            for (int i = 0; i < nReads; ++i) {
                int truth = r.nextBoolean() ? g1 : g2;
                for (int al = 0; al < nA; ++al) {
                    double v;
                    if (unsupported[al]) v = Double.NEGATIVE_INFINITY;
                    else if (al == truth) v = -r.nextDouble() * 0.5;
                    else v = r.nextInt(6) == 0 ? -r.nextDouble() * 0.5 : -0.5 - r.nextDouble() * 4.0;
                    if (r.nextInt(50) == 0) v = -4.5; // normalization cap value
                    m.set(al, i, v);
                    out.append("M ").append(al).append(' ').append(i).append(' ').append(h(v)).append('\n');
                }
            }
            HaplotypeCallerGenotypingEngine engine = new HaplotypeCallerGenotypingEngine(hcArgs, samples, false, false);
            VariantContext merged = new VariantContextBuilder("HC0", "1", 1000, 1000 + refLen - 1, alleles).make();
            GenotypesContext gts = engine.calculateGLsForThisEvent(lk, merged, GATKVariantContextUtils.noCallAlleles(2), new byte[] {'A'}, 0, null);
            Genotype g = gts.get(0);
            out.append("GL");
            for (double v : g.getLikelihoods().getAsVector()) out.append(' ').append(h(v));
            out.append("\nPL");
            for (int p : g.getPL()) out.append(' ').append(p);
            out.append('\n');
            VariantContext call = engine.calculateGenotypes(new VariantContextBuilder(merged).genotypes(gts).make(), gpc, Collections.emptyList());
            if (call == null) {
                out.append("NULL\n");
            } else {
                Genotype cg = call.getGenotype(0);
                StringBuilder kept = new StringBuilder();
                for (Allele al : call.getAlleles()) kept.append(kept.length() > 0 ? "," : "").append(alleles.indexOf(al));
                StringBuilder gt = new StringBuilder();
                for (Allele al : cg.getAlleles()) gt.append(gt.length() > 0 ? "/" : "").append(call.getAlleles().indexOf(al));
                out.append("CALL ").append(h(call.getLog10PError())).append(' ').append(kept).append(' ').append(gt).append(' ')
                   .append(Arrays.toString(cg.getPL()).replace(" ", "")).append(' ').append(cg.getGQ()).append(' ')
                   .append(String.valueOf(call.getAttribute("MLEAC")).replace(" ", "")).append(' ')
                   .append(String.valueOf(call.getAttribute("MLEAF")).replace(" ", "")).append('\n');
            }
        }
        System.out.print(out);
    }
}
