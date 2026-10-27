// Oracle for nomics-hc stage 4: GATK 4.6.2.0's AssemblyBasedCallerUtils.finalizeRegion
// on given assembly regions and read memberships (from nomics' traversal).
//   java -cp gatk-package-4.6.2.0-local.jar GatkFinalize.java <recal.bam> <regions.tsv> > finalized.tsv
// regions.tsv: "contig<TAB>start<TAB>end<TAB>name:flag..." (active regions, members in region order).
// Output per region: "R contig start end n" then per finalized read, in region order after
// finalization: "name flag start cigar bases quals(phred, comma-separated)".
import htsjdk.samtools.*;
import java.io.*;
import java.util.*;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.read.*;

public class GatkFinalize {
    public static void main(String[] a) throws Exception {
        SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(a[0]));
        SAMFileHeader header = sr.getFileHeader();
        String sample = header.getReadGroups().get(0).getSample();
        IndexedSampleList samples = new IndexedSampleList(sample);
        PrintWriter out = new PrintWriter(new BufferedWriter(new OutputStreamWriter(System.out), 1 << 20));
        try (BufferedReader br = new BufferedReader(new FileReader(a[1]))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] f = line.split("\t");
                String contig = f[0];
                int start = Integer.parseInt(f[1]), end = Integer.parseInt(f[2]);
                AssemblyRegion region = new AssemblyRegion(new SimpleInterval(contig, start, end), true, 100, header);
                SimpleInterval padded = region.getPaddedSpan();
                Map<String, Deque<SAMRecord>> byKey = new HashMap<>();
                try (SAMRecordIterator it = sr.queryOverlapping(contig, padded.getStart(), padded.getEnd())) {
                    while (it.hasNext()) {
                        SAMRecord r = it.next();
                        byKey.computeIfAbsent(r.getReadName() + ":" + r.getFlags(), k -> new ArrayDeque<>()).add(r);
                    }
                }
                List<GATKRead> members = new ArrayList<>();
                for (int i = 3; i < f.length; ++i) {
                    Deque<SAMRecord> q = byKey.get(f[i]);
                    if (q == null || q.isEmpty()) throw new IllegalStateException("read " + f[i] + " not found for " + line.substring(0, 40));
                    members.add(new SAMRecordToGATKReadAdapter(q.poll()));
                }
                region.addAll(members);
                AssemblyBasedCallerUtils.finalizeRegion(region, false, false, (byte) 9, header, samples, true, false, false, true);
                out.println("R " + contig + " " + start + " " + end + " " + region.getReads().size());
                for (GATKRead r : region.getReads()) {
                    StringBuilder q = new StringBuilder();
                    for (byte b : r.getBaseQualities()) { if (q.length() > 0) q.append(','); q.append(b); }
                    out.println(r.getName() + " " + r.convertToSAMRecord(header).getFlags() + " " + r.getStart() + " " + r.getCigar() + " " + new String(r.getBases()) + " " + q);
                }
            }
        }
        out.flush();
    }
}
