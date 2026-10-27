// Oracle for nomics-hc stage 5: GATK 4.6.2.0's AssemblyBasedCallerUtils.assembleReads
// (finalizeRegion + ReadThreadingAssembler.runLocalAssembly) on given regions/memberships.
//   java -cp gatk-package-4.6.2.0-local.jar GatkAssemble.java <recal.bam> <ref.fasta> <regions.tsv> > assembled.tsv
// regions.tsv as for GatkFinalize. Output per region: "R contig start end nHaplotypes" then
// "REF paddedRefLoc bases" (reference with 500 bp padding), then
// one line per haplotype in AssemblyResultSet.getHaplotypeList() order:
// "H isRef alignmentStartHapwrtRef cigar genomeLoc score kmerSize bases".
import htsjdk.samtools.*;
import java.io.*;
import java.nio.file.Paths;
import java.util.*;
import org.apache.logging.log4j.LogManager;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.*;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;

public class GatkAssemble {
    public static void main(String[] a) throws Exception {
        SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(a[0]));
        SAMFileHeader header = sr.getFileHeader();
        IndexedSampleList samples = new IndexedSampleList(header.getReadGroups().get(0).getSample());
        CachingIndexedFastaSequenceFile ref = new CachingIndexedFastaSequenceFile(Paths.get(a[1]));
        HaplotypeCallerArgumentCollection args = new HaplotypeCallerArgumentCollection();
        ReadThreadingAssembler engine = args.createReadThreadingAssembler();
        SmithWatermanAligner aligner = SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.FASTEST_AVAILABLE);
        PrintWriter out = new PrintWriter(new BufferedWriter(new OutputStreamWriter(System.out), 1 << 20));
        try (BufferedReader br = new BufferedReader(new FileReader(a[2]))) {
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
                for (int i = 3; i < f.length; ++i) members.add(new SAMRecordToGATKReadAdapter(byKey.get(f[i]).poll()));
                region.addAll(members);
                AssemblyResultSet rs = AssemblyBasedCallerUtils.assembleReads(region, args, header, samples,
                        LogManager.getLogger(GatkAssemble.class), ref, engine, aligner, true, new FlowBasedArgumentCollection(), false);
                List<Haplotype> haps = rs.getHaplotypeList();
                out.println("R " + contig + " " + start + " " + end + " " + haps.size());
                // The reference the assembler saw: region reference with 500 bp of padding.
                byte[] fullRef = region.getAssemblyRegionReference(ref, 500);
                SimpleInterval refLoc = AssemblyBasedCallerUtils.getPaddedReferenceLoc(region, 500, ref);
                out.println("REF " + refLoc + " " + new String(fullRef));
                for (Haplotype h : haps) {
                    out.println("H " + h.isReference() + " " + h.getAlignmentStartHapwrtRef() + " " + h.getCigar() + " "
                            + h.getGenomeLocation() + " " + h.getScore() + " " + h.getKmerSize() + " " + h.getBaseString());
                }
            }
        }
        out.flush();
    }
}
