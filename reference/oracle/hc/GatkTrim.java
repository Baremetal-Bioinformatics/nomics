// Oracle for nomics-hc stage 6: GATK 4.6.2.0 region trimming after assembly, as callRegion does it
// (getVariationEvents(0) -> AssemblyRegionTrimmer.trim -> AssemblyResultSet.trimTo -> read stubs and
// filterNonPassingReads), on given regions/memberships.
//   java -cp gatk.jar GatkTrim.java <recal.bam> <ref.fasta> <regions.tsv> > trimmed.tsv
// Per region: "R contig start end"; "NOVAR" if no variation; else
// "SPAN <variantSpan> <paddedSpan>", "H isRef alignStart cigar loc bases" per genotyping haplotype (order),
// "READ name flag start cigar" per read of regionForGenotyping after stub removal and filtering, and
// "FILTERED name flag" per read removed by filterNonPassingReads.
import htsjdk.samtools.*;
import java.io.*;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.*;
import org.apache.logging.log4j.LogManager;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.haplotype.*;
import org.broadinstitute.hellbender.utils.read.*;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;

public class GatkTrim {
    public static void main(String[] a) throws Exception {
        SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(a[0]));
        SAMFileHeader header = sr.getFileHeader();
        IndexedSampleList samples = new IndexedSampleList(header.getReadGroups().get(0).getSample());
        CachingIndexedFastaSequenceFile ref = new CachingIndexedFastaSequenceFile(Paths.get(a[1]));
        ReferenceDataSource refSource = ReferenceDataSource.of(Paths.get(a[1]));
        HaplotypeCallerArgumentCollection args = new HaplotypeCallerArgumentCollection();
        AssemblyRegionArgumentCollection regionArgs = new AssemblyRegionArgumentCollection();
        System.err.println("snp padding " + regionArgs.snpPaddingForGenotyping + " indel " + regionArgs.indelPaddingForGenotyping
                + " str " + regionArgs.strPaddingForGenotyping + " maxExt " + regionArgs.maxExtensionIntoRegionPadding + " legacy " + regionArgs.enableLegacyAssemblyRegionTrimming);
        AssemblyRegionTrimmer trimmer = new AssemblyRegionTrimmer(regionArgs, header.getSequenceDictionary());
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
                out.println("R " + contig + " " + start + " " + end);
                ReferenceContext refContext = new ReferenceContext(refSource, padded);
                AssemblyResultSet untrimmed = AssemblyBasedCallerUtils.assembleReads(region, args, header, samples,
                        LogManager.getLogger(GatkTrim.class), ref, engine, aligner, true, new FlowBasedArgumentCollection(), false);
                SortedSet<org.broadinstitute.hellbender.utils.haplotype.Event> events = untrimmed.getVariationEvents(0).stream()
                        .collect(Collectors.toCollection(() -> new TreeSet<>(AssemblyResultSet.HAPLOTYPE_EVENT_COMPARATOR)));
                AssemblyRegionTrimmer.Result tr = trimmer.trim(region, events, refContext);
                if (!tr.isVariationPresent()) { out.println("NOVAR"); continue; }
                AssemblyResultSet rs = untrimmed.trimTo(tr.getVariantRegion());
                AssemblyRegion g = rs.getRegionForGenotyping();
                out.println("SPAN " + g.getSpan() + " " + g.getPaddedSpan() + " " + rs.isVariationPresent());
                for (Haplotype h : rs.getHaplotypeList())
                    out.println("H " + h.isReference() + " " + h.getAlignmentStartHapwrtRef() + " " + h.getCigar() + " " + h.getGenomeLocation() + " " + h.getBaseString());
                List<GATKRead> stubs = g.getReads().stream().filter(r -> AlignmentUtils.unclippedReadLength(r) < 10).collect(Collectors.toList());
                g.removeAll(stubs);
                Set<GATKRead> remove = new LinkedHashSet<>();
                for (GATKRead r : g.getReads())
                    if (AlignmentUtils.unclippedReadLength(r) < 10 || r.getMappingQuality() < 20 || !ReadFilterLibrary.MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE.test(r))
                        remove.add(r);
                g.removeAll(remove);
                for (GATKRead r : g.getReads()) out.println("READ " + r.getName() + " " + r.convertToSAMRecord(header).getFlags() + " " + r.getStart() + " " + r.getCigar());
                for (GATKRead r : remove) out.println("FILTERED " + r.getName() + " " + r.convertToSAMRecord(header).getFlags());
            }
        }
        out.flush();
    }
}
