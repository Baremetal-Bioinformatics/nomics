// Oracle for nomics-hc stages 7-8 (GatkCall): as GatkLikelihoods through realignment, then
// HaplotypeCallerGenotypingEngine.assignGenotypeLikelihoods with the HC annotation set (no dbSNP),
// on given regions/memberships.
//   java -cp gatk.jar GatkCall.java <recal.bam> <ref.fasta> <regions.tsv> > calls.tsv
// Per region: "R contig start end"; then for every annotateContext call:
//   "VC <VCF line of the unannotated call>" (QUAL/GT/PL/GQ/MLEAC/MLEAF),
//   "A <allele>..." (annotation matrix alleles), per evidence
//   "E name flag mapq reverse start cigar bases quals hex..." (realigned read; 0.0 rows are the
//   appended filtered reads). The VC line carries the genotype annotations (AD, DP), which run first.
// Finally "CALL <VCF line>" per call returned to callRegion (after allele trimming and finalizeRawMQ).
// The shared RNG is reset when the call is dumped (before QD), so QD jitter is the first nextGaussian of
// java.util.Random(47382911).
import htsjdk.samtools.*;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import java.io.*;
import java.nio.file.Paths;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.*;
import org.apache.logging.log4j.LogManager;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.haplotype.*;
import org.broadinstitute.hellbender.utils.read.*;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.variant.*;

public class GatkCall {
    static PrintWriter out;
    static VCFEncoder enc;
    static SAMFileHeader header;

    static String h(double d) { return Long.toHexString(Double.doubleToRawLongBits(d)); }

    // First info annotation: sees the call (with genotype annotations) and the annotation matrix.
    public static class Dump implements InfoFieldAnnotation {
        @Override
        public Map<String, Object> annotate(ReferenceContext ref, VariantContext vc, AlleleLikelihoods<GATKRead, Allele> lk) {
            out.println("VC " + enc.encode(vc));
            StringBuilder sb = new StringBuilder("A");
            for (Allele a : lk.alleles()) sb.append(' ').append(a.getDisplayString());
            out.println(sb);
            LikelihoodMatrix<GATKRead, Allele> m = lk.sampleMatrix(0);
            for (int r = 0; r < m.evidenceCount(); ++r) {
                GATKRead read = m.evidence().get(r);
                StringBuilder q = new StringBuilder();
                for (byte b : read.getBaseQualities()) { if (q.length() > 0) q.append(','); q.append(b); }
                StringBuilder e = new StringBuilder("E " + read.getName() + " " + read.convertToSAMRecord(header).getFlags() + " "
                        + read.getMappingQuality() + " " + read.isReverseStrand() + " " + read.getStart() + " " + read.getCigar() + " "
                        + new String(read.getBases()) + " " + q);
                for (int a = 0; a < m.numberOfAlleles(); ++a) e.append(' ').append(h(m.get(a, r)));
                out.println(e);
            }
            Utils.resetRandomGenerator();
            return Collections.emptyMap();
        }

        @Override
        public List<String> getKeyNames() { return Collections.emptyList(); }

        @Override
        public List<VCFCompoundHeaderLine> getDescriptions() { return Collections.emptyList(); }
    }

    public static void main(String[] a) throws Exception {
        SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(a[0]));
        header = sr.getFileHeader();
        String sample = header.getReadGroups().get(0).getSample();
        IndexedSampleList samples = new IndexedSampleList(sample);
        CachingIndexedFastaSequenceFile ref = new CachingIndexedFastaSequenceFile(Paths.get(a[1]));
        ReferenceDataSource refSource = ReferenceDataSource.of(Paths.get(a[1]));
        HaplotypeCallerArgumentCollection args = new HaplotypeCallerArgumentCollection();
        AssemblyRegionArgumentCollection regionArgs = new AssemblyRegionArgumentCollection();
        AssemblyRegionTrimmer trimmer = new AssemblyRegionTrimmer(regionArgs, header.getSequenceDictionary());
        ReadThreadingAssembler engine = args.createReadThreadingAssembler();
        ReadLikelihoodCalculationEngine lkEngine = AssemblyBasedCallerUtils.createLikelihoodCalculationEngine(args.likelihoodArgs, args.fbargs, !args.softClipLowQualityEnds, args.likelihoodArgs.likelihoodEngineImplementation);
        SmithWatermanAligner aligner = SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.FASTEST_AVAILABLE);

        List<Annotation> anns = Arrays.asList(new Dump(), new Coverage(), new ChromosomeCounts(), new BaseQualityRankSumTest(),
                new ExcessHet(), new FisherStrand(), new InbreedingCoeff(), new RMSMappingQuality(),
                new MappingQualityRankSumTest(), new QualByDepth(), new ReadPosRankSumTest(), new StrandOddsRatio(),
                new DepthPerAlleleBySample(), new DepthPerSampleHC());
        VariantAnnotatorEngine annEngine = new VariantAnnotatorEngine(anns, null, Collections.emptyList(), false, false);
        Set<VCFHeaderLine> lines = new LinkedHashSet<>(annEngine.getVCFAnnotationDescriptions(false));
        lines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        lines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
        VCFStandardHeaderLines.addStandardFormatLines(lines, true, VCFConstants.GENOTYPE_KEY, VCFConstants.GENOTYPE_QUALITY_KEY,
                VCFConstants.DEPTH_KEY, VCFConstants.GENOTYPE_PL_KEY);
        lines.add(new VCFFilterHeaderLine("LowQual", "Low quality"));
        VCFHeader vh = new VCFHeader(lines, Collections.singletonList(sample));
        vh.setSequenceDictionary(header.getSequenceDictionary());
        enc = new VCFEncoder(vh, true, false);

        HaplotypeCallerGenotypingEngine gt = new HaplotypeCallerGenotypingEngine(args, samples, false, args.applyBQD);
        gt.setAnnotationEngine(annEngine);

        out = new PrintWriter(new BufferedWriter(new OutputStreamWriter(System.out), 1 << 20));
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
                        LogManager.getLogger(GatkCall.class), ref, engine, aligner, true, new FlowBasedArgumentCollection(), false);
                SortedSet<Event> events = untrimmed.getVariationEvents(0).stream()
                        .collect(Collectors.toCollection(() -> new TreeSet<>(AssemblyResultSet.HAPLOTYPE_EVENT_COMPARATOR)));
                AssemblyRegionTrimmer.Result tr = trimmer.trim(region, events, refContext);
                if (!tr.isVariationPresent()) continue;
                AssemblyResultSet rs = untrimmed.trimTo(tr.getVariantRegion());
                AssemblyRegion g = rs.getRegionForGenotyping();
                List<GATKRead> stubs = g.getReads().stream().filter(r -> AlignmentUtils.unclippedReadLength(r) < 10).collect(Collectors.toList());
                g.removeAll(stubs);
                List<GATKRead> remove = new ArrayList<>();
                for (GATKRead r : g.getReads())
                    if (AlignmentUtils.unclippedReadLength(r) < 10 || r.getMappingQuality() < 20 || !ReadFilterLibrary.MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE.test(r))
                        remove.add(r);
                g.removeAll(remove);
                if (!rs.isVariationPresent() || g.size() == 0) continue;
                Map<String, List<GATKRead>> byS = new HashMap<>();
                byS.put(sample, new ArrayList<>(g.getReads()));
                Map<String, List<GATKRead>> filteredByS = new HashMap<>();
                filteredByS.put(sample, remove);
                AlleleLikelihoods<GATKRead, Haplotype> lk = lkEngine.computeReadLikelihoods(rs, samples, byS, true);
                Map<GATKRead, GATKRead> realigned = AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype(lk, rs.getReferenceHaplotype(), rs.getPaddedReferenceLoc(), aligner, args.getReadToHaplotypeSWParameters());
                lk.changeEvidence(realigned);
                CalledHaplotypes called = gt.assignGenotypeLikelihoods(lk.alleles(), lk, filteredByS,
                        rs.getFullReferenceWithPadding(), rs.getPaddedReferenceLoc(), g.getSpan(),
                        new FeatureContext(), Collections.emptyList(), false, args.maxMnpDistance, header, false, new HashSet<>(), lk);
                for (VariantContext vc : called.getCalls())
                    out.println("CALL " + enc.encode(RMSMappingQuality.getInstance().finalizeRawMQ(vc)));
            }
        }
        out.flush();
    }
}
