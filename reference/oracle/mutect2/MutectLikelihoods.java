// Test oracle only: unmodified jar traversal, finalizer, assembler and trimmer.
// Run only after run_oracle.sh qualifies this JVM on chr21 VCF/stats/regions.
import java.io.*;
import java.util.*;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMM;
import org.apache.logging.log4j.LogManager;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.utils.downsampling.*;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.read.*;

public class MutectLikelihoods extends AssemblyRegionWalker {
  @ArgumentCollection public M2ArgumentCollection mtac = new M2ArgumentCollection();
  @Argument(fullName="dump-prefix", doc="Dump prefix") public String prefix;
  Mutect2Engine engine;
  PrintWriter assembly;
  CachingIndexedFastaSequenceFile reference;
  ReadThreadingAssembler assembler;
  SmithWatermanAligner aligner;
  AssemblyRegionTrimmer trimmer;
  int ordinal;
  public boolean shouldTrackPileupsForAssemblyRegions() { return false; }
  public List<ReadFilter> getDefaultReadFilters() { return Mutect2Engine.makeStandardMutect2ReadFilters(); }
  public ReadTransformer makePostReadFilterTransformer() {
    return Mutect2Engine.makeStandardMutect2PostFilterReadTransformer(referenceArguments.getReferencePath(), true);
  }
  protected ReadsDownsampler createDownsampler() { return new MutectDownsampler(50, 0, 1); }
  public void onTraversalStart() {
    try {
      assembly = new PrintWriter(new BufferedWriter(new FileWriter(prefix+".likelihoods.tsv"), 1<<20));
      reference = new CachingIndexedFastaSequenceFile(referenceArguments.getReferencePath());
      assembler = mtac.createReadThreadingAssembler();
      aligner = SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.FASTEST_AVAILABLE);
      trimmer = new AssemblyRegionTrimmer(assemblyRegionArgs, getHeaderForReads().getSequenceDictionary());
    } catch(IOException e) { throw new RuntimeException(e); }
    engine = new Mutect2Engine(mtac, assemblyRegionArgs, false, false, getHeaderForReads(), getBestAvailableSequenceDictionary(),
       referenceArguments.getReferenceSpecifier(), new VariantAnnotatorEngine(Collections.emptyList(), null, Collections.emptyList(), false, false));
  }
  public AssemblyRegionEvaluator assemblyRegionEvaluator() {
    return engine;
  }
  public void apply(AssemblyRegion region, ReferenceContext ref, FeatureContext features) {
    ++ordinal;
    AssemblyBasedCallerUtils.cleanOverlappingReadPairs(region.getReads(), new IndexedSampleList(ReadUtils.getSamplesFromHeader(getHeaderForReads())),
        getHeaderForReads(), false, OptionalInt.of(20), OptionalInt.of(20));
    if (!region.isActive() || region.size() == 0) return;
    try {
      var remove = Mutect2Engine.class.getDeclaredMethod("removeUnmarkedDuplicates", AssemblyRegion.class);
      remove.setAccessible(true); remove.invoke(engine, region);
    } catch (Exception e) { throw new RuntimeException(e); }
    var rs = AssemblyBasedCallerUtils.assembleReads(region, mtac, getHeaderForReads(),
       new IndexedSampleList(ReadUtils.getSamplesFromHeader(getHeaderForReads())), LogManager.getLogger(MutectLikelihoods.class),
       reference, assembler, aligner, false, mtac.fbargs, false);
    var result = trimmer.trim(region, rs.getVariationEvents(mtac.maxMnpDistance), ref);
    if (!result.isVariationPresent()) return;
    var trimmed = rs.trimTo(result.getVariantRegion());
    if (!trimmed.isVariationPresent()) return;
    var genotyping = trimmed.getRegionForGenotyping();
    int beforeStubs = genotyping.size();
    genotyping.removeAll(genotyping.getReads().stream().filter(r -> r.getLength() < 10).toList());
    assembly.println("R " + ordinal + " " + genotyping.getPaddedSpan() + " " + (beforeStubs-genotyping.size()));
    dump("H", trimmed);
    try {
      var field = Mutect2Engine.class.getDeclaredField("likelihoodCalculationEngine"); field.setAccessible(true);
      var calc = (PairHMMLikelihoodCalculationEngine) field.get(engine);
      var samples = new IndexedSampleList(ReadUtils.getSamplesFromHeader(getHeaderForReads()));
      Map<String,List<GATKRead>> bySample = new HashMap<>();
      for (GATKRead r : genotyping.getReads()) bySample.computeIfAbsent(ReadUtils.getSampleName(r,getHeaderForReads()), k -> new ArrayList<>()).add(r);
      var lk = new AlleleLikelihoods<GATKRead,Haplotype>(samples,new IndexedAlleleList<>(trimmed.getHaplotypeList()),bySample);
      var init = PairHMMLikelihoodCalculationEngine.class.getDeclaredMethod("initializePairHMM", List.class, Map.class); init.setAccessible(true);
      init.invoke(calc, trimmed.getHaplotypeList(), bySample);
      var prepare = PairHMMLikelihoodCalculationEngine.class.getDeclaredMethod("modifyReadQualities",List.class); prepare.setAccessible(true);
      var compute = PairHMMLikelihoodCalculationEngine.class.getDeclaredMethod("computeReadLikelihoods",LikelihoodMatrix.class); compute.setAccessible(true);
      for (int n=0;n<samples.numberOfSamples();n++) {
        assembly.println("S " + samples.getSample(n) + " " + lk.sampleEvidenceCount(n));
        @SuppressWarnings("unchecked") var processed = (List<GATKRead>)prepare.invoke(calc,lk.sampleEvidence(n));
        for (int i=0;i<processed.size();i++) {
          GATKRead r=processed.get(i), orig=lk.sampleEvidence(n).get(i);
          assembly.println("D " + orig.getName() + " " + orig.convertToSAMRecord(getHeaderForReads()).getFlags() + " " + orig.getStart() + " " + orig.getCigar() + " " + new String(r.getBases()) + " " + HexFormat.of().formatHex(r.getBaseQualities()) + " " + HexFormat.of().formatHex(ReadUtils.getBaseInsertionQualities(r)) + " " + HexFormat.of().formatHex(ReadUtils.getBaseDeletionQualities(r)) + " " + "0a".repeat(r.getLength()));
        }
        compute.invoke(calc,lk.sampleMatrix(n));
      }
      matrix("RAW",lk);
      lk.normalizeLikelihoods(AssemblyBasedCallerUtils.getGlobalMismatchingRateFromArgs(mtac.likelihoodArgs),!mtac.likelihoodArgs.disableSymmetricallyNormalizeAllelesToReference);
      matrix("CAP",lk);
      calc.filterPoorlyModeledEvidence(lk,mtac.likelihoodArgs.enableDynamicReadDisqualification,mtac.likelihoodArgs.expectedErrorRatePerBase,mtac.likelihoodArgs.readDisqualificationThresholdConstant);
      var publicResult=calc.computeReadLikelihoods(trimmed,samples,bySample,true);
      for(int n=0;n<samples.numberOfSamples();n++) {
        if(!lk.sampleEvidence(n).equals(publicResult.sampleEvidence(n))) throw new IllegalStateException("public engine read order differs");
        for(int h=0;h<lk.numberOfAlleles();h++) for(int r=0;r<lk.sampleEvidenceCount(n);r++)
          if(Double.doubleToRawLongBits(lk.sampleMatrix(n).get(h,r))!=Double.doubleToRawLongBits(publicResult.sampleMatrix(n).get(h,r))) throw new IllegalStateException("public engine bits differ");
      }
      matrix("KEEP",lk);
      lk.switchToNaturalLog(); matrix("LN",lk);
      var realigned=AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype(lk,trimmed.getReferenceHaplotype(),trimmed.getPaddedReferenceLoc(),aligner,mtac.getReadToHaplotypeSWParameters());
      lk.changeEvidence(realigned);
      for (int n=0;n<samples.numberOfSamples();n++) for(GATKRead r:lk.sampleEvidence(n))
        assembly.println("A " + samples.getSample(n) + " " + r.getName() + " " + r.convertToSAMRecord(getHeaderForReads()).getFlags() + " " + r.getStart() + " " + r.getCigar() + " " + r.getAttributeAsString("HC"));
    } catch(Exception e) {throw new RuntimeException(e);}
  }
  void matrix(String mark, AlleleLikelihoods<GATKRead,Haplotype> lk) {
    for(int n=0;n<lk.numberOfSamples();n++) {
      assembly.println(mark + " " + lk.getSample(n) + " " + lk.sampleEvidenceCount(n));
      for (int i=0;i<lk.sampleEvidenceCount(n);i++) assembly.println("E " + lk.sampleEvidence(n).get(i).getName());
      var m=lk.sampleMatrix(n);
      for(int h=0;h<lk.numberOfAlleles();h++) {
        assembly.print("L " + h);
        for(int i=0;i<lk.sampleEvidenceCount(n);i++) assembly.printf(" %016x",Double.doubleToRawLongBits(m.get(h,i)));
        assembly.println();
      }
    }
  }
  void dump(String marker, AssemblyResultSet rs) {
    assembly.println(marker + "S " + rs.getHaplotypeList().size());
    for (Haplotype h : rs.getHaplotypeList()) assembly.println(marker + " " + h.isReference() + " " + h.getAlignmentStartHapwrtRef() + " " + h.getCigar() + " "
       + h.getGenomeLocation() + " " + h.getScore() + " " + h.getKmerSize() + " " + h.getBaseString());
  }
  public Object onTraversalSuccess() { engine.writeExtraOutputs(new File(prefix+".stats")); return "SUCCESS"; }
  public void closeTool() {
    if (assembly != null) assembly.close();
    if (aligner != null) aligner.close();
    if (reference != null) reference.close();
    if (engine != null) engine.close();
  }
  public static void main(String[] args) { new MutectLikelihoods().instanceMain(args); }
}
