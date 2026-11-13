// Test oracle only: unmodified jar traversal, finalizer, assembler and trimmer.
// Run only after run_oracle.sh qualifies this JVM on chr21 VCF/stats/regions.
import java.io.*;
import java.util.*;
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

public class MutectAssembly extends AssemblyRegionWalker {
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
      assembly = new PrintWriter(new BufferedWriter(new FileWriter(prefix+".assembly.tsv"), 1<<20));
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
       new IndexedSampleList(ReadUtils.getSamplesFromHeader(getHeaderForReads())), LogManager.getLogger(MutectAssembly.class),
       reference, assembler, aligner, false, mtac.fbargs, false);
    assembly.println("R " + ordinal + " " + region.getSpan() + " " + region.getPaddedSpan() + " " + region.size());
    assembly.println("REF " + rs.getPaddedReferenceLoc() + " " + new String(rs.getFullReferenceWithPadding()));
    for (GATKRead read : region.getReads()) {
      assembly.println("D " + ReadUtils.getSampleName(read,getHeaderForReads()) + " " + read.getName() + " "
        + read.convertToSAMRecord(getHeaderForReads()).getFlags() + " " + read.getStart() + " " + read.getCigar()
        + " " + new String(read.getBases()) + " " + HexFormat.of().formatHex(read.getBaseQualities())
        + " " + HexFormat.of().formatHex(ReadUtils.getBaseInsertionQualities(read))
        + " " + HexFormat.of().formatHex(ReadUtils.getBaseDeletionQualities(read)));
    }
    dump("H", rs);
    var result = trimmer.trim(region, rs.getVariationEvents(mtac.maxMnpDistance), ref);
    if (result.isVariationPresent()) {
      var trimmed = rs.trimTo(result.getVariantRegion());
      assembly.println("T " + trimmed.getRegionForGenotyping().getSpan() + " " + trimmed.getRegionForGenotyping().getPaddedSpan() + " " + trimmed.isVariationPresent());
      dump("h", trimmed);
    } else assembly.println("T none");
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
  public static void main(String[] args) { new MutectAssembly().instanceMain(args); }
}
