// Test oracle only: unmodified jar engine, traversal, filters and downsampler.
import java.io.*;
import java.util.*;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.downsampling.*;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.read.*;

public class MutectRegions extends AssemblyRegionWalker {
  @ArgumentCollection public M2ArgumentCollection mtac = new M2ArgumentCollection();
  @Argument(fullName="dump-prefix", doc="Dump prefix") public String prefix;
  Mutect2Engine engine;
  PrintWriter activities, members;
  DataOutputStream prepared;
  int ordinal;
  public boolean shouldTrackPileupsForAssemblyRegions() { return false; }
  public List<ReadFilter> getDefaultReadFilters() { return Mutect2Engine.makeStandardMutect2ReadFilters(); }
  public ReadTransformer makePostReadFilterTransformer() {
    return Mutect2Engine.makeStandardMutect2PostFilterReadTransformer(referenceArguments.getReferencePath(), true);
  }
  protected ReadsDownsampler createDownsampler() { return new MutectDownsampler(50, 0, 1); }
  public void onTraversalStart() {
    try {
      prepared = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(prefix+".prepared.bin")));
      activities = new PrintWriter(prefix + ".activity.tsv"); members = new PrintWriter(prefix + ".reads.tsv");
    } catch(IOException e) { throw new RuntimeException(e); }
    engine = new Mutect2Engine(mtac, assemblyRegionArgs, false, false, getHeaderForReads(), getBestAvailableSequenceDictionary(),
       referenceArguments.getReferenceSpecifier(), new VariantAnnotatorEngine(Collections.emptyList(), null, Collections.emptyList(), false, false));
  }
  public AssemblyRegionEvaluator assemblyRegionEvaluator() {
    return (context, ref, features) -> {
      ActivityProfileState state = engine.isActive(context, ref, features);
      activities.println(ref.getInterval().getContig()+"\t"+ref.getInterval().getStart()+"\t"+(int)state.isActiveProb());
      return state;
    };
  }
  public void apply(AssemblyRegion region, ReferenceContext ref, FeatureContext features) {
    ++ordinal;
    for (GATKRead read : region.getReads()) {
      members.println(ordinal+"\t"+ReadUtils.getSampleName(read,getHeaderForReads())+"\t"+read.getName()+"\t"+read.convertToSAMRecord(getHeaderForReads()).getFlags()+"\t"+read.getStart()+"\t"+read.getCigar());
    }
    // callRegion's only mutation before assembly; preserve shared-read timing.
    AssemblyBasedCallerUtils.cleanOverlappingReadPairs(region.getReads(), new IndexedSampleList(ReadUtils.getSamplesFromHeader(getHeaderForReads())),
        getHeaderForReads(), false, OptionalInt.of(20), OptionalInt.of(20));
    try { for (GATKRead read : region.getReads()) {
      prepared.writeInt(read.getLength()); prepared.write(read.getBases()); prepared.write(read.getBaseQualities());
      prepared.write(ReadUtils.getBaseInsertionQualities(read)); prepared.write(ReadUtils.getBaseDeletionQualities(read));
    } } catch(IOException e) { throw new RuntimeException(e); }
  }
  public Object onTraversalSuccess() { engine.writeExtraOutputs(new File(prefix+".stats")); return "SUCCESS"; }
  public void closeTool() { try { if(prepared!=null)prepared.close(); } catch(IOException e){throw new RuntimeException(e); } if(activities!=null) activities.close(); if(members!=null) members.close(); if(engine!=null)engine.close(); }
  public static void main(String[] args) { new MutectRegions().instanceMain(args); }
}
