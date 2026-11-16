// Test-only read-only instrumentation around the unmodified GATK genotyper.
// The simple class name deliberately preserves the standard Mutect2 header.
package nomics.oracle;
import java.io.*;
import java.lang.reflect.*;
import java.util.*;
import htsjdk.samtools.*;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.downsampling.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.*;
import org.broadinstitute.hellbender.tools.walkers.mutect.*;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.haplotype.*;
import org.broadinstitute.hellbender.utils.read.*;

public class Mutect2 extends AssemblyRegionWalker {
    @ArgumentCollection public M2ArgumentCollection MTAC = new M2ArgumentCollection();
    @Argument(fullName="output",shortName="O",doc="VCF output") public String output;
    @Argument(fullName="calls-dump", doc="Test-only pre-annotation dump") public String dumpPath;
    PrintWriter dump;
    Mutect2Engine engine;
    VariantContextWriter writer;
    static Object field(Object obj, Class<?> cls, String name) throws Exception {
        var f=cls.getDeclaredField(name); f.setAccessible(true); return f.get(obj);
    }
    public String getVersion() { return new org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2().getVersion(); }
    public boolean useVariantAnnotations() {return true;}
    public boolean shouldTrackPileupsForAssemblyRegions() {return false;}
    public List<ReadFilter> getDefaultReadFilters() {return Mutect2Engine.makeStandardMutect2ReadFilters();}
    public ReadTransformer makePostReadFilterTransformer() {return Mutect2Engine.makeStandardMutect2PostFilterReadTransformer(referenceArguments.getReferencePath(),true);}
    public List<Class<? extends Annotation>> getDefaultVariantAnnotationGroups() {return Mutect2Engine.getStandardMutect2AnnotationGroups();}
    protected ReadsDownsampler createDownsampler() {return new MutectDownsampler(50,0,1);}
    public AssemblyRegionEvaluator assemblyRegionEvaluator() {return engine;}
    public void onTraversalStart() {
        var annotation=new VariantAnnotatorEngine(makeVariantAnnotations(),null,Collections.emptyList(),false,false);
        engine=new Mutect2Engine(MTAC,assemblyRegionArgs,false,false,getHeaderForReads(),getBestAvailableSequenceDictionary(),referenceArguments.getReferenceSpecifier(),annotation);
        writer=createVCFWriter(new File(output));
        engine.writeHeader(writer,getDefaultToolVCFHeaderLines());
        try {
            dump=new PrintWriter(new BufferedWriter(new FileWriter(dumpPath),1<<20));
            var annotations=new PrintWriter(new BufferedWriter(new FileWriter(dumpPath+".annotations.tsv"),1<<20));
            @SuppressWarnings("unchecked") var normals=(Set<String>)field(engine,Mutect2Engine.class,"normalSamples");
            var f=Mutect2Engine.class.getDeclaredField("genotypingEngine"); f.setAccessible(true);
            ((SomaticGenotypingEngine)f.get(engine)).close();
            f.set(engine,new LoggingEngine(MTAC,normals,annotation,getHeaderForReads(),getBestAvailableSequenceDictionary(),dump,annotations));
        } catch(Exception e) { throw new RuntimeException(e); }
    }
    public void apply(AssemblyRegion r,ReferenceContext rc,FeatureContext fc) {engine.callRegion(r,rc,fc).forEach(writer::add);}
    public Object onTraversalSuccess() {engine.writeExtraOutputs(new File(output+".stats"));return "SUCCESS";}
    public void closeTool() {if(writer!=null)writer.close();if(engine!=null)engine.close();if(dump!=null)dump.close();}
    public static void main(String[] args) { new Mutect2().instanceMain(args); }
}
class LoggingEngine extends SomaticGenotypingEngine {
    final PrintWriter out;
    final PrintWriter annotations;
    final Set<String> normals;
    final M2ArgumentCollection args;
    final List<long[]> predicted=new ArrayList<>();
    int actual=0;
    LoggingEngine(M2ArgumentCollection a,Set<String> n,VariantAnnotatorEngine e,SAMFileHeader h,SAMSequenceDictionary d,PrintWriter o,PrintWriter annotationOutput) {
        super(a,n,e,h,d); args=a;normals=n;out=o;annotations=annotationOutput;
    }
    void bits(String name, double[] v) {
        out.print(name); for(double x:v) out.printf(" %016x",Double.doubleToRawLongBits(x)); out.println();
    }
    <E extends Locatable,A extends Allele> void matrix(String mark,LikelihoodMatrix<E,A> m) {
        out.println(mark+" "+m.numberOfAlleles()+" "+m.evidenceCount());
        for(int a=0;a<m.numberOfAlleles();a++) {
            out.print("L"); for(int r=0;r<m.evidenceCount();r++) out.printf(" %016x",Double.doubleToRawLongBits(m.get(a,r))); out.println();
        }
    }
    static String indices(List<Integer> v) { return String.join(",",v.stream().map(Object::toString).toList()); }
    static double[] alts(PerAlleleCollection<Double> x,LikelihoodMatrix<?,Allele> m) {
        return m.alleles().stream().skip(1).mapToDouble(x::getAlt).toArray();
    }
    static String raw(Object value) {
        if(value==null) return "null";
        if(value instanceof Double x) return String.format("d:%016x",Double.doubleToRawLongBits(x));
        if(value instanceof Float x) return String.format("f:%08x",Float.floatToRawIntBits(x));
        if(value.getClass().isArray()) {
            List<String> elements=new ArrayList<>();
            for(int i=0;i<java.lang.reflect.Array.getLength(value);i++) elements.add(raw(java.lang.reflect.Array.get(value,i)));
            return "["+String.join(",",elements)+"]";
        }
        if(value instanceof Collection<?> xs) return "["+String.join(",",xs.stream().map(LoggingEngine::raw).toList())+"]";
        return value.getClass().getSimpleName()+":"+value.toString();
    }
    void annotationValues(VariantContext vc) {
        annotations.println("V "+vc.getContig()+":"+vc.getStart()+" "+String.join(",",vc.getAlleles().stream().map(Allele::getBaseString).toList()));
        for(String key:new TreeSet<>(vc.getAttributes().keySet())) annotations.println("I "+key+" "+raw(vc.getAttribute(key)));
        for(Genotype gt:vc.getGenotypesOrderedByName()) {
            annotations.println("S "+gt.getSampleName()+" "+gt.getGenotypeString());
            if(gt.hasAD()) annotations.println("AD "+raw(gt.getAD()));
            if(gt.hasDP()) annotations.println("DP "+raw(gt.getDP()));
            for(String key:new TreeSet<>(gt.getExtendedAttributes().keySet())) annotations.println("G "+key+" "+raw(gt.getExtendedAttribute(key)));
        }
    }
    @Override public void close() {annotations.close();super.close();}
    @Override protected <E extends Locatable> PerAlleleCollection<Double> somaticLogOdds(LikelihoodMatrix<E,Allele> m) {
        var result=super.somaticLogOdds(m);
        long[] bits=Arrays.stream(alts(result,m)).mapToLong(Double::doubleToRawLongBits).toArray();
        if(actual>=predicted.size() || !Arrays.equals(bits,predicted.get(actual++))) throw new IllegalStateException("actual genotyper differs from independently projected event "+actual);
        return result;
    }
    @Override public CalledHaplotypes callMutations(AlleleLikelihoods<GATKRead,Haplotype> original, AssemblyResultSet rs, ReferenceContext rc, SimpleInterval window, FeatureContext fc, List<Event> given, SAMFileHeader header, boolean bam, boolean gvcf, Set<Integer> suspicious) {
        if(args.genotypeGermlineSites || args.independentMates || gvcf || !given.isEmpty()) throw new IllegalStateException("unsupported oracle mode");
        actual=0;predicted.clear();
        out.println("R "+window);
        Map<String,List<GATKRead>> bySample=new LinkedHashMap<>();
        for(int sn=0;sn<original.numberOfSamples();sn++) bySample.put(original.getSample(sn),original.sampleEvidence(sn));
        var reads=new AlleleLikelihoods<GATKRead,Haplotype>(original,new IndexedAlleleList<>(original.alleles()),bySample);
        for(int sn=0;sn<original.numberOfSamples();sn++) for(int h=0;h<original.numberOfAlleles();h++) for(int r=0;r<original.sampleEvidenceCount(sn);r++) reads.sampleMatrix(sn).set(h,r,original.sampleMatrix(sn).get(h,r));
        reads.normalizeLikelihoods(NaturalLogUtils.qualToLogErrorProb(45),true);
        var fragments=reads.groupEvidence(GATKRead::getName,Fragment::createAndAvoidFailure);
        for(int s=0;s<reads.numberOfSamples();s++) {
            out.println("S "+reads.getSample(s)); matrix("C",reads.sampleMatrix(s));
            var sample=reads.sampleEvidence(s);
            var groups=new ArrayList<>(sample.stream().collect(java.util.stream.Collectors.groupingBy(GATKRead::getName)).values());
            for(int f=0;f<groups.size();f++) {
                var fragment=fragments.sampleEvidence(s).get(f);
                var all=groups.get(f).stream().map(sample::indexOf).toList();
                var chosen=fragment.getReads().stream().map(sample::indexOf).toList();
                out.println("F "+fragment.getStart()+" "+fragment.getEnd()+" "+indices(all)+" "+indices(chosen));
            }
            matrix("G",fragments.sampleMatrix(s));
        }
        EventMap.buildEventMapsForHaplotypes(reads.alleles(),rs.getFullReferenceWithPadding(),rs.getPaddedReferenceLoc(),false,1);
        for(int loc:EventMap.getEventStartPositions(reads.alleles())) {
            if(loc<window.getStart()||loc>window.getEnd()) continue;
            var merged=AssemblyBasedCallerUtils.makeMergedVariantContext(AssemblyBasedCallerUtils.getVariantsFromActiveHaplotypes(loc,reads.alleles(),false));
            if(merged==null) continue;
            var map=AssemblyBasedCallerUtils.createAlleleMapper(merged,loc,reads.alleles(),true);
            var alleles=fragments.marginalize(map);
            var overlap=new SimpleInterval(merged).expandWithinContig(2,header.getSequenceDictionary());
            alleles.retainEvidence(overlap::overlaps);
            out.println("V "+loc+" "+String.join(",",alleles.alleles().stream().map(Allele::getBaseString).toList()));
            LikelihoodMatrix<Fragment,Allele> tumor=null,normal=null;
            for(int s=0;s<alleles.numberOfSamples();s++) {
                final int sn=s;
                out.println("M "+alleles.getSample(s)+" "+indices(alleles.sampleEvidence(s).stream().map(f->fragments.sampleEvidence(sn).indexOf(f)).toList()));
                matrix("A",alleles.sampleMatrix(s));
                if(normals.contains(alleles.getSample(s))) normal=alleles.sampleMatrix(s);else tumor=alleles.sampleMatrix(s);
            }
            double[] tlod=alts(super.somaticLogOdds(tumor),tumor);
            double[] nalod=alts(super.somaticLogOdds(normal),normal);
            predicted.add(Arrays.stream(tlod).mapToLong(Double::doubleToRawLongBits).toArray());
            predicted.add(Arrays.stream(nalod).mapToLong(Double::doubleToRawLongBits).toArray());
            try {
                var m=SomaticGenotypingEngine.class.getDeclaredMethod("diploidAltLogOdds",LikelihoodMatrix.class);m.setAccessible(true);
                @SuppressWarnings("unchecked") var nl=(PerAlleleCollection<Double>)m.invoke(this,normal);
                double[] nlod=alts(nl,normal);
                bits("TLOD",tlod);bits("NLOD",nlod);bits("NALOD",nalod);
                boolean pass=false;
                for(int i=0;i<tlod.length;i++) if(!alleles.getAllele(i+1).equals(Allele.SPAN_DEL)&&tlod[i]>args.getEmissionLogOdds()&&nlod[i]>MathUtils.log10ToLog(args.normalLog10Odds)) pass=true;
                out.println("PASS "+pass);
            } catch(Exception e) {throw new RuntimeException(e);}
        }
        var result=super.callMutations(original,rs,rc,window,fc,given,header,bam,gvcf,suspicious);
        if(actual!=predicted.size()) throw new IllegalStateException("genotyper event count differs");
        annotations.println("R "+window);
        result.getCalls().forEach(this::annotationValues);
        return result;
    }
}
