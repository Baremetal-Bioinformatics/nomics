// Independently generated quality-array and binary64 fixtures from GATK 4.6.2.0.
import htsjdk.samtools.*;
import java.util.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.*;
import org.broadinstitute.hellbender.utils.read.*;
import org.broadinstitute.hellbender.utils.pairhmm.*;
public class MutectLikelihoodNumerics {
 public static void main(String[] args) throws Exception {
  var config=new LikelihoodEngineArgumentCollection(); config.pairHMM=PairHMM.Implementation.LOGLESS_CACHING;
  var engine=(PairHMMLikelihoodCalculationEngine)AssemblyBasedCallerUtils.createLikelihoodCalculationEngine(config,true);
  var method=PairHMMLikelihoodCalculationEngine.class.getDeclaredMethod("modifyReadQualities",List.class); method.setAccessible(true);
  var hex=HexFormat.of();
  String[] bases={"ACGTACGTACGT", "AAAAAAAAAAAA", "ATATATATATAT", "ACGTNNACGTAC", "GATTACAGATTA"};
  for(int k=0;k<bases.length;k++) {
   byte[] b=bases[k].getBytes(),q=new byte[b.length],i=new byte[b.length],d=new byte[b.length];
   for(int n=0;n<b.length;n++) {q[n]=(byte)(new int[]{0,5,6,17,18,19,30,40,60,127,128,254}[n]); i[n]=(byte)(new int[]{45,20,5,6,40,90,93,10,15,35,45,2}[n]); d[n]=(byte)(new int[]{20,45,6,5,10,93,90,40,35,15,2,45}[n]);}
   int mq=new int[]{60,20,30,127,255}[k];
   SAMRecord sam=new SAMRecord(new SAMFileHeader()); sam.setReadName("fixture");sam.setReadBases(b);sam.setBaseQualities(q);sam.setMappingQuality(mq);sam.setCigarString(b.length+"M");
   GATKRead r=new SAMRecordToGATKReadAdapter(sam); ReadUtils.setInsertionBaseQualities(r,i); ReadUtils.setDeletionBaseQualities(r,d);
   @SuppressWarnings("unchecked") var prepared=(List<GATKRead>)method.invoke(engine,List.of(r));
   var p=prepared.get(0);
   byte[] pi=ReadUtils.getBaseInsertionQualities(p),pd=ReadUtils.getBaseDeletionQualities(p),g=new byte[b.length];Arrays.fill(g,(byte)10);
   var hmm=new LoglessPairHMM();hmm.initialize(b.length,32);
   String hap="TT"+bases[k]+"GG";
   var cm=PairHMM.class.getDeclaredMethod("computeReadLikelihoodGivenHaplotypeLog10",byte[].class,byte[].class,byte[].class,byte[].class,byte[].class,byte[].class,boolean.class,byte[].class);cm.setAccessible(true);
   double v=(Double)cm.invoke(hmm,hap.getBytes(),b,p.getBaseQualities(),pi,pd,g,true,null);
   System.out.println(bases[k]+" "+mq+" "+hex.formatHex(q)+" "+hex.formatHex(i)+" "+hex.formatHex(d)+" "+hex.formatHex(p.getBaseQualities())+" "+hex.formatHex(pi)+" "+hex.formatHex(pd)+" "+hap+" "+String.format("%016x",Double.doubleToRawLongBits(v))+" "+String.format("%016x",Double.doubleToRawLongBits(v*Math.log(10))));
  }
  engine.close();
 }
}
