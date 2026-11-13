// Test fixture generator; invokes the unmodified GATK jar.
import java.util.*;
import htsjdk.samtools.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.activityprofile.*;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
public class MutectNumerics {
  public static void main(String[] args) {
    for(int r : new int[]{0,1,10,50,100,1000,1500}) {
      for(String raw : new String[]{"", "25", "30,30", "7,20,25,30,40", "127,127", "254,10"}) {
        List<Byte> q = new ArrayList<>();
        if(!raw.isEmpty())for(String s:raw.split(","))q.add((byte)Integer.parseInt(s));
        System.out.printf("LLR\t%d\t%s\t%016x%n",r,raw,Double.doubleToRawLongBits(Mutect2Engine.logLikelihoodRatio(r,q)));
      }
    }
    org.broadinstitute.hellbender.utils.Utils.resetRandomGenerator();
    org.broadinstitute.hellbender.utils.downsampling.MutectDownsampler ds = new org.broadinstitute.hellbender.utils.downsampling.MutectDownsampler(50,0,1);
    SAMFileHeader header = new SAMFileHeader();
    header.addSequence(new SAMSequenceRecord("test",1000));
    for(int id=0;id<143;id++) {
      SAMRecord sam=new SAMRecord(header); sam.setReadName(Integer.toString(id)); sam.setReferenceName("test");
      sam.setAlignmentStart(id<80?10:id<140?20:30); sam.setCigarString("1M"); sam.setReadString("A"); sam.setBaseQualityString("I");
      sam.setMappingQuality(id<80?(id%7==0?50:60):id<140?(id<90?60:50):20);
      ds.submit(new org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter(sam));
    }
    ds.signalEndOfInput();
    System.out.print("DS"); for(var read:ds.consumeFinalizedItems())System.out.print("\t"+read.getName());System.out.println();
    BandPassActivityProfile profile = new BandPassActivityProfile(50,0.002,50,17.0,header);
    for(int pos=1;pos<=1000;pos++) {
      for(AssemblyRegion region:profile.popReadyAssemblyRegions(100,50,300,false))dump(region);
      profile.add(new ActivityProfileState(new SimpleInterval("test",pos,pos),(pos==400||pos==600)?1.0:0.0));
    }
    for(AssemblyRegion region:profile.popReadyAssemblyRegions(100,50,300,true))dump(region);
  }
  static void dump(AssemblyRegion r) { System.out.printf("REGION\t%d\t%d\t%s%n",r.getStart(),r.getEnd(),r.isActive()); }
}
