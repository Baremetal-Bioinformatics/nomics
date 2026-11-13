// Independent oracle for the ln tie window and literal-0.2 informative test.
import java.util.*;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
public class MutectLikelihoodTies {
 public static void main(String[] args) {
  var alt=new Haplotype("AAAAAAAAAAAA".getBytes(),false);alt.setCigar(TextCigarCodec.decode("12M"));
  var ref=new Haplotype("AAAAAAAAAAAC".getBytes(),true);ref.setCigar(TextCigarCodec.decode("12M"));
  double t=0.2*Math.log(10);
  for(double d:new double[]{0,0.199,0.2,0.201,Math.nextDown(t),t,Math.nextUp(t),1}) {
   var lk=new AlleleLikelihoods<org.broadinstitute.hellbender.utils.SimpleInterval,Haplotype>(new IndexedSampleList("s"),new IndexedAlleleList<>(List.of(alt,ref)),Map.of("s",List.of(new org.broadinstitute.hellbender.utils.SimpleInterval("chr1",1,12))));
   lk.switchToNaturalLog();lk.sampleMatrix(0).set(0,0,0);lk.sampleMatrix(0).set(1,0,-d);
   var best=lk.bestAllelesBreakingTies(h->(h.isReference()?1.0:0.0)+1-h.getCigar().numCigarElements()).iterator().next();
   System.out.printf("%016x %d %s %016x%n",Double.doubleToRawLongBits(d),best.allele==alt?0:1,best.isInformative(),Double.doubleToRawLongBits(best.confidence));
  }
 }
}
