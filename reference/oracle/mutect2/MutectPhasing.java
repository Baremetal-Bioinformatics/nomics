// Test-only independent phase-set fixtures from the unmodified jar.
import java.util.*;
import java.lang.reflect.*;
import htsjdk.variant.variantcontext.*;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
public class MutectPhasing {
 public static void main(String[] args) throws Exception {
  var method=AssemblyBasedCallerUtils.class.getDeclaredMethod("constructPhaseSetMapping",List.class,Map.class);method.setAccessible(true);
  List<Haplotype> haps=new ArrayList<>();for(int i=0;i<8;i++)haps.add(new Haplotype(("A"+"C".repeat(i)).getBytes(),false));
  Random random=new Random(74519);
  for(int k=0;k<120;k++) {
   int n=2+k%9;List<VariantContext> calls=new ArrayList<>();Map<VariantContext,Set<Haplotype>> sets=new HashMap<>();List<String> masks=new ArrayList<>();
   for(int i=0;i<n;i++) {int mask=k<8?(i+k)%4:random.nextInt(16);masks.add(""+mask);var vc=new VariantContextBuilder("test","chr1",100+i,100+i,List.of(Allele.create("A",true),Allele.create("C"))).make();calls.add(vc);Set<Haplotype> s=new HashSet<>();for(int h=0;h<8;h++)if((mask&(1<<h))!=0)s.add(haps.get(h));sets.put(vc,s);}
   @SuppressWarnings("unchecked") Map<VariantContext,Pair<Integer,?>> result=(Map<VariantContext,Pair<Integer,?>>)method.invoke(null,calls,sets);
   List<String> expected=new ArrayList<>();for(var vc:calls){var p=result.get(vc);expected.add(p==null?".":p.getLeft()+":"+(p.getRight().toString().equals("PHASE_10")?1:0));}
   System.out.println(String.join(",",masks)+" "+String.join(",",expected));
  }
 }
}
