// Test-only independent annotation fixtures from the unmodified 4.6.2.0 jar.
import java.nio.file.*;
import java.util.*;
import java.util.stream.*;
import htsjdk.samtools.*;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.tools.walkers.mutect.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.*;
public class MutectAnnotationNumerics {
 static void bits(double x) {System.out.printf(" %016x",Double.doubleToRawLongBits(x));}
 public static void main(String[] args) throws Exception {
  var engine=new SomaticGenotypingEngine(new M2ArgumentCollection(),Set.of(),null,null,null);
  var effective=SomaticGenotypingEngine.class.getDeclaredMethod("getEffectiveCounts",LikelihoodMatrix.class);effective.setAccessible(true);
  var add=SomaticGenotypingEngine.class.getDeclaredMethod("addGenotypes",AlleleLikelihoods.class,List.class,VariantContextBuilder.class);add.setAccessible(true);
  var lines=new ArrayList<>(Files.readAllLines(Path.of(args[0])));
  // Empty evidence is a separate caller branch, not a zero-column RealMatrix.
  for(int a=1;a<=5;a++) lines.add(a+" 0");
  double window=Math.log(10)*0.2;
  double[] deltas={0,0.2,Math.nextUp(0.2),0.3,Math.nextDown(window),window,Math.nextUp(window)};
  String weak="2 "+deltas.length;
  for(double x:deltas)weak+=" "+String.format("%016x",Double.doubleToRawLongBits(-x));
  for(double x:deltas)weak+=" 0000000000000000";
  lines.add(weak);
  for(int a:new int[]{13,20,49}) {
   String wide=a+" 2";
   for(int i=0;i<a;i++)for(int j=0;j<2;j++)wide+=" "+String.format("%016x",Double.doubleToRawLongBits(i==j?0:-10.0-i));
   lines.add(wide);
  }
  for(String line:lines) {
   var f=line.split(" ");int a=Integer.parseInt(f[0]),n=Integer.parseInt(f[1]);
   var alleles=IntStream.range(0,a).mapToObj(i->Allele.create("A"+"C".repeat(i),i==0)).toList();
   List<GATKRead> reads=new ArrayList<>();
   for(int j=0;j<n;j++) {var r=new SAMRecordToGATKReadAdapter(new SAMRecord(null));r.setName("r"+j);r.setPosition("chr1",100);r.setBases("A".getBytes());r.setBaseQualities(new byte[]{30});r.setCigar("1M");reads.add(r);}
   var lk=new AlleleLikelihoods<GATKRead,Allele>(new IndexedSampleList("S"),new IndexedAlleleList<>(alleles),Map.of("S",reads));
   lk.switchToNaturalLog();
   System.out.print("E "+a+" "+n);
   for(int i=0;i<a;i++)for(int j=0;j<n;j++){double x=Double.longBitsToDouble(Long.parseUnsignedLong(f[2+i*n+j],16));lk.sampleMatrix(0).set(i,j,x);bits(x);}
   for(double x:(double[])effective.invoke(null,lk.sampleMatrix(0)))bits(x);
   var vb=new VariantContextBuilder("test","chr1",100,100,alleles);add.invoke(engine,lk,alleles,vb);
   var af=(double[])vb.make().getGenotype("S").getExtendedAttribute("AF");for(double x:af)bits(x);System.out.println();
   var map=alleles.stream().collect(Collectors.toMap(x->x,x->List.of(x)));
   var subset=lk.marginalize(map);
   System.out.print("O "+a);for(var x:subset.alleles())System.out.print(" "+alleles.indexOf(x));System.out.println();
   System.out.print("B "+n);for(var b:subset.bestAllelesBreakingTies()) {System.out.print(" "+alleles.indexOf(b.allele));bits(b.confidence);}System.out.println();
  }
  for(int[] v:new int[][]{{},{0},{0,1},{1,2},{-2,-1},{-1,0},{Integer.MAX_VALUE,Integer.MAX_VALUE},{Integer.MIN_VALUE,Integer.MAX_VALUE},{3,8,1,4,9}}) {
   System.out.print("M "+v.length);for(int x:v)System.out.print(" "+x);System.out.println(" "+(v.length==0?17:MathUtils.median(v)));
  }
  for(String cigar:List.of("10M","3S7M","2H10M3H","3M2D7M","3M2I7M","2I10M","3M2N7M"))for(int pos=98;pos<=114;pos++) {
   var r=new SAMRecordToGATKReadAdapter(new SAMRecord(null));r.setName("r");r.setPosition("chr1",100);r.setCigar(cigar);int n=r.getCigar().getReadLength();r.setBases("A".repeat(n).getBytes());byte[] q=new byte[n];for(int i=0;i<n;i++)q[i]=(byte)(i+15);r.setBaseQualities(q);
   var vc=new VariantContextBuilder("test","chr1",pos,pos,List.of(Allele.create("A",true),Allele.create("C"))).make();
   var b=BaseQuality.getBaseQuality(r,vc);var p=ReadPosition.getPosition(r,vc);
   System.out.println("R "+cigar+" "+pos+" "+(b.isPresent()?b.getAsInt():"null")+" "+(p.isPresent()?p.getAsInt():"null"));
  }
  engine.close();
 }
}
