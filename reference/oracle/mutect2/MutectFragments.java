// Independent groupEvidence / Fragment.createAndAvoidFailure fixtures.
import java.util.*;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.*;
import org.broadinstitute.hellbender.utils.NaturalLogUtils;
public class MutectFragments {
    public static void main(String[] args) {
        List<GATKRead> rs=new ArrayList<>(); List<double[]> values=new ArrayList<>();
        // All names here collide under String.hashCode, exercising tree bins.
        for(int i=0;i<32;i++) {
            String name="";
            for(int b=0;b<5;b++) name+=(i&(1<<b))==0?"Aa":"BB";
            add(rs,values,name,0,100+i,20,-0.1*i,-70+i*0.25);
        }
        for(int i=0;i<100;i++) add(rs,values,"read"+i,0,200+i,10+i%7,-0.25*i,-0.4*i);
        // <=2 reads are not filtered for flags. >2 has every fallback case.
        add(rs,values,"paired",0x800,700,25,-0.1,-100);
        add(rs,values,"paired",0x100,800,10,-100,-0.2);
        add(rs,values,"many",0,600,10,-0.25,-40);
        add(rs,values,"many",0x800,550,50,-0.5,-50);
        add(rs,values,"many",0,800,10,-10,-1);
        add(rs,values,"many",0,900,10,-1,-10);
        for(int i=0;i<4;i++) add(rs,values,"allflagged",0x800,1000+i*10,10,-0.5,-50);
        add(rs,values,"oneclean",0x800,1300,10,-0.1,-40);
        add(rs,values,"oneclean",0,1350,10,-10,-1);
        add(rs,values,"oneclean",0x100,1400,10,-1,-10);
        var lk=new AlleleLikelihoods<GATKRead,Allele>(new IndexedSampleList("S"),new IndexedAlleleList<>(Allele.REF_A,Allele.ALT_C),Map.of("S",rs));
        for(int r=0;r<rs.size();r++) for(int a=0;a<2;a++) lk.sampleMatrix(0).set(a,r,values.get(r)[a]);
        lk.normalizeLikelihoods(NaturalLogUtils.qualToLogErrorProb(45),true);
        matrix("C",lk.sampleMatrix(0));
        var grouped=lk.groupEvidence(GATKRead::getName,Fragment::createAndAvoidFailure);
        var groups=new ArrayList<>(rs.stream().collect(java.util.stream.Collectors.groupingBy(GATKRead::getName)).values());
        for(int i=0;i<groups.size();i++) {
            var fragment=grouped.sampleEvidence(0).get(i);
            System.out.println("F "+fragment.getStart()+" "+fragment.getEnd()+" "+join(groups.get(i).stream().map(rs::indexOf).toList())+" "+join(fragment.getReads().stream().map(rs::indexOf).toList()));
        }
        matrix("G",grouped.sampleMatrix(0));
    }
    static String join(List<Integer> xs) {return String.join(",",xs.stream().map(Object::toString).toList());}
    static void matrix(String mark,LikelihoodMatrix<?,Allele> m) {
        System.out.println(mark+" "+m.numberOfAlleles()+" "+m.evidenceCount());
        for(int a=0;a<m.numberOfAlleles();a++) {
            System.out.print("L");for(int r=0;r<m.evidenceCount();r++) System.out.printf(" %016x",Double.doubleToRawLongBits(m.get(a,r)));System.out.println();
        }
    }
    static void add(List<GATKRead> rs,List<double[]> values,String name,int flags,int start,int len,double r,double a) {
        var record=new SAMRecord(null);record.setFlags(flags);
        var read=new SAMRecordToGATKReadAdapter(record);read.setName(name);read.setPosition("chr1",start);read.setCigar(len+"M");read.setBases("A".repeat(len).getBytes());read.setBaseQualities(new byte[len]);
        rs.add(read);values.add(new double[]{r,a});
        System.out.printf("R %s %d %d %d %016x %016x%n",name,flags,start,len,Double.doubleToRawLongBits(r),Double.doubleToRawLongBits(a));
    }
}
