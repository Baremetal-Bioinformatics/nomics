// Independent unit fixtures from the unmodified GATK pruner, including exact LLR bits.
import java.util.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.*;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
public class AdaptivePruner {
  public static void main(String[] args) throws Exception {
    var find = ChainPruner.class.getDeclaredMethod("findAllChains", BaseGraph.class); find.setAccessible(true);
    var remove = AdaptiveChainPruner.class.getDeclaredMethod("chainsToRemove", List.class); remove.setAccessible(true);
    for (int n : new int[]{1,2,3,10,100,1000}) for (int alt : new int[]{1,2,5,20}) {
      var g = new SeqGraph(1);
      var vertices = new ArrayList<SeqVertex>();
      for (String s : new String[]{"A","C","G","T","AC","CC","GC"}) { var v = new SeqVertex(s); g.addVertex(v); vertices.add(v); }
      int[][] edges = {{0,1,n,1},{1,2,n,1},{2,5,n,1},{1,3,alt,0},{3,5,alt,0},{1,4,1,0},{4,5,1,0},{5,6,n+alt,1}};
      for (int[] e : edges) g.addEdge(vertices.get(e[0]), vertices.get(e[1]), new BaseEdge(e[3]==1,e[2]));
      dump(g, vertices, find, remove, ""+n+"\t"+alt);
    }
    // More than 100 seeded variants: terminal counts are small relative to
    // internal weights, so the second-pass error estimate still retains seeds.
    var cap = new SeqGraph(1); var cv = new ArrayList<SeqVertex>();
    for (String b : new String[]{"A","C","G","T"}) { var v = new SeqVertex(b); cap.addVertex(v); cv.add(v); }
    cap.addEdge(cv.get(0),cv.get(1),new BaseEdge(true,1000000));
    cap.addEdge(cv.get(2),cv.get(3),new BaseEdge(true,1000000));
    for (int i=0;i<110;i++) {
      var x = new SeqVertex(String.format("A%03d",i)); var y = new SeqVertex(String.format("C%03d",i));
      cap.addVertex(x); cap.addVertex(y); cv.add(x); cv.add(y);
      cap.addEdge(cv.get(1),x,new BaseEdge(i==0,1000));
      cap.addEdge(x,y,new BaseEdge(i==0,100000));
      cap.addEdge(y,cv.get(2),new BaseEdge(i==0,1000));
    }
    dump(cap,cv,find,remove,"cap");

    for (int ref : new int[]{0,1,10,1000}) for(int alt : new int[]{0,1,2,100}) for (double e : new double[]{0,0.001,0.01,0.1,0.5,1})
      System.out.println("LLR\t"+ref+"\t"+alt+"\t"+e+"\t"+Long.toHexString(Double.doubleToLongBits(Mutect2Engine.logLikelihoodRatio(ref,alt,e))));
  }
  static void dump(SeqGraph g, ArrayList<SeqVertex> vertices, java.lang.reflect.Method find, java.lang.reflect.Method remove, String label) throws Exception {
      var pruner = new AdaptiveChainPruner<SeqVertex,BaseEdge>(0.001,Math.log(10),Math.log(10)*4,100);
      @SuppressWarnings("unchecked") var chains = (List<Path<SeqVertex,BaseEdge>>)find.invoke(pruner,g);
      @SuppressWarnings("unchecked") var removed = (Collection<Path<SeqVertex,BaseEdge>>)remove.invoke(pruner,chains);
      System.out.println("CASE\t"+label);
      for (var c : chains) {
        int l = g.outgoingEdgesOf(c.getFirstVertex()).stream().mapToInt(BaseEdge::getMultiplicity).sum();
        int r = g.incomingEdgesOf(c.getLastVertex()).stream().mapToInt(BaseEdge::getMultiplicity).sum();
        String mult = c.getEdges().stream().map(e->""+e.getMultiplicity()).collect(java.util.stream.Collectors.joining(","));
        System.out.println("CHAIN\t"+vertices.indexOf(c.getFirstVertex())+"\t"+vertices.indexOf(c.getLastVertex())+"\t"+new String(c.getFirstVertex().getSequence())+"\t"+new String(c.getBases())+"\t"+mult+"\t"+c.getEdges().get(0).isRef()+"\t"+c.getEdges().stream().anyMatch(BaseEdge::isRef)+"\t"+l+"\t"+r+"\t"+g.isSource(c.getFirstVertex())+"\t"+g.isSink(c.getLastVertex())+"\t"+removed.contains(c));
      }
  }
}
