import org.apache.commons.math3.util.*;
import org.apache.commons.math3.stat.descriptive.moment.*;
import org.apache.commons.math3.stat.descriptive.rank.*;
import java.util.*;
public class MoreNumerics {
 static void out(double x){System.out.printf("%016x%n",Double.doubleToRawLongBits(x));}
 public static void main(String[] args)throws Exception {
  for(double x:new double[]{Double.MIN_VALUE,1e-300,.001,.1,.995,1,1.005,2,100,1e300})out(Math.log(x));
  for(double[] x:new double[][]{{1,2,3,4},{1e16,1,1,-1e16},{.1,.2,.3,.4,.5},{}}){out(new Mean().evaluate(x));out(new Percentile(50).evaluate(x));out(new Percentile(90).evaluate(x));out(Arrays.stream(x).sum());}
  Random r=new Random(1216);for(int n:new int[]{2,100,1000,1073741825})for(int i=0;i<8;i++)System.out.println(r.nextInt(n));
  for(double[] data:new double[][]{{1,0,1,0,1},{3,1,2,1,3,0,4},{0,0,0},{1}})System.out.println(Class.forName("org.broadinstitute.hellbender.tools.copynumber.utils.optimization.PersistenceOptimizer").getMethod("getMinimaIndices").invoke(Class.forName("org.broadinstitute.hellbender.tools.copynumber.utils.optimization.PersistenceOptimizer").getConstructor(double[].class).newInstance((Object)data)));
 }
}
