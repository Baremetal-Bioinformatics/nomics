import org.apache.commons.math3.util.*;
import org.apache.commons.math3.distribution.*;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.optim.*;
import org.apache.commons.math3.optim.nonlinear.scalar.*;
import org.apache.commons.math3.optim.univariate.*;
public class Numerics{
 static void out(double x){System.out.printf("%016x%n",Double.doubleToRawLongBits(x));}
 public static void main(String[] args)throws Exception{
 for(double x:new double[]{-745,-710,-709,-100,-4,-.1,0,.01,1,100,709})out(FastMath.exp(x));
 for(double x:new double[]{Double.MIN_VALUE,1e-300,.001,.1,.995,1,1.005,2,100,1e300})out(FastMath.log(x));
 for(double[] xy:new double[][]{{3,4},{1e-300,1e-300},{1e300,1e300},{0,0},{1,.123}})out(FastMath.hypot(xy[0],xy[1]));
 for(int[] nk:new int[][]{{0,0},{10,0},{10,10},{10,3},{150,60},{300,4}})for(double p:new double[]{0,.01,.1,.4,1})out(new BinomialDistribution(null,nk[0],p).probability(nk[1]));
 for(double guess:new double[]{.02,.05,.1,.2,.4}){var r=new BrentOptimizer(.0001,.0001).optimize(new UnivariateObjectiveFunction(x->-(x-.137)*(x-.137)),GoalType.MAXIMIZE,new SearchInterval(0,.5,guess),new MaxEval(30));out(r.getPoint());out(r.getValue());}
 for(double[][] m:new double[][][]{{{1,2},{2,4}},{{1,.2,.3},{.2,1,.7},{.3,.7,1}}}){var s=new SingularValueDecomposition(new Array2DRowRealMatrix(m));for(double x:s.getSingularValues())out(x);for(double[] row:s.getU().getData())for(double x:row)out(x);}
 var cls=Class.forName("org.apache.commons.math3.distribution.SaddlePointExpansion");
 var stirling=cls.getDeclaredMethod("getStirlingError",double.class);stirling.setAccessible(true);
 for(double z:new double[]{0,1,14,15,20,100})out((double)stirling.invoke(null,z));
 var deviance=cls.getDeclaredMethod("getDeviancePart",double.class,double.class);deviance.setAccessible(true);
 for(double[] xy:new double[][]{{10,9.9},{10,4},{1,1},{100,110}})out((double)deviance.invoke(null,xy[0],xy[1]));

 }
}
