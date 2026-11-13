import org.apache.commons.math3.special.Beta;
import org.apache.commons.math3.util.*;
public class Oracle {
 static double bb(int n,int k,double a,double b) { return CombinatoricsUtils.binomialCoefficientLog(n,k)+Beta.logBeta(k+a,n-k+b)-Beta.logBeta(a,b); }
 static double[] resp(int ref,int alt,int n,int k,int x,double[] p) {
  double[] a={1,1,1,1,1,1,1,1,3,5,2,10000}, b={9,9,9,9,9,9,9,9,10000,5,5,3};
  double[] r=new double[12]; int max=0;
  for(int s=0;s<12;s++) { r[s]=s<8&&(s%4==ref||s%4!=alt)?Double.NEGATIVE_INFINITY:Math.log(p[s])+bb(n,k,a[s],b[s])+bb(k,x,s<4?100:s<8?1:10,s<4?1:s<8?100:10); if(r[s]>r[max])max=s; }
  double sum=1; for(int s=0;s<12;s++)if(s!=max&&r[s]!=Double.NEGATIVE_INFINITY)sum+=Math.exp(r[s]-r[max]);
  double logsum=r[max]+(sum!=1?Math.log(sum):0); for(int s=0;s<12;s++)r[s]=Math.exp(r[s]-logsum);return r;
 }
 public static void main(String[] args) {
  for(double[] ab:new double[][]{{1,9},{3,10000},{5,5},{2,5},{10000,3},{100,1},{1,100},{10,10}})for(int[] nk:new int[][]{{0,0},{1,0},{20,2},{200,1},{66,33},{1029,513},{1200,17}})System.out.printf("bb %d %d %.0f %.0f %016x%n",nk[0],nk[1],ab[0],ab[1],Double.doubleToLongBits(bb(nk[0],nk[1],ab[0],ab[1])));
  for(double x:new double[]{Double.MIN_VALUE,0.1,0.99,0.99999,1,1.00001,1.009,2,Math.PI,10000,Double.MAX_VALUE}) System.out.printf("log %016x %016x%n",Double.doubleToLongBits(x),Double.doubleToLongBits(FastMath.log(x)));
  for(double x:new double[]{-0.99999999,-1e-7,0,1e-7,1e-6,1.00001e-6,0.01,10,1e15}) System.out.printf("log1p %016x %016x%n",Double.doubleToLongBits(x),Double.doubleToLongBits(FastMath.log1p(x)));
  for(double a:new double[]{1,2})for(double b:new double[]{1,2})System.out.printf("beta %.0f %.0f %016x%n",a,b,Double.doubleToLongBits(Beta.logBeta(a,b)));
  double[] p=new double[12];java.util.Arrays.fill(p,.1);p[1]=p[5]=0;
  for(int[] nkx:new int[][]{{30,0,0},{30,2,2},{30,2,0},{100,50,23}}) { double[] r=resp(1,nkx[1]==0?1:0,nkx[0],nkx[1],nkx[2],p);System.out.print("resp "+java.util.Arrays.toString(nkx));for(double v:r)System.out.printf(" %016x",Double.doubleToLongBits(v));System.out.println(); }
 }
}
