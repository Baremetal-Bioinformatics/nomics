// Validation only: GATK 4.6.2.0 StrandArtifactFilter expressions, Commons Math 3.5.
// Neither this oracle nor its dependencies are called by nomics.
import org.apache.commons.math3.special.Beta;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.*;
import java.util.*;

public class FilteringLearning {
    static void emit(double x) { System.out.printf("%016x%n", Double.doubleToLongBits(x)); }
    static double bb(int n, int k, double a, double b) {
        return CombinatoricsUtils.binomialCoefficientLog(n, k)
            + Beta.logBeta(k + a, n - k + b) - Beta.logBeta(a, b);
    }
    // nf,nr,kf,kr,indel length; first two include observed platform-pow failures.
    static final int[][] COUNTS = {
        {1,37,0,15,0}, {6,23,1,1,0}, {70,45,28,0,0}, {50,50,0,21,1},
        {90,70,35,2,2}, {100,100,1,40,3}, {20,30,10,0,4}
    };
    static double[] step(int[] c, double a, double b, double prior) {
        int nf=c[0], nr=c[1], kf=c[2], kr=c[3], len=c[4];
        double beta=len==0?1000:len<3?5000:50000;
        double f=bb(nf,kf,a,b)+bb(nr,kr,1,beta);
        double r=bb(nr,kr,a,b)+bb(nf,kf,1,beta);
        double none=CombinatoricsUtils.binomialCoefficientLog(nf,kf)
            + CombinatoricsUtils.binomialCoefficientLog(nr,kr)
            - CombinatoricsUtils.binomialCoefficientLog(nf+nr,kf+kr)+bb(nf+nr,kf+kr,1,1);
        double[] x={(f+Math.log(prior/2))*Math.log10(Math.E),
            (r+Math.log(prior/2))*Math.log10(Math.E),
            (none+Math.log(1-prior))*Math.log10(Math.E)};
        int max=0; for(int i=1;i<3;i++) if(x[i]>x[max]) max=i;
        double sum=1; for(int i=0;i<3;i++) if(i!=max && x[i]!=Double.NEGATIVE_INFINITY)
            sum+=Math.pow(10,x[i]-x[max]);
        double ls=x[max]+(sum!=1?Math.log10(sum):0);
        return new double[]{Math.pow(10,x[0]-ls),Math.pow(10,x[1]-ls),nf,nr,kf,kr};
    }
    public static void main(String[] args) {
        emit(-Math.log(3)); emit(Math.log(1.0/3));
        for(double x:new double[]{-20,-.7,-.1,-1e-8,0,.01,.9}) {
            emit(Math.log1p(x)); emit(Math.expm1(x));
        }
        for(double x:new double[]{-300,-17.61,-5.186226408608968,-.1,0,2,3,8})
            emit(Math.pow(10,x));
        for(double x:new double[]{1e-300,.001,.1,1,1.5,3,1000}) emit(Math.log10(x));
        double[] probs={0.04163411606541967,0.006318578253251994,0.03865375732824039};
        for(int[] order:new int[][]{{0,1,2},{0,2,1}}) {
            double truth=1;for(int i:order)truth*=1-probs[i];double combined=1-truth;
            emit(combined);emit(Math.round((combined+Math.ulp(combined))*1000)/1000.0);
        }
        double a=1,b=20,prior=.001;
        for(int pass=0;pass<3;pass++) {
            List<double[]> steps=new ArrayList<>();
            for(int[] c:COUNTS) { double[] s=step(c,a,b,prior);steps.add(s);emit(s[0]);emit(s[1]); }
            if(pass==2)break;
            List<double[]> potential=new ArrayList<>();
            for(double[] s:steps) if(s[0]+s[1]>.1) potential.add(s);
            double artifacts=potential.stream().mapToDouble(s->s[0]+s[1]).sum();
            double non=steps.stream().mapToDouble(s->1-(s[0]+s[1])).sum();
            prior=(artifacts+1)/(artifacts+1+non+1000);
            double alt=potential.stream().mapToDouble(s->s[0]*s[4]+s[1]*s[5]).sum();
            double depth=potential.stream().mapToDouble(s->s[0]*s[2]+s[1]*s[3]).sum();
            double mean=(alt+1)/(depth+1+20);
            UnivariateObjectiveFunction objective=new UnivariateObjectiveFunction(alpha -> {
                double beta=(1/mean-1)*alpha;
                return potential.stream().mapToDouble(s->s[0]*bb((int)s[2],(int)s[4],alpha,beta)
                    +s[1]*bb((int)s[3],(int)s[5],alpha,beta)).sum();
            });
            a=new BrentOptimizer(.01,.01).optimize(new MaxEval(100),objective,
                GoalType.MAXIMIZE,new SearchInterval(.01,100,1)).getPoint();
            b=(1/mean-1)*a;
            emit(a);emit(b);emit(prior);
        }
    }
}
