import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.special.Beta;
public class FilteringNumerics {
    static void emit(double x) { System.out.printf("%016x%n", Double.doubleToRawLongBits(x)); }
    public static void main(String[] args) {
        for (double x : new double[]{.01,.1,.5,.7,1,1.2,1.7,2,2.5,3,7.9,8,10,49,1000}) {
            emit(Gamma.logGamma(x)); emit(Gamma.digamma(x));
        }
        for (double[] ab : new double[][]{{.01,.2},{.3,20},{1.79,3.17},{10.02,.5},{99,10001},{1,50000},{.5,2},{2,2}}) {
            emit(Beta.logBeta(ab[0],ab[1]));
            for(double x : new double[]{.001,.1,.5,.9}) emit(Beta.regularizedBeta(x,ab[0],ab[1]));
        }
    }
}
