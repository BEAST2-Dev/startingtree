
package beast.math;


/**
 * Interface for a function of one variable.
 *
 * @author Korbinian Strimmer
 */
public interface UnivariateFunction {
    /**
     * compute function value
     *
     * @param argument
     * @return function value
     */
    double evaluate(double argument);

    /**
     * @return lower bound of argument
     */
    double getLowerBound();

    /**
     * @return upper bound of argument
     */
    double getUpperBound();

    public abstract class AbstractLogEvaluatableUnivariateFunction implements UnivariateFunction {

        public double logEvaluate(double argument) {
            return Math.log(evaluate(argument));
        }

    }

}
