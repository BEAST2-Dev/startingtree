package beast.math.statistic;

import org.apache.commons.math3.stat.regression.SimpleRegression;

/**
 * @author Walter Xie
 */
public class Regression extends SimpleRegression {

    public Regression(double[] xData, double[] yData) {
        super();

        if (xData.length != yData.length)
            throw new IllegalArgumentException("Array xData and yData should have the same length !");

        for (int i = 0; i < xData.length; i++)
            addData(xData[i], yData[i]);
    }

    public double getYIntercept() {
        return getIntercept();
    }

    public double getXIntercept() {
        return -getIntercept() / getSlope();
    }
}
