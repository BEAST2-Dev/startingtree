package beast.evolution.tree;


import beast.core.Citation;
import beast.core.Description;
import beast.evolution.alignment.TaxonSet;
import beast.math.UnivariateFunction;
import beast.math.UnivariateMinimum;
import beast.math.statistic.DiscreteStatistics;
import beast.math.statistic.Regression;

import java.util.*;

/**
 *  Usage:
 *  <code>
 *  Tree rootedTree = tree;
 *  if (!keepRoot)
 *    rootedTree = temporalRooting.findRoot(tree, TemporalRooting.RootingFunction.CORRELATION);
 *  regressions.add(temporalRooting.getRootToTipRegression(rootedTree));
 *  </code>
 *
 *  @author Andrew Rambaut
 */
@Citation("Rambaut et al 2016, 'Exploring the temporal structure of heterochronous sequences using TempEst (formerly Path-O-Gen)', Virus Evolution")
@Description("Regression method to analyse temporally sampled sequence data. Imported from BEAST 1 TempEst.")
public class TemporalRooting {

    public enum RootingFunction {
        HEURISTIC_RESIDUAL_MEAN_SQUARED("heuristic residual mean squared"),
        RESIDUAL_MEAN_SQUARED("residual mean squared"),
        //        SUM_RESIDUAL_SQUARED("sum squared residuals"),
        CORRELATION("correlation"),
        R_SQUARED("R squared");

        RootingFunction(final String name) {
            this.name = name;
        }

        @Override
        public String toString() {
            return name;
        }

        private final String name;
    }

    private boolean contemporaneous = false;
    private final TaxonSet taxa;
    private final Map<String, Double> dates;
//    private final Map<String, Double> precisions;
    private boolean useTargetRate = false;
    private double targetRate = 0.0;
    private double dateMin;
    private double dateMax;

    private int currentRootBranch = 0;
    private int totalRootBranches = 0;

    private boolean forcePositiveRate = false;

    // todo only suit for forward time represented by double at the moment
    public TemporalRooting(TraitSet timeTraitSet) {
        this.taxa = timeTraitSet.taxaInput.get();

        dates = new HashMap<String, Double>();
//        precisions = new HashMap<String, Double>();

        dateMin = timeTraitSet.minValue;
        dateMax = timeTraitSet.maxValue;

        for (int i = 0; i < taxa.getTaxonCount(); i++) {
            String taxon = taxa.getTaxonId(i);
            // note: timeTraitSet.getValue(taxon) is not the original value
            double date = Double.parseDouble(timeTraitSet.getStringValue(taxon));
            dates.put(taxon, date);
        }

        if (Math.abs(dateMax - dateMin) < 1.0E-8) {
            // probably contemporaneous tips
            contemporaneous = true;
        }
    }

    public void setForcePositiveRate(boolean forcePositiveRate) {
        this.forcePositiveRate = forcePositiveRate;
    }

    public void setTargetRate(double targetRate) {
        this.targetRate = targetRate;
    }

    public boolean isContemporaneous() {
        return contemporaneous;
    }

    public double getDateRange() {
        return dateMax - dateMin;
    }

    /**
     * Find the optimised root looping through each branch and resetting a new root.
     *
     * @param tree binary tree
     * @param rootingFunction
     * @return
     */
    public Tree findRoot(Tree tree, RootingFunction rootingFunction) {

        double[] dates = getTipDates(tree);
        return findGlobalRoot(tree, dates, rootingFunction, forcePositiveRate);
    }

    /**
     * Find the optimised root between its two children without changing all other child nodes.
     *
     * @param tree binary tree
     * @param rootingFunction
     * @return
     */
    public Tree findLocalRoot(Tree tree, RootingFunction rootingFunction) {

        double[] dates = getTipDates(tree);
        FlexibleTree bestTree = new FlexibleTree(tree.getRoot());

        double score = findLocalRoot(bestTree, dates, rootingFunction, forcePositiveRate);
        System.out.println("score = " + score);
        return bestTree;
    }

    /**
     * Root to tip distance vs. time of sampling.
     * Time is given in <code>TraitSet</code>
     *
     * @param tree
     * @return <code>Regression</code> result.
     */
    public Regression getRootToTipRegression(Tree tree) {

        if (contemporaneous) {
            throw new IllegalArgumentException("Cannot do a root to tip regression on contemporaneous tips");
        }

        double[] dates = getTipDates(tree);
        double[] distances = getRootToTipDistances(tree);

        return new Regression(dates, distances);
    }

    public Regression getNodeDensityRegression(Tree tree) {

        if (contemporaneous) {
            throw new IllegalArgumentException("Cannot do a node density regression on contemporaneous tips");
        }

        double[] dates = getTipDates(tree);
//        double[] distances = getRootToTipDistances(tree);
        double[] density = getNodeDensity(tree);

        return new Regression(dates, density);
    }

    public Regression getAncestorRootToTipRegression(Tree tree, Regression regression) {

        if (contemporaneous) {
            throw new IllegalArgumentException("Cannot do a root to tip regression on contemporaneous tips");
        }

        List<Node> externalNodes = tree.getExternalNodes();
        double[] dates = new double[externalNodes.size()];
        double[] distances = new double[externalNodes.size()];
        for (int i = 0; i < externalNodes.size(); i++) {
            Node tip = externalNodes.get(i);
            Node parent = tip.getParent();
            distances[i] = getRootToTipDistance(tree, parent);

            dates[i] = regression.getXIntercept() + (distances[i] / regression.getGradient());
        }

        return new Regression(dates, distances);
    }

    public double[] getRootToTipDistances(Tree tree) {

        List<Node> externalNodes = tree.getExternalNodes();
        double[] d = new double[externalNodes.size()];
        for (int i = 0; i < externalNodes.size(); i++) {
            Node tip = externalNodes.get(i);
            d[i] = getRootToTipDistance(tree, tip);
        }
        return d;
    }

    public double[] getParentRootToTipDistances(Tree tree) {

        List<Node> externalNodes = tree.getExternalNodes();
        double[] d = new double[externalNodes.size()];
        for (int i = 0; i < externalNodes.size(); i++) {
            Node tip = externalNodes.get(i);
            Node parent = tip.getParent();
            d[i] = getRootToTipDistance(tree, parent);
        }
        return d;
    }

    public double[] getRootToTipResiduals(Tree tree, Regression regression) {

        List<Node> externalNodes = tree.getExternalNodes();
        double[] r = new double[externalNodes.size()];
        for (int i = 0; i < externalNodes.size(); i++) {
            Node tip = externalNodes.get(i);
            Double date = dates.get(tree.getTaxonId(tip));
            double d = getRootToTipDistance(tree, tip);

            r[i] = regression.getResidual(date, d);
        }
        return r;
    }

    public double[] getNodeDensity(Tree tree) {

        List<Node> externalNodes = tree.getExternalNodes();
        double[] d = new double[externalNodes.size()];
        for (int i = 0; i < externalNodes.size(); i++) {
            Node tip = externalNodes.get(i);
            d[i] = getNodeDensity(tree, tip);
        }
        return d;
    }

    public double[] getTipDates(Tree tree) {
        List<Node> externalNodes = tree.getExternalNodes();
        double[] d = new double[externalNodes.size()];
        for (int i = 0; i < externalNodes.size(); i++) {
            Node tip = externalNodes.get(i);
            Double date = dates.get(tree.getTaxonId(tip));
            if (date == null) {
                throw new IllegalArgumentException("Taxon, " + tree.getTaxonId(tip) + ", not found in taxon list");
            }
            d[i] = date;
        }
        return d;
    }

//    public double[] getTipDatePrecisions(Tree tree) {
//        List<Node> externalNodes = tree.getExternalNodes();
//        double[] p = new double[externalNodes.size()];
//        for (int i = 0; i < externalNodes.size(); i++) {
//            Node tip = externalNodes.get(i);
//            Double precision = precisions.get(tree.getNodeTaxon(tip).getId());
//            if (precision == null) {
//                precision = 0.0;
//            }
//            p[i] = precision;
//        }
//        return p;
//    }

    public String[] getTipLabels(Tree tree) {
        List<Node> externalNodes = tree.getExternalNodes();
        String[] labels = new String[externalNodes.size()];
        for (int i = 0; i < externalNodes.size(); i++) {
            Node tip = externalNodes.get(i);
            labels[i] = tip.getID();
        }
        return labels;
    }

    private Tree findGlobalRoot(final Tree source, final double[] dates, RootingFunction rootingFunction, boolean forcePositiveRate) {

        FlexibleTree bestTree = new FlexibleTree(source.getRoot());
        double minF = findLocalRoot(bestTree, dates, rootingFunction, forcePositiveRate);
        double minDiff = Double.MAX_VALUE;

        totalRootBranches = source.getNodeCount();
        for (currentRootBranch = 0; currentRootBranch < source.getNodeCount(); currentRootBranch++) {
            FlexibleTree tmpTree = new FlexibleTree(source.getRoot());
            Node node = tmpTree.getNode(currentRootBranch);
            if (!tmpTree.isRoot(node)) {
                double length = node.getLength();
                tmpTree.changeRootTo(node, 0.5);

                double f = findLocalRoot(tmpTree, dates, rootingFunction, forcePositiveRate);
                if (useTargetRate) {
                    Regression r = getRootToTipRegression(tmpTree);
                    if (Math.abs(r.getGradient() - targetRate) < minDiff) {
                        minDiff = Math.abs(r.getGradient() - targetRate);
                        bestTree = tmpTree;
                    }
                } else {
                    if (f < minF) {
                        minF = f;
                        bestTree = tmpTree;
                    }
                }
            }
        }

        return bestTree;
    }

    private double findLocalRoot(final FlexibleTree tree,
                                 final double[] dates,
                                 final RootingFunction rootingFunction,
                                 final boolean forcePositiveRate) {

        if (rootingFunction == RootingFunction.RESIDUAL_MEAN_SQUARED) {
            return findAnalyticalLocalRoot(tree, dates, rootingFunction);
        }

        Node node1 = tree.getRoot().getChild(0);
        Node node2 = tree.getRoot().getChild(1);

        final double length1 = node1.getLength();
        final double length2 = node2.getLength();

        final double sumLength = length1 + length2;

        final List<Node> tipSet1 = node1.getAllLeafNodes();
        final List<Node> tipSet2 = node2.getAllLeafNodes();

        final double[] y = new double[tree.getLeafNodeCount()];

        UnivariateFunction f = new UnivariateFunction() {
            //        MultivariateFunction f = new MultivariateFunction() {
            public double evaluate(final double argument) {
                double l1 = argument * sumLength;

                for (Node tip : tipSet1) {
                    y[tip.getNr()] = getRootToTipDistance(tree, tip) - length1 + l1;
                }

                double l2 = (1.0 - argument) * sumLength;

                for (Node tip : tipSet2) {
                    y[tip.getNr()] = getRootToTipDistance(tree, tip) - length2 + l2;
                }

                double score;

                if (!contemporaneous) {
                    Regression r = new Regression(dates, y);

                    switch (rootingFunction) {

                        case CORRELATION:
                            score = -r.getCorrelationCoefficient();
                            break;
                        case R_SQUARED:
                            score = -r.getRSquared();
                            break;
                        case HEURISTIC_RESIDUAL_MEAN_SQUARED:
                        case RESIDUAL_MEAN_SQUARED:
                            score = r.getResidualMeanSquared();
                            break;
                        default:
                            throw new RuntimeException("Unknown enum value");
                    }

                    if (forcePositiveRate) {
                        score = (r.getGradient() < 0.0 ? -score : score);
                    }

                } else {
                    score = DiscreteStatistics.variance(y);
                }

                return score;
            }

            public int getNumArguments() {
                return 1;
            }

            public double getLowerBound() {
                return 0;
            }

            public double getUpperBound() {
                return 1.0;
            }
        };

//        DifferentialEvolution minimum = new DifferentialEvolution(1);
//        ConjugateDirectionSearch minimum = new ConjugateDirectionSearch();
//        double[] minx = new double[] { 0.5 };
//
//        double fminx = minimum.findMinimum(f, minx);
//        double x = minx[0];

        UnivariateMinimum minimum = new UnivariateMinimum();
        double x = minimum.findMinimum(f);
// todo x == 0 ?
        double fminx = minimum.fminx;
        double l1 = x * sumLength;
        double l2 = (1.0 - x) * sumLength;

        tree.setBranchLength(node1, l1);
        tree.setBranchLength(node2, l2);

        return fminx;
    }

    private double findAnalyticalLocalRoot(final FlexibleTree tree,
                                           final double[] t,
                                           final RootingFunction rootingFunction) {

        if (rootingFunction != RootingFunction.RESIDUAL_MEAN_SQUARED) {
            throw new UnsupportedOperationException("Analytical local root solution only for residual mean squared");
        }

        Node node1 = tree.getRoot().getChild(0);
        Node node2 = tree.getRoot().getChild(1);

        final double length1 = node1.getLength();
        final double length2 = node2.getLength();

        final double sumLength = length1 + length2;

        final List<Node> tipSet1 = node1.getAllLeafNodes();
        final List<Node> tipSet2 = node2.getAllLeafNodes();

        int N = tipSet1.size() + tipSet2.size();
        int n = tipSet2.size();
//todo if n==0, then N < y.length
        final double[] c = new double[N];
        for (Node tip : tipSet2) {
            int i = tip.getNr();
            c[i] = 1;
        }

        final double[] y = getRootToTipDistances(tree);
        for (int j = 0; j < y.length; j++) { // little fiddling with the root-to-tip divergences to get the right input vector
        	y[j] = y[j] + (1-c[j])*(sumLength-length1) - c[j]*(sumLength-length1);
        }

        double sum_tt = 0.0;
        double sum_t = 0.0;
        double sum_y = 0.0;
        double sum_ty = 0.0;
        double sum_tc = 0.0;
        double Nd = N;
        double nd = n;  // need to set these naughty guys to doubles

        for (int i = 0; i < N; i++) {
            sum_tt += t[i] * t[i];
            sum_t += t[i];
            sum_y += y[i];
            sum_ty += t[i] * y[i];
            sum_tc += t[i] * c[i];
        }
        double y_bar = sum_y / Nd;
        double t_bar = sum_t / Nd;

        double C = sum_tt - (sum_t * sum_t / Nd);
        double sumAB = 0.0;
        double sumAA = 0.0;

        for (int i = 0; i < N; i++) {
            double Ai = 2*c[i] -
            		    ((2*nd-Nd)/Nd) +
            		(2*(t_bar-t[i])/(C*Nd))*(Nd*sum_tc - nd*sum_t) - 1;
            double Bi = (y[i] - y_bar)
                    + ((t_bar - t[i]) / (C * Nd)) * ((Nd * sum_ty) - (sum_t * sum_y));

            sumAB += Ai * Bi;
            sumAA += Ai * Ai;
        }
        double x = -sumAB / (sumLength * sumAA);
        x = Math.min(Math.max(x, 0.0), 1.0);

        double l1 = (1.0 - x) * sumLength;
        double l2 = x * sumLength;


        tree.setBranchLength(node1, l1);
        tree.setBranchLength(node2, l2);

        Regression r = new Regression(t, getRootToTipDistances(tree));

        return r.getResidualMeanSquared();
    }

    public double getRootToTipDistance(Tree tree, Node node) {
        double distance = 0;
        while (node != null) {
            distance += node.getLength();
            node = node.getParent();
        }
        return distance;
    }

    public double getNodeDensity(Tree tree, Node node) {
        double density = 0;
        while (node != null) {
            density ++;
            node = node.getParent();
        }
        return density;
    }

    public Tree adjustTreeToConstraints(Tree source, Map<Set<String>, double[]> cladeHeights) {

        FlexibleTree tree = new FlexibleTree(source.getRoot());
        setHeightsFromDates(tree);

        adjustTreeToConstraints(tree, tree.getRoot(), null, cladeHeights);

        return tree;
    }

    public int getCurrentRootBranch() {
        return currentRootBranch;
    }

    public int getTotalRootBranches() {
        return totalRootBranches;
    }

    private double adjustTreeToConstraints(FlexibleTree tree, Node node,
                                           Set<String> leaves,
                                           Map<Set<String>, double[]> cladeHeights) {

        if (!node.isLeaf()) {
            Set<String> l = new HashSet<String>();
            double maxChildHeight = 0.0;

            for (Node child : node.getChildren()) {
                double h = adjustTreeToConstraints(tree, child, l, cladeHeights);
                if (h > maxChildHeight) {
                    maxChildHeight = h;
                }
            }

            double height = node.getHeight();

            double lower = maxChildHeight;
            double upper = Double.POSITIVE_INFINITY;

            if (cladeHeights != null) {
                for (Set<String> clade : cladeHeights.keySet()) {
                    if (clade.equals(l)) {
                        double[] bounds = cladeHeights.get(clade);
                        lower = Math.max(bounds[0], maxChildHeight);
                        upper = bounds[1];
                    }
                }
            }

            if (lower > upper) {
                throw new IllegalArgumentException("incompatible constraints");
            }

            if (height < lower) {
                height = lower + 1E-6;
            } else if (height > upper) {
                height = (upper + lower) / 2;
            }
            node.setHeight(height);

            if (leaves != null) {
                leaves.addAll(l);
            }
        } else {
            leaves.add(tree.getTaxonId(node));
        }
        return node.getHeight();
    }

    private void setHeightsFromDates(FlexibleTree tree) {
        throw new UnsupportedOperationException("In development");
//        beast.evolution.util.Date mostRecent = null;
//        for (int i = 0; i < taxa.getTaxonCount(); i++) {
//            Date date = taxa.getTaxon(i).getDate();
//            if ((date != null) && (mostRecent == null || date.after(mostRecent))) {
//                mostRecent = date;
//            }
//        }
//
//        if (mostRecent != null) {
//            TimeScale timeScale = new TimeScale(mostRecent.getUnits(), true, mostRecent.getAbsoluteTimeValue());
//            double time0 = timeScale.convertTime(mostRecent.getTimeValue(), mostRecent);
//
//            double rootHeight = tree.getRoot().getHeight();
//            for (Node tip : tree.getExternalNodes()) {
//
//                Date date = tree.getNodeTaxon(tip).getDate();
//                if (date != null) {
//                    tree.setNodeHeight(tip, timeScale.convertTime(date.getTimeValue(), date) - time0);
//                } else {
//                    tree.setNodeHeight(tip, 0.0);
//                }
//            }
//        }
    }

}

