package beast.evolution.tree;

import beast.core.Citation;
import beast.core.Description;
import beast.evolution.alignment.TaxonSet;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Linear Dating algorithm without constraints
 */
@Citation("To et al 2016, 'Fast dating using least-squares criteria and algorithms', Systematic Biology")
@Description("Linear dating without constraints")
public class LinearDating {

    private final Map<String, Double> dates;
//    private final Map<String, Double> precisions;
    private double dateMin;
    private double dateMax;

    private final FlexibleTree flexibleTree; // do not change the tree in this class

    private final Map<Node, LS> lsMap; // root or internal nodes

    private final int s; // sequence length
    private double c = 10; // recommended in paper

    // output
    private double omega;

    public LinearDating(FlexibleTree flexibleTree, TraitSet timeTraitSet, final int s, final double minOmega) {
        this.flexibleTree = flexibleTree;
        this.s = s;

        dates = new HashMap<>();
//        precisions = new HashMap<>();
        setDates(timeTraitSet);

        lsMap = new HashMap<>();

        analyseRootedTree(c, s, minOmega);
    }

    private void setDates(TraitSet timeTraitSet) {
        TaxonSet taxa = timeTraitSet.taxaInput.get();

        dateMin = timeTraitSet.minValue;
        dateMax = timeTraitSet.maxValue;

        for (int i = 0; i < taxa.getTaxonCount(); i++) {
            String taxon = taxa.getTaxonId(i);
            // note: timeTraitSet.getValue(taxon) is not the original value
            double date = Double.parseDouble(timeTraitSet.getStringValue(taxon));
            dates.put(taxon, date);
        }

        if (Math.abs(dateMax - dateMin) < 1.0E-8) // probably contemporaneous tips
            throw new IllegalArgumentException("Tips are contemporaneous ! ");
    }

    /**
     * A constant in eq.4, 10 is recommended in paper
     * @param c
     */
    public void setC(double c) {
        this.c = c;
    }

    /**
     * Start the Linear Dating algorithm.
     * <ol>
     * <li>post-order of internal nodes to deduce x, y, z;</li>
     * <li>pre-order of internal nodes to deduce u, v;</li>
     * <li>compute omega that minimize WLSRooted;</li>
     * <li>if &omega; < &omega;_min then &omega; = &omega;_min;</li>
     * <li>compute t(1, ..., n-1).</li>
     * </ol>
     *
     * @param c a constant in eq.4, 10 is recommended in paper
     * @param s sequence length
     * @param omegaMin &omega;_min where to fix the min of the estimated rate &omega;>= &omega;_min > 0
     */
    public void analyseRootedTree(final double c, final int s, final double omegaMin) {
        // deduce x, y, z
        Node root = flexibleTree.getRoot();
        List<Node> children = root.getChildren();
        double[] bs = flexibleTree.getBranchLengths(children);
        WLSRooted rootWLSRooted = new WLSRooted(root, bs, c, s);
        deduceXYZ(root, rootWLSRooted);

        // deduce u, v
        rootWLSRooted = (WLSRooted) lsMap.get(root);
        // u(root) = x * 0 + y(root), where set w(root) = b(root) = 0
        double u = rootWLSRooted.getY();
        double v = rootWLSRooted.getZ();
        for (Node child : root.getChildren()) {
            deduceUV(child, u, v);
        }

        // compute omega that minimize WLSRooted
        omega = computeOmega();

        // omega not lower than omegaMin
        if (omega < omegaMin)
            omega = omegaMin;

        // compute t(1, ..., n-1)
        computeT(omega);
    }

    public double getOmega() {
        return omega;
    }

    /**
     * get time given the root or an internal node in the final result
     *
     * @param node the root or an internal node
     * @return
     */
    public double getTime(Node node) {
        LS nodeWLSRooted = lsMap.get(node);
        if (nodeWLSRooted == null)
            throw new RuntimeException("WLSRooted for internal node " + node.getID() + " is not calculated !");

        return nodeWLSRooted.getT();
    }

    private double getTipDate(Node tip) {
        String taxon = tip.getID();
        Double date = dates.get(taxon);
        if (date == null)
            throw new IllegalArgumentException("Tip " + taxon + " does not have date ! ");
        return date;
    }

    private void deduceXYZ(Node node, WLSRooted nodeWLSRooted) {
        List<Node> children = node.getChildren();
        for (int n = 0; n < children.size(); n++) {
            Node child = children.get(n);
            if (child.isLeaf()) {
                double ts = getTipDate(node);
                nodeWLSRooted.setTS(n, ts);
            } else {
                List<Node> grandChildren = child.getChildren();
                double b = flexibleTree.getBranchLength(child);
                double[] bs = flexibleTree.getBranchLengths(grandChildren);
                WLSRooted childWLSRooted = new WLSRooted(child, b, bs, c, s);
                deduceXYZ(child, childWLSRooted);
            }
        }
        nodeWLSRooted.calculateXYZ();
        lsMap.put(node, nodeWLSRooted);
    }

    private void deduceUV(Node node, double ua, double va) {
        if (!node.isLeaf()) {
            LS nodeWLSRooted = lsMap.get(node);
            if (nodeWLSRooted == null)
                throw new IllegalArgumentException("Cannot find x, y, z for node : " + node.getID());

            nodeWLSRooted.calculateUV(ua, va);

            for (Node child : node.getChildren()) {
                deduceUV(child, nodeWLSRooted.getU(), nodeWLSRooted.getV());
            }
        }
    }

    private double computeOmega() {
        double omega = Double.MAX_VALUE;
        for (Node node : flexibleTree.getInternalNodes()) { // include root
            LS nodeWLSRooted = lsMap.get(node);
            if (nodeWLSRooted == null)
                throw new RuntimeException("WLSRooted for internal node " + node.getID() + " is not calculated !");
            double omegaN = nodeWLSRooted.getOmegaNew();
            if (omega > omegaN)
                omega = omegaN;
        }
        return omega;
    }

    // nodeWLS.setT(t)
    private void computeT(final double omega) {
        for (Node node : flexibleTree.getInternalNodes()) { // include root
            LS nodeWLSRooted = lsMap.get(node);
            if (nodeWLSRooted == null)
                throw new RuntimeException("WLSRooted for internal node " + node.getID() + " is not calculated !");
            double t = nodeWLSRooted.getTime(omega);
            nodeWLSRooted.setT(t);
        }
    }

    // Note: do not use weights (variances) in the objective function eq.9,
    // since weights depend on their associated branch length,
    // which are unknown for the two branches containing the assumed root.
    public void estimateRoot(final double c, final int s, final double omegaMin) {

    }
}
