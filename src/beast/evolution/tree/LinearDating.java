package beast.evolution.tree;

import beast.core.Citation;
import beast.core.Description;
import beast.evolution.alignment.TaxonSet;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * binary tree
 */
@Citation("To et al 2015, 'Fast dating using least-squares criteria and algorithms', Systematic Biology")
@Description("Linear dating without constraints")
public class LinearDating {

    private final Map<String, Double> dates;
//    private final Map<String, Double> precisions;
    private double dateMin;
    private double dateMax;

    private final FlexibleTree flexibleTree; // do not change the tree in this class

    private final Map<Node, WLS> wlsMap;

    private final int s; // sequence length
    private double c = 10; // recommended in paper

    // output
    private double omega;
    private double[] time;

    public LinearDating(FlexibleTree flexibleTree, TraitSet timeTraitSet, final int s, final double minOmega) {
        this.flexibleTree = flexibleTree;
        this.s = s;

        dates = new HashMap<>();
//        precisions = new HashMap<>();
        setDates(timeTraitSet);

        wlsMap = new HashMap<>();

        analyse(c, s, minOmega);
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

    public void analyse(final double c, final int s, final double minOmega) {
        // deduce x, y, z
        Node root = flexibleTree.getRoot();
        List<Node> children = root.getChildren();
        double[] bs = flexibleTree.getBranchLengths(children);
        WLS rootWLS = new WLS(root, bs, c, s);
        deduceXYZ(root, rootWLS);

        // deduce u, v
        rootWLS = wlsMap.get(root);
        // u(root) = x * 0 + y(root), where set w(root) = b(root) = 0
        double u = rootWLS.getY();
        double v = rootWLS.getZ();
        for (Node child : root.getChildren()) {
            deduceUV(child, u, v);
        }

        // compute omega that minimize WLS
        omega = computeOmega();

        // omega not lower than minOmega
        if (omega < minOmega)
            omega = minOmega;

        // compute t(1, ..., n-1)
        time = computeT(omega);
    }

    private double getTipDate(Node tip) {
        String taxon = tip.getID();
        Double date = dates.get(taxon);
        if (date == null)
            throw new IllegalArgumentException("Tip " + taxon + " does not have date ! ");
        return date;
    }

    private void deduceXYZ(Node node, WLS nodeWLS) {
        List<Node> children = node.getChildren();
        for (int n = 0; n < children.size(); n++) {
            Node child = children.get(n);
            if (child.isLeaf()) {
                double ts = getTipDate(node);
                nodeWLS.setTS(n, ts);
            } else {
                List<Node> grandChildren = child.getChildren();
                double b = flexibleTree.getBranchLength(child);
                double[] bs = flexibleTree.getBranchLengths(grandChildren);
                WLS childWLS = new WLS(child, b, bs, c, s);
                deduceXYZ(child, childWLS);
            }
        }
        nodeWLS.calculateXYZ();
        wlsMap.put(node, nodeWLS);
    }

    private void deduceUV(Node node, double ua, double va) {
        if (!node.isLeaf()) {
            WLS nodeWLS = wlsMap.get(node);
            if (nodeWLS == null)
                throw new IllegalArgumentException("Cannot find x, y, z for node : " + node.getID());

            nodeWLS.calculateUV(ua, va);

            for (Node child : node.getChildren()) {
                deduceUV(child, nodeWLS.getU(), nodeWLS.getV());
            }
        }
    }

    private double computeOmega() {
        double omega = Double.MAX_VALUE;
        for (Node node : flexibleTree.getInternalNodes()) { // include root
            WLS nodeWLS = wlsMap.get(node);
            if (nodeWLS == null)
                throw new RuntimeException("WLS for internal node " + node.getID() + " is not calculated !");
            double residual = nodeWLS.getResidual();
            if (omega > residual)
                omega = residual;
        }
        return omega;
    }

    private double[] computeT(final double omega) {
        List<Node> internalNodes = flexibleTree.getInternalNodes();
        double[] time = new double[internalNodes.size()];
        for (int i = 0; i < internalNodes.size(); i++) {
            Node node = internalNodes.get(i); // include root
            WLS nodeWLS = wlsMap.get(node);
            if (nodeWLS == null)
                throw new RuntimeException("WLS for internal node " + node.getID() + " is not calculated !");
            double t = nodeWLS.getT(omega);
            time[i] =t;
        }
        return time;
    }
}
