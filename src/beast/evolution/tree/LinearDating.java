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

    private final TaxonSet taxa;
    private final Map<String, Double> dates;
//    private final Map<String, Double> precisions;
    private double dateMin;
    private double dateMax;

    private final FlexibleTree flexibleTree;
    private final int s; // sequence length
    private double c = 10; // recommended in paper

    private final Map<Node, XYZ> xyzMap;
//    private final Map<Node, UV> uvMap;

    public LinearDating(FlexibleTree flexibleTree, TraitSet timeTraitSet, final int s) {
        this.flexibleTree = flexibleTree;
        this.s = s;
        this.taxa = timeTraitSet.taxaInput.get();

        dates = new HashMap<>();
//        precisions = new HashMap<>();

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

        xyzMap = new HashMap<>();
//        uvMap = new HashMap<>();
    }

    public double getDate(Node tip) {
        String taxon = tip.getID();
        Double date = dates.get(taxon);
        if (date == null)
            throw new IllegalArgumentException("Tip " + taxon + " does not have date ! ");
        return date;
    }



    private void deduceUV() {
        // root 1; internal nodes 2, ..., n-1;
        Node root = flexibleTree.getRoot();
        List<Node> children = root.getChildren();
        double[] bs = flexibleTree.getBranchLengths(children);
        XYZ rootXYZ = new XYZ(bs);
        deduceXYZ(root, rootXYZ);

        rootXYZ = xyzMap.get(root);
        // u(root) = x * 0 + y(root), where set w(root) = b(root) = 0
        double u = rootXYZ.getY();
        double v = rootXYZ.getZ();
        for (Node child : root.getChildren()) {
            deduceUV(child, u, v);
        }
    }

    private void deduceXYZ(Node node, XYZ nodeXYZ) {
        List<Node> children = node.getChildren();
        for (int n = 0; n < children.size(); n++) {
            Node child = children.get(n);
            if (child.isLeaf()) {
                double ts = getDate(node);
                nodeXYZ.setTS(n, ts);
            } else {
                List<Node> grandChildren = child.getChildren();
                double b = flexibleTree.getBranchLength(child);
                double[] bs = flexibleTree.getBranchLengths(grandChildren);
                XYZ childXYZ = new XYZ(b, bs);
                deduceXYZ(child, childXYZ);
            }
        }
        nodeXYZ.calculateXYZ();
        xyzMap.put(node, nodeXYZ);
    }

    private void deduceUV(Node node, double uA, double vA) {
        if (!node.isLeaf()) {
            XYZ xyz = xyzMap.get(node);
            if (xyz == null)
                throw new IllegalArgumentException("Cannot find x, y, z for node : " + node.getID());

            // substitute eq. 7 [t_a(i) = u_a(i) + v_a(i) / omega]
            // into eq. 6 [t(i) = x(i) * t_a(i) + y(i) + z(i) / omega]
            double u = xyz.getX() * uA + xyz.getY();
            double v = xyz.getX() * vA + xyz.getZ();

            for (Node child : node.getChildren()) {
                deduceUV(child, u, v);
            }
        }
    }



    // u(i) = x(i) * t_a(i)
    private double getU(double x, double ta, double y) {
        return x*ta + y;
    }
    // v(i) = z(i)
    private double getV(double z) {
        return z;
    }


    public class XYZ {
//        public final int i; // root 1; internal nodes 2, ..., n-1;
        private final double b; // branch length
        private final double w; // weight from eq. 4
        private final double[] bs;
        private final double[] ws;
        private final double[] ts;

        private double x;
        private double y;
        private double z;

        private double u;
        private double v;

        // i = 1, root
        XYZ (double[] bs) {
//            this.i = 1;
            this.b = 0;
            this.w = 0;
            this.bs = bs;
            this.ws = new double[bs.length]; // bs.length == 0 is tip
            this.ts = new double[bs.length];
            setWSFromBS();
        }

        // i > 1
        XYZ (double b, double[] bs) {
//            if (i < 2)
//                throw new IllegalArgumentException("Internal nodes index i should be 2 ... n-1 ! But i = " + i);
//            this.i = i;
            this.b = b;
            this.w = getW(b);
            this.bs = bs;
            this.ws = new double[bs.length]; // bs.length == 0 is tip
            this.ts = new double[bs.length];
            setWSFromBS();
        }

        public void setTS(int n, double tsi) {
            if (n >= ts.length || Double.isNaN(tsi))
                throw new IllegalArgumentException("Improper child index or date is NaN ! " +
                        "n = " + n + ", tsi = " + tsi);
            this.ts[n] = tsi;
        }

        public void calculateXYZ() {
            if (sum(ts) == 0 || sum(ws) == 0)
                throw new IllegalArgumentException("ts or ws are not initialised correctly, sum == 0 ! ");
            final double wsw = getWSW(w, ws);
            this.x = getX(w, wsw);
            this.y = getY(ws, ts, wsw);
            this.z = getZ(w, b, ws, bs, wsw);
        }

        public double getX() {
            return x;
        }

        public double getY() {
            return y;
        }

        public double getZ() {
            return z;
        }

        private double sum(double[] values) {
            double s = 0;
            for (double v : values)
                s += v;
            return s;
        }

        private void setWSFromBS() {
            if (bs.length < 1)
                throw new IllegalArgumentException("x, y, z is only working on internal nodes ! But child nodes = " + bs.length);
            for (int n = 0; n < bs.length; n++) {
                ws[n] = getW(bs[n]); // i+1 make sure it is not the root
            }
        }

        // w(i) / wsw(i)
        private double getX(double w, double wsw) {
            return w / wsw;
        }

        // (w_s1(i) * t_s1(i) + w_s2(i) * t_s2(i) + ...) / wsw(i)
        private double getY(double[] ws, double[] ts, double wsw) {
            if (ws.length != ts.length)
                throw new IllegalArgumentException("ws(i) should have the same length ts(i) !");
            double y = 0;
            for (int i=0; i<ws.length; i++)
                y += ws[i] * ts[i];
            return y / wsw;
        }

        // (w(i) * b(i) - w_s1(i) * b_s1(i) - w_s2(i) * b_s2(i) - ...) / wsw(i)
        private double getZ(double w, double b, double[] ws, double[] bs, double wsw) {
            if (ws.length != bs.length)
                throw new IllegalArgumentException("ws(i) should have the same length bs(i) !");
            double z = w*b;
            for (int i=0; i<ws.length; i++)
                z -= ws[i] * bs[i];
            return z / wsw;
        }

        // wsw(i), denominator of eq 5.i : w_s1(i) + w_s2(i) + ... + w(i), where w(1) = 0
        private double getWSW(double w, double... ws) {
            double wsw = w;
            for (double wsi : ws) {
                wsw += wsi;
            }
            return wsw;
        }

        // w(i), the weight that is calculated by eq 4
        private double getW(double b) {
            return s / (b + c/s);
        }

    }



}
