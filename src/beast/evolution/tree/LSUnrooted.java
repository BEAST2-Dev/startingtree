package beast.evolution.tree;

import beast.core.Citation;
import beast.core.Description;

/**
 * Unweighted least squares in eq.9. for root estimation
 * i = 1 root; internal nodes 2, ..., n-1; tips n ... 2n - 1
 */
@Citation("To et al 2016, 'Fast dating using least-squares criteria and algorithms', Systematic Biology")
@Description("Unweighted least squares used in modified Linear Dating for an unrooted tree")
public class LSUnrooted extends LS {

    private double p = 0.5; // root position [0, 1]

    /**
     * Constructor for root node, where w=0, b=0.
     * Detail at {@link #LSUnrooted(Node, double, double[]) LSUnrooted}
     */
    LSUnrooted(Node root, double[] bs) {
        this(root, 0, bs);
    }

    /**
     * Unweighted least squares in eq.9. for root estimation.
     * Note: do not use weights (variances) in the objective function eq.9,
     * since weights depend on their associated branch length,
     * which are unknown for the two branches containing the assumed root.
     *
     * @see beast.evolution.tree.LS#LS(Node, double, double[])
     */
    LSUnrooted(Node node, double b, double[] bs) {
        super(node, b, bs);

        // set ts[] later
    }


    public void calculateXYZ() {

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

    // w(i), also 1/(&sigma(i))^2, the weight that is calculated by eq 4
    private double getW(double b, double c, double s) {
        return s / (b + c/s);
    }

}
