package beast.evolution.tree;

import beast.core.Citation;
import beast.core.Description;

/**
 * Weighted least squares used in Linear Dating eq.2.
 * i = 1 root; internal nodes 2, ..., n-1; tips n ... 2n - 1
 */
@Citation("To et al 2016, 'Fast dating using least-squares criteria and algorithms', Systematic Biology")
@Description("Weighted least squares used in Linear Dating for a rooted tree")
public class WLSRooted extends LS {
    private final double w; // weight from eq. 4
    private final double[] ws;

    /**
     * Constructor for root node, where w=0, b=0.
     * Detail at {@link #WLSRooted(Node, double, double[], double, int) WLSRooted}
     */
    WLSRooted(Node root, double[] bs, final double c, final int s) {
        this(root, 0, bs, c, s);
    }

    /**
     * Weighted least squares used in Linear Dating for a rooted tree.
     * Set <i>sl</i>=0 to be the unweighted version, where <i>w_i</i>=1.
     *
     * @param node root or internal node <i>i</i>.
     * @param b    the length of the branch (<i>i</i>, <i>a(i)</i>),
     *             where <i>a(i)</i> is the ancestor of <i>i</i>.
     *             If <i>b</i>=0 (root), then make <i>w</i>=0
     *             from <code>getW(...)</code>, which satisfies eq.5.1.
     * @param bs   the array of the length of the branch (<i>i</i>, <i>s(i)</i>),
     *             where <i>s(i)</i> is the child of <i>i</i>.
     * @param sl    sequence length to calculate &sigma;.
     *             If <i>s</i>=0, then <i>w_i</i> from <code>getW(...)</code>
     *             always equals to 1, which is the unweighted version.
     */
    WLSRooted(Node node, double b, double[] bs, final double c, final int sl) {
        super(node, b, bs);
        w = getW(b, c, sl);
        ws = new double[bs.length]; // bs.length == 0 is tip
        setWSFromBS(c, sl);
        // set ts[] later
    }

    private void setWSFromBS(final double c, final int s) {
        if (bs.length < 1)
            throw new IllegalArgumentException("x, y, z is only working on internal nodes ! But child nodes = " + bs.length);
        for (int n = 0; n < bs.length; n++) {
            ws[n] = getW(bs[n], c, s); // i+1 make sure it is not the root
        }
    }

    public void calculateXYZ() {
        if (sum(ts) == 0 || sum(ws) == 0)
            throw new IllegalArgumentException("ts or ws are not initialised correctly, sum == 0 ! ");
        final double wsw = getWSW(w, ws);
        x = getX(w, wsw);
        y = getY(ws, ts, wsw);
        z = getZ(w, b, ws, bs, wsw);
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
