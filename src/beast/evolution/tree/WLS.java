package beast.evolution.tree;

/**
 * Weighted least squares used in Linear Dating eq.2.
 * i = 1 root; internal nodes 2, ..., n-1; tips n ... 2n - 1
 */
public class WLS {
    private final Node node; // do not work on the node in this class
    private final double b;
    private final double w; // weight from eq. 4
    private final double[] bs;
    private final double[] ws;
    private final double[] ts;

    private double x;
    private double y;
    private double z;

    private double u;
    private double v;
    private double ua;
    private double va;

    private double t; // the final result, time of internal node

    // root, w=0, b=0
    WLS(Node root, double[] bs, final double c, final int s) {
        this(root, 0, bs, c, s);
    }

    /**
     *
     * @param node root or internal node <i>i</i>.
     * @param b    the length of the branch (<i>i</i>, <i>a(i)</i>),
     *             where <i>a(i)</i> is the ancestor of <i>i</i>.
     * @param bs   the array of the length of the branch (<i>i</i>, <i>s(i)</i>),
     *             where <i>s(i)</i> is the child of <i>i</i>.
     * @param s    sequence length to calculate &sigma;
     */
    WLS(Node node, double b, double[] bs, final double c, final int s) {
        if (node.isLeaf())
            throw new IllegalArgumentException("Linear Dating requires internal nodes !");
        this.node = node;
        this.b = b;
        this.w = getW(b, c, s);
        this.bs = bs;
        this.ws = new double[bs.length]; // bs.length == 0 is tip
        this.ts = new double[bs.length];
        setWSFromBS(c, s);
        // set ts[] later
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

    public void calculateUV(double ua, double va) {
        this.ua = ua;
        this.va = va;
        this.u = getU(ua);
        this.v = getV(va);
    }

    // residual to minimise eq.2, (d)/(dx)(w (b - t x)^2) = 2 t w (t x - b) = 0
    public double getResidual() {
        if (node.isRoot())
            throw new RuntimeException("Least squares in eq.2 do not include root !");
        if (ua == 0 && va == 0)
            throw new RuntimeException("Ancestor u and v are not calculated previously !");
        return (b + va - v) / (u - ua);
    }

    // date of internal node calculated by eq.7
    public double getTime(double omega) {
        return u + v / omega;
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
    public double getU() {
        return u;
    }
    public double getV() {
        return v;
    }

    private double sum(double[] values) {
        double s = 0;
        for (double v : values)
            s += v;
        return s;
    }

    private void setWSFromBS(final double c, final int s) {
        if (bs.length < 1)
            throw new IllegalArgumentException("x, y, z is only working on internal nodes ! But child nodes = " + bs.length);
        for (int n = 0; n < bs.length; n++) {
            ws[n] = getW(bs[n], c, s); // i+1 make sure it is not the root
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

    // w(i), also 1/(&sigma(i))^2, the weight that is calculated by eq 4
    private double getW(double b, double c, double s) {
        return s / (b + c/s);
    }

    // substitute eq. 7 [t_a(i) = u_a(i) + v_a(i) / omega]
    // into eq. 6 [t(i) = x(i) * t_a(i) + y(i) + z(i) / omega]
    // u(i) = x(i) * u_a(i) + y(i)
    private double getU(double ua) {
        return getX() * ua + getY();
    }
    // v(i) = x(i) * v_a(i) + z(i)
    private double getV(double va) {
        return getX() * va + getZ();
    }

    /**
     * the time of internal node
     *
     * @return
     */
    public double getT() {
        return t;
    }

    /**
     * set time in the end of {@link LinearDating#analyse(double, int, double) analyse}
     *
     * @param t
     */
    public void setT(double t) {
        this.t = t;
    }
}
