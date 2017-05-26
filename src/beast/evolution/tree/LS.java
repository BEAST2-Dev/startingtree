package beast.evolution.tree;

import beast.core.Citation;
import beast.core.Description;

/**
 * Objective function using least squares introduced in Linear Dating eq.2.
 * i = 1 root; internal nodes 2, ..., n-1; tips n ... 2n - 1
 */
@Citation("To et al 2016, 'Fast dating using least-squares criteria and algorithms', Systematic Biology")
@Description("Objective function using least squares introduced in Linear Dating")
public abstract class LS {
    protected final Node node; // do not work on the node in this class
    protected final double b;
    protected final double[] bs;
    protected final double[] ts;

    protected double x;
    protected double y;
    protected double z;

    protected double u;
    protected double v;
    protected double ua;
    protected double va;

    protected double t; // the final result, time of internal node

    /**
     * Constructor for root node, where w=0, b=0.
     * Detail at {@link #LS(Node, double, double[]) LS}
     */
    LS(Node root, double[] bs) {
        this(root, 0, bs);
    }

    /**
     * Objective function using least squares introduced in Linear Dating.
     *
     * @param node root or internal node <i>i</i>.
     * @param b    the length of the branch (<i>i</i>, <i>a(i)</i>),
     *             where <i>a(i)</i> is the ancestor of <i>i</i>.
     *             If <i>b</i>=0 (root), then make <i>w</i>=0
     *             from <code>getW(...)</code>, which satisfies eq.5.1.
     * @param bs   the array of the length of the branch (<i>i</i>, <i>s(i)</i>),
     *             where <i>s(i)</i> is the child of <i>i</i>.
     */
    LS(Node node, double b, double[] bs) {
        if (node.isLeaf())
            throw new IllegalArgumentException("Linear Dating requires internal nodes !");
        this.node = node;
        this.b = b;
        this.bs = bs;
        this.ts = new double[bs.length];
        // set ts[] later
    }

    /**
     * set time of <i>n<i/>th child node,
     * @param n
     * @param ts_n
     */
    public void setTS(int n, double ts_n) {
        if (n >= ts.length || Double.isNaN(ts_n))
            throw new IllegalArgumentException("Incorrect child index or date is NaN ! " +
                    "n = " + n + ", ts_n = " + ts_n);
        this.ts[n] = ts_n;
    }

    /**
     * calculate x_i, y_i, z_i from eq.6.i
     */
    public abstract void calculateXYZ();

    /**
     * calculate u_i, v_i from eq.7.i
     * @param ua u_a(i) u of ancestor of node i
     * @param va v_a(i)
     */
    public void calculateUV(double ua, double va) {
        this.ua = ua;
        this.va = va;
        this.u = getU(ua);
        this.v = getV(va);
    }

    /**
     * &omega;^ to minimise eq.2, (d)/(dx)(w (b - t x)^2) = 2 t w (t x - b) = 0,
     * where x is &omega;^, t is (t_i - t_ai).
     * @return &omega;^ to minimise eq.2
     */
    public double getOmegaNew() {
        if (node.isRoot())
            throw new RuntimeException("Least squares in eq.2 do not include root !");
        if (ua == 0 && va == 0)
            throw new RuntimeException("Ancestor u and v are not calculated previously !");
        return (b + va - v) / (u - ua);
    }

    /**
     * date of internal node calculated by eq.7
     * @param omega &omega;^
     * @return
     */
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

    protected double sum(double[] values) {
        double sum = 0;
        for (double v : values)
            sum += v;
        return sum;
    }

    // substitute eq. 7 [t_a(i) = u_a(i) + v_a(i) / omega]
    // into eq. 6 [t(i) = x(i) * t_a(i) + y(i) + z(i) / omega]
    // u(i) = x(i) * u_a(i) + y(i)
    protected double getU(double ua) {
        return getX() * ua + getY();
    }
    // v(i) = x(i) * v_a(i) + z(i)
    protected double getV(double va) {
        return getX() * va + getZ();
    }

    /**
     * get the time of internal node
     *
     * @return
     */
    public double getT() {
        return t;
    }

    /**
     * set time in the end of {@link LinearDating#analyseRootedTree(double, int, double) analyseRootedTree}
     *
     * @param t
     */
    public void setT(double t) {
        this.t = t;
    }
}
