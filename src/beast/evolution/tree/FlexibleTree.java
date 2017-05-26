package beast.evolution.tree;

import beast.core.Description;
import beast.evolution.alignment.TaxonSet;
import beast.math.statistic.Regression;
import beast.util.TreeParser;

import java.util.*;

@Description("Tree can be changed, such as re-root. Imported from BEAST 1 FlexibleTree.")
public class FlexibleTree extends Tree {

    boolean heightsKnown = false;
    boolean lengthsKnown = false;

    // Tree class does not support setBranchLength(), use double[];
    // index is Nr, the length is from node Nr to its parent, the length of the root is 0.
    // For unrooted tree, ignore root, and the actual branch length
    // between its two children r_1 and r_2 is allBranchLengths[Nr(r_1)] + allBranchLengths[Nr(r_2)].
    protected double[] allBranchLengths = new double[getNodeCount()];


    public FlexibleTree(final String newick) {
        super(new TreeParser(newick, false, true, true, 1, false).getRoot());
//        heightsKnown = true;
    }

    /**
     * This constructor wraps a new root node in a <code>FlexibleTree</code> object,
     * which does not copy the tree.
     * To copy a <code>FlexibleTree</code>, use {@link FlexibleTree#copy() FlexibleTree#copy} only,
     * which performs true deep copy.
     * For example, <code>oldTree.copy()</code>.
     * Note: <code>new FlexibleTree(oldTree.getRoot().copy())</code>
     * does not copy flags and <code>allBranchLengths</code>
     *
     * @param rootNode root <code>Node</code>
     */
    public FlexibleTree(final Node rootNode) {
        super(rootNode); // todo multifurcating tree
//        heightsKnown = true;
    }

    /**
     * deep copy, returns a completely new tree
     * @return FlexibleTree
     */
    public FlexibleTree copy() {
        FlexibleTree flexibleTree = new FlexibleTree(getRoot().copy());
        flexibleTree.setID(getID());
        flexibleTree.index = index;
        flexibleTree.heightsKnown = this.heightsKnown;
        flexibleTree.lengthsKnown = this.lengthsKnown;
        flexibleTree.allBranchLengths = Arrays.copyOf(this.allBranchLengths, this.allBranchLengths.length);
        if (hasDateTrait())
            flexibleTree.setDateTrait(getDateTrait());
        return flexibleTree;
    }

    public boolean hasNodeHeights() {
        return heightsKnown;
    }

    public boolean hasBranchLengths() {
        return lengthsKnown;
    }

    public double getBranchLength(Node node) {
        if (!lengthsKnown)
            calculateBranchLengths();
        int nodeNr = node.getNr();
        return allBranchLengths[nodeNr];
    }

    public double[] getBranchLengths(List<Node> nodes) {
        double[] bl = new double[nodes.size()];
        for (int i = 0; i < nodes.size(); i++) {
            Node node = nodes.get(i);
            bl[i] = getBranchLength(node);
        }
        return bl;
    }

    public void setBranchLength(Node node, double length) {
        if (!lengthsKnown)
            calculateBranchLengths();

        int nodeNr = node.getNr();
        allBranchLengths[nodeNr] = length;

        heightsKnown = false;

//        fireTreeChanged();
    }

    public void setNodeHeight(Node node, double height) {
        if (!heightsKnown)
            calculateNodeHeights();

        node.setHeight(height);

        lengthsKnown = false;

//        fireTreeChanged();
    }

    /**
     * Set the node heights from the current branch lengths.
     */
    protected void calculateNodeHeights() {
        if (!lengthsKnown)
            throw new IllegalArgumentException("Branch lengths not known");

        nodeLengthsToHeights(getRoot(), 0.0);

        double maxHeight = 0.0;
        Node node;
        for (Node tip : getExternalNodes()) {
            if (tip.getHeight() > maxHeight)
                maxHeight = tip.getHeight();
        }

        for (int i = 0; i < getNodeCount(); i++) {
            node = getNode(i);
            node.setHeight(maxHeight - node.getHeight());
        }

        heightsKnown = true;
    }

    /**
     * Set the node heights from the current node branch lengths.
     * Actually sets distance from root so the heights then need to be reversed.
     */
    private void nodeLengthsToHeights(Node node, double height) {

        // getBranchLength call setAllBranchLengths() in the first time
        double branchLength = getBranchLength(node);
        if (branchLength > 0.0)
            height += branchLength;

        node.setHeight(height);

        for (Node child : node.getChildren()) {
            nodeLengthsToHeights(child, height);
        }

    }

    /**
     * Calculate branch lengths from the current node heights.
     */
    protected void calculateBranchLengths() {
        List<Node> allChildNodes = getRoot().getAllChildNodes();
        for (Node child : allChildNodes) {
            int nodeNr = child.getNr();
            double branchLength = child.getLength();
//            if (lengths[nodeNr] > 0)
//                throw new IllegalArgumentException("Duplicate node Nr is invalid !");
            allBranchLengths[nodeNr] = branchLength;
        }

        lengthsKnown = true;
    }

    /**
     * Calculate branch lengths from the current node heights.
     */
    private void nodeHeightsToLengths(Node node, double height) {

        setBranchLength(node, height - node.getHeight());

        for (Node child : node.getChildren())
            nodeHeightsToLengths(child, node.getHeight());

    }

    /**
     * Re-root the tree on the branch above the given <code>node</code>
     * with the given new root.
     * <code>len(node, new_root) = len(node, parent) * propLen </code>
     *
     * @param node the new root
     * @param propLen the proportion of the branch length between <code>node</code>
     *                and its parent node to define the new root, such as 0.5.
     */
    public void changeRootTo(Node node, double propLen) {
        // todo non-binary tree re-rooting incorrectly
        if (!TreeUtils.isBinary(this))
            throw new IllegalArgumentException("changeRootTo is only available to binary tree !");

        Node node1 = node;
        Node parent = node1.getParent();
        if (parent == null || parent == root) {
            // the node is already the root so nothing to do...
            return;
        }

        hasStartedEditing = true;
        // todo m_tree.getState() == null
//        startEditing(null); // called in rm / add

        if (!lengthsKnown)
            calculateBranchLengths();

        Node parent2 = parent.getParent();

        // only change topology
        swapParentNode(parent, parent2, null);

        // the root is now free so use it as the root again
        parent.removeChild(node1);
        getRoot().addChild(node1);
        getRoot().addChild(parent);
        // adjust lengths for children of new root
        double nodeToParent = getBranchLength(node1);
        // setBranchLength change getBranchLength(node1)
        setBranchLength(node1, nodeToParent * propLen);
        setBranchLength(parent, nodeToParent * (1 - propLen));

        heightsKnown = false;

        hasStartedEditing = false; // todo is it correct to use restore()? no proposal
    }


    /**
     * Work up through the tree putting the parent into the child.
     */
    private void swapParentNode(Node node, Node parent, Node child) {

        if (parent != null) {
            Node parent2 = parent.getParent();

            swapParentNode(parent, parent2, node);

            if (child != null) {
                node.removeChild(child);
                child.addChild(node);
                setBranchLength(node, getBranchLength(child));
            }

        } else {
            // First remove child from the root
            node.removeChild(child);

            // can't remove from list if browsing it with "for each" loop
            List<Node> children = new ArrayList<>(node.getChildren());

            int numChild = children.size();
            if (numChild > 1) {
                // todo insert new internal node in the same position of old root for > 2 children
//                Node newInternalNode = new Node();
                for (int i=0; i<numChild; i++) {
                    Node tmp = children.get(i);
                    node.removeChild(tmp);
                    child.addChild(tmp);
                    setBranchLength(tmp, getBranchLength(tmp) + getBranchLength(child));
                }
            } else {
                Node tmp = children.get(0);
                node.removeChild(tmp);
                child.addChild(tmp);
                setBranchLength(tmp, getBranchLength(tmp) + getBranchLength(child));
            }
        }

    }

    public String toNewick() {
        if (lengthsKnown)
            return toNewickLengthsKnown(getRoot(), false) + ";";
        return this.getRoot().toNewick() + ";";
    }

    // lengthsKnown == true
    private String toNewickLengthsKnown(Node node, boolean onlyTopology) {
        final StringBuilder buf = new StringBuilder();
        if (!node.isLeaf()) {
            buf.append("(");
            boolean isFirst = true;
            for (Node child : node.getChildren()) {
                if (isFirst)
                    isFirst = false;
                else
                    buf.append(",");
                buf.append(toNewickLengthsKnown(child, onlyTopology));
            }
            buf.append(")");

            if (node.getID() != null)
                buf.append(node.getID());
        } else {
            if (node.getID() != null)
                buf.append(node.getID());
            else
                buf.append(node.labelNr);
        }

        if (!onlyTopology) {
            buf.append(node.getNewickMetaData());
            buf.append(":").append(node.getNewickLengthMetaData()).append(getBranchLength(node));
        }
        return buf.toString();
    }

    public boolean isRoot(Node node) {
        return (node == getRoot());
    }




    //++++++++ Time tree ++++++++ // todo incorrect ?
    private Map<String, Double> dates;
    private double dateMin;
    private double dateMax;

    public void setDates(TraitSet timeTraitSet) {
        TaxonSet taxa = timeTraitSet.taxaInput.get();

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
    }


    /**
     * Calculate the sum of squared distances of branch length (distance)
     * given a tree or subtree <code>node</code>.
     *
     * @return the sum of squared residuals
     */
    public double getRSS(double mu) {

        if (!lengthsKnown)
            calculateBranchLengths();
        return this.getRSS(getRoot(), mu);
    }

    /**
     * Post order traversal to calculate the residual sum of squares
     * given a tree or subtree <code>node</code>.
     *
     * @param node the root of the given tree or subtree
     * @return the sum of squared residuals
     */
    private double getRSS(Node node, double mu) {
        double rss = 0;
        for (Node child : node.getChildren()) {
            // b_i
            double b = getBranchLength(child);
            // t_i - t_a(i)
            double t = child.getDate() - node.getDate(); // todo wrong
//            rss += getRSS(child, mu) + (b - mu*t)*(b - mu*t) / sigma * sigma;
//            System.out.println(child.getNr() + " : " + ss + " , " + d);
        }
        return rss;
    }

    /**
     *
     * @return
     */
    public FlexibleTree getMinRSSTree(double mu) {

        final FlexibleTree source = this.copy();
        if (mu <= 0) {
            TemporalRooting temporalRooting = new TemporalRooting(getDateTrait());
            Regression r = temporalRooting.getRootToTipRegression(source);
            // todo correct ?
            mu = r.getGradient();
        }

        double minRSS = source.getRSS(mu);
        System.out.println("Init minRSS = " + minRSS + ", tree = " + source.toNewick());
        // all child nodes including this node
        FlexibleTree bestTree = source.copy();
        for (Node node : source.getRoot().getAllChildNodes()) {
            if (!node.isRoot() && !node.getParent().isRoot()) {
                source.changeRootTo(node, 0.5);
                double rss = source.getRSS(mu);
                System.out.println("minRSS = " + minRSS + ", rss = " + rss + ", tree = " + source.toNewick());
                if (rss < minRSS) {
                    minRSS = rss;
                    bestTree = source.copy();
                }
            }
        }
        System.out.println("min residual sum of squares = " + minRSS);
        System.out.println("bestTree = " + bestTree.toNewick());
        return bestTree;
    }

}
