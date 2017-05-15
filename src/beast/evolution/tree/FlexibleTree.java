package beast.evolution.tree;

import beast.core.Description;
import beast.util.TreeParser;

import java.util.ArrayList;
import java.util.List;

@Description("Tree can be changed, such as re-root. Imported from BEAST 1 FlexibleTree.")
public class FlexibleTree extends Tree {

    boolean heightsKnown = false;
    boolean lengthsKnown = false;

    // Tree class does not support setBranchLength(), use double[];
    // index is Nr, the length is from node Nr to its parent, the length of the root is 0.
    protected double[] allBranchLengths = new double[getNodeCount()];


    public FlexibleTree(final String newick) {
        super(new TreeParser(newick, false, true, true, 1, false).getRoot());

        heightsKnown = hasNodeHeights();
        lengthsKnown = hasBranchLengths();
    }

    public FlexibleTree(final Node rootNode) {
        super(rootNode); // todo multifurcating tree ?

        heightsKnown = hasNodeHeights();
        lengthsKnown = hasBranchLengths();
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




    //++++++++ Time tree ++++++++

    /**
     * Calculate the sum of squared distances of branch length (distance)
     * given a tree or subtree <code>node</code>.
     *
     * @return the sum of squared residuals
     */
    public double getSumOfSquaredDistance() {
        return this.getSumOfSquaredDistances(getRoot());
    }

    /**
     * Post order traversal to calculate the sum of squared distances of branch length (distance)
     * given a tree or subtree <code>node</code>.
     *
     * @param node the root of the given tree or subtree
     * @return the sum of squared residuals
     */
    private double getSumOfSquaredDistances(Node node) {
        double ss = 0;
        for (Node child : node.getChildren()) {
            double d = node.getHeight() - child.getHeight();
            ss += getSumOfSquaredDistances(child) + d * d;
//            System.out.println(child.getNr() + " : " + ss + " , " + d);
        }
        return ss;
    }

    /**
     *
     * @param rootNode
     * @return
     */
    public FlexibleTree getMinSSDTree(final Node rootNode) {
        FlexibleTree tree = new FlexibleTree(rootNode);
        double minSSD = tree.getSumOfSquaredDistance();
        System.out.println("ssd = " + minSSD + ", tree = " + tree.toNewick());
        // all child nodes including this node
        for (Node node : rootNode.getAllChildNodes()) {
            if (!node.isRoot() && !node.getParent().isRoot()) {
                tree.changeRootTo(node, 0.5);
                double ssd = tree.getSumOfSquaredDistance();
                System.out.println("ssd = " + ssd + ", tree = " + tree.toNewick());
                if (ssd < minSSD) {
                    minSSD = ssd;
                }
            }
        }
        System.out.println("min sum of squared distances = " + minSSD);
        return tree;
    }


    public FlexibleTree getMinSSDTree() {
        return this.getMinSSDTree(getRoot());
    }
}
