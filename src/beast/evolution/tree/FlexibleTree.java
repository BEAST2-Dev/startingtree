package beast.evolution.tree;

import beast.core.Description;
import beast.util.TreeParser;

import java.util.ArrayList;
import java.util.List;

@Description("Tree can be changed, such as re-root. Imported from BEAST 1 FlexibleTree.")
public class FlexibleTree extends Tree {

    // Tree class does not support setBranchLength(), use double[];
    // index is Nr, the length is from node Nr to its parent, the length of the root is 0.
    protected double[] allBranchLengths;


    public FlexibleTree(final String newick) {
        super(new TreeParser(newick, false, true, true, 1, false).getRoot());
    }

    public FlexibleTree(final Node rootNode) {
        super(rootNode); // todo multifurcating tree ?
    }

    public void setAllBranchLengths() {
        allBranchLengths = getAllBranchLengths(getRoot(), getNodeCount());
    }

    public double[] getAllBranchLengths() {
        if (allBranchLengths == null)
            setAllBranchLengths();
        return allBranchLengths;
    }

    public double getBranchLength(Node node) {
        if (allBranchLengths == null)
            setAllBranchLengths();
        int nodeNr = node.getNr();
        return allBranchLengths[nodeNr];
    }

    public void setBranchLength(Node node, double branchLength) {
        int nodeNr = node.getNr();
        allBranchLengths[nodeNr] = branchLength;
    }

    /**
     * Get all branch lengths of sub/tree.
     * If lengths[nodeNr] == 0, then either is root or not child node
     *
     * @param node a given node, such as <code>getRoot()</code>
     * @param maxNr the max index of all nodes, such as <code>getNodeCount()</code>
     * @return
     */
    private double[] getAllBranchLengths(final Node node, int maxNr) {
        double[] lengths = new double[maxNr];
        List<Node> allChildNodes = node.getAllChildNodes();
        for (Node child : allChildNodes) {
            int nodeNr = child.getNr();
            double branchLength = child.getLength();
//            if (lengths[nodeNr] > 0)
//                throw new IllegalArgumentException("Duplicate node Nr is invalid !");
            lengths[nodeNr] = branchLength;
        }
        return lengths;
    }


    /**
     * Get the maximum node height of the sub/tree including given <code>node</code>
     *
     * @param node the root of the given tree or subtree
     * @return
     */
    public static double getMaxNodeHeight(Node node) {
        if (!node.isLeaf()) {
            double maxNodeHeight = 0;
            for (Node child : node.getAllChildNodes()) {
                double childHeight = child.getHeight();
                if (maxNodeHeight < childHeight)
                    maxNodeHeight = childHeight;
            }
            return maxNodeHeight;
        } else return node.getHeight();
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

        Node parent = node.getParent();
        if (parent == null || parent == root) {
            // the node is already the root so nothing to do...
            return;
        }

        hasStartedEditing = true;
        // todo m_tree.getState() == null
//        startEditing(null); // called in rm / add

        setAllBranchLengths();

        Node parent2 = parent.getParent();

        // only change topology
        swapParentNode(parent, parent2, null);

        // the root is now free so use it as the root again
        parent.removeChild(node);
        getRoot().addChild(node);
        getRoot().addChild(parent);
        // adjust lengths for children of new root
        double nodeToParent = getBranchLength(node);
        // setBranchLength change getBranchLength(node)
        setBranchLength(node, nodeToParent * propLen);
        setBranchLength(parent, nodeToParent * (1 - propLen));

        setNodeHeightsByLengths(node, parent, propLen);
        // update lengths after set heights
        setAllBranchLengths();

        hasStartedEditing = false; // todo is it correct to use restore()? no proposal
    }

    /**
     * Set the node heights from the given branch lengths.
     */
    private void setNodeHeightsByLengths(Node child1, Node child2, double propLen) {

        nodeLengthsToHeights(getRoot(), 0.0);

        double maxHeight = FlexibleTree.getMaxNodeHeight(getRoot());

        for (int i = 0; i < getNodeCount(); i++) {
            Node node = getNode(i);
            // Set the node heights to the reversed heights
            node.setHeight(maxHeight - node.getHeight());
        }

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
        return this.getRoot().toNewick() + ";";
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
