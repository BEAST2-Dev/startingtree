package beast.evolution.tree;

/**
 * @author Walter Xie
 */
public class StartingTree extends FlexibleTree {


    public StartingTree(final Node rootNode) {
        super(rootNode);
    }

    public FlexibleTree optimizeRoot() {
        FlexibleTree thisTree = new FlexibleTree(this.getRoot());
        for (int nr = 0; nr < this.getNodeCount(); nr++) {
            Node node = getNode(nr);
            // 2n-2
            if (!node.isRoot()) {
                thisTree.changeRootTo(node, 0.5);


            }
        }

        return thisTree;
    }

    public double getDistanceScore(FlexibleTree flexibleTree) {
        double score = 0;

        for (int i = 0; i < flexibleTree.getNodeCount(); i++) {
            for (int j = i+1; j < flexibleTree.getNodeCount(); j++) {

                System.out.println("i=" + i + ", j=" + j);

            }
        }


        return score;
    }

}
