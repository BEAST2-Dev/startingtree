package beast.evolution.tree;

import beast.util.ClusterTree;

/**
 * Starting point for BEAST using neighbor-joining tree topology and a treedater time-tree.
 * @see <a href="https://github.com/emvolz/treedater">treedater</a>
 */
public class TreeDater extends ClusterTree {

    public void initAndValidate() {
        clusterTypeInput.setValue(Type.neighborjoining, this);

        super.initAndValidate();



    }



}
