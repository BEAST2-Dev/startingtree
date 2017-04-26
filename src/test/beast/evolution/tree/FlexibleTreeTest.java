package test.beast.evolution.tree;

import beast.evolution.tree.FlexibleTree;
import beast.evolution.tree.Node;
import junit.framework.TestCase;

/**
 * @author Walter Xie
 */
public class FlexibleTreeTest extends TestCase {
    String[] trees = new String[]{
            "((((A:1.0,B:1.0):1.0,C:2.0):2.0,D:3.0):3.0,E:5.0);", // binary tree

    };
    String[] newTrees = new String[]{
            "(C:1.0,((A:1.0,B:1.0):1.0,(D:3.0,E:8.0):2.0):1.0):0.0;", // binary tree

    };

    public void testChangeRootTo() throws Exception {

        for (int i = 0; i < trees.length; i++) {
            String tree = trees[i];

            System.out.println("Change root of tree " + i );

            FlexibleTree flexibleTree = new FlexibleTree(tree);

            System.out.println(flexibleTree.getRoot().toNewick());

            // set new root between C and its parent with half length each side
            Node newRoot = flexibleTree.getNode(2);
            flexibleTree.changeRootTo(newRoot, 0.5);
            String newTree = flexibleTree.getRoot().toNewick();

            System.out.println(newTree);

            if (!newTree.endsWith(";"))
                newTree += ";";

            assertEquals(newTrees[i], newTree);
        }

    }
}
