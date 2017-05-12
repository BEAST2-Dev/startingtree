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
//            "(A:1.0,B:2.0,(C:3.0,D:4.0):5.0);" // multifurcating tree
    };
    String[] newTrees = new String[]{
            "(C:1.0,((A:1.0,B:1.0):1.0,(D:3.0,E:8.0):2.0):1.0):0.0;", // binary tree
            ""
    };

    public void testChangeRootTo() throws Exception {

        for (int i = 0; i < trees.length; i++) {
            String tree = trees[i];

            System.out.println("Change root of tree " + i );

            FlexibleTree flexibleTree = new FlexibleTree(tree);

            System.out.println(flexibleTree.toNewick());

            // set new root between C and its parent with half length each side
            Node newRoot = flexibleTree.getNode(2);
            flexibleTree.changeRootTo(newRoot, 0.5);
            String newTree = flexibleTree.toNewick();

            System.out.println(newTree);

            if (!newTree.endsWith(";"))
                newTree += ";";

            assertEquals(newTrees[i], newTree);
        }

    }

    public void testSumOfSquaredDistance() throws Exception {

        FlexibleTree flexibleTree = new FlexibleTree(trees[0]);
        System.out.println(flexibleTree.toNewick());

        double ss = flexibleTree.getSumOfSquaredDistance();
        System.out.println("sum of squared distances = " + ss);

        assertEquals(ss, 54.0);
    }

    public void testMinSSDTree() throws Exception {

        String minSSDTreeString = "(E:4.0,(((A:1.0,B:1.0):1.0,C:2.0):2.0,D:3.0):4.0):0.0;";

        FlexibleTree flexibleTree = new FlexibleTree(trees[0]);
        System.out.println(flexibleTree.toNewick());

        FlexibleTree minSSDTree = flexibleTree.getMinSSDTree();
        assertEquals(minSSDTree.toNewick(), minSSDTreeString);

        double ss = minSSDTree.getSumOfSquaredDistance();
        assertEquals(ss, 52.0);
    }

}
