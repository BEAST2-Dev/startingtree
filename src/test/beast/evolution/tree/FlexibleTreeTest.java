package test.beast.evolution.tree;

import beast.evolution.tree.FlexibleTree;
import beast.evolution.tree.Node;
import junit.framework.TestCase;

/**
 * @author Walter Xie
 */
public class FlexibleTreeTest extends TestCase {
    public static String[] binaryTrees = new String[]{ // binary tree
            "((((A:1.0,B:1.0):1.0,C:2.0):2.0,D:3.0):3.0,E:5.0);",
            "(A:0.5,(B:1.0,(C:2.0,(D:3.0,E:8.0):2.0):1.0):0.5):0.0;",
            "(B:0.5,((C:2.0,(D:3.0,E:8.0):2.0):1.0,A:1.0):0.5):0.0;",
            "(C:1.0,((D:3.0,E:8.0):2.0,(A:1.0,B:1.0):1.0):1.0):0.0;",
            "(D:1.5,(E:8.0,((A:1.0,B:1.0):1.0,C:2.0):2.0):1.5):0.0;",
            "(E:4.0,(((A:1.0,B:1.0):1.0,C:2.0):2.0,D:3.0):4.0):0.0;",
            "((A:1.0,B:1.0):0.5,(C:2.0,(D:3.0,E:8.0):2.0):0.5):0.0;"
    };
    // todo "(A:1.0,B:2.0,(C:3.0,D:4.0):5.0);" // multifurcating tree

    double decimal = 100000;
    double mu = 0.4;


    public void testChangeRootTo() throws Exception {
        String tree = binaryTrees[0];
        FlexibleTree flexibleTree = new FlexibleTree(tree);

        System.out.println(flexibleTree.toNewick() + "\n");

        for (int i = 1; i < binaryTrees.length; i++) {
            // set new root between i and its parent with half length each side
            int nr = i-1;
            Node newRoot = flexibleTree.getNode(nr);

            System.out.println("Change the root at the lineage " + nr + " ascended from " + newRoot.getID());

            flexibleTree.changeRootTo(newRoot, 0.5);
            String newTree = flexibleTree.toNewick();

            System.out.println(newTree);

            if (!newTree.endsWith(";"))
                newTree += ";";

            assertEquals(binaryTrees[i], newTree);
        }

    }

    public void testRSS() throws Exception {

        FlexibleTree flexibleTree = new FlexibleTree(binaryTrees[0]);
        flexibleTree.setDateTrait(TemporalRootingTest.setUpTimeTrait());
        System.out.println(flexibleTree.toNewick() + "\n");

        double ss = flexibleTree.getRSS(mu);
        System.out.println("residual sum of squares = " + ss + " given mu = " + mu);

        assertEquals(19.44, Math.round(ss * decimal) / decimal);
    }

    public void testMinRSSTree() throws Exception {

        String minRSSTreeString = "(E:4.0,(((A:1.0,B:1.0):1.0,C:2.0):2.0,D:3.0):4.0):0.0;";

        FlexibleTree flexibleTree = new FlexibleTree(binaryTrees[0]);
        flexibleTree.setDateTrait(TemporalRootingTest.setUpTimeTrait());
        System.out.println(flexibleTree.toNewick() + "\n");

        FlexibleTree bestTree = flexibleTree.getMinRSSTree(mu);
        assertEquals(minRSSTreeString, bestTree.toNewick());

        double ss = bestTree.getRSS(mu);
        assertEquals(19.04, ss);
    }

}
