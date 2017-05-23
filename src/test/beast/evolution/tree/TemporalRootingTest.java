package test.beast.evolution.tree;

import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.FlexibleTree;
import beast.evolution.tree.TemporalRooting;
import beast.evolution.tree.TraitSet;
import beast.math.statistic.Regression;
import beast.math.statistic.Variate;
import junit.framework.TestCase;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Walter Xie
 */
public class TemporalRootingTest extends TestCase {
//    String[] trees = new String[]{
//            "((((A:1.0,B:1.0):1.0,C:2.0):2.0,D:3.0):3.0,E:5.0);",
//            "((D:3.0,E:8.0):1.0,((A:1.0,B:1.0):1.0,C:2.0):1.0);",
//            "((((A:1.0,B:1.0):1.0,C:2.0):2.0,E:8.0):1.5,D:1.5);"
//    };

    double decimal = 100000;
    TemporalRooting temporalRooting;

    @Override
    public void setUp() throws Exception {
        TraitSet timeTraitSet = setUpTimeTrait();
        temporalRooting = new TemporalRooting(timeTraitSet);
    }

    public static TraitSet setUpTimeTrait() {
        String[] taxa = new String[]{ "A","B","C","D","E" };
        List<Sequence> seqList = new ArrayList<Sequence>();
        for (String taxonID : taxa)
            seqList.add(new Sequence(taxonID, "?"));

        Alignment alignment = new Alignment(seqList, "nucleotide");
        TaxonSet taxonSet = new TaxonSet(alignment);
        TraitSet timeTraitSet = new TraitSet();
        timeTraitSet.initByName(
                "traitname", "date",
                "taxa", taxonSet,
                "value", "A=2017, B=2017, C=2016, D=2015, E=2012");
        return timeTraitSet;
    }

    public void testRootToTipRegression() {

//        for (int i = 0; i < trees.length; i++) {
        FlexibleTree rootedTree = new FlexibleTree(FlexibleTreeTest.binaryTrees[0]);
        System.out.println(rootedTree.getRoot().toNewick());

        Regression r = temporalRooting.getRootToTipRegression(rootedTree);
        Variate dates = r.getXData();
        Variate distances = r.getYData();

        System.out.println("date\tdistance");
        for (int i = 0; i < dates.getCount(); i++)
            System.out.println(dates.get(i) + "\t" + distances.get(i));

        assertEquals(0.41860, r.getGradient(), 1e-5);
        assertEquals(2000.11111, r.getXIntercept(), 1e-5);
        assertEquals(-837.25581, r.getYIntercept(), 1e-5);
        assertEquals(0.06202, r.getResidualMeanSquared(), 1e-5);
        assertEquals(0.94186, r.getRSquared(), 1e-5);
        assertEquals(0.97049, r.getCorrelationCoefficient(), 1e-5);
//        }

        double[] residuals = temporalRooting.getRootToTipResiduals(rootedTree, r);
        String[] tipLabels = temporalRooting.getTipLabels(rootedTree);
        for (int i = 0; i < residuals.length; i++)
            System.out.println(i + " " + tipLabels[i] + " " + residuals[i]);

        assertEquals(-0.06977, residuals[0], 1e-5);

    }

    // todo limited to set new root between internal nodes
    public void testFindLocalRoot() {
        FlexibleTree rootedTree = new FlexibleTree(FlexibleTreeTest.binaryTrees[1]);
        System.out.println(rootedTree.toNewick() + "\n");

        FlexibleTree bestTree = temporalRooting.findLocalRoot(rootedTree, TemporalRooting.RootingFunction.HEURISTIC_RESIDUAL_MEAN_SQUARED);
        System.out.println(bestTree.toNewick());

        assertEquals("", bestTree.toNewick());

//        FlexibleTree bestTree = temporalRooting.findLocalRoot(rootedTree, TemporalRooting.RootingFunction.RESIDUAL_MEAN_SQUARED);
//        System.out.println(bestTree.toNewick());
//
//        assertEquals("", bestTree.toNewick());

    }

    public void testFindRoot() {
        FlexibleTree rootedTree = new FlexibleTree(FlexibleTreeTest.binaryTrees[0]);
        System.out.println(rootedTree.toNewick() + "\n");

        FlexibleTree bestTree = temporalRooting.findRoot(rootedTree, TemporalRooting.RootingFunction.HEURISTIC_RESIDUAL_MEAN_SQUARED);
        System.out.println(bestTree.toNewick());

    }

}
