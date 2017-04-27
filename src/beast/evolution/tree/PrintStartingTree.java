package beast.evolution.tree;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.util.ClusterTree;

import java.io.FileNotFoundException;
import java.io.PrintStream;

/**
 * Create a starting tree using {@see beast.util.ClusterTree}, and print it.
 *
 * examples/nj.xml
 *
 * @author Walter Xie
 */
public class PrintStartingTree extends beast.core.Runnable {

    final public Input<ClusterTree> startingTreeInput = new Input<>("init", "starting tree from ClusterTree", Validate.REQUIRED);

    final public Input<String> outputFileNameInput = new Input<>( "outputFileName",
            "If provided, starting tree is written to this file rather than to standard out.");

    /**
     * starting tree
     */
    protected ClusterTree startingTree;
    /**
     * name of output file *
     */
    String outputFileName;

    public void initAndValidate() {

        startingTree = startingTreeInput.get();

        outputFileName = outputFileNameInput.get();

    }


    public void run() throws Exception {
        FlexibleTree flexibleTree = new FlexibleTree(startingTree.getRoot());


        printTree(flexibleTree);


    }

    private void printTree(FlexibleTree flexibleTree) throws FileNotFoundException {
        // Write output to stdout or file
        PrintStream pstream;
        if (outputFileName == null)
            pstream = System.out;
        else
            pstream = new PrintStream(outputFileName);

//        pstream.println(startingTree.getRoot().toNewick());
        pstream.println(flexibleTree.toNewick());
    }


}
