package epictree;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.MigrationModel;
import beast.evolution.tree.MultiTypeNode;
import beast.evolution.tree.Node;
import multitypetree.distributions.MultiTypeTreeDistribution;

import java.util.*;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class BDMMTreeDensity extends MultiTypeTreeDistribution {

    public Input<MigrationModel> migrationModelInput = new Input<>(
            "migModel", "Migration model", Input.Validate.REQUIRED);

    public Input<RealParameter> originInput = new Input<>(
            "origin", "Origin time", Input.Validate.REQUIRED);

    private MigrationModel migModel;
    private RealParameter origin;

    class TreeInterval {
        double startTime, endTime;
        MultiTypeNode node;

        int[] k;

        int terminalSrcType, terminalDestType;
    }

    private List<TreeInterval> treeIntervalList;


    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();

        migModel = migrationModelInput.get();
        origin = originInput.get();
    }

    @Override
    public double calculateLogP() throws Exception {
        logP = 0.0;

        List<MultiTypeNode> nodeList = new ArrayList<>();
        for (Node node : mtTree.getNodesAsArray())
            nodeList.add((MultiTypeNode)node);

        nodeList.sort((o1, o2) -> {
            if (o1.getHeight() > o2.getHeight())
                return -1;

            if (o1.getHeight() < o2.getHeight())
                return 1;

            return 0;
        });

        MultiTypeNode root = nodeList.get(0);

        // Zero probability of events prior to origin.
        if (root.getFinalChangeTime()>origin.getValue())
            return Double.NEGATIVE_INFINITY;

        updateIntervalList();

        return logP;
    }

    /**
     * Determines the sequence of migration, coalescence and sampling events
     * which make up the coloured tree.
     */
    protected void updateIntervalList() {

        // Clean up previous list:
        treeIntervalList.clear();
        MultiTypeNode rootNode = (MultiTypeNode)mtTree.getRoot();

        Map<MultiTypeNode, Integer> activeChangeIndices = new HashMap<>();

        // Initialise map of active nodes to active change indices:
        activeChangeIndices.put(rootNode, rootNode.getChangeCount()-1);

        // Initialise lineage count per colour array:
        int[] k = new int[mtTree.getNTypes()];
        k[rootNode.getFinalType()] = 1;

        double t = 0.0;

        // Calculate event sequence:
        while (!activeChangeIndices.isEmpty()) {

            TreeInterval nextTreeInterval = new TreeInterval();
            nextTreeInterval.startTime = t;
            nextTreeInterval.endTime = Double.POSITIVE_INFINITY;
            nextTreeInterval.k = Arrays.copyOf(k, k.length);

            // Determine next event
            for (MultiTypeNode node : activeChangeIndices.keySet())
                if (activeChangeIndices.get(node)<0) {
                    if (origin.getValue() - node.getHeight() < nextTreeInterval.endTime) {
                        nextTreeInterval.endTime = origin.getValue() - node.getHeight();
                        nextTreeInterval.node = node;
                        nextTreeInterval.terminalSrcType = node.getNodeType();
                        nextTreeInterval.terminalDestType = -1;
                        nextTreeInterval.node = node;
                    }
                } else {
                    // Next event is a migration
                    int changeIdx = activeChangeIndices.get(node);
                    double thisChangeTime = origin.getValue() - node.getChangeTime(changeIdx);
                    if (thisChangeTime < nextTreeInterval.endTime) {
                        nextTreeInterval.endTime = thisChangeTime;
                        nextTreeInterval.terminalSrcType = node.getChangeType(changeIdx);
                        if (changeIdx>0)
                            nextTreeInterval.terminalDestType = node.getChangeType(changeIdx-1);
                        else
                            nextTreeInterval.terminalDestType = node.getNodeType();
                        nextTreeInterval.node = null;
                    }
                }

            // Update active node list (changeIdx) and lineage count appropriately:
            if (nextTreeInterval.terminalDestType<0) {
                if (nextTreeInterval.node.isLeaf()) {
                    // Sample
                    activeChangeIndices.remove(nextTreeInterval.node);
                    k[nextTreeInterval.terminalSrcType] -= 1;

                } else {
                    // Coalesce
                    MultiTypeNode leftChild = (MultiTypeNode) nextTreeInterval.node.getLeft();
                    MultiTypeNode rightChild = (MultiTypeNode) nextTreeInterval.node.getRight();

                    activeChangeIndices.remove(nextTreeInterval.node);
                    activeChangeIndices.put(leftChild, leftChild.getChangeCount() - 1);
                    activeChangeIndices.put(rightChild, rightChild.getChangeCount() - 1);

                    k[nextTreeInterval.terminalSrcType] += 1;
                }
            } else {
                // Migration
                int changeIdx = activeChangeIndices.get(nextTreeInterval.node);
                activeChangeIndices.put(nextTreeInterval.node, changeIdx - 1);

                k[nextTreeInterval.terminalSrcType] -= 1;
                k[nextTreeInterval.terminalDestType] += 1;
            }

            // Add interval to list:
            treeIntervalList.add(nextTreeInterval);

            // Update interval start time:
            t = nextTreeInterval.endTime;
        }
    }

}
