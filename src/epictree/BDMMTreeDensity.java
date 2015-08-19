package epictree;

import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.MigrationModel;
import beast.evolution.tree.SCMigrationModel;
import beast.evolution.tree.MultiTypeNode;
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

    public Input<IntegerParameter> initialPopSizesInput = new Input<>(
            "initialPopSizes", "Initial population sizes",
            Input.Validate.REQUIRED);

    public Input<Integer> nParticlesInput =new Input<>(
            "nParticles", "Number of particles to use in algorithm",
            Input.Validate.REQUIRED);

    protected MigrationModel migModel;
    protected RealParameter origin;
    protected  int nParticles;

    class TreeInterval {
        double startTime, endTime;
        MultiTypeNode node;

        int[] k;

        int terminalSrcType, terminalDestType;
    }

    class Particle {

        int[] n;
        double weight = 0.0;

        public Particle(int[] n) {
            this.n = Arrays.copyOf(n, n.length);
        }

        public Particle(Integer[] nList) {
            this.n = new int[nList.length];
            for (int i=0; i<nList.length; i++)
                this.n[i] = nList[i];
        }

    }

    private List<TreeInterval> treeIntervalList = new ArrayList<>();


    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();

        migModel = migrationModelInput.get();
        origin = originInput.get();
    }

    @Override
    public double calculateLogP() throws Exception {
        logP = 0.0;

        MultiTypeNode root = (MultiTypeNode) mtTree.getRoot();

        // Zero probability of events prior to origin.
        if (root.getFinalChangeTime()>origin.getValue())
            return Double.NEGATIVE_INFINITY;

        updateIntervalList();

        List<Particle> particles = new ArrayList<>();

        // Initialise particles:
        for (int i=0; i<nParticlesInput.get(); i++) {
            particles.add(new Particle(initialPopSizesInput.get().getValues()));
        }

        return logP;
    }

    protected void updateParticle(Particle particle, TreeInterval interval) {

        double t = interval.startTime;

        while (true) {

            // Compute propensities

            // Update time

            // Implement reaction and update weight
        }
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
                        nextTreeInterval.node = node;
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
