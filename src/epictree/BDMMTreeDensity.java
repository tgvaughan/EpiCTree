package epictree;

import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.MultiTypeNode;
import beast.evolution.tree.Node;
import beast.math.Binomial;
import beast.util.Randomizer;
import multitypetree.distributions.MultiTypeTreeDistribution;
import org.apache.commons.math.distribution.BinomialDistribution;
import org.apache.commons.math.distribution.BinomialDistributionImpl;

import java.util.*;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class BDMMTreeDensity extends MultiTypeTreeDistribution {

    public Input<BDMigrationModel> migrationModelInput = new Input<>(
            "migModel", "Migration model", Input.Validate.REQUIRED);

    public Input<RealParameter> originInput = new Input<>(
            "origin", "Origin time", Input.Validate.REQUIRED);

    public Input<IntegerParameter> initialPopSizesInput = new Input<>(
            "initialPopSizes", "Initial population sizes",
            Input.Validate.REQUIRED);

    public Input<RealParameter> contempSamplingProbInput = new Input<>(
            "contempSamplingProb",
            "Probability of sampling lineages that survive until present.",
            Input.Validate.REQUIRED);

    public Input<Integer> nParticlesInput =new Input<>(
            "nParticles", "Number of particles to use in algorithm",
            Input.Validate.REQUIRED);

    protected BDMigrationModel migModel;
    protected RealParameter origin;
    protected RealParameter contempSamplingProb;
    protected  int nParticles;

    class TreeInterval {
        double startTime, endTime;
        MultiTypeNode node;

        int[] k;

        int terminalSrcType, terminalDestType;

        public double getLength() {
            return endTime-startTime;
        }
    }

    class Particle {

        int[] n;

        double[] birthProp = new double[migModel.getNTypes()];
        double[] deathProp = new double[migModel.getNTypes()];
        double[] samplingProp = new double[migModel.getNTypes()];
        double[] migProp = new double[migModel.getNTypes()*migModel.getNTypes()];
        double totalBirthProp, totalDeathProp, totalSamplingProp, totalMigProp;
        double totalProp;

        public Particle(int[] n) {
            this.n = Arrays.copyOf(n, n.length);
        }

        public Particle(Integer[] nList) {
            this.n = new int[nList.length];
            for (int i=0; i<nList.length; i++)
                this.n[i] = nList[i];
        }

        public Particle(Particle other) {
            this.n = Arrays.copyOf(other.n, other.n.length);
        }

        public void updatePropensities() {
            totalBirthProp = 0;
            totalDeathProp = 0;
            totalSamplingProp = 0;
            totalMigProp = 0;

            for (int i=0; i<migModel.getNTypes(); i++) {
                birthProp[i] = n[i]*migModel.getBirthRate(i);
                totalBirthProp += birthProp[i];

                deathProp[i] = n[i]*migModel.getDeathRate(i);
                totalDeathProp += deathProp[i];

                samplingProp[i] = n[i]*migModel.getSamplingRate(i);
                totalSamplingProp += samplingProp[i];

                for (int j=0; j<migModel.getNTypes(); j++) {
                    if (i == j)
                        continue;
                    migProp[i*migModel.getNTypes() + j] = n[i]*migModel.getForwardRate(i,j);
                    totalMigProp += migProp[i*migModel.getNTypes() + j];
                }
            }

            totalProp = totalBirthProp + totalDeathProp
                    + totalSamplingProp + totalMigProp;
        }

        public double getMigProp(int i, int j) {
            return migProp[i*migModel.getNTypes() + j];
        }
    }

    private List<TreeInterval> treeIntervalList = new ArrayList<>();


    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();

        migModel = migrationModelInput.get();
        origin = originInput.get();
        nParticles = nParticlesInput.get();
        contempSamplingProb = contempSamplingProbInput.get();
    }

    @Override
    public double calculateLogP() throws Exception {
        logP = 0.0;

        MultiTypeNode root = (MultiTypeNode) mtTree.getRoot();

        // Zero probability of events prior to origin.
        if (root.getFinalChangeTime()>origin.getValue())
            return Double.NEGATIVE_INFINITY;

        updateIntervalList();

        Particle[] particles = new Particle[nParticles];
        Particle[] particlesPrime = new Particle[nParticles];
        double[] weights = new double[nParticles];

        // Initialise particles and weights:
        for (int i=0; i<nParticles; i++) {
            particles[i] = new Particle(initialPopSizesInput.get().getValues());
            weights[i] = 1.0;
        }

        for (TreeInterval interval : treeIntervalList) {

            // Skip empty intervals
            // (These should only occur when sampling contemporaneously.)
            if (interval.getLength()==0.0)
                continue;

            // Update Particles
            double totalWeight = 0.0;
            for (int i=0; i<nParticles; i++) {
                weights[i] = updateParticle(particles[i], interval);
                totalWeight += weights[i];
            }

            // Include partial tree density
            logP += Math.log(totalWeight/nParticles);

            // Resample Particles
            for (int i=0; i<nParticles; i++) {
                particlesPrime[i] = new Particle(
                        particles[Randomizer.randomChoicePDF(weights)]);
            }

            Particle[] tmp = particles;
            particles = particlesPrime;
            particlesPrime = tmp;
        }

        // Include probability of contemporaneous sampling:

        if (contempSamplingProb.getValue()>0.0) {
            int nContemp = 0;
            for (Node node : mtTree.getExternalNodes()) {
                if (node.getHeight()==0.0)
                    nContemp += 1;
            }

            for (int pidx=0; pidx<nParticles; pidx++) {
                int nTotal = 0;
                for (int i=0; i<migModel.getNTypes(); i++) {
                    nTotal += particles[pidx].n[i];
                }
                if (nTotal<nContemp) {
                    weights[pidx] = 0.0;
                    continue;
                }

                if (contempSamplingProb.getValue()<1.0) {
                    double p = Math.log(contempSamplingProb.getValue());
                    weights[pidx] = Math.exp(Binomial.logChoose(nTotal, nContemp)
                            + nContemp * Math.log(p)
                            + (nTotal - nContemp) * Math.log(1.0 - p));
                    continue;
                }

                if (nContemp<nTotal)
                    weights[pidx] = 0.0;
                else
                    weights[pidx] = 1.0;
            }

            // Resample Particles
            for (int i=0; i<nParticles; i++) {
                particlesPrime[i] = new Particle(
                        particles[Randomizer.randomChoicePDF(weights)]);
            }

            Particle[] tmp = particles;
            particles = particlesPrime;
            particlesPrime = tmp;
        }

        return logP;
    }

    protected double updateParticle(Particle particle, TreeInterval interval) {

        double weight = 1.0;

        double t = interval.startTime;

        double[] birthProp = new double[migModel.getNTypes()];
        double[] deathProp = new double[migModel.getNTypes()];

        while (true) {

            // Compute propensities
            particle.updatePropensities();

            // Update next reaction time
            if (particle.totalProp > 0.0)
                t += Randomizer.nextExponential(particle.totalProp);
            else
                t = Double.POSITIVE_INFINITY;

            if (t>interval.endTime)
                break;

            // Implement reaction and update weight

            double u = Randomizer.nextDouble()*particle.totalProp;

            if (u<particle.totalSamplingProp) {
                // Sampling

                return 0.0; // Sampling ALWAYS results in weight 0
            }

            if (u<particle.totalBirthProp) {
                // Birth

                for (int i=0; i<migModel.getNTypes(); i++) {
                    if (u<particle.birthProp[i]) {
                        if (Randomizer.nextInt(particle.n[i]) < interval.k[i])
                            weight *= 2.0;

                        particle.n[i] += 1;
                        break;
                    }
                    u -= particle.birthProp[i];
                }

                continue;
            }
            u -= particle.totalBirthProp;

            if (u<particle.totalDeathProp) {
                for (int i=0; i<migModel.getNTypes(); i++) {
                    if (u<particle.deathProp[i]) {
                        if (Randomizer.nextInt(particle.n[i]) < interval.k[i])
                            return 0.0;

                        particle.n[i] -= 1;
                        break;
                    }
                    u -= particle.deathProp[i];
                }

                continue;
            }
            u -= particle.totalDeathProp;

            // Migration
            outer: {
                for (int i = 0; i < migModel.getNTypes(); i++) {
                    for (int j = 0; j < migModel.getNTypes(); j++) {
                        if (i == j)
                            continue;

                        if (u < particle.getMigProp(i, j)) {
                            if (Randomizer.nextInt(particle.n[i]) < interval.k[i])
                                return 0.0;

                            particle.n[i] -= 1;
                            particle.n[j] += 1;
                            break outer;
                        }

                        u -= particle.getMigProp(i, j);
                    }
                }
            }
        }

        if (interval.terminalDestType<0) {
            if (!interval.node.isLeaf()) {
                weight *= migModel.getBirthRate(interval.terminalSrcType);
            } else {
                if (interval.endTime < origin.getValue()) {
                    weight *= migModel.getSamplingRate(interval.terminalSrcType);
                }
            }
        } else {
            weight *= migModel.getForwardRate(
                    interval.terminalSrcType, interval.terminalDestType);
        }

        return weight;
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
