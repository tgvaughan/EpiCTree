package epictree.util;

import beast.core.*;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Utility for mapping probability density surfaces.")
public class DensityMapper extends beast.core.Runnable {

    public Input<List<RealParameter>> realParamsInput = new Input<>(
            "realParam",
            "Real parameter to vary.",
            new ArrayList<>());

    public Input<List<IntegerParameter>> stepsInput = new Input<>(
            "steps",
            "Number of steps to take between parameter bounds (minimum 2).",
            new ArrayList<>());

    public Input<List<Distribution>> distribsInput = new Input<>(
            "distribution",
            "Distribution to map.",
            new ArrayList<>());


    public Input<List<Logger>> loggersInput = new Input<>(
            "logger",
            "Logger used to store output.",
            new ArrayList<>());

    int nValues;
    int sample;

    State dummyState;

    public DensityMapper() { }

    @Override
    public void initAndValidate() throws Exception {

        if (realParamsInput.get().size() != stepsInput.get().size())
            throw new IllegalArgumentException("Number of step sizes " +
                    "must match number of real params.");

        // Calculate total number of params to vary
        nValues = 0;
        for (int i=0; i<realParamsInput.get().size(); i++) {
            int thisN = realParamsInput.get().get(i).getDimension();
            if (stepsInput.get().get(i).getDimension() != thisN)
                throw new IllegalArgumentException("Dimension of step sizes " +
                        "must match dimension of real params.");

            nValues += thisN;
        }

        // Add RealParameters to dummy state:
        dummyState = new State();
        for (RealParameter param : realParamsInput.get())
            dummyState.stateNodeInput.setValue(param, dummyState);
        dummyState.initAndValidate();
        dummyState.initialise();

        sample = 0;

    }

    @Override
    public void run() throws Exception {
        // Initialise loggers:
        for (Logger logger : loggersInput.get())
            logger.init();

        // Loop over parameter values
        nestedLoop(0);

        // Finalize loggers
        for (Logger logger : loggersInput.get())
            logger.close();
    }

    /**
     * Perform nested loop over combinations of parameter values
     *
     * @param depth initialize to zero,
     */
    void nestedLoop(int depth) {

        if (depth<nValues) {
            int paramIdx = 0;
            int elIdx = 0;
            int count = 0;
            for (paramIdx = 0; paramIdx<realParamsInput.get().size(); paramIdx++) {
                RealParameter thisParam = realParamsInput.get().get(paramIdx);
                if (count + thisParam.getDimension() > depth) {
                    elIdx = depth - count;
                    break;
                }
                count += thisParam.getDimension();
            }

            int nSteps = stepsInput.get().get(paramIdx).getValue(elIdx);
            RealParameter param = realParamsInput.get().get(paramIdx);

            double delta = (param.getUpper()-param.getLower())/(nSteps-1);

            for (int i=0; i<stepsInput.get().get(paramIdx).getValue(elIdx); i++) {
                param.setValue(elIdx, param.getLower() + i*delta);

                nestedLoop(depth + 1);
            }

        } else {
            for (Distribution distrib : distribsInput.get())
                try {
                    distrib.calculateLogP();
                } catch (Exception e) {
                    System.out.println("Error computing density.");
                    e.printStackTrace();
                }
            for (Logger logger : loggersInput.get()) {
                logger.log(sample);
            }

            sample += 1;
        }

    }

}
