/*
 * Copyright (C) 2012 Tim Vaughan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package epictree;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.MigrationModel;
import org.jblas.DoubleMatrix;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Tim Vaughan
 */
@Description("Basic plugin describing a simple Markovian migration model.")
public class BDMigrationModel extends CalculationNode implements MigrationModel {

    public Input<RealParameter> rateMatrixInput = new Input<>(
            "rateMatrix",
            "Migration rate matrix (forward time)",
            Validate.REQUIRED);

    public Input<RealParameter> birthRatesInput = new Input<>(
            "birthRates",
            "Birth rates", Validate.REQUIRED);

    public Input<RealParameter> deathRatesInput = new Input<>(
            "deathRates",
            "Death rates", Validate.REQUIRED);

    protected RealParameter rateMatrix, birthRates, deathRates;
    protected double mu, muSym;
    protected int nTypes;
    protected DoubleMatrix Q, R;
    protected DoubleMatrix Qsym, Rsym;
    protected List<DoubleMatrix> RpowN, RsymPowN;
    protected DoubleMatrix RpowMax, RsymPowMax;
    protected boolean RpowSteady, RsymPowSteady;

    protected boolean rateMatrixIsSquare, symmetricRateMatrix;

    // Flag to indicate whether EV decompositions need updating.
    protected boolean dirty;

    public BDMigrationModel() {
        // Initialise caching array for powers of uniformized
        // transition matrix:
        RpowN = new ArrayList<>();
        RsymPowN = new ArrayList<>();
    }

    @Override
    public void initAndValidate() throws Exception {

        birthRates = birthRatesInput.get();
        deathRates = deathRatesInput.get();
        rateMatrix = rateMatrixInput.get();
        nTypes = birthRatesInput.get().getDimension();

        if (deathRates.getDimension() != birthRates.getDimension())
            throw new IllegalArgumentException("Birth and death rates do not " +
                    "have the same dimension.");

        rateMatrix.setLower(Math.max(rateMatrix.getLower(), 0.0));

        if (rateMatrix.getDimension() == nTypes*nTypes) {
            rateMatrixIsSquare = true;
            symmetricRateMatrix = false;
        } else {
            if (rateMatrix.getDimension() != nTypes*(nTypes-1)) {
                if (rateMatrix.getDimension() == nTypes*(nTypes-1)/2) {
                    symmetricRateMatrix = true;
                } else 
                    throw new IllegalArgumentException("Migration matrix has "
                            + "incorrect number of elements for given deme count.");
            } else {
                rateMatrixIsSquare = false;
                symmetricRateMatrix = false;
            }
        }
        

        dirty = true;
        updateMatrices();
    }

    /**
     * Ensure all local fields including matrices and eigenvalue decomposition
     * objects are consistent with current values held by inputs.
     */
    public void updateMatrices()  {

        if (!dirty)
            return;

        mu = 0.0;
        muSym = 0.0;
        Q = new DoubleMatrix(nTypes, nTypes);
        Qsym = new DoubleMatrix(nTypes, nTypes);

        // Set up backward transition rate matrix Q and symmetrized backward
        // transition rate matrix Qsym:
        for (int i = 0; i < nTypes; i++) {
            Q.put(i,i, 0.0);
            Qsym.put(i,i, 0.0);
            for (int j = 0; j < nTypes; j++) {
                if (i != j) {
                    Q.put(i, j, getBackwardRate(i, j));
                    Q.put(i, i, Q.get(i, i) - Q.get(i, j));
                    
                    Qsym.put(i, j, 0.5*(getBackwardRate(i, j) + getBackwardRate(j, i)));
                    Qsym.put(i, i, Qsym.get(i, i) - Qsym.get(i, j));
                }
            }

            if (-Q.get(i, i) > mu)
                mu = -Q.get(i, i);
            
            if (-Qsym.get(i,i) > muSym)
                muSym = -Qsym.get(i, i);
        }

        // Set up uniformized backward transition rate matrices R and Rsym:
        R = Q.mul(1.0/mu).add(DoubleMatrix.eye(nTypes));
        Rsym = Qsym.mul(1.0/muSym).add(DoubleMatrix.eye(nTypes));
        
        // Clear cached powers of R and Rsym and steady state flag:
        RpowN.clear();
        RsymPowN.clear();
        
        RpowSteady = false;
        RsymPowSteady = false;
        
        // Power sequences initially contain R^0 = I
        RpowN.add(DoubleMatrix.eye(nTypes));
        RsymPowN.add(DoubleMatrix.eye(nTypes));
        
        RpowMax = DoubleMatrix.eye(nTypes);
        RsymPowMax = DoubleMatrix.eye(nTypes);

        dirty = false;
    }

    /**
     * @return number of demes in the migration model.
     */
    @Override
    public int getNTypes() {
        return nTypes;
    }

    /**
     * Obtain element of rate matrix for migration model for use in likelihood
     * calculation.  (May be switched to zero in BSSVS calculation.)
     *
     * @param i
     * @param j
     * @return Rate matrix element.
     */
    @Override
    public double getBackwardRate(int i, int j) {
        if (i==j)
            return 0;

        return getForwardRate(j, i);
    }

    /**
     * Retrieve rate of migration from i to j forward in time.
     *
     * @param i source deme
     * @param j dest deme
     * @return migration rate
     */
    @Override
    public double getForwardRate(int i, int j) {
        if (i==j)
            return 0.0;

        return rateMatrix.getValue(getArrayOffset(i, j));
    }

    /**
     * Obtain offset into "rate matrix" and associated flag arrays.
     * 
     * @param i
     * @param j
     * @return Offset (or -1 if i==j)
     */
    protected int getArrayOffset(int i, int j) {
        
        if (i==j)
            throw new RuntimeException("Programmer error: requested migration "
                    + "rate array offset for diagonal element of "
                    + "migration rate matrix.");
        
        if (rateMatrixIsSquare) {
            return i*nTypes+j;
        } else {
            if (symmetricRateMatrix) {
             if (j<i)
                 return i*(i-1)/2 + j;
             else
                 return j*(j-1)/2 + i;
            } else {
                if (j>i)
                    j -= 1;
                return i*(nTypes-1)+j;
            }
        }
    }

    @Override
    public double getMu(boolean symmetric) {
        updateMatrices();
        if (symmetric)
            return muSym;
        else
            return mu;
    }
    
    @Override
    public DoubleMatrix getR(boolean symmetric) {
        updateMatrices();
        if (symmetric)
            return Rsym;
        else
            return R;
    }
    
    @Override
    public DoubleMatrix getQ(boolean symmetric) {
        updateMatrices();
        if (symmetric)
            return Qsym;
        else
            return Q;
    }
    
    @Override
    public DoubleMatrix getRpowN(int n, boolean symmetric) {
        updateMatrices();
        
        List <DoubleMatrix> matPowerList;
        DoubleMatrix mat, matPowerMax;
        if (symmetric) {
            matPowerList = RsymPowN;
            mat = Rsym;
            matPowerMax = RsymPowMax;
        } else {
            matPowerList = RpowN;
            mat = R;
            matPowerMax = RpowMax;
        }
        
        if (n>=matPowerList.size()) {
                
            // Steady state of matrix iteration already reached
            if ((symmetric && RsymPowSteady) || (!symmetric && RpowSteady)) {
                //System.out.println("Assuming R SS.");
                return matPowerList.get(matPowerList.size()-1);
            }
                
            int startN = matPowerList.size();
            for (int i=startN; i<=n; i++) {
                matPowerList.add(matPowerList.get(i-1).mmul(mat));
                
                matPowerMax.maxi(matPowerList.get(i));
                    
                // Occasionally check whether matrix iteration has reached steady state
                if (i%10 == 0) {
                    double maxDiff = 0.0;
                    for (double el : matPowerList.get(i).sub(matPowerList.get(i-1)).toArray())
                        maxDiff = Math.max(maxDiff, Math.abs(el));
                        
                    if (!(maxDiff>0)) {
                        if (symmetric)
                            RsymPowSteady = true;
                        else
                            RpowSteady = true;
                        
                        return matPowerList.get(i);
                    }
                }
            }
        }
        return matPowerList.get(n);
    }
    
    /**
     * Return matrix containing upper bounds on elements from the powers
     * of R if known.  Returns a matrix of ones if steady state has not yet
     * been reached.
     * 
     * @param symmetric
     * @return Matrix of upper bounds.
     */
    @Override
    public DoubleMatrix getRpowMax(boolean symmetric) {
        
        if (symmetric) {
            if (RsymPowSteady)
                return RsymPowMax;
            else
                return DoubleMatrix.ones(nTypes, nTypes);
        } else {
            if (RpowSteady)
                return RpowMax;
            else
                return DoubleMatrix.ones(nTypes, nTypes);
        }
    }
    
    
    /**
     * Power above which R is known to be steady.
     * 
     * @param symmetric
     * @return index of first known steady element.
     */
    @Override
    public int RpowSteadyN(boolean symmetric) {
        if (symmetric) {
            if (RsymPowSteady)
                return RsymPowN.size();
            else
                return -1;
        } else {
            if (RpowSteady)
                return RpowN.size();
            else
                return -1;
        }
    }

    /*
     * CalculationNode implementations.
     */
    
    @Override
    protected boolean requiresRecalculation() {
        // we only get here if something is dirty
        dirty = true;
        return true;
    }

    @Override
    protected void restore() {
        dirty = true;
        super.restore();
    }
}
