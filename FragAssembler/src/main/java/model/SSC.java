/*
 * The MIT License
 *
 * Copyright 2018 Michael Wenk [https://github.com/michaelwenk].
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package model;

import casekit.NMR.Utils;
import casekit.NMR.model.Assignment;
import casekit.NMR.model.Signal;
import casekit.NMR.model.Spectrum;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.CDKValencyChecker;
import org.openscience.cdk.tools.HOSECodeGenerator;

/**
 * Class for representing a subspectrum-substructure-correlation.
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class SSC {
    
    private final Spectrum subspectrum;
    private final Assignment assignment;
    private final IAtomContainer substructure;
    // spherical search limit
    private final int maxNoOfSpheres;
    // stores all shifts for each single HOSE code
    private final HashMap<String, ArrayList<Double>> hoseLookupShifts;
    // stores all atom indices for each single HOSE code
    private final HashMap<String, ArrayList<Integer>> hoseLookupIndices;
    // stores all atom indices for each occurring atom type in substructure
    private final HashMap<String, ArrayList<Integer>> atomTypeIndices;
    // for pre-search: map of multiplicities (no. of attached protons) consisting 
    // of map of shift values and its atom indices
    private final HashMap<String, HashMap<Integer, ArrayList<Integer>>> presenceMultiplicities;
    // atom type of subspectrum for which atoms should have an assigned shift value
    private final String atomType;
    // min/max shift range to consider 
    private int minShift, maxShift, index;
    // indices of open-sphere (unsaturated) atoms of substructure
    private final ArrayList<Integer> unsaturatedAtomIndices;

    public SSC(final Spectrum subspectrum, final Assignment assignment, final IAtomContainer substructure, final int maxNoOfSpheres) throws CDKException {
        this.subspectrum = subspectrum;
        this.assignment = assignment;
        this.substructure = substructure; 
        this.maxNoOfSpheres = maxNoOfSpheres;
        this.atomTypeIndices = Utils.getAtomTypeIndices(this.substructure);
        this.atomType = Utils.getAtomTypeFromSpectrum(this.subspectrum, 0);
        this.hoseLookupShifts = new HashMap<>();
        this.hoseLookupIndices = new HashMap<>();
        this.initHOSECodes();
        this.minShift = 0;
        this.maxShift = 220;
        this.index = -1;
        this.unsaturatedAtomIndices = new ArrayList<>();
        this.initUnsaturatedAtomIndices();
        this.presenceMultiplicities = new HashMap<>();
        this.initPresenceMultiplicities();        
    }
    
    private void initUnsaturatedAtomIndices() throws CDKException{
        final IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
        for (int i = 0; i < this.substructure.getAtomCount(); i++) {
            // set the indices of unsaturated atoms in substructure
            if (!CDKValencyChecker.getInstance(builder).isSaturated(this.substructure.getAtom(i), this.substructure)) {
                this.unsaturatedAtomIndices.add(i);
            }            
        }
    }
    
    
    /**
     * Specified for carbons only -> not generic.
     *
     * @return
     */
    private void initPresenceMultiplicities(){
        final String[] mults = new String[]{"S", "D", "T", "Q"};   
        // init
        for (final String mult : mults) {            
            this.presenceMultiplicities.put(mult, new HashMap<>());
            for (int i = this.minShift; i < this.maxShift + 1; i++) {
                this.presenceMultiplicities.get(mult).put(i, new ArrayList<>());
            }
        }
        if(this.getAtomTypeIndices().get(this.atomType) == null){
            return;
        }
        // pre-search and settings
        IAtom atom;
        int shift;
        for (final int i : this.getAtomTypeIndices().get(this.atomType)) {
            atom = this.substructure.getAtom(i);                                            
            if((atom.getProperty(Utils.getNMRShiftConstant(this.atomType)) != null)                     
                    && (atom.getImplicitHydrogenCount() != null)
                    && (Utils.getMultiplicityFromHydrogenCount(atom.getImplicitHydrogenCount()) != null)){                     
                shift = ((Double) atom.getProperty(Utils.getNMRShiftConstant(this.atomType))).intValue();
                if(Utils.checkMinMaxValue(this.minShift, this.maxShift, shift)){
                    this.presenceMultiplicities.get(Utils.getMultiplicityFromHydrogenCount(atom.getImplicitHydrogenCount())).get(shift).add(i);
                }                
            }
        }        
    }    
    

    /**
     * Stores the shifts and atom indices for each HOSE code.
     *
     * @throws CDKException
     */
    private void initHOSECodes() throws CDKException{
        if(this.getAtomTypeIndices().get(this.atomType) == null){
            return;
        }
        final HOSECodeGenerator hcg = new HOSECodeGenerator();
        String hose;
        for (final int i : this.getAtomTypeIndices().get(this.atomType)) {
            hose = hcg.getHOSECode(this.substructure, this.substructure.getAtom(i), this.maxNoOfSpheres);
            if (!this.hoseLookupShifts.containsKey(hose)) {
                this.hoseLookupShifts.put(hose, new ArrayList<>());
                this.hoseLookupIndices.put(hose, new ArrayList<>());
            }
            if (this.substructure.getAtom(i).getProperty(Utils.getNMRShiftConstant(this.atomType)) != null) {
                this.hoseLookupShifts.get(hose).add(this.substructure.getAtom(i).getProperty(Utils.getNMRShiftConstant(this.atomType)));
                this.hoseLookupIndices.get(hose).add(i);
            }
        }
    }
    
    /**
     * Sets the lower bound for shift matching. This bound is set to 0 by default.
     *
     * @param minShift
     */
    public void setMinShift(final int minShift){
        this.minShift = minShift;
    }
    
    public int getMinShift(){
        return this.minShift;
    }
    
    /**
     * Sets the upper bound for shift matching. This bound is set to 220 by default.
     *
     * @param maxShift
     */
    public void setMaxShift(final int maxShift){
        this.maxShift = maxShift;
    }
    
    public int getMaxShift(){
        return this.maxShift;
    }
    
    public void setIndex(final int index){
        this.index = index;
    }
    
    public int getIndex(){
        return this.index;
    }
    
    public Spectrum getSubspectrum(){
        return this.subspectrum;
    }
    
    public Assignment getAssignments(){
        return this.assignment;
    }
    
    public IAtomContainer getSubstructure(){
        return this.substructure;
    }
    
    public int getAtomCount(){
        return this.getSubstructure().getAtomCount();
    }
    
    public HashMap<String, ArrayList<Double>> getHOSELookupShifts(){
        return this.hoseLookupShifts;
    }
    
    public HashMap<String, ArrayList<Integer>> getHOSELookupIndices(){
        return this.hoseLookupIndices;
    }
    
    public HashMap<String, ArrayList<Integer>> getAtomTypeIndices(){
        return this.atomTypeIndices;
    }
    
    public HashMap<String, HashMap<Integer, ArrayList<Integer>>> getPresenceMultiplicities(){
        return this.presenceMultiplicities;
    }
    
    /**
     * Returns the indices of open-sphere atoms in substructure.
     *
     * @return
     */
    public ArrayList<Integer> getUnsaturatedAtomIndices(){
        return this.unsaturatedAtomIndices;
    }
    
    /**
     * Returns the full list of matched shifts as atom indices between a SSC 
     * and a query spectrum.
     * Intensities are still not considered here.
     *
     * @param multiplicity
     * @param shift
     * @param tol
     * @return
     */
    public ArrayList<Integer> findMatches(final String multiplicity, final double shift, final double tol) {

        if (this.presenceMultiplicities.get(multiplicity) == null) {
            return new ArrayList<>();
        }
        final HashSet<Integer> hs = new HashSet<>();
        for (int i = (int) (shift - tol); i <= (int) (shift + tol); i++) {
            if (Utils.checkMinMaxValue(this.minShift, this.maxShift, i)) {
                hs.addAll(this.presenceMultiplicities.get(multiplicity).get(i));
            }
        }

        return new ArrayList<>(hs);
    }

    /**
     * Returns the closest shift matches between a SSC and a
     * query spectrum as an Assignment object.
     * Intensities are still not considered here.
     *
     * @param querySpectrum Query spectrum
     * @param tol Tolerance value [ppm] while shift matching
     * @return
     */
    public Assignment findMatches(final Spectrum querySpectrum, final double tol) {
        final Assignment matchAssignment = new Assignment(querySpectrum);
        if (!Utils.getElementIdentifier(querySpectrum.getNuclei()[0]).equals(this.atomType)) {
            System.err.println("Wrong nucleus in query spectrum!!!");
            return matchAssignment;
        }
        Signal signal;
        for (int i = 0; i < querySpectrum.getSignalCount(); i++) {
            signal = querySpectrum.getSignal(i);
            matchAssignment.setAssignment(0, i, this.getClosestMatch(this.findMatches(signal.getMultiplicity(), signal.getShift(0), tol), signal.getShift(0), tol));
        }

        return matchAssignment;
    }

    private int getClosestMatch(final ArrayList<Integer> matchAtomIndices, final double queryShift, final double tol) {
        int closestMatchIndex = -1;
        double diff = tol, shiftMatchIndex;
        for (final int matchIndex : matchAtomIndices) {
            shiftMatchIndex = this.substructure.getAtom(matchIndex).getProperty(Utils.getNMRShiftConstant(this.atomType));
            if (Math.abs(shiftMatchIndex - queryShift) < diff) {
                diff = Math.abs(shiftMatchIndex - queryShift);
                closestMatchIndex = matchIndex;
            }
        }

        return closestMatchIndex;
    }

    /**
     * Returns deviatons between matched shifts in SSC and query query spectrum.
     * The matching procedure is already included here.
     *
     * @param querySpectrum
     * @param tol
     * @return
     *
     */
    public Double[] getDeviations(final Spectrum querySpectrum, final double tol) {
        final Double[] deviations = new Double[querySpectrum.getSignalCount()];
        final Assignment matches = this.findMatches(querySpectrum, tol);

        double shiftMatchIndex;
        for (int i = 0; i < querySpectrum.getSignalCount(); i++) {
            if (matches.getAssignment(0, i) == -1) {
                deviations[i] = null;
            } else {
                shiftMatchIndex = this.substructure.getAtom(matches.getAssignment(0, i)).getProperty(Utils.getNMRShiftConstant(this.atomType));
                deviations[i] = Math.abs(shiftMatchIndex - querySpectrum.getShift(i, 0));
            }
        }

        return deviations;
    }

    /**
     * Returns the average of all deviations of matched shifts between a SSC 
     * and a query spectrum. 
     * The calculation of deviations is already included here.
     *
     * @param querySpectrum Query spectrum
     * @param tol Tolerance value [ppm] while shift matching
     * @param minOverlap Minimum overlap threshold
     * @return
     * 
     * @see model.SSC#getDeviations(casekit.NMR.model.Spectrum, double) 
     */
    public Double getMatchFactor(final Spectrum querySpectrum, final double tol, final int minOverlap) {
        return this.getMatchFactor(this.getDeviations(querySpectrum, tol), minOverlap);
    }
    
    /**
     * Returns the average of all deviations of matched shifts between a SSC 
     * and a query spectrum. 
     * If the minimum overlap threshold is not reached, a null value will be 
     * returned.
     *
     * @param deviations Deviations between a SSC subspectrum and query spectrum.
     * @param minOverlap Minimum overlap threshold
     * @return
     * 
     * @see model.SSC#getDeviations(casekit.NMR.model.Spectrum, double) 
     */
    public Double getMatchFactor(final Double[] deviations, final int minOverlap) {
        int hitCounter = 0;
        for (final Double deviation : deviations) {
            if(deviation != null){
                hitCounter++;
            }
        }
        
        return (hitCounter >= minOverlap) ? Utils.getMean(deviations) : null;
    }
}
