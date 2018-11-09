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
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.CDKValencyChecker;

/**
 * Class for representing a subspectrum-substructure-correlation.
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public final class SSC implements Cloneable {
    
    private final Spectrum subspectrum;
    private final Assignment assignment;
    private final IAtomContainer substructure;
    // spherical search limit
    private final int rootAtomIndex, maxSphere;
    // stores all shifts for each single HOSE code
    private final HashMap<String, ArrayList<Double>> HOSELookupShifts;
    // stores all atom indices for each single HOSE code
    private final HashMap<String, ArrayList<Integer>> HOSELookupIndices;
    
    private final HashMap<Integer, HashMap<Integer, List<IAtom>>> atomsInHOSECodeSpheres;
    private final HashMap<Integer, HashMap<Integer, ArrayList<Integer>>> atomIndicesInHOSECodeSpheres;
    private final HashMap<Integer, HashMap<Integer, ArrayList<Integer>>> parentAtomIndicesInHOSECodeSpheres;
            
    // stores all atom indices for each occurring atom type in substructure    
    private final HashMap<String, ArrayList<Integer>> atomTypeIndices;
    // for pre-search: map of multiplicities as keys consisting 
    // of maps of shift values as keys and its atom indices
    private final HashMap<String, HashMap<Integer, ArrayList<Integer>>> presenceMultiplicities;
    // atom type of subspectrum for which atoms should have an assigned shift value
    private final String atomType;
    // min/max shift range to consider
    private int minShift, maxShift;
    // index to use in SSC library
    private int index;
    // indices of open-sphere (unsaturated) atoms of substructure
    private final ArrayList<Integer> unsaturatedAtomIndices;

    /**
     *
     * @param subspectrum
     * @param assignment
     * @param substructure
     * @param rootAtomIndex
     * @param maxSphere
     * @throws CDKException
     */
    public SSC(final Spectrum subspectrum, final Assignment assignment, final IAtomContainer substructure, final int rootAtomIndex, final int maxSphere) throws CDKException {
        this.subspectrum = subspectrum;
        this.assignment = assignment;
        this.substructure = substructure;
        this.rootAtomIndex = rootAtomIndex;
        this.maxSphere = maxSphere;
        this.atomTypeIndices = Utils.getAtomTypeIndices(this.substructure);
        this.atomType = Utils.getAtomTypeFromSpectrum(this.subspectrum, 0);
        this.HOSELookupShifts = new HashMap<>();
        this.HOSELookupIndices = new HashMap<>();
        this.atomsInHOSECodeSpheres = new HashMap<>();
        this.atomIndicesInHOSECodeSpheres = new HashMap<>();
        this.parentAtomIndicesInHOSECodeSpheres = new HashMap<>();
        this.minShift = 0;
        this.maxShift = 220;
        this.index = -1;
        this.unsaturatedAtomIndices = new ArrayList<>();
        this.updateUnsaturatedAtomIndices();
        this.presenceMultiplicities = new HashMap<>();
        this.updatePresenceMultiplicities();        
    }   
    
    @Override
    public SSC clone() throws CloneNotSupportedException{
      return (SSC) super.clone();
    }
    
    public void updateUnsaturatedAtomIndices() throws CDKException {
        final IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
        this.unsaturatedAtomIndices.clear();
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
     */
    public void updatePresenceMultiplicities(){
        final String[] mults = new String[]{"S", "D", "T", "Q"};   
        // init
        for (final String mult : mults) {            
            this.presenceMultiplicities.put(mult, new HashMap<>());
            for (int i = this.minShift; i <= this.maxShift; i++) {
                this.presenceMultiplicities.get(mult).put(i, new ArrayList<>());
            }
        }
        if(this.getAtomTypeIndices().get(this.atomType) == null){
            return;
        }
        // pre-search and settings
        IAtom atom;
        int shift;
        for (final int i : this.getAtomTypeIndices().get(this.atomType)) { // for all atoms of that atom type
            atom = this.substructure.getAtom(i);                                            
            if((atom.getProperty(Utils.getNMRShiftConstant(this.atomType)) != null) // that atom has a set shift value         
                    && (atom.getImplicitHydrogenCount() != null)                    // hydrogen count must not be null
                    && (Utils.getMultiplicityFromHydrogenCount(atom.getImplicitHydrogenCount()) != null)){  // multiplicity obtained by attached hydrogen count must not be null                     
                shift = ((Double) atom.getProperty(Utils.getNMRShiftConstant(this.atomType))).intValue();
                if(Utils.checkMinMaxValue(this.minShift, this.maxShift, shift)){
                    this.presenceMultiplicities.get(Utils.getMultiplicityFromHydrogenCount(atom.getImplicitHydrogenCount())).get(shift).add(i);
                }                
            }
        }        
    }        
    
    public boolean setHOSECode(final int atomIndexInSubstructure, final String HOSECode){
        if(!Utils.checkIndexInAtomContainer(this.substructure, atomIndexInSubstructure) || (HOSECode == null)){
            return false;
        }
        if (!this.HOSELookupShifts.containsKey(HOSECode)) {
            this.HOSELookupShifts.put(HOSECode, new ArrayList<>());
            this.HOSELookupIndices.put(HOSECode, new ArrayList<>());
        }
        if (this.substructure.getAtom(atomIndexInSubstructure).getProperty(Utils.getNMRShiftConstant(this.atomType)) != null) {
            this.HOSELookupShifts.get(HOSECode).add(this.substructure.getAtom(atomIndexInSubstructure).getProperty(Utils.getNMRShiftConstant(this.atomType)));
        }
        this.HOSELookupIndices.get(HOSECode).add(atomIndexInSubstructure);
        
        return true;
    }
    
    public String getHOSECode(final int atomIndexInSubstructure){
        if(!Utils.checkIndexInAtomContainer(this.substructure, atomIndexInSubstructure)){
            return null;
        }
        for (final String HOSECode : this.getHOSELookupIndices().keySet()) {
            if (this.getHOSELookupIndices().get(HOSECode).contains(atomIndexInSubstructure)) {
                return HOSECode;
            }
        }
        
        return null;
    }
    
    public boolean setAtomsInHOSECodeSpheres(final int atomIndexInSubstructure, final List<IAtom> atomsInSpheres, final int sphere) {
        if (!Utils.checkIndexInAtomContainer(this.substructure, atomIndexInSubstructure)
                || sphere > this.getMaxSphere()) {
            return false;
        }
        if (!this.atomsInHOSECodeSpheres.containsKey(atomIndexInSubstructure)) {
            this.atomsInHOSECodeSpheres.put(atomIndexInSubstructure, new HashMap<>());
        }
        this.atomsInHOSECodeSpheres.get(atomIndexInSubstructure).put(sphere, atomsInSpheres);
        this.updateAtomIndicesInHOSECodeSpheres(atomIndexInSubstructure, sphere);       
        
        return true;
    }
    
    public List<IAtom> getAtomsInHOSECodeSpheres(final int atomIndexInSubstructure, final int sphere) {
        if (!Utils.checkIndexInAtomContainer(this.substructure, atomIndexInSubstructure)
                || !this.atomsInHOSECodeSpheres.containsKey(atomIndexInSubstructure)
                || !this.atomsInHOSECodeSpheres.get(atomIndexInSubstructure).containsKey(sphere)) {
            return null;
        }               

        return this.atomsInHOSECodeSpheres.get(atomIndexInSubstructure).get(sphere);
    }
    
    public boolean updateAtomIndicesInHOSECodeSpheres(final int atomIndexInSubstructure, final int sphere) {
        // atom index and sphere have to exist in atom HOSE code sphere list      
        if(!this.atomsInHOSECodeSpheres.containsKey(atomIndexInSubstructure) || !this.atomsInHOSECodeSpheres.get(atomIndexInSubstructure).containsKey(sphere)){
            return false;
        }
        // create new and delete old instances if needed
        if(!this.atomIndicesInHOSECodeSpheres.containsKey(atomIndexInSubstructure)){
            this.atomIndicesInHOSECodeSpheres.put(atomIndexInSubstructure, new HashMap<>());
        }
        this.atomIndicesInHOSECodeSpheres.get(atomIndexInSubstructure).put(sphere, new ArrayList<>());
        // search for each atom of origin structure in this SSC substructure and set matched indices
        final List<IAtom> atomsInSphere = this.atomsInHOSECodeSpheres.get(atomIndexInSubstructure).get(sphere);
        for (int a = 0; a < atomsInSphere.size(); a++) { 
            if(atomsInSphere.get(a) == null){
                this.atomIndicesInHOSECodeSpheres.get(atomIndexInSubstructure).get(sphere).add(null);                
            } else {
                this.atomIndicesInHOSECodeSpheres.get(atomIndexInSubstructure).get(sphere).add(this.getSubstructure().indexOf(atomsInSphere.get(a)));
            }            
        }  
        this.updateParentAtomIndicesInHOSECodeSpheres(atomIndexInSubstructure, sphere);
        
        return true;
    }
    
    public ArrayList<Integer> getAtomIndicesInHOSECodeSpheres(final int atomIndexInSubstructure, final int sphere) {

        if (!Utils.checkIndexInAtomContainer(this.substructure, atomIndexInSubstructure)
                || !this.atomIndicesInHOSECodeSpheres.containsKey(atomIndexInSubstructure)
                || !this.atomIndicesInHOSECodeSpheres.get(atomIndexInSubstructure).containsKey(sphere)) {
            return new ArrayList<>();
        }               

        return this.atomIndicesInHOSECodeSpheres.get(atomIndexInSubstructure).get(sphere);
    }
    
    private boolean updateParentAtomIndicesInHOSECodeSpheres(final int atomIndexInSubstructure, final int sphere) {
        // atom index and sphere have to exist in atom HOSE code sphere list      
        if ((sphere <= 0) || !this.atomIndicesInHOSECodeSpheres.containsKey(atomIndexInSubstructure) || !this.atomsInHOSECodeSpheres.get(atomIndexInSubstructure).containsKey(sphere)) {
            return false;
        }
        // create new and delete old instances if needed
        if (!this.parentAtomIndicesInHOSECodeSpheres.containsKey(atomIndexInSubstructure)) {
            this.parentAtomIndicesInHOSECodeSpheres.put(atomIndexInSubstructure, new HashMap<>());
        }
        this.parentAtomIndicesInHOSECodeSpheres.get(atomIndexInSubstructure).put(sphere, new ArrayList<>());
        // sets for parents of atom in a sphere        
        this.parentAtomIndicesInHOSECodeSpheres.get(atomIndexInSubstructure).get(sphere).addAll(this.findParentAtomIndicesInHOSECodeSphere(atomIndexInSubstructure, sphere));

        return true;
    }

    public ArrayList<Integer> getParentAtomIndicesInHOSECodeSpheres(final int atomIndexInSubstructure, final int sphere) {

        if (!Utils.checkIndexInAtomContainer(this.substructure, atomIndexInSubstructure)
                || !this.parentAtomIndicesInHOSECodeSpheres.containsKey(atomIndexInSubstructure)
                || !this.parentAtomIndicesInHOSECodeSpheres.get(atomIndexInSubstructure).containsKey(sphere)) {
            return new ArrayList<>();
        }

        return this.parentAtomIndicesInHOSECodeSpheres.get(atomIndexInSubstructure).get(sphere);
    }
    
    private ArrayList<Integer> findParentAtomIndicesInHOSECodeSphere(final int atomIndexInSubstructure, final int sphere){
        // atom index and sphere have to exist in atom HOSE code sphere list      
        if ((sphere <= 0) 
                || !this.parentAtomIndicesInHOSECodeSpheres.containsKey(atomIndexInSubstructure) 
                || !this.parentAtomIndicesInHOSECodeSpheres.get(atomIndexInSubstructure).containsKey(sphere)) {
            return new ArrayList<>();
        }
                        
        final ArrayList<Integer> parentAtomIndicesInSphere = new ArrayList<>();
        for (final Integer atomIndexInSphere: this.atomIndicesInHOSECodeSpheres.get(atomIndexInSubstructure).get(sphere)) {
            if (atomIndexInSphere == null) {
                parentAtomIndicesInSphere.add(null);
            } else if(atomIndexInSphere == 0){
                parentAtomIndicesInSphere.add(-3);
            } else if(atomIndexInSphere == -1){
                parentAtomIndicesInSphere.add(-1);
            } else if(sphere == 1){
                parentAtomIndicesInSphere.add(atomIndexInSubstructure);
            } else {
                int atomIndexInSphereFrequency = Collections.frequency(this.atomIndicesInHOSECodeSpheres.get(atomIndexInSubstructure).get(sphere), atomIndexInSphere);
                if(atomIndexInSphereFrequency == 1){ // node in sphere exists only once, so unique assignment to a parent possible
                    int atomIndexInPrevSphereToAdd = -1;
                    for (final Integer atomIndexInPrevSphere : this.atomIndicesInHOSECodeSpheres.get(atomIndexInSubstructure).get(sphere-1)) {
                        if((atomIndexInPrevSphere != null) && (atomIndexInPrevSphere != -1)){
                            if (this.getSubstructure().getAtom(atomIndexInSphere).getBond(this.getSubstructure().getAtom(atomIndexInPrevSphere)) != null) {
                                atomIndexInPrevSphereToAdd = atomIndexInPrevSphere;                                
                                break;
                            }
                        }                            
                    }
                    parentAtomIndicesInSphere.add(atomIndexInPrevSphereToAdd);
                } else { // node exists multiple times in sphere
                    ArrayList<Integer> parentsInPrevSphere = new ArrayList<>();
                    for (final Integer atomIndexInPrevSphere : this.atomIndicesInHOSECodeSpheres.get(atomIndexInSubstructure).get(sphere-1)) {
                        if ((atomIndexInPrevSphere != null) && (atomIndexInPrevSphere != -1)) {
                            if (this.getSubstructure().getAtom(atomIndexInSphere).getBond(this.getSubstructure().getAtom(atomIndexInPrevSphere)) != null) {
                                parentsInPrevSphere.add(atomIndexInPrevSphere);
                            }
                        }
                    }
                    if(parentsInPrevSphere.size() == 1){ // node exists multiple times but all of them have the same parent node in previous sphere
                        parentAtomIndicesInSphere.add(parentsInPrevSphere.get(0));
                    } else { 
                        // node in sphere multiple times and parents in previous sphere multiple times too 
                        // so no unique assignment between node and parent node possible
                        // -> still open
                        parentAtomIndicesInSphere.add(-2);
                    }
                }
            }
        }
        
        return parentAtomIndicesInSphere;
    }
    
    /**
     * Returns the root atom index of this substructure.
     *
     * @return
     */
    public int getRootAtomIndex(){
        return this.rootAtomIndex;
    }
    
    public Boolean isUnsaturated(final int atomIndex){
        if(!Utils.checkIndexInAtomContainer(this.substructure, atomIndex)){
            return null;
        }                
        
        return this.getUnsaturatedAtomIndices().contains(atomIndex);
    }
    
    /**
     * Returns the used maximum number of spheres for building the SSC's 
     * HOSE codes.
     *
     * @return
     */
    public int getMaxSphere(){
        return this.maxSphere;
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
    
    public int getBondCount(){
        return this.getSubstructure().getBondCount();
    }
    
    public HashMap<String, ArrayList<Double>> getHOSELookupShifts(){
        return this.HOSELookupShifts;
    }
    
    public HashMap<String, ArrayList<Integer>> getHOSELookupIndices(){
        return this.HOSELookupIndices;
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
     * Despite intensities are given, they are still not considered here.
     *
     * @param multiplicity Multiplicity as multiplet string, e.g. "S" (Singlet)
     * @param shift Shift value [ppm]
     * @param intensity Intensity value
     * @param tol Tolerance value [ppm] used during shift matching
     * @return
     */
    public ArrayList<Integer> findMatches(final String multiplicity, final double shift, final double intensity, final double tol) {
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
     * Despite intensities are given, they are still not considered here.
     *
     * @param querySpectrum Query spectrum
     * @param tol Tolerance value [ppm] used during shift matching
     * @return
     */
    public Assignment findMatches(final Spectrum querySpectrum, final double tol) {
        final Assignment matchAssignments = new Assignment(querySpectrum);
        // wrong nucleus in query spectrum or this SSC contains no atoms of requested atom type
        if (!Utils.getElementIdentifier(querySpectrum.getNuclei()[0]).equals(this.atomType) 
                || (this.getAtomTypeIndices().get(this.atomType) == null)) {
            return matchAssignments;
        }        
        Signal signal;
        for (int i = 0; i < querySpectrum.getSignalCount(); i++) {
            signal = querySpectrum.getSignal(i);
            matchAssignments.setAssignment(0, i, this.getClosestMatch(this.findMatches(signal.getMultiplicity(), signal.getShift(0), signal.getIntensity(), tol), signal.getShift(0), tol));
        }
        
        // too many atoms of that atom type in SSC, but check for possible symmetry
        if (this.getAtomTypeIndices().get(this.atomType).size() > querySpectrum.getSignalCount()
                || matchAssignments.getSetAssignmentsCount(0) > this.getAtomTypeIndices().get(this.atomType).size()) {
//            System.out.println("Check for multiple assignments/symmetry!!!");
            return new Assignment(querySpectrum);
        }

        return matchAssignments;
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
     * @param matchAssignments 
     * @return
     *
     * @see #findMatches(casekit.NMR.model.Spectrum, double) 
     */
    public Double[] getDeviations(final Spectrum querySpectrum, final Assignment matchAssignments) {
        final Double[] deviations = new Double[matchAssignments.getAssignmentsCount()];

        double shiftMatchIndex;
        for (int i = 0; i < matchAssignments.getAssignmentsCount(); i++) {
            if (matchAssignments.getAtomIndex(0, i) == -1) {
                deviations[i] = null;
            } else {
                shiftMatchIndex = this.substructure.getAtom(matchAssignments.getAtomIndex(0, i)).getProperty(Utils.getNMRShiftConstant(this.atomType));
                deviations[i] = Math.abs(shiftMatchIndex - querySpectrum.getShift(i, 0));
            }
        }

        return deviations;
    }
    
    /**
     * Returns deviatons between matched shifts in SSC and query query spectrum.
     * The matching procedure is already included here.
     *
     * @param querySpectrum
     * @param tol
     * @return
     *
     * @see #findMatches(casekit.NMR.model.Spectrum, double) 
     */
    public Double[] getDeviations(final Spectrum querySpectrum, final double tol) {
        final Double[] deviations = new Double[querySpectrum.getSignalCount()];
        final Assignment matchAssignments = this.findMatches(querySpectrum, tol);

        double shiftMatchIndex;
        for (int i = 0; i < querySpectrum.getSignalCount(); i++) {
            if (matchAssignments.getAtomIndex(0, i) == -1) {
                deviations[i] = null;
            } else {
                shiftMatchIndex = this.substructure.getAtom(matchAssignments.getAtomIndex(0, i)).getProperty(Utils.getNMRShiftConstant(this.atomType));
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
     * @param tol Tolerance value [ppm] used during shift matching
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
