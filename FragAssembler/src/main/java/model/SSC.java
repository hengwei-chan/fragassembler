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
    private int minShift, maxShift;
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
        this.unsaturatedAtomIndices = new ArrayList<>();
        this.initUnsaturatedAtomIndices();
        this.presenceMultiplicities = new HashMap<>();
        this.initPresenceMultiplicities();        
    }
    
    private void initUnsaturatedAtomIndices() throws CDKException{
        IAtom atom;
        final IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
        for (int i = 0; i < this.substructure.getAtomCount(); i++) {
            atom = this.substructure.getAtom(i);
            // set the indices of unsaturated atoms in substructure
            if (!CDKValencyChecker.getInstance(builder).isSaturated(atom, this.substructure)) {
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
        // pre-search and settings
        IAtom atom;
        int shift;
        for (final int i : this.getAtomTypeIndices().get(this.atomType)) {
            atom = this.substructure.getAtom(i);                                
            shift = ((Double) atom.getProperty(Utils.getNMRShiftConstant(this.atomType))).intValue();
            if((atom.getProperty(Utils.getNMRShiftConstant(this.atomType)) != null) 
                    && Utils.checkMinMaxValue(this.minShift, this.maxShift, shift)
                    && (atom.getImplicitHydrogenCount() != null)
                    && (Utils.getMultiplicityFromHydrogenCount(atom.getImplicitHydrogenCount()) != null)){                     
                this.presenceMultiplicities.get(Utils.getMultiplicityFromHydrogenCount(atom.getImplicitHydrogenCount())).get(shift).add(i);
            }
        }        
    }    
    
    // stores the shifts and atom indices for each HOSE code
    private void initHOSECodes() throws CDKException{
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
    
    public ArrayList<Integer> getUnsaturatedAtomIndices(){
        return this.unsaturatedAtomIndices;
    }
    
    public ArrayList<Integer> findMatches(final String multiplicity, final double shift, final int tol){
                        
        if(this.presenceMultiplicities.get(multiplicity) == null){
            return new ArrayList<>();
        }
        final int ShiftInteger = ((Double) shift).intValue();
        final HashSet<Integer> hs = new HashSet<>();         
        for (int i = ShiftInteger - tol ; i <= ShiftInteger + tol; i++) {
            if(Utils.checkMinMaxValue(this.minShift, this.maxShift, i)){
                hs.addAll(this.presenceMultiplicities.get(multiplicity).get(i));
            }            
        }
        
        return new ArrayList<>(hs);
    }
}
