/*
 * The MIT License
 *
 * Copyright (c) 2019 Michael Wenk [https://github.com/michaelwenk]
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
package model;

import analysis.MultiplicitySectionsBuilder;
import casekit.NMR.Utils;
import casekit.NMR.model.Assignment;
import casekit.NMR.model.Signal;
import casekit.NMR.model.Spectrum;
import hose.HOSECodeBuilder;
import hose.model.ConnectionTree;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Class for representing a subspectrum-substructure-correlation.
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public final class SSC {
    
    private final Spectrum subspectrum;
    private final Assignment assignment;
    private final IAtomContainer substructure;
    // spherical search limit
    private final int rootAtomIndex;
    private final int maxSphere;
    // stores all shifts for each single HOSE code
    private final HashMap<String, ArrayList<Double>> HOSECodeLookupShifts;
    // stores all atom indices for each single HOSE code
    private final HashMap<String, ArrayList<Integer>> HOSECodeLookupIndices;
    // connection tree holding the spherical information about the substructure
    private final HashMap<Integer, ConnectionTree> connectionTrees;
    // stores all atom indices for each occurring atom type in substructure    
    private HashMap<String, ArrayList<Integer>> atomTypeIndices;
    // for pre-search: map of multiplicities as keys consisting 
    // of section indices with occurrency in this SSC
    private HashMap<String, ArrayList<Integer>> multiplicitySections;
    private final MultiplicitySectionsBuilder multiplicitySectionsBuilder;
    // index to use in SSC library
    private long index;
    // indices of open-sphere (unsaturated) atoms of substructure
    private final ArrayList<Integer> unsaturatedAtomIndices;    
    public final static int MIN_LIMIT = -20, MAX_LIMIT = 260, STEP_SIZE = 5, STEPS = (MAX_LIMIT - MIN_LIMIT) / STEP_SIZE; // ppm range from -20 to 260 in 5 ppm steps

    /**
     *
     * @param subspectrum
     * @param assignment
     * @param substructure  
     * @param rootAtomIndex
     * @param maxSphere
     * @throws CDKException
     * @throws java.lang.CloneNotSupportedException
     */
    public SSC(final Spectrum subspectrum, final Assignment assignment, final IAtomContainer substructure, final int rootAtomIndex, final int maxSphere) throws Exception {
        this.subspectrum = subspectrum.getClone();
        this.assignment = assignment.clone();
        this.substructure = substructure.clone();     
        this.rootAtomIndex = rootAtomIndex;
        this.maxSphere = maxSphere;
        this.HOSECodeLookupShifts = new HashMap<>();
        this.HOSECodeLookupIndices = new HashMap<>();
        this.connectionTrees = new HashMap<>();
        this.index = -1;
        this.unsaturatedAtomIndices = new ArrayList<>();
        this.multiplicitySections = new HashMap<>();      
        this.multiplicitySectionsBuilder = new MultiplicitySectionsBuilder();
        this.update();
    }   
    
    /**
     * Returns a full clone of that SSC, with one exception: The index of the 
     * SSC clone is set to default value (-1).
     *
     * @return
     * @throws CDKException
     * @throws java.lang.CloneNotSupportedException
     */
    public SSC getClone() throws Exception {
      return new SSC(this.subspectrum, this.assignment, this.substructure, this.rootAtomIndex, this.maxSphere);
    }
    
    /**
     * Updates all features of that SSC object.
     *
     * @throws CDKException
     * @see #updateAtomTypeIndices() 
     * @see #updateUnsaturatedAtomIndices() 
     * @see #updateMultiplicitySections() 
     * @see #updateHOSECodes() 
     */
    public void update() throws CDKException {
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(this.substructure);
        this.updateAtomTypeIndices();
        this.updateUnsaturatedAtomIndices();
        this.updateHOSECodes();
        this.updateMultiplicitySections(); 
    }
    
    @Override
    public String toString(){
        String output = "\natom types and indices:";
        for (final String atomType : this.atomTypeIndices.keySet()) {
            output += "\n-> " + atomType + ": " +  this.atomTypeIndices.get(atomType);
        }
        output += "\nunsaturated atoms (indices): \n-> " + this.unsaturatedAtomIndices;        
        output += "\nmultiplicities and shift section:";
        for (final String mult : this.multiplicitySections.keySet()) {
            output += "\n-> " + mult + ":\t" + this.multiplicitySections.get(mult);
            
        }        
        output += "\nmax sphere: " + this.maxSphere;
        output += "\natoms (indices), HOSE codes and connection trees: ";
        for (int i = 0; i < this.getAtomCount(); i++) {
            output += "\n-> " + i + ":\t" + this.getHOSECode(i) + " -> " + this.getConnectionTree(i).toString();
        }
        
        return output;
    }
    
    /**
     * Updates all connection trees in this SSC. 
     *
     * @return 
     * 
     * @throws CDKException
     *
     * @deprecated
     */
    private boolean updateConnectionTrees() throws CDKException {
        for (int i = 0; i < this.getAtomCount(); i++) {
            if(!this.updateConnectionTree(i)){
                return false;
            }            
        }
        
        return true;
    }
    
    private boolean updateConnectionTree(final int atomIndexInSubstructure) throws CDKException {
        if (!Utils.checkIndexInAtomContainer(this.substructure, atomIndexInSubstructure)) {
            return false;
        }
        this.connectionTrees.put(atomIndexInSubstructure, HOSECodeBuilder.buildConnectionTree(this.substructure, atomIndexInSubstructure, this.maxSphere));     
        
        return true;
    }    
    
    public void updateUnsaturatedAtomIndices() throws CDKException  {
        this.unsaturatedAtomIndices.clear();
        for (int i = 0; i < this.substructure.getAtomCount(); i++) {
            // set the indices of unsaturated atoms in substructure
            if (!Utils.isSaturated(this.substructure, i)) {
                this.unsaturatedAtomIndices.add(i);
            }            
        }
    }   
    
    /**
     * Returns the indices of open-sphere atoms in substructure.
     *
     * @return
     */
    public ArrayList<Integer> getUnsaturatedAtomIndices(){
        return this.unsaturatedAtomIndices;
    }
    
    public void updateAtomTypeIndices(){
        this.atomTypeIndices = Utils.getAtomTypeIndices(this.substructure);
    } 
    
    /**
     * Specified for carbons only -> still open.
     *
     * @throws org.openscience.cdk.exception.CDKException
     */
    public void updateMultiplicitySections() throws CDKException {            
        this.multiplicitySections = this.multiplicitySectionsBuilder.getMultiplicitySections(this.subspectrum);
    }    

    public boolean updateHOSECodes() throws CDKException {
        for (int i = 0; i < this.getAtomCount(); i++) {
            if(!updateHOSECode(i)){
                return false;
            }
        }
        
        return true;
    }
    
    public boolean updateHOSECode(final int atomIndexInSubstructure) throws CDKException {
        if(!Utils.checkIndexInAtomContainer(this.substructure, atomIndexInSubstructure)){
            return false;
        }
        if(!this.updateConnectionTree(atomIndexInSubstructure)){
            return false;
        }        
        final String HOSECodePrev = this.getHOSECode(atomIndexInSubstructure);
        final String HOSECode = HOSECodeBuilder.buildHOSECode(this.connectionTrees.get(atomIndexInSubstructure), false);
        if((HOSECodePrev != null) && !HOSECode.equals(HOSECodePrev)){
            // remove old HOSE code entries of that atom
            final int positionInLookupLists = this.HOSECodeLookupIndices.get(HOSECodePrev).indexOf(atomIndexInSubstructure);
            this.HOSECodeLookupIndices.get(HOSECodePrev).remove(positionInLookupLists);
            if(this.HOSECodeLookupIndices.get(HOSECodePrev).isEmpty()){
                this.HOSECodeLookupIndices.remove(HOSECodePrev);
                this.HOSECodeLookupShifts.remove(HOSECodePrev);
            }            
        }
        if (!this.HOSECodeLookupShifts.containsKey(HOSECode)) {
            this.HOSECodeLookupShifts.put(HOSECode, new ArrayList<>());
            this.HOSECodeLookupIndices.put(HOSECode, new ArrayList<>());
        }
        if(!this.HOSECodeLookupIndices.get(HOSECode).contains(atomIndexInSubstructure)){
            this.HOSECodeLookupIndices.get(HOSECode).add(atomIndexInSubstructure);
            final Integer signalIndex = this.assignment.getIndex(0, atomIndexInSubstructure);
            if (signalIndex != null) {
                final Signal signal = this.subspectrum.getSignal(signalIndex);
                if ((signal != null) && (signal.getShift(0) != null)) {
                    this.HOSECodeLookupShifts.get(HOSECode).add(signal.getShift(0));
                }
            }
        }                
                
        return true;
    }
    
    public String getHOSECode(final int atomIndexInSubstructure){
        if(!Utils.checkIndexInAtomContainer(this.substructure, atomIndexInSubstructure)){
            return null;
        }
        for (final String HOSECode : this.getHOSECodeLookupIndices().keySet()) {
            if (this.getHOSECodeLookupIndices().get(HOSECode).contains(atomIndexInSubstructure)) {
                return HOSECode;
            }
        }        
        
        return null;
    }

    /**
     * @param atomIndexInSubstructure
     * @param sphere
     * @return
     *
     * @deprecated
     */
    public ArrayList<Integer> getAtomIndicesInHOSECodeSpheres(final int atomIndexInSubstructure, final int sphere) {
        if (!Utils.checkIndexInAtomContainer(this.substructure, atomIndexInSubstructure)) {
            return null;
        }

        return this.connectionTrees.get(atomIndexInSubstructure).getNodeKeysInSphere(sphere);
    }

    public ConnectionTree getConnectionTree(final int atomIndexInSubstructure){
        if(!Utils.checkIndexInAtomContainer(this.substructure, atomIndexInSubstructure)){
            return null;
        }
        
        return this.connectionTrees.get(atomIndexInSubstructure);
    }
    
    /**
     * Returns the root atom index of this substructure.
     *
     * @return
     */
    public int getRootAtomIndex(){
        return this.rootAtomIndex;
    }
    
    /**
     * Checks whether an atom is unsaturated or not.
     *
     * @param atomIndex index of atom to check
     * @return
     */
    public Boolean isUnsaturated(final int atomIndex){
        if(!Utils.checkIndexInAtomContainer(this.substructure, atomIndex)){
            return null;
        }                
        
        return this.getUnsaturatedAtomIndices().contains(atomIndex);
    }
    
    /**
     * Checks whether the whole substructure is unsaturated or not.
     *
     * @return
     */
    public boolean hasUnsaturatedAtoms(){
        return !this.getUnsaturatedAtomIndices().isEmpty();
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
    
    public void setIndex(final long index){
        this.index = index;
    }
    
    public long getIndex(){
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
    
    public HashMap<String, ArrayList<Double>> getHOSECodeLookupShifts(){
        return this.HOSECodeLookupShifts;
    }
    
    public HashMap<String, ArrayList<Integer>> getHOSECodeLookupIndices(){
        return this.HOSECodeLookupIndices;
    }        
    
    public HashMap<String, ArrayList<Integer>> getAtomTypeIndices(){
        return this.atomTypeIndices;
    }
    
    public HashMap<String, ArrayList<Integer>> getMultiplicitySections(){
        return this.multiplicitySections;
    }        
    
    public String getSubspectrumAtomType(){
        return Utils.getAtomTypeFromSpectrum(this.subspectrum, 0);
    }
}
