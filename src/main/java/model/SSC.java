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
import hose.model.ConnectionTreeNode;
import org.openscience.cdk.Bond;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.ArrayList;
import java.util.Collections;
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
    private int maxSphere;
    // stores all shifts for each single atom in SSC and its HOSE code
    private final HashMap<Integer, String> HOSECodes;
    private final HashMap<Integer, ArrayList<Double>> shifts;
    private final HashMap<Integer, Double[]> shiftsRanges;
    private Integer[] attachedHydrogensInOuterSphere;
    // connection tree holding the spherical information about the substructure
    private final HashMap<Integer, ConnectionTree> connectionTrees;
    /**
     * stores all atom indices for each occurring atom type in substructure
     *
     * @deprecated
     */
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
        this.HOSECodes = new HashMap<>();
        this.shifts = new HashMap<>();
        this.shiftsRanges = new HashMap<>();

        this.initSignalShifts();

        this.connectionTrees = new HashMap<>();
        this.index = -1;
        this.unsaturatedAtomIndices = new ArrayList<>();
        this.multiplicitySections = new HashMap<>();      
        this.multiplicitySectionsBuilder = new MultiplicitySectionsBuilder();
        this.update();
    }

    /**
     * Constructor to "restore" a SSC object from given inputs without recreating/calculating every class members.
     *
     * @param subspectrum
     * @param assignment
     * @param substructure
     * @param rootAtomIndex
     * @param maxSphere
     * @param multiplicitySections
     * @param shifts
     * @param shiftsRanges
     * @param unsaturatedAtomIndices
     * @throws Exception
     */
//    public SSC(final Spectrum subspectrum, final Assignment assignment, final IAtomContainer substructure, final int rootAtomIndex, final int maxSphere, final HashMap<Integer, String> HOSECodes, final HashMap<Integer, ConnectionTree> connectionTrees, final int[] attachedHydrogensInOuterSphere, final HashMap<Integer, ArrayList<Double>> shifts, final ArrayList<Integer> unsaturatedAtomIndices, final HashMap<String, ArrayList<Integer>> multiplicitySections) throws Exception {
    public SSC(final Spectrum subspectrum, final Assignment assignment, final IAtomContainer substructure, final int rootAtomIndex, final int maxSphere, final HashMap<Integer, ArrayList<Double>> shifts, final HashMap<Integer, Double[]> shiftsRanges, final ArrayList<Integer> unsaturatedAtomIndices, final HashMap<String, ArrayList<Integer>> multiplicitySections) throws Exception {
        this.subspectrum = subspectrum;
        this.assignment = assignment;
        this.substructure = substructure;
        this.rootAtomIndex = rootAtomIndex;
        this.maxSphere = maxSphere;
//        this.HOSECodes = HOSECodes;
//        this.connectionTrees = connectionTrees;
//        this.attachedHydrogensInOuterSphere = attachedHydrogensInOuterSphere;
        this.HOSECodes = new HashMap<>();
        this.connectionTrees = new HashMap<>();
        this.updateHOSECodes();

        this.shifts = shifts;
        this.shiftsRanges = shiftsRanges;
        this.unsaturatedAtomIndices = unsaturatedAtomIndices;
        this.multiplicitySections = multiplicitySections;

        this.index = -1;
        this.multiplicitySectionsBuilder = new MultiplicitySectionsBuilder();
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
      return new SSC(this.subspectrum, this.assignment, this.substructure, this.rootAtomIndex, this.maxSphere);//, this.shifts, this.unsaturatedAtomIndices, this.multiplicitySections);
    }

    /**
     * Returns the SSC substructure as HOSE code with the beginning at the SSC's root atom.
     *
     * @return HOSE code of SSC's root atom
     *
     */
    public String getAsHOSECode() {
        return this.getHOSECode(this.getRootAtomIndex());
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
//        this.updateAtomTypeIndices();
        this.updateUnsaturatedAtomIndices();
        this.updateHOSECodes();
        this.updateMultiplicitySections();
    }

    private void updateAttachedHydrogensInOuterSphere(){
        final ArrayList<Integer> atomIndicesInHOSECodeSphere = this.getAtomIndicesInHOSECodeSphere(this.getRootAtomIndex(), this.getMaxSphere());
        if((atomIndicesInHOSECodeSphere == null) || atomIndicesInHOSECodeSphere.isEmpty()){
            this.attachedHydrogensInOuterSphere = new Integer[0];

            return;
        }
        this.attachedHydrogensInOuterSphere = new Integer[atomIndicesInHOSECodeSphere.size()];

        for (int i = 0; i < this.attachedHydrogensInOuterSphere.length; i++) {
            this.attachedHydrogensInOuterSphere[i] = this.getSubstructure().getAtom(atomIndicesInHOSECodeSphere.get(i)).getImplicitHydrogenCount();
        }

    }

    public Integer[] getAttachedHydrogensInOuterSphere(){
        return this.attachedHydrogensInOuterSphere;
    }
    
    @Override
    public String toString(){
        String output = "\natom types and indices:";
//        for (final String atomType : this.atomTypeIndices.keySet()) {
//            output += "\n-> " + atomType + ": " +  this.atomTypeIndices.get(atomType);
//        }
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

    private void initSignalShifts(){
        for (int i = 0; i < this.subspectrum.getSignalCount(); i++) {
            this.shifts.put(this.assignment.getAssignment(0, i), new ArrayList<>());
            this.shifts.get(this.assignment.getAssignment(0, i)).add(this.subspectrum.getSignal(i).getShift(0));
            this.shiftsRanges.put(this.assignment.getAssignment(0, i), new Double[]{this.subspectrum.getSignal(i).getShift(0), this.subspectrum.getSignal(i).getShift(0)});
        }
    }

    public boolean addAtom(final IAtom atom, final IBond bond, final int parentAtomIndexinSubstructure){
        if (!Utils.checkIndexInAtomContainer(this.substructure, parentAtomIndexinSubstructure)) {
            return false;
        }
        this.substructure.addAtom(atom);

        if(!this.addBond(bond, parentAtomIndexinSubstructure, this.substructure.indexOf(atom))){
            this.substructure.removeAtom(atom);

            return false;
        }

        return true;
    }

    public boolean addBond(final IBond bond, final int parentAtomIndexInSubstructure1, final int parentAtomIndexInSubstructure2){
        if (!Utils.checkIndexInAtomContainer(this.substructure, parentAtomIndexInSubstructure1)
                || !Utils.checkIndexInAtomContainer(this.substructure, parentAtomIndexInSubstructure2)) {
            return false;
        }
        final IBond bondToAdd = new Bond(this.substructure.getAtom(parentAtomIndexInSubstructure1), this.substructure.getAtom(parentAtomIndexInSubstructure2), bond.getOrder());
        bondToAdd.setIsInRing(bond.isInRing());
        bondToAdd.setIsAromatic(bond.isAromatic());

        if(!Utils.isValidBondAddition(this.substructure, parentAtomIndexInSubstructure1, bondToAdd)
                || !Utils.isValidBondAddition(this.substructure, parentAtomIndexInSubstructure2, bondToAdd)){
            return false;
        }
        this.substructure.addBond(bondToAdd);

        return true;
    }

    public boolean addSignal(final Signal newSignal, final int atomIndexInSubstructure){
        if (!Utils.checkIndexInAtomContainer(this.substructure, atomIndexInSubstructure)) {
            return false;
        }
        if(newSignal == null){
            return false;
        }
        this.subspectrum.addSignal(newSignal);
        this.assignment.addAssignment(new int[]{atomIndexInSubstructure});

        return true;
    }

    public boolean addShiftsFromSubspectrum(final Spectrum subspectrum){
        // 1: check spectra dimensions and sizes
        if((subspectrum.getNDim() != this.subspectrum.getNDim())
                || (!subspectrum.getNuclei()[0].equals(this.subspectrum.getNuclei()[0]))
                || (subspectrum.getSignalCount() != this.subspectrum.getSignalCount())){
            return false;
        }
        // 2: check each input signal's multiplicity (in order)
        Signal signal, newSignal;
        for (int i = 0; i < this.subspectrum.getSignalCount(); i++) {
            signal = this.subspectrum.getSignal(i);
            newSignal = subspectrum.getSignal(i);
            if(!signal.getMultiplicity().equals(newSignal.getMultiplicity())){
                return false;
            }
        }
        // 3: add shift values of each input signal to the shifts list; for each signal in subspectrum of this SSC
        for (int i = 0; i < this.subspectrum.getSignalCount(); i++) {
            this.shifts.get(this.assignment.getAssignment(0, i)).add(subspectrum.getSignal(i).getShift(0));
        }
        // 4: calculate/update the median shift of all stored shifts for a signal in subspectrum of this SSC
        //    and update the shift ranges (minimum, maximum)
        for (int i = 0; i < this.subspectrum.getSignalCount(); i++) {
            this.subspectrum.setShift(Utils.getMedian(this.shifts.get(this.getAssignments().getAssignment(0, i))), 0, i);
            this.shiftsRanges.put(i, new Double[]{Collections.min(this.shifts.get(this.getAssignments().getAssignment(0, i))), Collections.max(this.shifts.get(this.getAssignments().getAssignment(0, i)))});
        }

        return true;
    }

    public HashMap<Integer, Double[]> getShiftsRanges(){
        return this.shiftsRanges;
    }

    /**
     * Removes shifts for each signal of this SSC subspectrum which are considered as outlier ({@link Utils#removeOutliers(ArrayList, double)}).
     */
    public void removeSignalShiftsOutlier(){
        for (int i = 0; i < this.subspectrum.getSignalCount(); i++) {
            this.removeSignalShiftsOutlier(i);
        }
    }

    private void removeSignalShiftsOutlier(final int atomIndexInSubstructure){
        this.subspectrum.setShift(Utils.getMedian(Utils.removeOutliers(this.shifts.get(this.getAssignments().getAssignment(0, atomIndexInSubstructure)), 1.5)), 0, atomIndexInSubstructure);
    }

    private boolean updateConnectionTree(final int atomIndexInSubstructure) throws CDKException {
        if (!Utils.checkIndexInAtomContainer(this.substructure, atomIndexInSubstructure)) {
            return false;
        }
        this.connectionTrees.put(atomIndexInSubstructure, HOSECodeBuilder.buildConnectionTree(this.substructure, atomIndexInSubstructure, null));//this.maxSphere));
        if(atomIndexInSubstructure == this.getRootAtomIndex()){
            boolean sphereContainsNonRingClosureNodes;
            int maxSphereBuilt = this.getConnectionTree(atomIndexInSubstructure).getMaxSphere();
            for (int s = 0; s <= maxSphereBuilt; s++) {
                sphereContainsNonRingClosureNodes = false;
                for (final ConnectionTreeNode nodeInLastSphere : this.getConnectionTree(atomIndexInSubstructure).getNodesInSphere(s)){
                    if(!nodeInLastSphere.isRingClosureNode()){
                        sphereContainsNonRingClosureNodes = true;
                        break;
                    }
                }
                if (sphereContainsNonRingClosureNodes){
                    this.maxSphere = s;
                } else {
                    break;
                }
            }
        }
        
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

    /**
     *
     * @deprecated
     */
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
            if(!this.updateHOSECode(i)){
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
        this.HOSECodes.put(atomIndexInSubstructure, HOSECodeBuilder.buildHOSECode(this.connectionTrees.get(atomIndexInSubstructure), false));

        if(atomIndexInSubstructure == this.getRootAtomIndex()){
            this.updateAttachedHydrogensInOuterSphere();
        }

        return true;
    }
    
    public String getHOSECode(final int atomIndexInSubstructure){
        if(!Utils.checkIndexInAtomContainer(this.substructure, atomIndexInSubstructure)){
            return null;
        }
        
        return this.HOSECodes.get(atomIndexInSubstructure);
    }

    /**
     * Returns the atom indices in an HOSE code sphere starting at a certain root atom.
     * The indices of ring closure nodes in connection trees are removed.
     *
     * @param atomIndexInSubstructure
     * @param sphere
     * @return
     *
     */
    public ArrayList<Integer> getAtomIndicesInHOSECodeSphere(final int atomIndexInSubstructure, final int sphere) {
        if (!Utils.checkIndexInAtomContainer(this.substructure, atomIndexInSubstructure)) {
            return null;
        }
        final ArrayList<ConnectionTreeNode> nodesInSphere = this.connectionTrees.get(atomIndexInSubstructure).getNodesInSphere(sphere);
        final ArrayList<Integer> nodeIndicesInSphere = new ArrayList<>();
        // ignore ring closure node indices
        for (final ConnectionTreeNode nodeInSphere : nodesInSphere){
            if(!nodeInSphere.isRingClosureNode()){
                nodeIndicesInSphere.add(nodeInSphere.getKey());
            }
        }

        return nodeIndicesInSphere;
    }

    public ConnectionTree getConnectionTree(final int atomIndexInSubstructure){
        if(!Utils.checkIndexInAtomContainer(this.substructure, atomIndexInSubstructure)){
            return null;
        }
        
        return this.connectionTrees.get(atomIndexInSubstructure);
    }

    public HashMap<Integer, ConnectionTree> getConnectionTrees(){
        return this.connectionTrees;
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

    public HashMap<Integer, String> getHOSECodes(){
        return this.HOSECodes;
    }
    public HashMap<Integer, ArrayList<Double>> getShifts(){
        return this.shifts;
    }

    /**
     * @return
     *
     * @deprecated
     */
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
