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
import com.google.gson.Gson;
import com.google.gson.reflect.TypeToken;
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

    private static final MultiplicitySectionsBuilder MULTIPLICITY_SECTIONS_BUILDER = new MultiplicitySectionsBuilder();;

    private final Spectrum subspectrum;
    private final Assignment assignment;
    private final IAtomContainer substructure;
    // spherical search limit
    private final int rootAtomIndex;
    private int maxSphere;
    // stores all shifts for each single atom in SSC and its HOSE code
    private final HashMap<Integer, String> HOSECodes;
    private final HashMap<Integer, ArrayList<Double>> shifts;
    private final HashMap<Integer, Double[]> shiftsMinMax;
    private Integer[] attachedHydrogensCountsInOuterSphere;
    private String[] multiplicitiesInOuterSphere;
    // connection tree holding the spherical information about the substructure
    private final HashMap<Integer, ConnectionTree> connectionTrees;
    // for pre-search: map of multiplicities as keys consisting 
    // of section indices with occurrence in this SSC
    private HashMap<String, ArrayList<Integer>> multiplicitySections;
    // index to use in SSC library
    private long index;
    // indices of open-sphere (unsaturated) atoms of substructure
    private final ArrayList<Integer> unsaturatedAtomIndices;    
    //public final static int MIN_LIMIT = -20, MAX_LIMIT = 260, STEP_SIZE = 5, STEPS = (MAX_LIMIT - MIN_LIMIT) / STEP_SIZE; // ppm range from -20 to 260 in 5 ppm steps


    /**
     * @param subspectrum
     * @param assignment
     * @param substructure
     * @param rootAtomIndex
     * @param maxSphere
     *
     * @throws Exception
     */
    public SSC(final Spectrum subspectrum, final Assignment assignment, final IAtomContainer substructure, final int rootAtomIndex, final int maxSphere) throws Exception {
        this.subspectrum = subspectrum.getClone();
        this.assignment = assignment.clone();
        this.substructure = substructure.clone();     
        this.rootAtomIndex = rootAtomIndex;
        this.maxSphere = maxSphere;
        this.HOSECodes = new HashMap<>();
        this.shifts = new HashMap<>();
        this.shiftsMinMax = new HashMap<>();
        this.initSignalShifts();
        this.connectionTrees = new HashMap<>();
        this.index = -1;
        this.unsaturatedAtomIndices = new ArrayList<>();
        this.multiplicitySections = new HashMap<>();      
        this.update();
    }


    /**
     * Constructor to "restore" a SSC object from given inputs without creating/calculating every class members again.
     *
     * @param subspectrum
     * @param assignment
     * @param substructure
     * @param rootAtomIndex
     * @param maxSphere
//     * @param HOSECodes
//     * @param connectionTrees
     * @param shifts
     * @param shiftsMinMax
     * @param unsaturatedAtomIndices
     * @param multiplicitySections
//     * @param attachedHydrogensCountsInOuterSphere
//     * @param multiplicitiesInOuterSphere
     *
     */
    public SSC(final Spectrum subspectrum, final Assignment assignment, final IAtomContainer substructure, final int rootAtomIndex, final int maxSphere,
//               final HashMap<Integer, String> HOSECodes, final HashMap<Integer, ConnectionTree> connectionTrees,
               final HashMap<Integer, ArrayList<Double>> shifts, final HashMap<Integer, Double[]> shiftsMinMax,
               final ArrayList<Integer> unsaturatedAtomIndices, final HashMap<String, ArrayList<Integer>> multiplicitySections
//               final Integer[] attachedHydrogensCountsInOuterSphere, final String[] multiplicitiesInOuterSphere
            ) throws Exception {
        this.subspectrum = subspectrum.getClone();
        this.assignment = assignment.clone();
        this.substructure = substructure.clone();
        this.rootAtomIndex = rootAtomIndex;
        this.maxSphere = maxSphere;

        // Deep clones
        final Gson gson = new Gson();
//        this.HOSECodes = gson.fromJson(gson.toJson(HOSECodes), HashMap.class);
//        this.connectionTrees = gson.fromJson(gson.toJson(connectionTrees), new TypeToken<HashMap<Integer, ConnectionTree>>(){}.getType());
        this.HOSECodes = new HashMap<>();
        this.connectionTrees = new HashMap<>();
        this.updateHOSECodes();

        this.shifts = gson.fromJson(gson.toJson(shifts), new TypeToken<HashMap<Integer, ArrayList<Double>>>(){}.getType());
        this.shiftsMinMax = gson.fromJson(gson.toJson(shiftsMinMax), new TypeToken<HashMap<Integer, Double[]>>(){}.getType());
        this.unsaturatedAtomIndices = gson.fromJson(gson.toJson(unsaturatedAtomIndices), new TypeToken<ArrayList<Integer>>(){}.getType());
        this.multiplicitySections = gson.fromJson(gson.toJson(multiplicitySections), new TypeToken<HashMap<String, ArrayList<Integer>>>(){}.getType());;

        this.updateAttachedHydrogensCountAndMultiplicitiesInOuterSphere();
//        this.attachedHydrogensCountsInOuterSphere = attachedHydrogensCountsInOuterSphere;
//        this.multiplicitiesInOuterSphere = multiplicitiesInOuterSphere;

        this.index = -1;
    }

    /**
     * Returns a full clone of that SSC. The index of the SSC clone is set to default value (-1). <br>
     * {@code reduceShifts == false} means that all class members are deep cloned, incl. all shifts in shifts lists.
     * The opposite ({@code reduceShifts == true}) means that the lists of shifts for each signal
     * will only contain one (!!!) shift: the previously set representative shift.
     *
     * @param reduceShifts indicates whether a all shifts in shift lists should be reduced
     *
     * @return
     * @throws CDKException
     * 
     * @see #SSC(Spectrum, Assignment, IAtomContainer, int, int, HashMap, HashMap, ArrayList, HashMap)
     */
    public SSC getClone(final boolean reduceShifts) throws Exception {
        if(reduceShifts){
            final HashMap<Integer, ArrayList<Double>> reducedShifts = new HashMap<>();
            for (int i = 0; i < this.getSubspectrum().getSignalCount(); i++) {
                reducedShifts.put(i, new ArrayList<>());
                reducedShifts.get(i).add(this.getSubspectrum().getShift(i, 0));
            }

            return new SSC(this.subspectrum, this.assignment, this.substructure, this.rootAtomIndex, this.maxSphere,
//                this.HOSECodes, this.connectionTrees,
                    reducedShifts, this.shiftsMinMax, this.unsaturatedAtomIndices, this.multiplicitySections
//                this.attachedHydrogensCountsInOuterSphere, this.multiplicitiesInOuterSphere
            );
        }

        return new SSC(this.subspectrum, this.assignment, this.substructure, this.rootAtomIndex, this.maxSphere,
//                    this.HOSECodes, this.connectionTrees,
                this.shifts, this.shiftsMinMax, this.unsaturatedAtomIndices, this.multiplicitySections
//                    this.attachedHydrogensCountsInOuterSphere, this.multiplicitiesInOuterSphere
        );
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
     * @see #updateUnsaturatedAtomIndices()
     * @see #updateMultiplicitySections() 
     * @see #updateHOSECodes()
     */
    public void update() throws CDKException {
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(this.substructure);
        this.updateUnsaturatedAtomIndices();
        this.updateHOSECodes();
        this.updateMultiplicitySections();
        this.updateAttachedHydrogensCountAndMultiplicitiesInOuterSphere();
    }

    private void updateAttachedHydrogensCountAndMultiplicitiesInOuterSphere(){
        final ArrayList<Integer> atomIndicesInOuterSphere = this.getAtomIndicesInHOSECodeSphere(this.getRootAtomIndex(), this.getMaxSphere());
        if((atomIndicesInOuterSphere == null) || atomIndicesInOuterSphere.isEmpty()){
            this.attachedHydrogensCountsInOuterSphere = new Integer[0];
            this.multiplicitiesInOuterSphere = new String[0];

            return;
        }
        // update attached hydrogen counts
        this.attachedHydrogensCountsInOuterSphere = new Integer[atomIndicesInOuterSphere.size()];
        for (int i = 0; i < this.attachedHydrogensCountsInOuterSphere.length; i++) {
            this.attachedHydrogensCountsInOuterSphere[i] = this.getSubstructure().getAtom(atomIndicesInOuterSphere.get(i)).getImplicitHydrogenCount();
        }
        // update multiplicities
        final ArrayList<Integer> atomIndicesInOuterSphereWithSignal = new ArrayList<>();
        for (final int atomIndexInOuterSphere : atomIndicesInOuterSphere){
            if(this.assignment.getIndex(0, atomIndexInOuterSphere) != null){
                atomIndicesInOuterSphereWithSignal.add(atomIndexInOuterSphere);
            }
        }
        this.multiplicitiesInOuterSphere = new String[atomIndicesInOuterSphereWithSignal.size()];
        int counter = 0;
        for (final int atomIndexInOuterSphere : atomIndicesInOuterSphereWithSignal){
            this.multiplicitiesInOuterSphere[counter] = this.getSubspectrum().getMultiplicity(this.assignment.getIndex(0, atomIndexInOuterSphere));
            counter++;
        }
    }

    public Integer[] getAttachedHydrogensCountsInOuterSphere(){
        return this.attachedHydrogensCountsInOuterSphere;
    }

    public String[] getMultiplicitiesInOuterSphere(){
        return this.multiplicitiesInOuterSphere;
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
            this.shifts.put(i, new ArrayList<>());
            this.shifts.get(i).add(this.subspectrum.getSignal(i).getShift(0));
            this.shiftsMinMax.put(i, new Double[]{this.subspectrum.getSignal(i).getShift(0), this.subspectrum.getSignal(i).getShift(0)});
        }
    }

//    public boolean addAtom(final IAtom atom, final IBond bond, final int parentAtomIndexinSubstructure){
//        return this.addAtom(atom, bond, parentAtomIndexinSubstructure, null);
//    }

    public boolean addAtom(final IAtom atom, final IBond bond, final int parentAtomIndexinSubstructure, final Signal signal){
        if (!this.hasAtom(parentAtomIndexinSubstructure)) {
            return false;
        }
        final IAtom atomClone;
        try {
            atomClone = atom.clone();
            this.getSubstructure().addAtom(atomClone);
        } catch (CloneNotSupportedException e) {
            return false;
        }

        boolean removeInsertedAtom = false;
        if(!this.addBond(bond, parentAtomIndexinSubstructure, this.getSubstructure().indexOf(atomClone))){
            removeInsertedAtom = true;
        }
        if(!removeInsertedAtom && (signal != null)){
            try {
                if(!this.getSubspectrum().addSignal(signal.getClone())
                        || !this.getAssignments().addAssignment(new int[]{this.getAtomCount() - 1})){
                    removeInsertedAtom = true;
                }
            } catch (Exception e) {
                removeInsertedAtom = true;
            }
        }
        if(removeInsertedAtom){
            this.getSubstructure().removeAtom(atomClone);

            return false;
        }

        return true;
    }

    public boolean addBond(final IBond bond, final int atomIndexInSubstructure1, final int atomIndexInSubstructure2){
        if(this.hasBond(atomIndexInSubstructure1, atomIndexInSubstructure2)){
            return false;
        }
        IAtom atom1 = this.getAtom(atomIndexInSubstructure1);
        IAtom atom2 = this.getAtom(atomIndexInSubstructure2);
        if((atom1 == null) || (atom2 == null)){
            return false;
        }
        final IBond bondToAdd = new Bond(atom1, atom2, bond.getOrder());
        bondToAdd.setIsInRing(bond.isInRing());
        bondToAdd.setIsAromatic(bond.isAromatic());

        if(!Utils.isValidBondAddition(this.getSubstructure(), atomIndexInSubstructure1, bondToAdd)
                || !Utils.isValidBondAddition(this.getSubstructure(), atomIndexInSubstructure2, bondToAdd)){
            return false;
        }
        this.getSubstructure().addBond(bondToAdd);

        return true;
    }

    public boolean removeBond(final int atomIndexInSubstructure1, final int atomIndexInSubstructure2){
        if (!this.hasBond(atomIndexInSubstructure1, atomIndexInSubstructure2)) {
            return false;
        }
        this.getSubstructure().removeBond(this.getBond(atomIndexInSubstructure1, atomIndexInSubstructure2));

        return true;
    }

    public boolean removeLastBond(){
        if(this.getBondCount() < 1){
            return false;
        }
        this.getSubstructure().removeBond((this.getBondCount() - 1));

        return true;
    }

    public IBond getBond(final int atomIndexInSubstructure1, final int atomIndexInSubstructure2){
        if (!Utils.checkIndexInAtomContainer(this.getSubstructure(), atomIndexInSubstructure1)
                || !Utils.checkIndexInAtomContainer(this.getSubstructure(), atomIndexInSubstructure2)) {
            return null;
        }
        return this.getSubstructure().getBond(this.getSubstructure().getAtom(atomIndexInSubstructure1), this.getSubstructure().getAtom(atomIndexInSubstructure2));
    }

    public boolean hasBond(final int atomIndexInSubstructure1, final int atomIndexInSubstructure2){
        return this.getBond(atomIndexInSubstructure1, atomIndexInSubstructure2) != null;
    }

    public boolean removeAtom(final int atomIndexInSubstructure){
        if (!this.hasAtom(atomIndexInSubstructure)) {
            return false;
        }
        this.getSubstructure().removeAtom(atomIndexInSubstructure);
        final Integer signalIndex = this.getAssignments().getIndex(0, atomIndexInSubstructure);
        if(signalIndex != null){
            this.getSubspectrum().removeSignal(signalIndex);
            this.getAssignments().removeAssignment(signalIndex);
        }

        return true;
    }

    public boolean removeLastAtom(){
        return this.removeAtom((this.getAtomCount() - 1));
    }

    public IAtom getAtom(final int atomIndexInSubstructure){
        if (!Utils.checkIndexInAtomContainer(this.getSubstructure(), atomIndexInSubstructure)) {
            return null;
        }
        return this.getSubstructure().getAtom(atomIndexInSubstructure);
    }

    public boolean hasAtom(final int atomIndexInSubstructure){
        return this.getAtom(atomIndexInSubstructure) != null;
    }

//    public boolean addSignal(final Signal newSignal, final int atomIndexInSubstructure){
//        if (!Utils.checkIndexInAtomContainer(this.getSubstructure(), atomIndexInSubstructure)) {
//            return false;
//        }
//        if(newSignal == null){
//            return false;
//        }
//        this.subspectrum.addSignal(newSignal);
//        this.assignment.addAssignment(new int[]{atomIndexInSubstructure});
//
//        return true;
//    }

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
            this.shifts.get(i).add(subspectrum.getSignal(i).getShift(0));
        }
        // 4: calculate/update the median shift of all stored shifts for a signal in subspectrum of this SSC
        //    and update the shift ranges (minimum, maximum)
        for (int i = 0; i < this.subspectrum.getSignalCount(); i++) {
            this.subspectrum.setShift(Utils.getMedian(this.shifts.get(i)), 0, i);
            this.shiftsMinMax.put(i, new Double[]{Collections.min(this.shifts.get(i)), Collections.max(this.shifts.get(i))});
        }

        return true;
    }

    public HashMap<Integer, Double[]> getShiftsMinMax(){
        return this.shiftsMinMax;
    }

    /**
     * Removes shifts for each signal of this SSC subspectrum which are considered as outlier ({@link Utils#removeOutliers(ArrayList, double)}).
     */
    public void removeSignalShiftsOutlier(){
        for (int i = 0; i < this.subspectrum.getSignalCount(); i++) {
            this.removeSignalShiftsOutlier(i);
        }
    }

    private void removeSignalShiftsOutlier(final int signalIndexInSubspectrum){
        this.shifts.put(signalIndexInSubspectrum, Utils.removeOutliers(this.shifts.get(signalIndexInSubspectrum), 1.5));
        this.subspectrum.setShift(Utils.getMedian(this.shifts.get(signalIndexInSubspectrum)), 0, signalIndexInSubspectrum);
    }

    private boolean updateConnectionTree(final int atomIndexInSubstructure) {
        if (!Utils.checkIndexInAtomContainer(this.substructure, atomIndexInSubstructure)) {
            return false;
        }
        this.connectionTrees.put(atomIndexInSubstructure, HOSECodeBuilder.buildConnectionTree(this.substructure, atomIndexInSubstructure, null));//this.maxSphere));
        if(atomIndexInSubstructure == this.getRootAtomIndex()){
            boolean sphereContainsNonRingClosureNodes;
            for (int s = 0; s <= this.getConnectionTree(atomIndexInSubstructure).getMaxSphere(); s++) {
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
    
    private void updateUnsaturatedAtomIndices() {
        this.unsaturatedAtomIndices.clear();
        for (int i = 0; i < this.substructure.getAtomCount(); i++) {
            // set the indices of unsaturated atoms in substructure
            if (!Utils.isSaturated(this.substructure, i)) {
                this.unsaturatedAtomIndices.add(i);
            }            
        }
    }   
    
    /**
     * Returns the indices of unsaturated atoms in substructure.
     *
     * @return
     */
    public ArrayList<Integer> getUnsaturatedAtomIndices(){
        return this.unsaturatedAtomIndices;
    }

    /**
     * Specified for carbons only -> still open.
     *
     * @throws org.openscience.cdk.exception.CDKException
     */
    private void updateMultiplicitySections() throws CDKException {
        this.multiplicitySections = SSC.MULTIPLICITY_SECTIONS_BUILDER.buildMultiplicitySections(this.subspectrum);
    }    

    private boolean updateHOSECodes() throws CDKException {
        for (int i = 0; i < this.getAtomCount(); i++) {
            if(!this.updateHOSECode(i)){
                return false;
            }
        }
        
        return true;
    }
    
    private boolean updateHOSECode(final int atomIndexInSubstructure) throws CDKException {
        if(!Utils.checkIndexInAtomContainer(this.substructure, atomIndexInSubstructure)){
            return false;
        }
        if(!this.updateConnectionTree(atomIndexInSubstructure)){
            return false;
        }
        this.HOSECodes.put(atomIndexInSubstructure, HOSECodeBuilder.buildHOSECode(this.connectionTrees.get(atomIndexInSubstructure), false));

        if(atomIndexInSubstructure == this.getRootAtomIndex()){
            this.updateAttachedHydrogensCountAndMultiplicitiesInOuterSphere();
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
    
    public HashMap<String, ArrayList<Integer>> getMultiplicitySections(){
        return this.multiplicitySections;
    }        
    
    public String getSubspectrumAtomType(){
        return Utils.getAtomTypeFromSpectrum(this.subspectrum, 0);
    }
}
