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
package assembly;

import fragmentation.HOSECodeBuilder;
import casekit.NMR.Utils;
import casekit.NMR.model.Assignment;
import casekit.NMR.model.Signal;
import casekit.NMR.model.Spectrum;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import model.ConnectionTree;
import model.ConnectionTreeNode;
import model.SSC;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import rank.SSCRanker;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Assembly {  
    
    /**
     * Returns pairwise structural identities between two atoms in both 
     * substructures, including the maximum matching sphere count. 
     * 
     * @param ssc1
     * @param ssc2
     * @param minSphereMatchCount
     * @param maxSphereMatchCount
     * @return
     */
    public static HashMap<Integer, ArrayList<Integer[]>> getStructuralOverlapsHOSECode(final SSC ssc1, final SSC ssc2, final int minSphereMatchCount, final int maxSphereMatchCount) {
        final HashMap<Integer, ArrayList<Integer[]>> overlaps = new HashMap<>();
        
        int maxMatchingSphere;
        for (int i = 0; i < ssc1.getAtomCount(); i++) {
            for (int j = 0; j < ssc2.getAtomCount(); j++) {
                if(i == j){
                    continue;
                }
                maxMatchingSphere = Assembly.getMaximumMatchingSphereHOSECode(ssc1.getHOSECode(i), ssc2.getHOSECode(j));
                
                if((maxMatchingSphere >= minSphereMatchCount) && (maxMatchingSphere <= maxSphereMatchCount)){
                    if (!overlaps.containsKey(maxMatchingSphere)) {
                        overlaps.put(maxMatchingSphere, new ArrayList<>());
                    }
                    overlaps.get(maxMatchingSphere).add(new Integer[]{i, j});
                }    
            }
        }
       
        return overlaps;                                 
    }    
    
    public static int getMaximumMatchingSphereHOSECode(final String HOSECode1, final String HOSECode2) {
        final ArrayList<String> HOSECodeIntoSpheresSSC1 = HOSECodeBuilder.splitHOSECodeIntoSpheres(HOSECode1), HOSECodeIntoSpheresSSC2 = HOSECodeBuilder.splitHOSECodeIntoSpheres(HOSECode2);
        int maxMatchingSphere = -1;
        for (int s = 0; s < Integer.min(HOSECodeIntoSpheresSSC1.size(), HOSECodeIntoSpheresSSC2.size()); s++) {
            if (!HOSECodeIntoSpheresSSC1.get(s).equals(HOSECodeIntoSpheresSSC2.get(s))) {
                break;
            }
            maxMatchingSphere = s;
        }

        return maxMatchingSphere;
    }
    
    /**
     * Returns pairwise atom indices of two structural overlapping atoms in both
     * substructures, including the maximum matching sphere counts as keys.
     *
     * @param ssc1 first SSC
     * @param ssc2 second SSC
     * @param shiftTol
     * @param minSphereMatchCount number of minimum matching spheres
     * @param atomType
     * @return pairs of atom indices; first in SSC1, second in SSC2
     * @throws java.lang.CloneNotSupportedException
     */
    public static HashMap<Integer, ArrayList<Integer[]>> getStructuralOverlaps(final SSC ssc1, final SSC ssc2, final double shiftTol, final int minSphereMatchCount, final String atomType) throws CloneNotSupportedException {
        final HashMap<Integer, ArrayList<Integer[]>> overlaps = new HashMap<>();
        int maxMatchingSphere;
        for (int i = 0; i < ssc1.getAtomCount(); i++) {    
            if(!ssc1.getSubstructure().getAtom(i).getSymbol().equals(atomType)){
                continue;
            }
            for (int j = 0; j < ssc2.getAtomCount(); j++) {
                if (!ssc2.getSubstructure().getAtom(j).getSymbol().equals(atomType)) {
                    continue;
                }
                
                maxMatchingSphere = -1;
//                if(Assembly.isEqualNode(ssc1.getConnectionTree(i).getRootNode(), ssc2.getConnectionTree(j).getRootNode(), 
//                            ssc1.getSubspectrum().getSignal(ssc1.getAssignments().getSignalIndex(0, i)), 
//                            ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getSignalIndex(0, j)), shiftTol)){
//                    maxMatchingSphere = Assembly.getMaximumMatchingSphereHOSECode(ssc1.getHOSECode(i), ssc2.getHOSECode(j));
                    maxMatchingSphere = Assembly.getMaximumMatchingSphere(ssc1, ssc2, i, j, shiftTol);
////                    System.out.println("i: " + i + ", " + j + ": " + maxMatchingSphere);
//                }
                
                if ((maxMatchingSphere >= minSphereMatchCount)){// && (maxMatchingSphere <= maxSphereMatchCount)) {
                    if(!overlaps.containsKey(maxMatchingSphere)){
                        overlaps.put(maxMatchingSphere, new ArrayList<>());
                    }
                    overlaps.get(maxMatchingSphere).add(new Integer[]{i, j});
                }
            }
        }

        return overlaps;
    }
    
    public static int getMaximumMatchingSphere(final SSC ssc1, final SSC ssc2, final int atomIndexInSubstructure1, final int atomIndexInSubstructure2, final double shiftTol) throws CloneNotSupportedException {
        int maxMatchingSphere = -1;
        for (int s = 0; s <= Integer.min(ssc1.getConnectionTree(atomIndexInSubstructure1).getMaxSphere(), ssc2.getConnectionTree(atomIndexInSubstructure2).getMaxSphere()); s++) {
//            if(!Assembly.hasStructuralIdentity(ssc1, ssc2, atomIndexInSubstructure1, atomIndexInSubstructure2, s, shiftTol)){
            if(!Assembly.hasEqualNodesInSphere(ssc1, ssc2, atomIndexInSubstructure1, atomIndexInSubstructure2, s, shiftTol)){
                break;
            }
            maxMatchingSphere = s;
        }
        
        return maxMatchingSphere;
    }
    
//    public static boolean hasStructuralIdentity(final SSC ssc1, final SSC ssc2, final int atomIndexInSubstructure1, final int atomIndexInSubstructure2, final int sphere, final double shiftTol) throws CloneNotSupportedException {
//        if(!Assembly.hasEqualNodesInSphere(ssc1, ssc2, atomIndexInSubstructure1, atomIndexInSubstructure2, sphere, shiftTol)){
//            return false;
//        }
//        for (int i = 0; i < nodesInSphere1.size(); i++) {
//            if (!Assembly.hasEqualChildNodesInNextSphere(connectionTree1, connectionTree2, nodesInSphere1.get(i), nodesInSphere2.get(i))
//                    || !Assembly.hasEqualChildNodesInNextSphere(connectionTree2, connectionTree1, nodesInSphere2.get(i), nodesInSphere1.get(i))
//                    ) {
//                return false;
//            }
//        }
//
//        return true;
//    }
    
    /**
     * Checks whether connectionTree1 contains the exactly same node or bond 
     * (at ring closure node) properties as in connectionTree2 in a specific 
     * sphere.  
     *
     * @param ssc1
     * @param ssc2
     * @param atomIndexInSubstructure1
     * @param atomIndexInSubstructure2
     * @param sphere
     * @param shiftTol
     * @return
     */
    public static boolean hasEqualNodesInSphere(final SSC ssc1, final SSC ssc2, final int atomIndexInSubstructure1, final int atomIndexInSubstructure2, final int sphere, final double shiftTol) {
        final ArrayList<ConnectionTreeNode> nodesInSphere1 = ssc1.getConnectionTree(atomIndexInSubstructure1).getNodesInSphere(sphere);
        final ArrayList<ConnectionTreeNode> nodesInSphere2 = ssc2.getConnectionTree(atomIndexInSubstructure2).getNodesInSphere(sphere);
        if (nodesInSphere1.size() != nodesInSphere2.size()) {
            return false;
        }
        Signal signal1, signal2;
        for (int j = 0; j < nodesInSphere1.size(); j++) {
            signal1 = ssc1.getSubspectrum().getSignal(ssc1.getAssignments().getSignalIndex(0, nodesInSphere1.get(j).getKey()));
            signal2 = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getSignalIndex(0, nodesInSphere2.get(j).getKey()));
            if(!isEqualNode(nodesInSphere1.get(j), nodesInSphere2.get(j), signal1, signal2, shiftTol)){
                return false;                
            }            
        }

        return true;
    }
    
    
    public static boolean isValidBondAddition(final IAtomContainer ac, final int atomIndex, final IBond bondToAdd){
        double existingBondsOrderSum = 0;
        final IAtom atom = ac.getAtom(atomIndex);
        for (final IBond bond : ac.getConnectedBondsList(atom)) {
            existingBondsOrderSum += bond.getOrder().numeric();
        }
        int implicitHydrogenCount = (atom.getImplicitHydrogenCount() != null) ? atom.getImplicitHydrogenCount() : 0;
        
        System.out.println(" --> " + existingBondsOrderSum + " + " + implicitHydrogenCount + " + " + bondToAdd.getOrder().numeric() + " = " + (existingBondsOrderSum + implicitHydrogenCount + bondToAdd.getOrder().numeric()) + " <= " + atom.getValency() + " ? ");
        
        return (existingBondsOrderSum + implicitHydrogenCount + bondToAdd.getOrder().numeric()) <= atom.getValency();
    }
    
    
//    /**
//     * Checks whether the exact same child nodes of a node1 in connectionTree1 
//     * can be found at node2 in connectionTree2 in same next sphere and with 
//     * same stored bond properties.
//     *
//     * @param connectionTree1
//     * @param connectionTree2
//     * @param node1
//     * @param node2
//     * @return
//     * @throws CloneNotSupportedException
//     */
//    public static boolean hasEqualChildNodesInNextSphere(final ConnectionTree connectionTree1, final ConnectionTree connectionTree2, final ConnectionTreeNode node1, final ConnectionTreeNode node2) throws CloneNotSupportedException{
//        // check whether all atoms from sphere in SSC1 were covered in SSC2
//        return Assembly.getRetainedChildNodesInNextSphere(connectionTree1, connectionTree2, node1, node2).isEmpty();
//    }
//    
//    /**
//     * Returns the child nodes of node1 in connectionTree1 which can not be 
//     * found as child nodes from node2 in connectionTree2.
//     *
//     * @param connectionTree1
//     * @param connectionTree2
//     * @param node1
//     * @param node2
//     * @return 
//     * @throws CloneNotSupportedException
//     */
//    public static ArrayList<ConnectionTreeNode> getRetainedChildNodesInNextSphere(final ConnectionTree connectionTree1, final ConnectionTree connectionTree2, final ConnectionTreeNode node1, final ConnectionTreeNode node2) throws CloneNotSupportedException{
//        final ArrayList<ConnectionTreeNode> childNodesInNextSphere1 = new ArrayList<>(node1.getChildNodes());
//        final ArrayList<ConnectionTreeNode> childNodesInNextSphere2 = node2.getChildNodes();
//        final ArrayList<ConnectionTreeNode> childNodesInNextSphere1ToRemove = new ArrayList<>();
//        // remove all occurring child nodes in SSC2 from child nodes list in SSC1
//        for (final ConnectionTreeNode childNodeInNextSphere2 : childNodesInNextSphere2) {
//            if (childNodesInNextSphere1.isEmpty()) {
//                break;
//            }
//            childNodesInNextSphere1ToRemove.clear();
//            for (final ConnectionTreeNode childNodeInNextSphere1 : childNodesInNextSphere1) {
//                if(Assembly.isEqualNode(childNodeInNextSphere1, childNodeInNextSphere2)){
//                    childNodesInNextSphere1ToRemove.add(childNodeInNextSphere1);
//                }
//            }
//            childNodesInNextSphere1.removeAll(childNodesInNextSphere1ToRemove);
//        } 
//           
//        return childNodesInNextSphere1;
//    }
//    
//    
//    /**
//     * Returns the nodes in connectionTree1 which can not be
//     * found in connectionTree2 in the same sphere.
//     *
//     * @param connectionTree1
//     * @param connectionTree2
//     * @param sphere
//     * @return
//     * @throws CloneNotSupportedException
//     */
//    public static ArrayList<ConnectionTreeNode> getRetainedNodesInSphere(final ConnectionTree connectionTree1, final ConnectionTree connectionTree2, final int sphere) throws CloneNotSupportedException {
//        if(sphere > Integer.min(connectionTree1.getMaxSphere(), connectionTree2.getMaxSphere())){
//            return null;
//        }
//        final ArrayList<ConnectionTreeNode> retainedNodesInSphere = new ArrayList<>(connectionTree1.getNodesInSphere(sphere));
//        final ArrayList<ConnectionTreeNode> equalNodesInSphere = Assembly.getEqualNodesInSphere(connectionTree1, connectionTree2, sphere);
//        for (final ConnectionTreeNode equalNodeInSphere : equalNodesInSphere) {
//            retainedNodesInSphere.remove(equalNodeInSphere);
//        }
//        
//        return retainedNodesInSphere;
//    }
    
    /**
     * Returns the nodes in connectionTree1 which can not be
     * found in connectionTree2 in the same sphere.
     *
     * @param ssc1
     * @param ssc2
     * @param atomIndexInSubstructure1
     * @param atomIndexInSubstructure2
     * @param sphere
     * @param shiftTol
     * @return
     * @throws CloneNotSupportedException
     */
    public static HashMap<Integer, Integer> mapEqualNodesInSphere(final SSC ssc1, final SSC ssc2, final int atomIndexInSubstructure1, final int atomIndexInSubstructure2, final int sphere, final double shiftTol) throws CloneNotSupportedException {
        final ConnectionTree connectionTree1 = ssc1.getConnectionTree(atomIndexInSubstructure1);
        final ConnectionTree connectionTree2 = ssc2.getConnectionTree(atomIndexInSubstructure2);
        if (sphere > Integer.min(connectionTree1.getMaxSphere(), connectionTree2.getMaxSphere())) {
            return null;
        }
        final HashMap<Integer, Integer> map = new HashMap<>();
        Signal signal1, signal2;
        // add matched node indices in connectionTree2 from nodes list in connectionTree1 in given sphere
        for (final ConnectionTreeNode nodeInSphere1 : connectionTree1.getNodesInSphere(sphere)) {
            // search for nodes from conn. tree 2 in conn. tree 1; if found then add it to map
            for (final ConnectionTreeNode nodeInSphere2 : connectionTree2.getNodesInSphere(sphere)) {
                signal1 = ssc1.getSubspectrum().getSignal(ssc1.getAssignments().getSignalIndex(0, nodeInSphere1.getKey()));
                signal2 = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getSignalIndex(0, nodeInSphere2.getKey()));
                if (!map.containsValue(nodeInSphere2.getKey()) && Assembly.isEqualNode(nodeInSphere1, nodeInSphere2, signal1, signal2, shiftTol)) {
                    map.put(nodeInSphere1.getKey(), nodeInSphere2.getKey());
                    break;
                }
            }
            if(!map.containsKey(nodeInSphere1.getKey())){
                map.put(nodeInSphere1.getKey(), -1);
            }
        }

        return map;
    }
    
    /**
     * Returns whether two nodes have certain identical properties.
     *
     * @param nodeInSphere1
     * @param nodeInSphere2
     * @param signal1
     * @param signal2
     * @param shiftTol
     * @return
     */
    public static boolean isEqualNode(final ConnectionTreeNode nodeInSphere1, final ConnectionTreeNode nodeInSphere2, final Signal signal1, final Signal signal2, final double shiftTol){
        final IAtom atom1 = nodeInSphere1.getAtom(), atom2 = nodeInSphere2.getAtom();
        final IBond bondToParent1, bondToParent2;        
        final ConnectionTreeNode parentNode1, parentNode2;
        if ((atom1 != null) && (atom2 != null) && !nodeInSphere1.isRingClosureNode() && !nodeInSphere2.isRingClosureNode()) {
            if (atom1.getSymbol().equals(atom2.getSymbol())
                    && (Integer.compare(atom1.getImplicitHydrogenCount(), atom2.getImplicitHydrogenCount()) == 0)
                    && atom1.isAromatic() == atom2.isAromatic()
                    && atom1.isInRing() == atom2.isInRing()
                    ) {
                if((signal1 != null) && (signal2 != null)
                    && ((Math.abs(signal1.getShift(0) - signal2.getShift(0)) > shiftTol)
                        || !signal1.getMultiplicity().equals(signal2.getMultiplicity()))){
                    return false;
                }
                
//                if ((nodeInSphere1.getSphere() > 0) && (nodeInSphere2.getSphere() > 0)) {
//                    parentNode1 = nodeInSphere1.getParentNodes().get(0);
//                    parentNode2 = nodeInSphere2.getParentNodes().get(0);
//                    bondToParent1 = nodeInSphere1.getBondsToParents().get(0);
//                    bondToParent2 = nodeInSphere2.getBondsToParents().get(0);
//                    if (parentNode1.getAtom().getSymbol().equals(parentNode2.getAtom().getSymbol())
//                            && bondToParent1.getOrder() == bondToParent2.getOrder()
//                            && bondToParent1.isAromatic() == bondToParent2.isAromatic()
//                            && bondToParent1.isInRing() == bondToParent2.isInRing()) {
//                        return true;
//                    } 
//                } else if((nodeInSphere1.getSphere() == 0) && (nodeInSphere2.getSphere() == 0)){
                    return true;                
//                }                                
            }
        } else if (nodeInSphere1.isRingClosureNode() && nodeInSphere2.isRingClosureNode()) {
            parentNode1 = nodeInSphere1.getParentNodes().get(0);
            parentNode2 = nodeInSphere2.getParentNodes().get(0);
            bondToParent1 = nodeInSphere1.getBondsToParents().get(0);
            bondToParent2 = nodeInSphere2.getBondsToParents().get(0);
            if (parentNode1.getAtom().getSymbol().equals(parentNode2.getAtom().getSymbol())
                    && bondToParent1.getOrder() == bondToParent2.getOrder()
                    && bondToParent1.isAromatic() == bondToParent2.isAromatic()
                    && bondToParent1.isInRing() == bondToParent2.isInRing()) {
                return true;
            }
        }
        
        return false;
    }
    
    /**
     *
     * @param ssc
     * @param HOSECodeLookupTable
     * @return
     * @deprecated 
     */
    public static Spectrum predictSpectrum(final SSC ssc, final HashMap<String, ArrayList<Double>> HOSECodeLookupTable){
        final Spectrum predictedSpectrum = new Spectrum(ssc.getSubspectrum().getNuclei());
        int atomIndexCounter = 0;
        for (final int i: ssc.getAtomTypeIndices().get(ssc.getSubspectrumAtomType())) {
            if(HOSECodeLookupTable.containsKey(ssc.getHOSECode(i))){
                predictedSpectrum.addSignal(new Signal(
                        ssc.getSubspectrum().getNuclei(), 
                        new Double[]{Utils.getRMS(HOSECodeLookupTable.get(ssc.getHOSECode(i)))},
                        ssc.getSubspectrum().getIntensity(atomIndexCounter), // RMS or median of intensities too? (when used)
                        ssc.getSubspectrum().getMultiplicity(atomIndexCounter)
                ));
            } else {
                return null;
            }
            atomIndexCounter++;
        }
        
        return predictedSpectrum;
    } 
    
    public static boolean isValidSpectrum(final Spectrum spectrum, final Spectrum querySpectrum, final double tol, final double thrsMatchFactor){
        if(spectrum == null){
            return false;
        }
        final Assignment matchAssignments = Assembly.matchSpectra(spectrum, querySpectrum, tol);
//        System.out.println("predicted spectrum:\t" + spectrum.getShifts(0));
//        System.out.println("predicted spectrum equ.:" + spectrum.getEquivalences());
//        System.out.println("query spectrum   :\t" + querySpectrum.getShifts(0));
//        System.out.println("query spectrum equ.:\t" + querySpectrum.getEquivalences());
//        System.out.println("matchAssignments:\t" + Arrays.toString(matchAssignments.getAtomIndices(0)));        
        // filter for unset assignments
        if (matchAssignments.getSetAssignmentsCount(0) < matchAssignments.getAssignmentsCount()) {
//            System.out.println("-> set assigments not allowed!!!");
            return false;
        }
        // filter for multiple assignments
        // there might be multiple assignments to same signals, so check for possible symmetry (equivalences)
        for (final int matchedSignalIndexInQuerySpectrum : matchAssignments.getAtomIndices(0)) {
            if (Collections.frequency(Utils.ArrayToArrayList(matchAssignments.getAtomIndices(0)), matchedSignalIndexInQuerySpectrum)
                    > querySpectrum.getEquivalentSignals(matchedSignalIndexInQuerySpectrum).size() + 1) { 
                // note: + 1 to query spectrum equ. signal count, because in frequencies count the requested (equ.) signal always occurs at least once and in query spectrum equ. signal count it could be zero
//                System.out.println("-> equivalences not allowed!!!");
                return false;
            }
        }
        // filter for intensities
        // -> still open
//        for (int i = 0; i < matchAssignments.getAssignmentsCount(); i++) {            
//            if(Math.abs(predictedSpectrum.getIntensity(i) - this.querySpectrum.getIntensity(matchAssignments.getAtomIndex(0, i))) > this.toleranceIntensityMatching){
//                return false;
//            }
//        }   
        // filter for match factor
        if (Assembly.getMatchFactor(spectrum, querySpectrum, tol)> thrsMatchFactor) {
//            System.out.println("-> match factor not allowed!!!");
            return false;
        }
        
        return true;
    }    
        
    /**
     * Returns the closest shift matches between two spectra as an 
     * Assignment object.
     * Despite intensities are expected, they are still not considered here.
     *
     * @param spectrum1
     * @param spectrum2
     * @param pickPrecision Tolerance value [ppm] used during peak peaking
     * @return Assignments with signal indices of spectrum1 and matched indices 
     * in spectrum2; assignments are unset (-1) if spectrum1 has a different nucleus
     * than spectrum2
     */
    public static Assignment matchSpectra(final Spectrum spectrum1, final Spectrum spectrum2, final double pickPrecision) {
        final Assignment matchAssignments = new Assignment(spectrum1);
        // first nuclei in both spectra are not the same
        if (!spectrum1.getNuclei()[0].equals(spectrum2.getNuclei()[0])) {
            return matchAssignments;
        }
        final HashSet<Integer> pickedSignalIndices = new HashSet<>();
        int pickedSignalIndexSpectrum2, pickedSignalIndexSpectrum2Prev;
        for (int i = 0; i < spectrum1.getSignalCount(); i++) {
            if(spectrum1.getShift(i, 0) == null){
                pickedSignalIndexSpectrum2 = -1;
            } else {           
                pickedSignalIndexSpectrum2 = spectrum2.pickClosestSignal(spectrum1.getShift(i, 0), 0, pickPrecision);
                // if matched signal is already assigned, then consider symmetries (equiv. signals)
                if (pickedSignalIndices.contains(pickedSignalIndexSpectrum2)) {
                    // symmetry exists
                    if (spectrum2.hasEquivalences(pickedSignalIndexSpectrum2)) {
                        pickedSignalIndexSpectrum2Prev = pickedSignalIndexSpectrum2;
                        // assign the next signal in equivalence list
                        for (final int equivalentSignalIndexSpectrum2 : spectrum2.getEquivalentSignals(pickedSignalIndexSpectrum2)) {
                            if(!pickedSignalIndices.contains(equivalentSignalIndexSpectrum2)){
                                pickedSignalIndexSpectrum2 = equivalentSignalIndexSpectrum2;
                                break;
                            }
                        }
                        // if no further equivalent signal exists then that match is not valid
                        if(pickedSignalIndexSpectrum2 == pickedSignalIndexSpectrum2Prev){
                            pickedSignalIndexSpectrum2 = -1;
                        }
                    } else {
                        // not symmetric signals but the same (predicted) or very similar shifts and multiple assignments to catch
                        // -> still open
                        pickedSignalIndexSpectrum2 = -1;
                    }
                }
                // check multiplicity 
                if((spectrum1.getMultiplicity(i) == null) 
                        || (spectrum2.getMultiplicity(pickedSignalIndexSpectrum2) == null)
                        || !spectrum1.getMultiplicity(i).equals(spectrum2.getMultiplicity(pickedSignalIndexSpectrum2))) {
                    pickedSignalIndexSpectrum2 = -1;
                }
            }     
            // add only truly assigned signal to list of already assigned signals
            if (pickedSignalIndexSpectrum2 != -1) {
                pickedSignalIndices.add(pickedSignalIndexSpectrum2);
            }
            // set picked signal index in assignment object
            matchAssignments.setAssignment(0, i, pickedSignalIndexSpectrum2);
        }   
        // try to assign the still unassigned shifts in spectrum1 to shifts in spectrum2
//        System.out.println("--> assignments before:\t" + Utils.ArrayToArrayList(matchAssignments.getAtomIndices(0)));
//        ArrayList<Integer> pickedSignalIndicesInSpectrum2;
//        for (int i = 0; i < matchAssignments.getAssignmentsCount(); i++) {
//            final Double queryShiftSpectrum1 = spectrum1.getShift(i, 0);
//            if ((matchAssignments.getAtomIndex(0, i) == -1) && (queryShiftSpectrum1 != null)) {
//                pickedSignalIndicesInSpectrum2 = spectrum2.pickSignals(queryShiftSpectrum1, 0, pickPrecision);
//                for (final int pickedSignalIndexInSpectrum2 : pickedSignalIndicesInSpectrum2) {
//                    if (!pickedSignalIndices.contains(pickedSignalIndexInSpectrum2)
//                            && (spectrum1.getMultiplicity(i) != null)
//                            && (spectrum2.getMultiplicity(pickedSignalIndexInSpectrum2) != null)
//                            && spectrum1.getMultiplicity(i).equals(spectrum2.getMultiplicity(pickedSignalIndexInSpectrum2))) {
//                        matchAssignments.setAssignment(0, i, pickedSignalIndexInSpectrum2);
//                        pickedSignalIndices.add(pickedSignalIndexInSpectrum2);
//                        break;
//                    }
//                }
//            }
//        }
//        System.out.println("--> assignments after:\t" + Utils.ArrayToArrayList(matchAssignments.getAtomIndices(0)));
        

        return matchAssignments;
    }    

    /**
     * Returns deviatons between matched shifts in SSC and query query spectrum.
     * The matching procedure is already included here.
     *
     * @param spectrum1 
     * @param spectrum2
     * @param pickPrecision
     * @return
     *
     * @see #matchSpectra(casekit.NMR.model.Spectrum, casekit.NMR.model.Spectrum, double) 
     */
    public static Double[] getDeviations(final Spectrum spectrum1, final Spectrum spectrum2, final double pickPrecision) {
        final Double[] deviations = new Double[spectrum1.getSignalCount()];
        final Assignment matchAssignments = Assembly.matchSpectra(spectrum1, spectrum2, pickPrecision);
        Signal matchedSignalInSpectrum2;
        for (int i = 0; i < spectrum1.getSignalCount(); i++) {
            if (matchAssignments.getAtomIndex(0, i) == -1) {
                deviations[i] = null;
            } else {
                matchedSignalInSpectrum2 = spectrum2.getSignal(matchAssignments.getAtomIndex(0, i));
                deviations[i] = Math.abs(spectrum1.getSignal(i).getShift(0) - matchedSignalInSpectrum2.getShift(0));
            }
        }

        return deviations;
    }
        
    /**
     * Returns the average of all deviations within a given input array.
     * If the minimum overlap threshold is not reached, a null value will be
     * returned.
     *
     * @param deviations Deviations
     * @return
     *
     * @see #getDeviations(model.SSC, casekit.NMR.model.Spectrum, double)
     */
    public static Double getMatchFactor(final Double[] deviations) {
        for (final Double deviation : deviations) {
            if (deviation == null) {
                return null;
            }
        }

        return Utils.getMean(deviations);
    }

    /**
     * Returns the average of all deviations of matched shifts between two 
     * spectra.
     * The calculation of deviations is already included here.
     *
     * @param spectrum1 
     * @param spectrum2
     * @param pickPrecision Tolerance value [ppm] used during peak picking in 
     * shift comparison
     * @return
     *
     * @see #getDeviations(casekit.NMR.model.Spectrum, casekit.NMR.model.Spectrum, double) 
     * @see casekit.NMR.Utils#getMatchFactor(java.lang.Double[], int) 
     */
    public static Double getMatchFactor(final Spectrum spectrum1, final Spectrum spectrum2, final double pickPrecision) {
        return Assembly.getMatchFactor(Assembly.getDeviations(spectrum1, spectrum2, pickPrecision));
    }
    
    /**
     * Combines two 1D spectra while considering possible equivalent signals 
     * via the pickPrecision parameter and multiplicity comparison.
     * 
     *
     * @param spectrum1 first spectrum
     * @param spectrum2 second spectrum
     * @param pickPrecision tolerance value used for signal shift matching to 
     * find equivalent signals
     * @return
     */
    public static Spectrum combineSpectra(final Spectrum spectrum1, final Spectrum spectrum2, final double pickPrecision){
        // check for same spectrum dimension count and same nuclei in both spectra
        if((spectrum1.getDimCount() != 1) 
                || (spectrum1.getDimCount() != spectrum2.getDimCount()) 
                || !spectrum1.getNuclei()[0].equals(spectrum2.getNuclei()[0])){
            return null;
        }           
//        // create new spectra which is to fill with signals of both spectra
//        final Spectrum combinedSpectrum = new Spectrum(spectrum1.getNuclei());        
//        // fill in signals from spectrum1
//        // consider the possibility of potential equivalent signals here
        int equivalentSignalIndex;
//        for (final Signal signalSpectrum1 : spectrum1.getSignals()) {
//            equivalentSignalIndex = -1;
//            for (final int closestSignalIndex : combinedSpectrum.pickSignals(signalSpectrum1.getShift(0), 0, pickPrecision)) {
//                if(signalSpectrum1.getMultiplicity().equals(combinedSpectrum.getSignal(closestSignalIndex).getMultiplicity())){
//                    equivalentSignalIndex = closestSignalIndex;
//                    break;
//                }
//            }
//            combinedSpectrum.addSignal(signalSpectrum1, equivalentSignalIndex);      
//        }
        
        // create new spectra which is to fill with signals of both spectra
        final Spectrum combinedSpectrum = spectrum1.getClone();
        // fill in signals from spectrum2
        // consider the possibility of potential equivalent signals here
        for (final Signal signalSpectrum2 : spectrum2.getSignals()) {
            equivalentSignalIndex = -1;
            for (final int closestSignalIndex : combinedSpectrum.pickSignals(signalSpectrum2.getShift(0), 0, pickPrecision)) {
                if (signalSpectrum2.getMultiplicity().equals(combinedSpectrum.getSignal(closestSignalIndex).getMultiplicity())) {
                    equivalentSignalIndex = closestSignalIndex;
                }
            }
            combinedSpectrum.addSignal(signalSpectrum2, equivalentSignalIndex);                        
        }        
        
        return combinedSpectrum;
    }
    
    public static HashMap<String, SSC> assemble(final int nStartSSCs, final int nThreads, final SSCRanker sscRanker, final int minMatchingSphereCount, 
            final Spectrum querySpectrum, final double shiftTol, final double thrsMatchFactor, final double pickPrecision) throws InterruptedException{
        final HashMap<String, SSC> solutions = new HashMap<>();
        // initialize an executor for parallelization
        final ExecutorService executor = Utils.initExecuter(nThreads);
        final ArrayList<Callable<HashMap<String, SSC>>> callables = new ArrayList<>();
        // add all task to do
        for (int i = 0; i < nStartSSCs; i++) {
            final int j = i;
            callables.add((Callable<HashMap<String, SSC>>) () -> {                                
                return Assembly.assemble(sscRanker, j, minMatchingSphereCount, querySpectrum, shiftTol, thrsMatchFactor, pickPrecision);
            });
        }
        // execute all task in parallel
        executor.invokeAll(callables)
                .stream()
                .map(future -> {
                    try {
                        return future.get();
                    } catch (InterruptedException | ExecutionException e) {
                        throw new IllegalStateException(e);
                    }
                })
                .forEach((tempHashMap) -> {                    
                    solutions.putAll(tempHashMap);
                });
        // shut down the executor service
        Utils.stopExecuter(executor, 5);
        
        return solutions;
    } 
    
    public static HashMap<String, SSC> assemble(final SSCRanker sscRanker, final int startSSCIndex, final int minMatchingSphereCount, 
            final Spectrum querySpectrum, final double shiftTol, final double thrsMatchFactor, final double pickPrecision) throws CloneNotSupportedException, CDKException, IOException {
        
        final SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Absolute);
        final SmilesParser smilesParser = new SmilesParser(new SilentChemObjectBuilder());
        String structureAsSMILES;
        final HashMap<String, SSC> solutions = new HashMap<>();
        final ArrayList<Integer> rankedSSCIndices = sscRanker.getRankedSSCIndices();
        SSC ssc2, ssc1Backup;
        ArrayList<Object> resultSSCExtension;
        // copy(!!!) the SSC1 contents only; don't use the object (reference) itself because of modifications
        SSC ssc1 = sscRanker.getSSCLibrary().getSSC(rankedSSCIndices.get(startSSCIndex)).getClone();
        
        int maxMatchingSphere, mappedAtomIndexSSC1, mappedAtomIndexSSC2;
        Signal signalToRemoveSSC2;
        Spectrum subspectrumSSC2Clone, combinedSpectrum;
        HashMap<Integer, Integer> mappedAtomIndices;
        HashSet<Integer> missingAtomIndicesSSC2;
        HashMap<Integer, ArrayList<Integer[]>> overlaps;
        
        
        for (int ssc2Index = 0; ssc2Index < rankedSSCIndices.size(); ssc2Index++) {
            if (ssc2Index == startSSCIndex) {
                continue;
            }
            ssc2 = sscRanker.getSSCLibrary().getSSC(rankedSSCIndices.get(ssc2Index));

            System.out.println("\n\n-------------------------------- " + startSSCIndex + ", " + ssc2Index + " --------------------------------");

            
            // 1. check for partial structural identity (overlaps); via HOSE code or connection tree comparison                                        
            overlaps = Assembly.getStructuralOverlaps(ssc2, ssc1, shiftTol, minMatchingSphereCount, ssc1.getSubspectrumAtomType());
            // if there is no structural identity then skip that SSC pair comparison            
            if (overlaps.isEmpty()) {
                continue;
            }
            
            // check whether the current SSC is already a final SSC
            if (Assembly.isFinalSSC(ssc1, querySpectrum)) {
                structureAsSMILES = smilesGenerator.create(ssc1.getSubstructure());
                if (!solutions.containsKey(structureAsSMILES)) {
                    solutions.put(structureAsSMILES, ssc1);
                    System.out.println("--> new solution found!!! -> " + solutions.size() + " -> " + structureAsSMILES);
                    System.out.println("-> atom count: " + ssc1.getAtomCount() + ", bond count: " + ssc1.getBondCount());
                    System.out.println("-> query spectrum:\t" + querySpectrum.getShifts(0));
                    System.out.println("-> equivalences:\t" + querySpectrum.getEquivalences());
                    System.out.println("-> pred. spectrum:\t" + ssc1.getSubspectrum().getShifts(0));
                    System.out.println("-> equivalences:\t" + ssc1.getSubspectrum().getEquivalences());
                    Utils.generatePicture(smilesParser.parseSmiles(structureAsSMILES), "/Users/mwenk/Downloads/outputs/out_final_" + structureAsSMILES + ".png");
                }
                break;
            }
            
            ssc1Backup = ssc1.getClone();
            subspectrumSSC2Clone = ssc2.getSubspectrum().getClone();                                      
            
            
            maxMatchingSphere = Collections.max(overlaps.keySet());
            System.out.println("-> maxMatchingSphere: " + maxMatchingSphere);
            
            // 1.1 map structural overlaps
            mappedAtomIndices = Assembly.mapStructuralOverlaps(ssc1, ssc2, overlaps, shiftTol);
            if (mappedAtomIndices.isEmpty()) {
                continue;
            }

            missingAtomIndicesSSC2 = Assembly.getMissingAtomIndices(ssc2, mappedAtomIndices.keySet());
            System.out.println(" -> mapped atom indices: \t" + mappedAtomIndices);
            System.out.println(" -> missing atom indices:\t" + missingAtomIndicesSSC2);
            if (missingAtomIndicesSSC2.isEmpty()) {
                continue;
            }

            // 1.2 remove signals from SSC2's subspectrum associated with mapping (overlapping) atoms; for subspectra combination
            for (final Map.Entry<Integer, Integer> entry : mappedAtomIndices.entrySet()) {
                mappedAtomIndexSSC2 = entry.getKey();
                if (ssc2.getSubstructure().getAtom(mappedAtomIndexSSC2).getSymbol().equals(ssc1.getSubspectrumAtomType())) {
                    signalToRemoveSSC2 = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getSignalIndex(0, mappedAtomIndexSSC2));
                    subspectrumSSC2Clone.removeSignal(signalToRemoveSSC2);
                }
            }

            // 2. combine subspectra of SSC to extend and matched SSC and validate it
            combinedSpectrum = Assembly.combineSpectra(ssc1.getSubspectrum(), subspectrumSSC2Clone, pickPrecision);
//            System.out.println("\nspectrum1:\t" + ssc1.getSubspectrum().getShifts(0));
//            System.out.println("spectrum2:\t" + ssc2.getSubspectrum().getShifts(0));
//            System.out.println("spectrum2a:\t" + subspectrumSSC2Clone.getShifts(0));
//            System.out.println("spectrum3:\t" + combinedSpectrum.getShifts(0));
//            System.out.println("spectrum3 equ:\t" + combinedSpectrum.getEquivalences());

            // if no valid spectrum could be built then go to next pairwise SSC comparison
            if ((combinedSpectrum == null) || !Assembly.isValidSpectrum(combinedSpectrum, querySpectrum, shiftTol, thrsMatchFactor)) {
//                System.out.println("-> no valid combined spectrum!");
                continue;
            }
            System.out.println("-> !!!valid combined spectrum!!!");

            // 3. substructure extension in SSC1
//            System.out.println("entering substructure extension!!!");
            resultSSCExtension = Assembly.extendSSC(ssc1, ssc2, mappedAtomIndices, pickPrecision);
            if (resultSSCExtension == null) {
                System.out.println("---> could not extend ssc1 -> set to backup!");
                ssc1 = ssc1Backup.getClone();
                continue;
            }
            ssc1 = (SSC) resultSSCExtension.get(0);
            mappedAtomIndices = (HashMap<Integer, Integer>) resultSSCExtension.get(1);

            if (mappedAtomIndices.isEmpty()) {
                ssc1 = ssc1Backup.getClone();
                continue;
            }
            
//            // update overlaps
//            overlaps = Assembly.getStructuralOverlaps(ssc2, ssc1, shiftTol, minMatchingSphereCount, maxMatchingSphereCount);
//            // map structural overlaps again
//            mappedAtomIndices = Assembly.mapStructuralOverlaps(ssc1, ssc2, overlaps, shiftTol);
                        
            missingAtomIndicesSSC2 = Assembly.getMissingAtomIndices(ssc2, mappedAtomIndices.keySet());
            System.out.println(" -> mapped atoms indices  at the end:\t" + mappedAtomIndices);
            System.out.println(" -> missing atom indices  at the end:\t" + missingAtomIndicesSSC2);
            if (!missingAtomIndicesSSC2.isEmpty()) {
//                System.out.println("---> could not add all missing atoms!");
                ssc1 = ssc1Backup.getClone();
                continue;
            }
            
            
//            Utils.generatePicture(smilesParser.parseSmiles(smilesGenerator.create(ssc1.getSubstructure())), "/Users/mwenk/Downloads/outputs/temp_" + startSSCIndex + "_" + ssc2Index + ".png");
            Utils.generatePicture(ssc1.getSubstructure(), "/Users/mwenk/Downloads/outputs/temp_" + startSSCIndex + "_" + ssc2Index + ".png");
            
            if(Assembly.isFinalSSC(ssc1, querySpectrum)){
                structureAsSMILES = smilesGenerator.create(ssc1.getSubstructure());
                if (!solutions.containsKey(structureAsSMILES)) {
                    solutions.put(structureAsSMILES, ssc1);
                    System.out.println("--> new solution found!!! -> " + solutions.size() + " -> " + structureAsSMILES);
                    System.out.println("-> atom count: " + ssc1.getAtomCount() + ", bond count: " + ssc1.getBondCount());
                    System.out.println("-> query spectrum:\t" + querySpectrum.getShifts(0));
                    System.out.println("-> equivalences:\t" + querySpectrum.getEquivalences());
                    System.out.println("-> pred. spectrum:\t" + ssc1.getSubspectrum().getShifts(0));
                    System.out.println("-> equivalences:\t" + ssc1.getSubspectrum().getEquivalences());
                    Utils.generatePicture(smilesParser.parseSmiles(structureAsSMILES), "/Users/mwenk/Downloads/outputs/out_final_" + structureAsSMILES + ".png");
                }
                break;
            }
            
            // in case there are two hetero atoms the same (added two times) and unsaturated, then one could remove the first atom and link it to the neighbor(s) atom(s) of the second one
            // this step should not be necessary anymore if one use the molecula formula of the unknown compound to elucidate
            if((ssc1.getSubspectrum().getSignalCount() == querySpectrum.getSignalCount())){
                System.out.println("query spectrum covered but unsaturated atoms are left: " + ssc1.getUnsaturatedAtomIndices());
                System.out.println(" -> unsaturated atoms: " + ssc1.getUnsaturatedAtomIndices());
                if(ssc1.getUnsaturatedAtomIndices().size() == 2){
                    int unsaturatedAtomIndex1 = ssc1.getUnsaturatedAtomIndices().get(0);
                    int unsaturatedAtomIndex2 = ssc1.getUnsaturatedAtomIndices().get(1);                    
                    ConnectionTreeNode node1 = ssc1.getConnectionTree(unsaturatedAtomIndex1).getRootNode();
                    ConnectionTreeNode node2 = ssc1.getConnectionTree(unsaturatedAtomIndex2).getRootNode();
                    Signal signal1 = ssc1.getSubspectrum().getSignal(ssc1.getAssignments().getSignalIndex(0, unsaturatedAtomIndex1));
                    Signal signal2 = ssc1.getSubspectrum().getSignal(ssc1.getAssignments().getSignalIndex(0, unsaturatedAtomIndex2));
                    if(Assembly.isEqualNode(node1, node2, signal1, signal2, shiftTol)){
                        System.out.println(" -> link!!!");
                        IAtom unsaturatedAtom1 = ssc1.getSubstructure().getAtom(unsaturatedAtomIndex1);
                        IAtom unsaturatedAtom2 = ssc1.getSubstructure().getAtom(unsaturatedAtomIndex2);
                        for (final IBond bond2 : ssc1.getSubstructure().getConnectedBondsList(unsaturatedAtom2)) {
                            if(isValidBondAddition(ssc1.getSubstructure(), unsaturatedAtomIndex1, bond2)){
                                bond2.setAtom(bond2.getOther(unsaturatedAtom2), 0);   
                                bond2.setAtom(unsaturatedAtom1, 1);
                            }
                        }
                        ssc1.getSubstructure().removeAtom(unsaturatedAtom2);
                        ssc1.update();
                        
                        if (Assembly.isFinalSSC(ssc1, querySpectrum)) {
                            structureAsSMILES = smilesGenerator.create(ssc1.getSubstructure());
                            if (!solutions.containsKey(structureAsSMILES)) {
                                solutions.put(structureAsSMILES, ssc1);
                                System.out.println("--> new solution found!!! -> " + solutions.size() + " -> " + structureAsSMILES);
                                System.out.println("-> atom count: " + ssc1.getAtomCount() + ", bond count: " + ssc1.getBondCount());
                                System.out.println("-> query spectrum:\t" + querySpectrum.getShifts(0));
                                System.out.println("-> equivalences:\t" + querySpectrum.getEquivalences());
                                System.out.println("-> pred. spectrum:\t" + ssc1.getSubspectrum().getShifts(0));
                                System.out.println("-> equivalences:\t" + ssc1.getSubspectrum().getEquivalences());
                                Utils.generatePicture(smilesParser.parseSmiles(structureAsSMILES), "/Users/mwenk/Downloads/outputs/out_final_" + structureAsSMILES + ".png");
                            }
                            break;
                        }
                    }
                }
            }
        }
        
        return solutions;
    }
    
    public static boolean isFinalSSC(final SSC ssc, final Spectrum querySpectrum){
//        System.out.println("\n\nno more unsaturated atoms left? -> " + !ssc.hasUnsaturatedAtoms());
//        System.out.println("query spectrum size reached? -> " + ssc.getSubspectrum().getSignalCount() + " == " + querySpectrum.getSignalCount() + " -> " + (ssc.getSubspectrum().getSignalCount() == querySpectrum.getSignalCount()));
        return !ssc.hasUnsaturatedAtoms() && (ssc.getSubspectrum().getSignalCount() == querySpectrum.getSignalCount());
    }
    
    public static HashSet<Integer> getMissingAtomIndices(final SSC ssc, final Set<Integer> mappedAtomIndices){
        final HashSet<Integer> missingAtomIndices = new HashSet<>();
        for (int k = 0; k < ssc.getAtomCount(); k++) {
            if (!mappedAtomIndices.contains(k)) {
                missingAtomIndices.add(k);
            }
        }
        
        return missingAtomIndices;
    }
    
    public static HashMap<Integer, Integer> mapStructuralOverlaps(final SSC ssc1, final SSC ssc2, final HashMap<Integer, ArrayList<Integer[]>> overlaps, final double shiftTol) throws CDKException, CloneNotSupportedException{
        
        final HashMap<Integer, Integer> mappedAtomIndices = new HashMap<>();
        HashMap<Integer, Integer> mappedEqualNodesKeys;
        int overlapAtomIndexSSC1, overlapAtomIndexSSC2, mappedAtomIndexSSC1, mappedAtomIndexSSC2;
        boolean containsUnsaturatedAtomsSSC1;
        
        for (int m = Collections.max(overlaps.keySet()); m >= Collections.min(overlaps.keySet()); m--) {
            if (!overlaps.containsKey(m)) {
                continue;
            }
                        
            System.out.println("\n -> m: " + m);
            for (final Integer[] overlap : overlaps.get(m)) {
                overlapAtomIndexSSC2 = overlap[0];
                overlapAtomIndexSSC1 = overlap[1];
                System.out.println("\n--> overlap: " + overlapAtomIndexSSC2 + ", " + overlapAtomIndexSSC1);                
                System.out.println("HOSE code ssc1:\t" + HOSECodeBuilder.buildHOSECode(ssc1.getConnectionTree(overlapAtomIndexSSC1), false));
                System.out.println("HOSE code ssc2:\t" + HOSECodeBuilder.buildHOSECode(ssc2.getConnectionTree(overlapAtomIndexSSC2), false));                

                containsUnsaturatedAtomsSSC1 = false;
                for (final int nodeKey : ssc1.getConnectionTree(overlapAtomIndexSSC1).getKeysInOrder(true)) {                   
                    if(ssc1.isUnsaturated(nodeKey)){
                        containsUnsaturatedAtomsSSC1 = true;
                        break;
                    }
                }
                
                if (!containsUnsaturatedAtomsSSC1) {
                    System.out.println("\nnothing unsaturated -> skip this overlap");
                    continue;
                }                                

                // add mappings of all atoms between SSC1 and SSC2 which have an identical structural overlap
                for (int s = 0; s <= m; s++) {
                    for (int i = 0; i < ssc2.getConnectionTree(overlapAtomIndexSSC2).getNodesCountInSphere(s); i++) {
                        mappedAtomIndexSSC2 = ssc2.getConnectionTree(overlapAtomIndexSSC2).getNodeKeysInSphere(s).get(i);
                        mappedAtomIndexSSC1 = ssc1.getConnectionTree(overlapAtomIndexSSC1).getNodeKeysInSphere(s).get(i);
                        if (ssc2.getConnectionTree(overlapAtomIndexSSC2).getNode(mappedAtomIndexSSC2).isRingClosureNode()) {
                            continue;
                        }
                        if (mappedAtomIndexSSC1 == -1) {
                            continue;
                        }
                        if (mappedAtomIndices.containsKey(mappedAtomIndexSSC2) && (mappedAtomIndices.get(mappedAtomIndexSSC2) != -1)) {
//                            if (mappedAtomIndices.get(mappedAtomIndexSSC2) != mappedAtomIndexSSC1) {
//                                System.out.println(" in s: " + s + " -> !!!tried to set mappedAtomIndexSSC1 more than one time!!! -> "
//                                        + mappedAtomIndexSSC2 + " : " + mappedAtomIndices.get(mappedAtomIndexSSC2) + " vs. " + mappedAtomIndexSSC1);
//                            }
                            continue;
                        }
                        if (mappedAtomIndices.containsValue(mappedAtomIndexSSC1)) {
//                            System.out.println(" in s: " + s + " -> !!!tried to set mappedAtomIndexSSC2 more than one time!!! -> "
//                                    + mappedAtomIndexSSC2 + " : " + mappedAtomIndexSSC1);
                            continue;
                        }
//                        System.out.println(" in s: " + s + " -> new map: " + mappedAtomIndexSSC2 + ", " + mappedAtomIndexSSC1);
                        mappedAtomIndices.put(mappedAtomIndexSSC2, mappedAtomIndexSSC1);
                    }
                }
                // add mappings of all atoms between SSC1 and SSC2 which have to be assigned to each other (not structurally identical)
                for (int s = m + 1; s <= Integer.min(ssc2.getConnectionTree(overlapAtomIndexSSC2).getMaxSphere(), ssc1.getConnectionTree(overlapAtomIndexSSC1).getMaxSphere()); s++) {
                    mappedEqualNodesKeys = Assembly.mapEqualNodesInSphere(ssc2, ssc1, overlapAtomIndexSSC2, overlapAtomIndexSSC1, s, shiftTol);
//                    System.out.println(" -> in s: " + s + " -> mapped equal node:  " + mappedEqualNodesKeys);
                    for (final Map.Entry<Integer, Integer> entry : mappedEqualNodesKeys.entrySet()) {
                        mappedAtomIndexSSC2 = entry.getKey();
                        if (ssc2.getConnectionTree(overlapAtomIndexSSC2).getNode(mappedAtomIndexSSC2).isRingClosureNode()) {
                            continue;
                        }
                        mappedAtomIndexSSC1 = entry.getValue();
                        if (mappedAtomIndexSSC1 == -1) {
                            continue;
                        }
                        if (mappedAtomIndices.containsKey(mappedAtomIndexSSC2) && (mappedAtomIndices.get(mappedAtomIndexSSC2) != -1)) {
//                            if (mappedAtomIndices.get(mappedAtomIndexSSC2) != mappedAtomIndexSSC1) {
//                                System.out.println(" in s: " + s + " -> !!!tried to set mappedAtomIndexSSC1 more than one time!!! -> "
//                                        + mappedAtomIndexSSC2 + " : " + mappedAtomIndices.get(mappedAtomIndexSSC2) + " vs. " + mappedAtomIndexSSC1);
//                            }
                            continue;
                        }
                        if (mappedAtomIndices.containsValue(mappedAtomIndexSSC1)) {
//                            System.out.println(" in s: " + s + " -> !!!tried to set mappedAtomIndexSSC2 more than one time!!! -> "
//                                    + mappedAtomIndexSSC2 + " : " + mappedAtomIndexSSC1);
                            continue;
                        }
//                        System.out.println(" in s: " + s + " -> new map: " + mappedAtomIndexSSC2 + ", " + mappedAtomIndexSSC1);
                        mappedAtomIndices.put(mappedAtomIndexSSC2, mappedAtomIndexSSC1);
                    }
                }
            }
            System.out.println(" -> mapped atom indices after m: " + m + ":\t" + mappedAtomIndices);
            if (!mappedAtomIndices.isEmpty()) {
                System.out.println(" -> not empty");
                break;
            }
        }
        
        return mappedAtomIndices;
    }
    
    private static ArrayList<Object> extendSSC(final SSC ssc1, final SSC ssc2, final HashMap<Integer, Integer> mappedAtomIndices, final double pickPrecision) throws CDKException, CloneNotSupportedException{
        int mappedAtomIndexSSC1, mappedAtomIndexSSC2, equivalentSignalIndex;
        IBond bondToAdd;
        IAtom atomToAdd, parentAtomSSC1;
        Signal signalToAddSSC2;
        ConnectionTree connectionTreeToAddSSC2;
        final HashMap<Integer, Integer> mappedAtomIndices2 = new HashMap<>(mappedAtomIndices);
        
        for (final Map.Entry<Integer, Integer> entry : mappedAtomIndices.entrySet()) {
            mappedAtomIndexSSC2 = entry.getKey();
            mappedAtomIndexSSC1 = entry.getValue();
            if (!ssc1.isUnsaturated(mappedAtomIndexSSC1)) {
                continue;
            }
            System.out.println("\nmapping: " + mappedAtomIndexSSC2 + ", " + mappedAtomIndexSSC1
                    + "\n-> " + ssc1.getHOSECode(mappedAtomIndexSSC1) + "\n-> " + ssc2.getHOSECode(mappedAtomIndexSSC2));

            // BFS to build connection tree which contains atoms in SSC2 to add to SSC1                                                
            connectionTreeToAddSSC2 = HOSECodeBuilder.buildConnectionTree(ssc2.getSubstructure(), mappedAtomIndexSSC2, null, new HashSet<>(mappedAtomIndices2.keySet()));
            System.out.println(" -> BFS: " + connectionTreeToAddSSC2.toString()
                    + "\n -> maxSphere: " + connectionTreeToAddSSC2.getMaxSphere());

            for (int s = 1; s <= connectionTreeToAddSSC2.getMaxSphere(); s++) {
                for (final ConnectionTreeNode connectedNodeInSphereToAddSSC2 : connectionTreeToAddSSC2.getNodesInSphere(s)) {
                    if (connectedNodeInSphereToAddSSC2.isRingClosureNode()) {
//                        System.out.println("ringClosureNode at: " + connectedNodeInSphereToAddSSC2.getKey());

                        ConnectionTreeNode parentNodeSSC2 = connectedNodeInSphereToAddSSC2.getParentNodes().get(0);
                        ConnectionTreeNode ringClosureParentNodeSSC2 = connectedNodeInSphereToAddSSC2.getParentNodes().get(0).getParentNodes().get(1);
                        if (!mappedAtomIndices2.containsKey(parentNodeSSC2.getKey())
                                || !mappedAtomIndices2.containsKey(ringClosureParentNodeSSC2.getKey())) {
                            continue;
                        }

                        bondToAdd = parentNodeSSC2.getBondsToParents().get(1).clone();
                        if (Assembly.isValidBondAddition(ssc1.getSubstructure(), mappedAtomIndices2.get(parentNodeSSC2.getKey()), bondToAdd)
                                && Assembly.isValidBondAddition(ssc1.getSubstructure(), mappedAtomIndices2.get(ringClosureParentNodeSSC2.getKey()), bondToAdd)) {
                            bondToAdd.setAtom(ssc1.getSubstructure().getAtom(mappedAtomIndices2.get(parentNodeSSC2.getKey())), 0);
                            bondToAdd.setAtom(ssc1.getSubstructure().getAtom(mappedAtomIndices2.get(ringClosureParentNodeSSC2.getKey())), 1);
                            ssc1.getSubstructure().addBond(bondToAdd);
                        }

                        continue;
                    }
                    bondToAdd = connectedNodeInSphereToAddSSC2.getBondsToParents().get(0).clone();
                    System.out.println(" -> key: " + connectedNodeInSphereToAddSSC2.getParentNodes().get(0).getKey() + " -> " + mappedAtomIndices2.get(connectedNodeInSphereToAddSSC2.getParentNodes().get(0).getKey()));
                    parentAtomSSC1 = ssc1.getSubstructure().getAtom(mappedAtomIndices2.get(connectedNodeInSphereToAddSSC2.getParentNodes().get(0).getKey()));
                    System.out.println("in s: " + s + " -> bond: (" + connectedNodeInSphereToAddSSC2.getKey() + ") to " + mappedAtomIndices2.get(connectedNodeInSphereToAddSSC2.getParentNodes().get(0).getKey()));

                    if (Assembly.isValidBondAddition(ssc1.getSubstructure(), parentAtomSSC1.getIndex(), bondToAdd)) {

                        atomToAdd = connectedNodeInSphereToAddSSC2.getAtom().clone();
                        ssc1.getSubstructure().addAtom(atomToAdd);

                        bondToAdd.setAtom(atomToAdd, 0);
                        bondToAdd.setAtom(parentAtomSSC1, 1);
                        ssc1.getSubstructure().addBond(bondToAdd);
                        
                        // add belonging signal from SSC2 to SSC1 
                        if (atomToAdd.getSymbol().equals(ssc1.getSubspectrumAtomType())) {
                            signalToAddSSC2 = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getSignalIndex(0, connectedNodeInSphereToAddSSC2.getKey()));
                            if (signalToAddSSC2 == null) {
                                ssc1.getSubstructure().removeBond(bondToAdd);
                                ssc1.getSubstructure().removeAtom(atomToAdd);
                               
                                return null;
                            }
                            equivalentSignalIndex = ssc1.getSubspectrum().pickClosestSignal(signalToAddSSC2.getShift(0), 0, pickPrecision);
                            ssc1.getSubspectrum().addSignal(signalToAddSSC2, equivalentSignalIndex);
                            ssc1.getAssignments().addAssignment(new int[]{ssc1.getAtomCount() - 1});
                        }

                        mappedAtomIndices2.put(connectedNodeInSphereToAddSSC2.getKey(), ssc1.getAtomCount() - 1);
//                        System.out.println("-> added atom " + connectedNodeInSphereToAddSSC2.getKey() + " from SSC2 to SSC1");

                        // ring closure
                        IAtom atomSSC2 = ssc2.getSubstructure().getAtom(connectedNodeInSphereToAddSSC2.getKey());
                        for (final IAtom connectedAtomSSC2 : ssc2.getSubstructure().getConnectedAtomsList(atomSSC2)) {
                            if ((connectedAtomSSC2.getIndex() != connectedNodeInSphereToAddSSC2.getParentNodes().get(0).getKey())
                                    && mappedAtomIndices2.containsKey(connectedAtomSSC2.getIndex())
                                    && !connectionTreeToAddSSC2.containsKey(connectedAtomSSC2.getIndex())) {
//                                System.out.println("connected but mapped atom found -> ring closure!!!"
//                                        + "\n -> " + connectedNodeInSphereToAddSSC2.getKey() + ", " + connectedAtomSSC2.getIndex());
                                bondToAdd = ssc2.getSubstructure().getBond(atomSSC2, connectedAtomSSC2).clone();
                                if (Assembly.isValidBondAddition(ssc1.getSubstructure(), ssc1.getAtomCount() - 1, bondToAdd)
                                        && Assembly.isValidBondAddition(ssc1.getSubstructure(), mappedAtomIndices2.get(connectedAtomSSC2.getIndex()), bondToAdd)) {
//                                    System.out.println("!!! ring closure !!!");
                                    bondToAdd.setAtom(ssc1.getSubstructure().getAtom(ssc1.getAtomCount() - 1), 0);
                                    bondToAdd.setAtom(ssc1.getSubstructure().getAtom(mappedAtomIndices2.get(connectedAtomSSC2.getIndex())), 1);
                                    ssc1.getSubstructure().addBond(bondToAdd);
                                }
                            }
                        }
                    } else {
                        System.out.println("not a valid bond addition!!!");
                        return null;
                    }
                }                
            }
            ssc1.update();
        }
        
//        System.out.println("\n\n -> mapped atom indices at the end: \t" + mappedAtomIndices2);
        
        final ArrayList<Object> toReturn = new ArrayList<>();
        toReturn.add(ssc1);
        toReturn.add(mappedAtomIndices2);
        return toReturn;
    }    
    
}
