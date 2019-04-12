/*
 * The MIT License
 *
 * Copyright 2019 Michael Wenk [https://github.com/michaelwenk].
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
package match;

import hose.HOSECodeBuilder;
import casekit.NMR.Utils;
import casekit.NMR.model.Assignment;
import casekit.NMR.model.Signal;
import casekit.NMR.model.Spectrum;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import model.ConnectionTree;
import model.ConnectionTreeNode;
import model.SSC;
import org.apache.commons.lang3.ArrayUtils;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.similarity.Tanimoto;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Match {
    
//    /**
//     * Checks whether connectionTree1 contains the exactly same node or bond
//     * (at ring closure node) properties as in connectionTree2 in a specific
//     * sphere.
//     *
//     * @param ssc1
//     * @param ssc2
//     * @param atomIndexInSubstructure1
//     * @param atomIndexInSubstructure2
//     * @param sphere
//     * @param shiftTol
//     * @return
//     */
//    public static boolean hasEqualNodesInSphere(final SSC ssc1, final SSC ssc2, final int atomIndexInSubstructure1, final int atomIndexInSubstructure2, final int sphere, final double shiftTol) {
//        final ArrayList<ConnectionTreeNode> nodesInSphere1 = ssc1.getConnectionTree(atomIndexInSubstructure1).getNodesInSphere(sphere);
//        final ArrayList<ConnectionTreeNode> nodesInSphere2 = ssc2.getConnectionTree(atomIndexInSubstructure2).getNodesInSphere(sphere);
//        if (nodesInSphere1.size() != nodesInSphere2.size()) {
//            return false;
//        }
//        Signal signal1;
//        Signal signal2;
//        for (int j = 0; j < nodesInSphere1.size(); j++) {
//            signal1 = ssc1.getSubspectrum().getSignal(ssc1.getAssignments().getSignalIndex(0, nodesInSphere1.get(j).getKey()));
//            signal2 = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getSignalIndex(0, nodesInSphere2.get(j).getKey()));
//            if (!isEqualNode(nodesInSphere1.get(j), nodesInSphere2.get(j), signal1, signal2, shiftTol)) {
//                return false;
//            }
//        }
//        return true;
//    }

    /**
     * Combines two 1D spectra while considering possible equivalent signals
     * via the {@code shiftTol} parameter and multiplicity comparison.
     * In {@code spectrum1}, the equivalent signals have to be set.
     *
     *
     * @param spectrum1 first spectrum
     * @param spectrum2 second spectrum
     * @param pickPrecision tolerance value used for signal shift matching to
     * find equivalent signals
     * @return null if both spectra have not the same dimension or nuclei
     */
    public static Spectrum combineSpectra(final Spectrum spectrum1, final Spectrum spectrum2, final double pickPrecision) {
        // check for same spectrum dimension count and same nuclei in both spectra
        if ((spectrum1.getDimCount() != 1) || (spectrum1.getDimCount() != spectrum2.getDimCount()) || !spectrum1.getNuclei()[0].equals(spectrum2.getNuclei()[0])) {
            return null;
        }
        int equivalentSignalIndex;
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
            combinedSpectrum.addSignal(signalSpectrum2.getClone(), equivalentSignalIndex);
        }
        return combinedSpectrum;
    }
    
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
     */
    public static HashMap<Integer, Integer> mapEqualNodesInSphere(final SSC ssc1, final SSC ssc2, final int atomIndexInSubstructure1, final int atomIndexInSubstructure2, final int sphere, final double shiftTol) {
        final ConnectionTree connectionTree1 = ssc1.getConnectionTree(atomIndexInSubstructure1);
        final ConnectionTree connectionTree2 = ssc2.getConnectionTree(atomIndexInSubstructure2);
        if (sphere > Integer.min(connectionTree1.getMaxSphere(), connectionTree2.getMaxSphere())) {
            return null;
        }
        final HashMap<Integer, Integer> map = new HashMap<>();
        Signal signal1;
        Signal signal2;
        // add matched node indices in connectionTree2 from nodes list in connectionTree1 in given sphere
        for (final ConnectionTreeNode nodeInSphere1 : connectionTree1.getNodesInSphere(sphere)) {
            // search for nodes from conn. tree 2 in conn. tree 1; if found then add it to map
            for (final ConnectionTreeNode nodeInSphere2 : connectionTree2.getNodesInSphere(sphere)) {
                signal1 = ssc1.getSubspectrum().getSignal(ssc1.getAssignments().getSignalIndex(0, nodeInSphere1.getKey()));
                signal2 = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getSignalIndex(0, nodeInSphere2.getKey()));
                if (!map.containsValue(nodeInSphere2.getKey()) && Match.isEqualNode(nodeInSphere1, nodeInSphere2, signal1, signal2, shiftTol)) {
                    map.put(nodeInSphere1.getKey(), nodeInSphere2.getKey());
                    break;
                }
            }
            if (!map.containsKey(nodeInSphere1.getKey())) {
                map.put(nodeInSphere1.getKey(), -1);
            }
        }
        return map;
    }


    /**
     * Returns the average of all deviations within a given input array.
     *
     * @param deviations Deviations
     * @return
     *
     * @see #getDeviations(casekit.NMR.model.Spectrum, casekit.NMR.model.Spectrum, double) 
     */
    public static Double calculateMatchFactor(final Double[] deviations) {
        for (final Double deviation : deviations) {
            if (deviation == null) {
                return null;
            }
        }
        
        return Utils.getMean(deviations);        
    }
    
    public static Float calculateTanimotoCoefficient(final Spectrum spectrum1, final Spectrum spectrum2, final int dim) throws CDKException {
        if(spectrum1.getSignalCount() != spectrum2.getSignalCount()){
            return null;
        }
        final double[] shiftsSpectrum1 = ArrayUtils.toPrimitive(spectrum1.getShifts(dim).toArray(new Double[spectrum1.getSignalCount()]));
        Arrays.parallelSort(shiftsSpectrum1);
        final double[] shiftsSpectrum2 = ArrayUtils.toPrimitive(spectrum2.getShifts(dim).toArray(new Double[spectrum2.getSignalCount()]));
        Arrays.parallelSort(shiftsSpectrum2);

        return Tanimoto.calculate(shiftsSpectrum1, shiftsSpectrum2);
    }

    /**
     * Returns the average of all deviations of matched shifts between two
     * spectra.
     * The calculation of deviations is already included here.
     *
     * @param spectrum1
     * @param spectrum2
     * @param shiftTol Tolerance value [ppm] used during peak picking in
     * shift comparison
     * @return
     *
     * @see #getDeviations(casekit.NMR.model.Spectrum, casekit.NMR.model.Spectrum, double)
     * @see #calculateMatchFactor(java.lang.Double[]) 
     */
    public static Double calculateMatchFactor(final Spectrum spectrum1, final Spectrum spectrum2, final double shiftTol) {
        return Match.calculateMatchFactor(Match.getDeviations(spectrum1, spectrum2, shiftTol));
    }

//    public static int getMaximumMatchingSphere(final SSC ssc1, final SSC ssc2, final int atomIndexInSubstructure1, final int atomIndexInSubstructure2, final double shiftTol) throws CloneNotSupportedException {
//        int maxMatchingSphere = -1;
//        for (int s = 0; s <= Integer.min(ssc1.getConnectionTree(atomIndexInSubstructure1).getMaxSphere(), ssc2.getConnectionTree(atomIndexInSubstructure2).getMaxSphere()); s++) {
//            //            if(!Assembly.hasStructuralIdentity(ssc1, ssc2, atomIndexInSubstructure1, atomIndexInSubstructure2, s, shiftTol)){
//            if (!Match.hasEqualNodesInSphere(ssc1, ssc2, atomIndexInSubstructure1, atomIndexInSubstructure2, s, shiftTol)) {
//                break;
//            }
//            maxMatchingSphere = s;
//        }
//        return maxMatchingSphere;
//    }

    /**
     * Returns deviatons between matched shifts in SSC and query query spectrum.
     * The matching procedure is already included here.
     *
     * @param spectrum1
     * @param spectrum2
     * @param shiftTol
     * @return
     *
     * @see #matchSpectra(casekit.NMR.model.Spectrum, casekit.NMR.model.Spectrum, double)
     */
    public static Double[] getDeviations(final Spectrum spectrum1, final Spectrum spectrum2, final double shiftTol) {
        final Double[] deviations = new Double[spectrum1.getSignalCount()];
        final Assignment matchAssignments = Match.matchSpectra(spectrum1, spectrum2, shiftTol);
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
     * Returns the closest shift matches between two spectra as an
     * Assignment object.
     * Despite intensities are expected, they are still not considered here.
     *
     * @param spectrum1
     * @param spectrum2
     * @param shiftTol Tolerance value [ppm] used during spectra shift 
     * comparison
     * @return Assignments with signal indices of spectrum1 and matched indices
     * in spectrum2; assignments are unset (-1) if spectrum1 has a different nucleus
     * than spectrum2
     */
    public static Assignment matchSpectra(final Spectrum spectrum1, final Spectrum spectrum2, final double shiftTol) {
        final Assignment matchAssignments = new Assignment(spectrum1);
        // first nuclei in both spectra are not the same
        if (!spectrum1.getNuclei()[0].equals(spectrum2.getNuclei()[0])) {
            return matchAssignments;
        }
        final HashSet<Integer> pickedSignalIndices = new HashSet<>();
        int pickedSignalIndexSpectrum2;
        int pickedSignalIndexSpectrum2Prev;
        for (int i = 0; i < spectrum1.getSignalCount(); i++) {
            if (spectrum1.getShift(i, 0) == null) {
                pickedSignalIndexSpectrum2 = -1;
            } else {
                pickedSignalIndexSpectrum2 = spectrum2.pickClosestSignal(spectrum1.getShift(i, 0), 0, shiftTol);
                // if matched signal is already assigned, then consider symmetries (equiv. signals)
                if (pickedSignalIndices.contains(pickedSignalIndexSpectrum2)) {
                    // symmetry exists
                    if (spectrum2.hasEquivalences(pickedSignalIndexSpectrum2)) {
                        pickedSignalIndexSpectrum2Prev = pickedSignalIndexSpectrum2;
                        // assign the next signal in equivalence list
                        for (final int equivalentSignalIndexSpectrum2 : spectrum2.getEquivalentSignals(pickedSignalIndexSpectrum2)) {
                            if (!pickedSignalIndices.contains(equivalentSignalIndexSpectrum2)) {
                                pickedSignalIndexSpectrum2 = equivalentSignalIndexSpectrum2;
                                break;
                            }
                        }
                        // if no further equivalent signal exists then that match is not valid
                        if (pickedSignalIndexSpectrum2 == pickedSignalIndexSpectrum2Prev) {
                            pickedSignalIndexSpectrum2 = -1;
                        }
                    } else {
                        // not symmetric signals but the same (predicted) or very similar shifts and multiple assignments to catch
                        // -> still open
                        pickedSignalIndexSpectrum2 = -1;
                    }
                }
                // check multiplicity
                if ((spectrum1.getMultiplicity(i) == null) || (spectrum2.getMultiplicity(pickedSignalIndexSpectrum2) == null) || !spectrum1.getMultiplicity(i).equals(spectrum2.getMultiplicity(pickedSignalIndexSpectrum2))) {
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
        //                pickedSignalIndicesInSpectrum2 = spectrum2.pickSignals(queryShiftSpectrum1, 0, shiftTol);
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
     * Returns whether two nodes have certain identical properties.
     *
     * @param nodeInSphere1
     * @param nodeInSphere2
     * @param signal1
     * @param signal2
     * @param shiftTol
     * @return
     */
    public static boolean isEqualNode(final ConnectionTreeNode nodeInSphere1, final ConnectionTreeNode nodeInSphere2, final Signal signal1, final Signal signal2, final double shiftTol) {
        final IAtom atom1 = nodeInSphere1.getAtom();
        final IAtom atom2 = nodeInSphere2.getAtom();
        final IBond bondToParent1;
        final IBond bondToParent2;
        final ConnectionTreeNode parentNode1;
        final ConnectionTreeNode parentNode2;
        if ((atom1 != null) && (atom2 != null) && !nodeInSphere1.isRingClosureNode() && !nodeInSphere2.isRingClosureNode()) {
            if (atom1.getSymbol().equals(atom2.getSymbol()) && (Integer.compare(atom1.getImplicitHydrogenCount(), atom2.getImplicitHydrogenCount()) == 0) 
                    && (atom1.isAromatic() == atom2.isAromatic()) && (atom1.isInRing() == atom2.isInRing())
                    ) {
                
                if ((signal1 != null) && (signal2 != null) && ((Math.abs(signal1.getShift(0) - signal2.getShift(0)) > shiftTol) || !signal1.getMultiplicity().equals(signal2.getMultiplicity()))) {                    
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
//            parentNode1 = nodeInSphere1.getParentNodes().get(0);
//            parentNode2 = nodeInSphere2.getParentNodes().get(0);
//            bondToParent1 = nodeInSphere1.getBondsToParents().get(0);
//            bondToParent2 = nodeInSphere2.getBondsToParents().get(0);
//            if (parentNode1.getAtom().getSymbol().equals(parentNode2.getAtom().getSymbol()) 
//                    && (bondToParent1.getOrder() == bondToParent2.getOrder()) 
//                    && (bondToParent1.isAromatic() == bondToParent2.isAromatic()) 
//                    && (bondToParent1.isInRing() == bondToParent2.isInRing())) {
                return true;
//            }
        }
        return false;
    }

    /**
     *
     * @param ssc1
     * @param ssc2
     * @param rootAtomIndexSSC1
     * @param rootAtomIndexSSC2
     * @return
     * @throws org.openscience.cdk.exception.CDKException
     */
    public static int getMaximumMatchingSphereHOSECode(final SSC ssc1, final SSC ssc2, final int rootAtomIndexSSC1, final int rootAtomIndexSSC2) throws CDKException {
        int maxMatchingSphere = -1;
        String HOSECodeSSC1;
        String HOSECodeSSC2;
        for (int s = 0; s <= Integer.min(ssc1.getMaxSphere(), ssc2.getMaxSphere()); s++) {
            HOSECodeSSC1 = HOSECodeBuilder.buildHOSECode(ssc1.getSubstructure(), rootAtomIndexSSC1, s, false);
            HOSECodeSSC2 = HOSECodeBuilder.buildHOSECode(ssc2.getSubstructure(), rootAtomIndexSSC2, s, false);
            System.out.println(" --> in s: " + s + " -> " + HOSECodeSSC1 + " vs. " + HOSECodeSSC2);
            if (!HOSECodeSSC1.equals(HOSECodeSSC2)) {
                break;
            }
            maxMatchingSphere = s;
            System.out.println(" --> in s: " + s + " -> " + HOSECodeSSC1 + " vs. " + HOSECodeSSC2);
        }
        return maxMatchingSphere;
    }

    
}
