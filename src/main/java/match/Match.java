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
import casekit.NMR.model.Signal;

import java.util.HashMap;

import hose.model.ConnectionTree;
import hose.model.ConnectionTreeNode;
import model.SSC;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;

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
                signal1 = ssc1.getSubspectrum().getSignal(ssc1.getAssignments().getIndex(0, nodeInSphere1.getKey()));
                signal2 = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getIndex(0, nodeInSphere2.getKey()));
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
//            System.out.println(" --> in s: " + s + " -> " + HOSECodeSSC1 + " vs. " + HOSECodeSSC2);
            if (!HOSECodeSSC1.equals(HOSECodeSSC2)) {
                break;
            }
            maxMatchingSphere = s;
//            System.out.println(" --> in s: " + s + " -> " + HOSECodeSSC1 + " vs. " + HOSECodeSSC2);
        }
        return maxMatchingSphere;
    }

    
}
