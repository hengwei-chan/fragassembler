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
package match;

import casekit.NMR.Utils;
import casekit.NMR.model.Signal;
import hose.HOSECodeBuilder;
import hose.model.ConnectionTree;
import hose.model.ConnectionTreeNode;
import model.SSC;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;

import java.util.ArrayList;
import java.util.HashMap;

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
     * @param atomIndexSSC1
     * @param atomIndexSSC2
     * @return
     */
    public static int getMaximumMatchingSphereHOSECode(final SSC ssc1, final SSC ssc2, final int atomIndexSSC1, final int atomIndexSSC2, final double shiftTol) throws CDKException {
        int maxMatchingSphere = -1;

        if(!Utils.checkIndexInAtomContainer(ssc1.getSubstructure(), atomIndexSSC1)
                || !Utils.checkIndexInAtomContainer(ssc2.getSubstructure(), atomIndexSSC2)){
            return maxMatchingSphere;
        }

        ConnectionTree connectionTreeSSC1, connectionTreeSSC2;
        String HOSECodeSSC1, HOSECodeSSC2;
        // iterates over each sphere until max. matching sphere and builds new conn. trees because of comparable order of the elements
        for (int s = 0; s <= Integer.min(ssc1.getMaxSphere(), ssc2.getMaxSphere()); s++) {
            connectionTreeSSC1 = HOSECodeBuilder.buildConnectionTree(ssc1.getSubstructure(), atomIndexSSC1, s);
            connectionTreeSSC2 = HOSECodeBuilder.buildConnectionTree(ssc2.getSubstructure(), atomIndexSSC2, s);
            HOSECodeSSC1 = HOSECodeBuilder.buildHOSECode(connectionTreeSSC1, false);
            HOSECodeSSC2 = HOSECodeBuilder.buildHOSECode(connectionTreeSSC2, false);
            if (!HOSECodeSSC1.equals(HOSECodeSSC2)) {
                break;
            }
            // new added part:
            // for each conn. tree node pair in conn. tree (which should at least same multiplicity until max. matching sphere)
            if(!Match.rearrangeInvalidNodePairsInSphere(ssc1, ssc2, connectionTreeSSC1, connectionTreeSSC2, s, shiftTol)){
                break;
            }

            maxMatchingSphere = s;
//            System.out.println("   --> maxMatchingSphere: " + maxMatchingSphere + " -> " + HOSECodeSSC1 + " vs. " + HOSECodeSSC2);
        }

        return maxMatchingSphere;
    }

    public static boolean rearrangeInvalidNodePairsInSphere(final SSC ssc1, final SSC ssc2, final ConnectionTree connectionTreeSSC1, final ConnectionTree connectionTreeSSC2, final int sphere, final double shiftTol){
        ArrayList<Integer> invalids = Match.getInvalidNodePairsInSphere(ssc1, ssc2, connectionTreeSSC1, connectionTreeSSC2, sphere, shiftTol);
        if(invalids == null){
            return false;
        }
        if (!invalids.isEmpty()){
            // @TODO add code to handle  bigger sizes than 2?
            if(invalids.size() == 2){
                final ArrayList<Integer> nodeKeysInSphereSSC2 = connectionTreeSSC2.getNodeKeysInSphere(sphere);
                if(connectionTreeSSC2.swapChildNodes(connectionTreeSSC2.getNode(nodeKeysInSphereSSC2.get(invalids.get(0))).getParentNodes().get(0).getKey(),
                        nodeKeysInSphereSSC2.get(invalids.get(0)),
                        nodeKeysInSphereSSC2.get(invalids.get(1)))){
                    invalids = Match.getInvalidNodePairsInSphere(ssc1, ssc2, connectionTreeSSC1, connectionTreeSSC2, sphere, shiftTol);
                }
            }

            else {
                if(invalids.size() > 2){
                    System.out.println("\n\n\n --> HUHUHUHUHUHU -> " + invalids + "\n\n\n");
                }
            }

        }

        return invalids.isEmpty();
    }

    public static ArrayList<Integer> getInvalidNodePairsInSphere(final SSC ssc1, final SSC ssc2, final ConnectionTree connectionTreeSSC1, final ConnectionTree connectionTreeSSC2, final int sphere, final double shiftTol){
        if(connectionTreeSSC1.getNodesCountInSphere(sphere) != connectionTreeSSC2.getNodesCountInSphere(sphere)){
            return null;
        }
        final ArrayList<ConnectionTreeNode> nodesInSphereSSC1 = connectionTreeSSC1.getNodesInSphere(sphere);
        final ArrayList<ConnectionTreeNode> nodesInSphereSSC2 = connectionTreeSSC2.getNodesInSphere(sphere);
        final ArrayList<Integer> invalids = new ArrayList<>();
        ConnectionTreeNode nodeInSphereSSC1, nodeInSphereSSC2;
        Signal signalSSC1, signalSSC2;
        for (int i = 0; i < connectionTreeSSC1.getNodesCountInSphere(sphere); i++) {
            nodeInSphereSSC1 = nodesInSphereSSC1.get(i);
            nodeInSphereSSC2 = nodesInSphereSSC2.get(i);
            if(nodeInSphereSSC1.isRingClosureNode()){
                continue;
            }

            if(((nodeInSphereSSC1.getAtom().getHybridization() != null) && (nodeInSphereSSC2.getAtom().getHybridization() != null))
                    && (nodeInSphereSSC1.getAtom().getHybridization().compareTo(nodeInSphereSSC2.getAtom().getHybridization()) != 0)){
                invalids.add(i);
                continue;
            }
            if(((nodeInSphereSSC1.getAtom().getCharge() != null) && (nodeInSphereSSC2.getAtom().getCharge() != null))
                    && (nodeInSphereSSC1.getAtom().getCharge().compareTo(nodeInSphereSSC2.getAtom().getCharge()) != 0)){
                invalids.add(i);
                continue;
            }

            signalSSC1 = ssc1.getSubspectrum().getSignal(ssc1.getAssignments().getIndex(0, nodeInSphereSSC1.getKey()));
            signalSSC2 = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getIndex(0, nodeInSphereSSC2.getKey()));
            if((signalSSC1 == null) || (signalSSC2 == null)){
                continue;
            }
            if(!signalSSC1.getMultiplicity().equals(signalSSC2.getMultiplicity()) // this is, actually, for the last matching sphere
                            || (Math.abs(signalSSC1.getShift(0) - signalSSC2.getShift(0)) > shiftTol) // @TODO is this shift comparison needed?
            ){
                invalids.add(i);
            }

        }

        return invalids;
    }

}
