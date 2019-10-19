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

import java.util.ArrayList;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Match {

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
            // @TODO add code to handle bigger sizes than 2?
//            final ArrayList<Integer> nodeKeysInSphereSSC2 = connectionTreeSSC2.getNodeKeysInSphere(sphere);
            if(invalids.size() == 2){
                final ArrayList<Integer> nodeKeysInSphereSSC2 = connectionTreeSSC2.getNodeKeysInSphere(sphere);
//                System.out.println("\n -> exactly two invalids -> between " + connectionTreeSSC1.getRootNode().getKey() + ", " + connectionTreeSSC2.getRootNode().getKey() + " -> in sphere: " + sphere + " -> " + invalids);
                if(connectionTreeSSC2.swapChildNodes(connectionTreeSSC2.getNode(nodeKeysInSphereSSC2.get(invalids.get(0))).getParentNodes().get(0).getKey(),
                        nodeKeysInSphereSSC2.get(invalids.get(0)),
                        nodeKeysInSphereSSC2.get(invalids.get(1)))){
                    invalids = Match.getInvalidNodePairsInSphere(ssc1, ssc2, connectionTreeSSC1, connectionTreeSSC2, sphere, shiftTol);
                }
            } else {
                if(invalids.size() > 2){
                    System.out.println("\n -> more than two invalids -> between " + connectionTreeSSC1.getRootNode().getKey() + ", " + connectionTreeSSC2.getRootNode().getKey() + " -> in sphere: " + sphere + " -> " + invalids);
//                    for (int i = 0; i < invalids.size(); i++) {
//                        for (int j = i + 1; j < invalids.size(); j++) {
//                            System.out.println("--> would position i: " + i + " (" + connectionTreeSSC1.getNodesInSphere(sphere).get(i).getKey() + ")" + ", j: " + j + " (" + connectionTreeSSC2.getNodesInSphere(sphere).get(j).getKey() + ")" + " be a valid node pair? -> " +
//                                    Match.isValidNodePair(ssc1, ssc2, connectionTreeSSC1.getNodesInSphere(sphere).get(i), connectionTreeSSC2.getNodesInSphere(sphere).get(j), shiftTol));
//                        }
//                    }
                }
            }
        }

        return (invalids != null) && invalids.isEmpty();
    }

    public static ArrayList<Integer> getInvalidNodePairsInSphere(final SSC ssc1, final SSC ssc2, final ConnectionTree connectionTreeSSC1, final ConnectionTree connectionTreeSSC2, final int sphere, final double shiftTol){
        if(connectionTreeSSC1.getNodesCountInSphere(sphere) != connectionTreeSSC2.getNodesCountInSphere(sphere)){
            return null;
        }
        final ArrayList<ConnectionTreeNode> nodesInSphereSSC1 = connectionTreeSSC1.getNodesInSphere(sphere);
        final ArrayList<ConnectionTreeNode> nodesInSphereSSC2 = connectionTreeSSC2.getNodesInSphere(sphere);
        final ArrayList<Integer> invalids = new ArrayList<>();
        Boolean isValidNodePair;
        for (int i = 0; i < connectionTreeSSC1.getNodesCountInSphere(sphere); i++) {
            isValidNodePair = Match.isValidNodePair(ssc1, ssc2, nodesInSphereSSC1.get(i), nodesInSphereSSC2.get(i), shiftTol);
            if((isValidNodePair != null) && !isValidNodePair){
                invalids.add(i);
            }
        }

        return invalids;
    }

    public static Boolean isValidNodePair(final SSC ssc1, final SSC ssc2, final ConnectionTreeNode node1, final ConnectionTreeNode node2, final double shiftTol){
        if(node1.isRingClosureNode() || node2.isRingClosureNode()){
            return null;
        }
        if(((node1.getAtom().getImplicitHydrogenCount() != null) && (node2.getAtom().getImplicitHydrogenCount() != null))
                && (node1.getAtom().getImplicitHydrogenCount().compareTo(node2.getAtom().getImplicitHydrogenCount()) != 0)){
            return false;
        }
//        if(((node1.getAtom().getHybridization() != null) && (node2.getAtom().getHybridization() != null))
//                && (node1.getAtom().getHybridization().compareTo(node2.getAtom().getHybridization()) != 0)){
//            return false;
//        }
        final Signal signalSSC1 = ssc1.getSubspectrum().getSignal(ssc1.getAssignments().getIndex(0, node1.getKey()));
        final Signal signalSSC2 = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getIndex(0, node2.getKey()));
        if((signalSSC1 == null) || (signalSSC2 == null)){
            return null;
        }
        if(!signalSSC1.getMultiplicity().equals(signalSSC2.getMultiplicity()) // this is, actually, for the last matching sphere
                || (Math.abs(signalSSC1.getShift(0) - signalSSC2.getShift(0)) > shiftTol) // @TODO is this shift comparison needed?
        ){
            return false;
        }

        return true;
    }
}
