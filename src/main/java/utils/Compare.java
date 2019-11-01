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
package utils;

import casekit.NMR.Utils;
import casekit.NMR.model.Signal;
import casekit.NMR.model.Spectrum;
import hose.model.ConnectionTree;
import hose.model.ConnectionTreeNode;
import model.SSC;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Compare {

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

//        ConnectionTree connectionTreeSSC1, connectionTreeSSC2;
//        String HOSECodeSSC1, HOSECodeSSC2;
//        // iterates over each sphere until max. built sphere and builds new conn. trees because of comparable order of the elements
//        for (int s = 0; s <= Integer.min(ssc1.getConnectionTree(atomIndexSSC1).getMaxSphere(), ssc2.getConnectionTree(atomIndexSSC2).getMaxSphere()); s++) {
//            connectionTreeSSC1 = HOSECodeBuilder.buildConnectionTree(ssc1.getSubstructure(), atomIndexSSC1, s);
//            connectionTreeSSC2 = HOSECodeBuilder.buildConnectionTree(ssc2.getSubstructure(), atomIndexSSC2, s);
//            HOSECodeSSC1 = HOSECodeBuilder.buildHOSECode(connectionTreeSSC1, false);
//            HOSECodeSSC2 = HOSECodeBuilder.buildHOSECode(connectionTreeSSC2, false);
//            if (!HOSECodeSSC1.equals(HOSECodeSSC2)) {
//                break;
//            }
////            // new added part:
////            // for each conn. tree node pair in conn. tree (which should at least same multiplicity until max. matching sphere)
////            if(!Compare.rearrangeInvalidNodePairsInSphere(ssc1, ssc2, connectionTreeSSC1, connectionTreeSSC2, s, shiftTol)){
////                break;
////            }
//
//            maxMatchingSphere = s;
////            System.out.println("   --> maxMatchingSphere: " + maxMatchingSphere + " -> " + HOSECodeSSC1 + " vs. " + HOSECodeSSC2);
//        }

        final ArrayList<String> HOSECodeSpheresSSC1 = hose.Utils.splitHOSECodeIntoSpheres(ssc1.getHOSECode(atomIndexSSC1));
        final ArrayList<String> HOSECodeSpheresSSC2 = hose.Utils.splitHOSECodeIntoSpheres(ssc2.getHOSECode(atomIndexSSC2));
//        System.out.println("\n for " + atomIndexSSC1 + ", " + atomIndexSSC2 + ": ");
        for (int s = 0; s < Integer.min(HOSECodeSpheresSSC1.size(), HOSECodeSpheresSSC2.size()); s++) {
//            System.out.println(" --> in s: " + s + "\n-> " + HOSECodeSpheresSSC1.get(s));
//            System.out.println("-> " + HOSECodeSpheresSSC2.get(s));
            if(!HOSECodeSpheresSSC1.get(s).equals(HOSECodeSpheresSSC2.get(s))){
                break;
            }
            maxMatchingSphere = s;
        }

        return maxMatchingSphere;
    }

    /**
     * @param ssc1
     * @param ssc2
     * @param connectionTreeSSC1
     * @param connectionTreeSSC2
     * @param sphere
     * @param shiftTol
     * @return
     *
     * @deprecated
     */
    public static boolean rearrangeInvalidNodePairsInSphere(final SSC ssc1, final SSC ssc2, final ConnectionTree connectionTreeSSC1, final ConnectionTree connectionTreeSSC2, final int sphere, final double shiftTol){
        ArrayList<Integer> invalids = Compare.getInvalidNodePairsInSphere(ssc1, ssc2, connectionTreeSSC1, connectionTreeSSC2, sphere, shiftTol);
        if(invalids == null){
            return false;
        }
        if (!invalids.isEmpty()){
            // @TODO add code to handle bigger sizes than 2?
//            final ArrayList<Integer> nodeKeysInSphereSSC2 = connectionTreeSSC2.getNodeKeysInSphere(sphere);
            if(invalids.size() == 2){
                final ArrayList<Integer> nodeKeysInSphereSSC2 = connectionTreeSSC2.getNodeKeysInSphere(sphere);
//                System.out.println("\n -> exactly two invalids -> between " + connectionTreeSSC1.getRootNode().getKey() + ", " + connectionTreeSSC2.getRootNode().getKey() + " -> in sphere: " + sphere + " -> " + invalids);
                if(connectionTreeSSC2.swapChildNodes(connectionTreeSSC2.getNode(nodeKeysInSphereSSC2.get(invalids.get(0))).getParent().getKey(),
                        nodeKeysInSphereSSC2.get(invalids.get(0)),
                        nodeKeysInSphereSSC2.get(invalids.get(1)))){
                    invalids = Compare.getInvalidNodePairsInSphere(ssc1, ssc2, connectionTreeSSC1, connectionTreeSSC2, sphere, shiftTol);
                }
            } else {
                if(invalids.size() > 2){
//                    System.out.println("\n -> more than two invalids -> between " + connectionTreeSSC1.getRootNode().getKey() + ", " + connectionTreeSSC2.getRootNode().getKey() + " -> in sphere: " + sphere + " -> " + invalids);
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
        if(connectionTreeSSC1.getNodesCountInSphere(sphere, true) != connectionTreeSSC2.getNodesCountInSphere(sphere, true)){
            return null;
        }
        final ArrayList<ConnectionTreeNode> nodesInSphereSSC1 = connectionTreeSSC1.getNodesInSphere(sphere, false);
        final ArrayList<ConnectionTreeNode> nodesInSphereSSC2 = connectionTreeSSC2.getNodesInSphere(sphere, false);
        final ArrayList<Integer> invalids = new ArrayList<>();
        Boolean isValidNodePair;
        for (int i = 0; i < nodesInSphereSSC1.size(); i++) {
            isValidNodePair = Compare.isValidNodePair(ssc1, ssc2, nodesInSphereSSC1.get(i), nodesInSphereSSC2.get(i), shiftTol);
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

    public static boolean compareSSC(final SSC ssc1, final SSC ssc2){
//        if(!ssc1.toHOSECode().equals(ssc2.toHOSECode())){
//            return false;
//        }
//        // pre-filter
//        if(!Arrays.equals(ssc1.getAttachedHydrogensCountsInOuterSphere(), ssc2.getAttachedHydrogensCountsInOuterSphere())
//                || !Arrays.equals(ssc1.getMultiplicitiesInOuterSphere(), ssc2.getMultiplicitiesInOuterSphere())) {
//            return false;
//        }

        // pre-filter
        if(!Compare.getExtendedHOSECode(ssc1).equals(Compare.getExtendedHOSECode(ssc2))){
            return false;
        }

        return Compare.compareSubspectra(ssc1, ssc2);
    }

    public static boolean compareSubspectra(final SSC ssc1, final SSC ssc2){
//        final Spectrum subspectrumSSC1 = Compare.buildHOSECodeOrderedSubspectrum(ssc1);
//        final Spectrum subspectrumSSC2 = Compare.buildHOSECodeOrderedSubspectrum(ssc2);

        if(ssc1.getSubspectrum().getSignalCount() != ssc2.getSubspectrum().getSignalCount()){
            return false;
        }
        for (int i = 0; i < ssc1.getSubspectrum().getSignalCount(); i++) {
            // if both signals exist and the multiplicities are not equal
            if(!ssc1.getSubspectrum().getSignal(i).getMultiplicity().equals(ssc2.getSubspectrum().getSignal(i).getMultiplicity())){
                return false;
            }
        }

        return true;
    }

    public static boolean compareWithMolecularFormula(final IAtomContainer structure, final IMolecularFormula molecularFormula){
        if(molecularFormula != null){
            final HashMap<String, Integer> atomTypeCounts = new HashMap<>();
            IAtom atom;
            for (int i = 0; i < structure.getAtomCount(); i++) {
                atom = structure.getAtom(i);
                if(!atomTypeCounts.containsKey(atom.getSymbol())){
                    atomTypeCounts.put(atom.getSymbol(), 1);
                } else {
                    atomTypeCounts.put(atom.getSymbol(), atomTypeCounts.get(atom.getSymbol()) + 1);
                }
            }
            for (final String atomType : atomTypeCounts.keySet()){
//                System.out.println(" ---> intermediate atom type count for \"" + atomType + "\" : " + atomTypeCounts.get(atomType) + " <= " + MolecularFormulaManipulator.getElementCount(molecularFormula, atomType)
//                        + " ? -> " + (atomTypeCounts.get(atomType) <= MolecularFormulaManipulator.getElementCount(molecularFormula, atomType)));
                if(atomTypeCounts.get(atomType) > MolecularFormulaManipulator.getElementCount(molecularFormula, atomType)){
                    return false;
                }
            }

//            final double maxMolecularWeight = MolecularFormulaManipulator.getMass(molecularFormula);
//            final double molecularWeight = Double.parseDouble(new WeightDescriptor().calculate(structure).getValue().toString());
//            final double allowedDeviation = 0.1;
//            System.out.println(" ---> intermediate weight: " + molecularWeight + " <= " + maxMolecularWeight + " ? -> " + (molecularWeight <= (maxMolecularWeight + allowedDeviation)));
//            return molecularWeight <= (maxMolecularWeight + allowedDeviation);
        }

        return true;
    }

    public static Spectrum buildHOSECodeOrderedSubspectrum(final SSC ssc){
        final Spectrum orderedSubspectrum = new Spectrum(ssc.getSubspectrum().getNuclei());
        final ArrayList<Integer> connectionTreeKeys = ssc.getConnectionTree(ssc.getRootAtomIndex()).getKeys();
        Signal signal;
        for (int k = 0; k < connectionTreeKeys.size(); k++) {
            signal = ssc.getSubspectrum().getSignal(ssc.getAssignments().getIndex(0, connectionTreeKeys.get(k)));
            if(signal != null){
                orderedSubspectrum.addSignal(signal);
            }
        }

        return orderedSubspectrum;
    }

    public static String getExtendedHOSECode(final SSC ssc){
        return ssc.toHOSECode() + "_" + Arrays.toString(ssc.getAttachedHydrogensCountsInOuterSphere()) + "_" + ssc.getSubspectrum().getMultiplicities();//Arrays.toString(ssc.getMultiplicitiesInOuterSphere());
    }

}
