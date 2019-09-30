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
package assembly;

import casekit.NMR.Utils;
import casekit.NMR.match.Matcher;
import casekit.NMR.model.Assignment;
import casekit.NMR.model.Signal;
import casekit.NMR.model.Spectrum;
import hose.HOSECodeBuilder;
import hose.model.ConnectionTree;
import hose.model.ConnectionTreeNode;
import match.Match;
import model.SSC;
import model.SSCLibrary;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.isomorphism.mcss.RMap;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import start.Start;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.*;
import java.util.Map.Entry;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Assembly {  


    private static HashMap<Integer, ArrayList<Integer[]>> getOverlapsHOSECodeNew(final SSC ssc1, final SSC ssc2, final int minMatchingSphereCount, final double shiftTol){
        final HashMap<Integer, ArrayList<Integer[]>> overlapsInSpheres = new HashMap<>();
        int maxMatchingSphere;
        for (int i = 0; i < ssc1.getAtomCount(); i++) {
            for (int j = 0; j < ssc2.getAtomCount(); j++) {
                try {
                    maxMatchingSphere = Match.getMaximumMatchingSphereHOSECode(ssc1, ssc2, i, j, shiftTol);
                } catch (CDKException e) {
                   maxMatchingSphere = -1;
                }
                if(maxMatchingSphere < minMatchingSphereCount){
                    continue;
                }

                // create new key for the found max. matching sphere if it's not existing
                if (!overlapsInSpheres.containsKey(maxMatchingSphere)) {
                    overlapsInSpheres.put(maxMatchingSphere, new ArrayList<>());
                }
                // insert matching atom pair into mappings hashmap with found max. matching sphere
                overlapsInSpheres.get(maxMatchingSphere).add(new Integer[]{i, j});
            }
        }

        return overlapsInSpheres;
    }

    private static HashMap<Integer, ArrayList<Double[]>> getOverlapsHOSECodeCore(final SSC ssc1, final SSC ssc2, final HashSet<Integer> atomIndicesSSC1, final HashSet<Integer> atomIndicesSSC2, final int minSphereMatchCount, final double shiftTol) throws CDKException {
        final HashMap<Integer, ArrayList<Double[]>> overlapsInSpheres = new HashMap<>();

        ArrayList<ConnectionTree> connectionTrees;
        int maxMatchingSphere, overlappingAtomsCount;
        ConnectionTree connectionTreeSSC1;
        ConnectionTree connectionTreeSSC2;
        Signal signalSSC1, signalSSC2;
        for (final int i : atomIndicesSSC1) {
            for (final int j : atomIndicesSSC2) {
                connectionTrees = Assembly.mapConnectionTrees(ssc1, ssc2, i, j, minSphereMatchCount, shiftTol);
                if (connectionTrees.isEmpty()) {
//                    System.out.println("-> no valid connection trees found!!!");
                    continue;
                }
                connectionTreeSSC1 = connectionTrees.get(0);
                connectionTreeSSC2 = connectionTrees.get(1);
                maxMatchingSphere = connectionTreeSSC1.getMaxSphere();
                System.out.println("-> maxMatchingSphere " + maxMatchingSphere + " is valid!!!");
                System.out.println("conn. tree 1: -> " + connectionTreeSSC1.getMaxSphere() + " -> " + connectionTreeSSC1);
                System.out.println("conn. tree 2: -> " + connectionTreeSSC2.getMaxSphere() + " -> " + connectionTreeSSC2);
                // count number and calculate deviations of overlapping atoms until max. matching sphere; for ranking later
                overlappingAtomsCount = 0;
                final ArrayList<Double> deviations = new ArrayList<>();
                ArrayList<ConnectionTreeNode> nodesInSphereSSC1, nodesInSphereSSC2;
                for (int s = 0; s <= maxMatchingSphere; s++) {
                    nodesInSphereSSC1 = connectionTreeSSC1.getNodesInSphere(s);
                    nodesInSphereSSC2 = connectionTreeSSC2.getNodesInSphere(s);
                    for (int k = 0; k < nodesInSphereSSC1.size(); k++) {
                        if (!nodesInSphereSSC1.get(k).isRingClosureNode()) {
                            overlappingAtomsCount++;
                        }
                        signalSSC1 = ssc1.getSubspectrum().getSignal(ssc1.getAssignments().getIndex(0, nodesInSphereSSC1.get(k).getKey()));
                        signalSSC2 = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getIndex(0, nodesInSphereSSC2.get(k).getKey()));
                        if((signalSSC1 != null) && (signalSSC2 != null)){
                            deviations.add(Math.abs(signalSSC1.getShift(0) - signalSSC2.getShift(0)));
                        }
                    }
                }
                // create new key for the found max. matching sphere if it's not existing
                if (!overlapsInSpheres.containsKey(maxMatchingSphere)) {
                    overlapsInSpheres.put(maxMatchingSphere, new ArrayList<>());
                }
                // insert matching atom pair into mappings hashmap with found max. matching sphere
                overlapsInSpheres.get(maxMatchingSphere).add(new Double[]{(double) i, (double) j, (double) overlappingAtomsCount, Utils.getMean(deviations)});
            }
        }

        return overlapsInSpheres;
    }

//    private static HashMap<Integer, ArrayList<Double[]>> getOverlapsHOSECodeCore(final SSC ssc1, final SSC ssc2, final HashMap<Integer, Integer> atomMappings, final int minSphereMatchCount, final double shiftTol) throws CDKException {
//        final HashMap<Integer, ArrayList<Double[]>> overlapsInSpheres = new HashMap<>();
//
//        ArrayList<ConnectionTree> connectionTrees;
//        int maxMatchingSphere, overlappingAtomsCount;
//        ConnectionTree connectionTreeSSC1;
//        ConnectionTree connectionTreeSSC2;
//        Signal signalSSC1, signalSSC2;
//            for (final Entry<Integer, Integer> entry : atomMappings.entrySet()) {
//                connectionTrees = Assembly.mapConnectionTrees(ssc1, ssc2, entry.getKey(), entry.getValue(),minSphereMatchCount, shiftTol);
//                if (connectionTrees.isEmpty()) {
//                    continue;
//                }
//                connectionTreeSSC1 = connectionTrees.get(0);
//                connectionTreeSSC2 = connectionTrees.get(1);
//                maxMatchingSphere = connectionTreeSSC1.getMaxSphere();
//                System.out.println("-> maxMatchingSphere " + maxMatchingSphere + " is valid!!!");
//                // count number and calculate deviations of overlapping atoms until max. matching sphere; for ranking later
//                overlappingAtomsCount = 0;
//                final ArrayList<Double> deviations = new ArrayList<>();
//                ArrayList<ConnectionTreeNode> nodesInSphereSSC1, nodesInSphereSSC2;
//                for (int s = 0; s <= maxMatchingSphere; s++) {
//                    nodesInSphereSSC1 = connectionTreeSSC1.getNodesInSphere(s);
//                    nodesInSphereSSC2 = connectionTreeSSC2.getNodesInSphere(s);
//                    for (int k = 0; k < nodesInSphereSSC1.size(); k++) {
//                        if (!nodesInSphereSSC1.get(k).isRingClosureNode()) {
//                            overlappingAtomsCount++;
//                        }
//                        signalSSC1 = ssc1.getSubspectrum().getSignal(ssc1.getAssignments().getIndex(0, nodesInSphereSSC1.get(k).getKey()));
//                        signalSSC2 = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getIndex(0, nodesInSphereSSC2.get(k).getKey()));
//                        if((signalSSC1 != null) && (signalSSC2 != null)){
//                            deviations.add(Math.abs(signalSSC1.getShift(0) - signalSSC2.getShift(0)));
//                        }
//                    }
//                }
//                // create new key for the found max. matching sphere if it's not existing
//                if (!overlapsInSpheres.containsKey(maxMatchingSphere)) {
//                    overlapsInSpheres.put(maxMatchingSphere, new ArrayList<>());
//                }
//                // insert matching atom pair into mappings hashmap with found max. matching sphere
//                overlapsInSpheres.get(maxMatchingSphere).add(new Double[]{(double) entry.getKey(), (double) entry.getValue(), (double) overlappingAtomsCount, Utils.getMean(deviations)});
//            }
//
//        return overlapsInSpheres;
//    }

    private static ArrayList<ConnectionTree> mapConnectionTrees(final SSC ssc1, final SSC ssc2, final int atomIndexSSC1, final int atomIndexSSC2, final int minSphereMatchCount, final double shiftTol) throws CDKException {
        final ArrayList<ConnectionTree> connectionTrees = new ArrayList<>();

        Signal signalSSC1, signalSSC2;
        ArrayList<Integer> nodeKeysInSphereToCheckSSC1, nodeKeysInSphereToCheckSSC2;
        HashMap<Integer, Integer> mappedEqualNodesInSphere;
        boolean stop;

        IAtom rootAtomSSC1 = ssc1.getSubstructure().getAtom(atomIndexSSC1);
        Signal signalRootAtomSSC1 = ssc1.getSubspectrum().getSignal(ssc1.getAssignments().getIndex(0, atomIndexSSC1));
        IAtom rootAtomSSC2 = ssc2.getSubstructure().getAtom(atomIndexSSC2);
        // check for same atom types
        if(!rootAtomSSC1.getSymbol().equals(rootAtomSSC2.getSymbol())){
            return connectionTrees;
        }
        // get signals of both atoms to compare (if available) and check for same signal properties
        // to set a starting point with higher probability to be correct
        Signal signalRootAtomSSC2 = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getIndex(0, atomIndexSSC2));
        if ((signalRootAtomSSC1 != null) && (signalRootAtomSSC2 != null)
                &&  (!signalRootAtomSSC1.getMultiplicity().equals(signalRootAtomSSC2.getMultiplicity())
                || Math.abs(signalRootAtomSSC1.getShift(0) - signalRootAtomSSC2.getShift(0)) > shiftTol)
        ){
            return connectionTrees;
        }
        System.out.println("\n-> i: " + atomIndexSSC1 + ", j: " + atomIndexSSC2);
        // check for pure structural identity via HOSE code for each further sphere
        int maxMatchingSphere = Match.getMaximumMatchingSphereHOSECode(ssc1, ssc2, atomIndexSSC1, atomIndexSSC2, shiftTol);
        System.out.println("-> maxMatchingSphere: " + maxMatchingSphere);
        // skip non-matching atom pairs
        if (maxMatchingSphere == -1) {
            return connectionTrees;
        }
        // skip atom pairs which do not match at least in given minimum number of spheres
        if(maxMatchingSphere < minSphereMatchCount) {
//                        // but atoms at open sites are allowed because those sites could enable
//                        // an extension if no other atom mappings exist
            if(!ssc1.isUnsaturated(atomIndexSSC1) && !ssc2.isUnsaturated(atomIndexSSC2)){
                return connectionTrees;
            }
        }
        final ConnectionTree connectionTreeSSC1 = HOSECodeBuilder.buildConnectionTree(ssc1.getSubstructure(), atomIndexSSC1, maxMatchingSphere);
        final ConnectionTree connectionTreeSSC2 = HOSECodeBuilder.buildConnectionTree(ssc2.getSubstructure(), atomIndexSSC2, maxMatchingSphere);
        // if a atom pair matched succesfully via HOSE codes (structurally) then check now for spectral identity/similarity
        boolean notValid = false;
        // check all spheres until max. matching sphere for spectral similarity/identity
        for (int s = 0; s <= maxMatchingSphere; s++) {
            nodeKeysInSphereToCheckSSC1 = connectionTreeSSC1.getNodeKeysInSphere(s);
            nodeKeysInSphereToCheckSSC2 = connectionTreeSSC2.getNodeKeysInSphere(s);
//                    System.out.println("in s: " + s + " -> node keys SSC1: " + nodeKeysInSphereToCheckSSC1);
//                    System.out.println("in s: " + s + " -> node keys SSC2: " + nodeKeysInSphereToCheckSSC2);
            stop = false;
            // check for same/similar signal properties in same order as the belonging HOSE codes until max. matching sphere
            for (int k = 0; k < nodeKeysInSphereToCheckSSC1.size(); k++) {
                signalSSC1 = ssc1.getSubspectrum().getSignal(ssc1.getAssignments().getIndex(0, nodeKeysInSphereToCheckSSC1.get(k)));
                signalSSC2 = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getIndex(0, nodeKeysInSphereToCheckSSC2.get(k)));
                if ((signalSSC1 != null) && (signalSSC2 != null)
                        && (!signalSSC1.getMultiplicity().equals(signalSSC2.getMultiplicity()) || Math.abs(signalSSC1.getShift(0) - signalSSC2.getShift(0)) > shiftTol)) {
                    stop = true;
                    break;
                }
            }
            // it could happen that a structural identity until a certain sphere is there but the actual order of atoms
            // within a HOSE code is different than in the second HOSE code;
            // then one could search for a possible atom mapping via spectral properties
            // within a sphere in different order than HOSE code order
            if (stop) {
                mappedEqualNodesInSphere = Match.mapEqualNodesInSphere(ssc1, ssc2, atomIndexSSC1, atomIndexSSC2, s, shiftTol);
//                        System.out.println("---> in s: " + s + " (spectral matching) --> maps in sphere: " + mappedEqualNodesInSphere);
                if (mappedEqualNodesInSphere.containsValue(-1)) {
                    notValid = true;
                    break;
                }
//                        notValid = true;
//                        break;
            }
        }
        // skip non-valid matching atom pairs
        if (notValid) {
//                    System.out.println("-> maxMatchingSphere is NOT valid!!!");
            return connectionTrees;
        }
//        System.out.println("-> maxMatchingSphere " + maxMatchingSphere + " is valid!!!");


        connectionTrees.add(connectionTreeSSC1);
        connectionTrees.add(connectionTreeSSC2);

        return connectionTrees;
    }

    /**
     * Returns pairwise structural identities between two atoms (only specific ones) in both
     * substructures, including the maximum matching sphere count.
     *
     * @param ssc1 SSC to map against {@code ssc2};
     * @param ssc2 SSC to map against {@code ssc1};
     * @param atomIndicesSSC1 specific indices of atoms in {@code ssc1} to use
     * @param atomIndicesSSC2 specific indices of atoms in {@code ssc2} to use
     * @param minSphereMatchCount minimum matching sphere count during
     * HOSE code (sphere) comparison
     * @param shiftTol shift tolerance value [ppm] used during atom and shift
     * comparisons in HOSE code spheres
     * @return
     *
     * @throws org.openscience.cdk.exception.CDKException
     *
     */
    public static HashMap<Integer, ArrayList<Double[]>> getOverlapsHOSECode(final SSC ssc1, final SSC ssc2, final HashSet<Integer> atomIndicesSSC1, final HashSet<Integer> atomIndicesSSC2, final int minSphereMatchCount, final double shiftTol) throws CDKException {
        return Assembly.getOverlapsHOSECodeCore(ssc1, ssc2, atomIndicesSSC1, atomIndicesSSC2, minSphereMatchCount, shiftTol);
    }

    /**
     * Returns pairwise structural identities between two atoms in both
     * substructures, including the maximum matching sphere count.
     *
     * @param ssc1 SSC to map against {@code ssc2};
     * @param ssc2 SSC to map against {@code ssc1};
     * @param minSphereMatchCount minimum matching sphere count during
     * HOSE code (sphere) comparison
     * @param shiftTol shift tolerance value [ppm] used during atom and shift
     * comparisons in HOSE code spheres
     * @return
     *
     * @throws org.openscience.cdk.exception.CDKException
     *
     */
    public static HashMap<Integer, ArrayList<Double[]>> getOverlapsHOSECode(final SSC ssc1, final SSC ssc2, final int minSphereMatchCount, final double shiftTol) throws CDKException {
        final HashSet<Integer> atomIndicesSSC1 = new HashSet<>();
        for (final IAtom atom : ssc1.getSubstructure().atoms()){
            atomIndicesSSC1.add(atom.getIndex());
        }
        final HashSet<Integer> atomIndicesSSC2 = new HashSet<>();
        for (final IAtom atom : ssc2.getSubstructure().atoms()){
            atomIndicesSSC2.add(atom.getIndex());
        }

        return Assembly.getOverlapsHOSECodeCore(ssc1, ssc2, atomIndicesSSC1, atomIndicesSSC2, minSphereMatchCount, shiftTol);
    }
    
    /**
     * Maps sets of given atom pairs as starting point of structural overlaps 
     * in different spheres into one single set of atom pairs. 
     * This method starts with the highest matching sphere to increase the 
     * probability of correctness.
     *
     * @param ssc1 SSC to map against {@code ssc2}; used as keys in returning 
     * object
     * @param ssc2 SSC to map against {@code ssc1}; used as values in returning 
     * @param overlapsInSpheres HashMap containing the maximum 
     * matching spheres as keys and the belonging mapped atom pairs as well as
     * overlapping atoms count and mean of deviations as values
     * @return 
     * @throws CDKException
     * 
     * @see #getOverlapsHOSECode(model.SSC, model.SSC, int, double) 
     */
    public static HashMap<Integer, Integer> getAtomMappingsHOSECode(final SSC ssc1, final SSC ssc2, final HashMap<Integer, ArrayList<Double[]>> overlapsInSpheres, final HashMap<Integer, Integer> preMappedAtomIndices) throws CDKException{
        // insertion of all valid mappings in different spheres into one map, 
        // starting with highest matching sphere (higher chance of correctness)
        if (overlapsInSpheres.isEmpty()) {
            return new HashMap<>();
        }        
        
        final HashMap<Integer, Integer> mappedAtomIndices = new HashMap<>(preMappedAtomIndices);
        ArrayList<Integer> nodeKeysSSC1, nodeKeysSSC2;
        int rootMappedAtomIndexSSC1, rootMappedAtomIndexSSC2, mappedAtomIndexSSC1, mappedAtomIndexSSC2;
        ConnectionTree connectionTreeSSC1, connectionTreeSSC2;
//        boolean containsUnsaturatedAtomsSSC2;
        ArrayList<Double[]> overlapsInSphere;
        // for all found spheres try to insert the mappings
        for (int m = Collections.max(overlapsInSpheres.keySet()); m >= Collections.min(overlapsInSpheres.keySet()); m--) {
            if (!overlapsInSpheres.containsKey(m)) {
                continue;
            }
            // findHits overlaps in current max. sphere by (1) mean of deviations and (2) total number of overlapping atoms
            overlapsInSphere = overlapsInSpheres.get(m);
            overlapsInSphere.sort((entry1, entry2) -> {
                int deviationsComp;// = ((entry1[3] == null) || (entry2[3] == null)) ? 0 : Double.compare(entry1[3], entry2[3]);
                if((entry1[3] == null) || (entry2[3] == null)){
                    deviationsComp = 0;
                } else {
                    deviationsComp = Double.compare(entry1[3], entry2[3]);
                }
                if(deviationsComp != 0){
                    return deviationsComp;
                }
                return -1 * Integer.compare(entry1[2].intValue(), entry2[2].intValue());
            });
            // for each atom pair as valid starting point insert its belonging mapped atom pairs until its max. matching sphere            
            for (final Double[] entry : overlapsInSphere) {
                rootMappedAtomIndexSSC1 = entry[0].intValue();
                rootMappedAtomIndexSSC2 = entry[1].intValue();
                connectionTreeSSC1 = HOSECodeBuilder.buildConnectionTree(ssc1.getSubstructure(), rootMappedAtomIndexSSC1, m);
                connectionTreeSSC2 = HOSECodeBuilder.buildConnectionTree(ssc2.getSubstructure(), rootMappedAtomIndexSSC2, m);
                System.out.println("rootMappedAtomIndexSSC1: " + rootMappedAtomIndexSSC1 + "\t-> " + connectionTreeSSC1.toString());
                System.out.println("rootMappedAtomIndexSSC2: " + rootMappedAtomIndexSSC2 + "\t-> " + connectionTreeSSC2.toString());
                nodeKeysSSC1 = new ArrayList<>(connectionTreeSSC1.getKeys(true));
                nodeKeysSSC2 = new ArrayList<>(connectionTreeSSC2.getKeys(true));
                for (int k = 0; k < nodeKeysSSC1.size(); k++) {
                    mappedAtomIndexSSC1 = nodeKeysSSC1.get(k);
                    mappedAtomIndexSSC2 = nodeKeysSSC2.get(k);
                    System.out.println(" -> trying to insert: " + mappedAtomIndexSSC1 + ", " + mappedAtomIndices.get(mappedAtomIndexSSC1));
                    if (!mappedAtomIndices.containsKey(mappedAtomIndexSSC1)
                            && !mappedAtomIndices.containsValue(mappedAtomIndexSSC2)) {

//                        // check whether at least one mapped atom in SSC2 (SSC to extend) is unsaturated
//                        containsUnsaturatedAtomsSSC2 = false;
//                        for (final int nodeKey : ssc2.getConnectionTree(mappedAtomIndexSSC2).getKeys(true)) {
//                            if (ssc2.isUnsaturated(nodeKey)) {
//                                containsUnsaturatedAtomsSSC2 = true;
//                                break;
//                            }
//                        }
//                        if (!containsUnsaturatedAtomsSSC2) {
//                            System.out.println("\nmap index ssc1: " + mappedAtomIndexSSC1
//                                    + ", map index ssc2: " + mappedAtomIndexSSC2
//                                    + " -> nothing unsaturated -> skip this overlap");
//                            continue;
//                        }

                        mappedAtomIndices.put(mappedAtomIndexSSC1, mappedAtomIndexSSC2);
                        System.out.println(" ---> inserted: " + mappedAtomIndexSSC1 + ", " + mappedAtomIndices.get(mappedAtomIndexSSC1));
                    }
                }
                System.out.println();
            }
        }

        return mappedAtomIndices;
    }
    
    public static boolean isValidBondAddition(final IAtomContainer ac, final int atomIndex, final IBond bondToAdd){        

        System.out.println(atomIndex + " --> " + Utils.getBondOrderSum(ac, atomIndex, true) + " + " + Utils.getBondOrderAsNumeric(bondToAdd) + " = " + (Utils.getBondOrderSum(ac, atomIndex, true) + Utils.getBondOrderAsNumeric(bondToAdd)) + " <= " + ac.getAtom(atomIndex).getValency() + " ? -> " + ((Utils.getBondOrderSum(ac, atomIndex, true) + Utils.getBondOrderAsNumeric(bondToAdd)) <= ac.getAtom(atomIndex).getValency()));
        
        return (Utils.getBondOrderSum(ac, atomIndex, true) + Utils.getBondOrderAsNumeric(bondToAdd)) <= ac.getAtom(atomIndex).getValency();
    }
    
    public static boolean isValidSubspectrum(final Spectrum subspectrum, final Spectrum querySpectrum, final double shiftTol, final double thrsMatchFactor){
        if((subspectrum == null) || (subspectrum.getSignalCount() > querySpectrum.getSignalCount())){
            return false;
        }
        final Assignment matchAssignments = Matcher.matchSpectra(subspectrum, querySpectrum, 0, 0, shiftTol);
        // filter for unset assignments
        if (!matchAssignments.isFullyAssigned(0)){//matchAssignments.getSetAssignmentsCount(0) < matchAssignments.getAssignmentsCount()) {
//            System.out.println("-> set assigments not allowed!!!");
            return false;
        }
        // filter for multiple assignments
        // there might be multiple assignments to same signals, so check for possible symmetry (equivalences)
        for (final int matchedSignalIndexInQuerySpectrum : matchAssignments.getAssignments(0)) {
            if (Collections.frequency(matchAssignments.getAssignments(0), matchedSignalIndexInQuerySpectrum)
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
        if (Utils.roundDouble(Matcher.calculateAverageDeviation(subspectrum, querySpectrum, 0, 0, shiftTol),  Start.DECIMAL_PLACES) > thrsMatchFactor) {
//            System.out.println("-> match factor not allowed!!!");
            return false;
        }
        
        return true;
    }
    
    public static HashMap<String, SSC> assemble(final long nStarts, final int nThreads, final SSCLibrary rankedSSCLibrary, final int minMatchingSphereCount,
            final Spectrum querySpectrum, final double thrsMatchFactor, final double shiftTol, final String pathToOutputsFolder, final long querySpectrumCounter) throws InterruptedException {

//        long counter = 0;
//        for (final long sscindex : rankedSSCLibrary.getSSCIndices()) {
//            if (counter <= 100) {
//                try {
//                    Utils.generatePicture(rankedSSCLibrary.getSSC(sscindex).getSubstructure(), "results/out_" + counter + ".png");
//                } catch (Exception e) {
//                    System.out.println("could not depict for ranked ssc index " + sscindex + ": " + e.getMessage());
//                }
//                counter++;
//            } else {
//                break;
//            }
//        }


        final HashMap<String, SSC> solutions = new HashMap<>();
        // initialize an executor for parallelization
        final ExecutorService executor = Utils.initExecuter(nThreads);
        final ArrayList<Callable<HashMap<String, SSC>>> callables = new ArrayList<>();
        // add all task to do
        for (int i = 0; i < nStarts; i++) {
            final int j = i;
            callables.add(() -> {
                return Assembly.assembleDFS(rankedSSCLibrary, j, minMatchingSphereCount, querySpectrum, thrsMatchFactor, shiftTol, pathToOutputsFolder, querySpectrumCounter);
//                return Assembly.assembleBFS(rankedSSCLibrary, j, minMatchingSphereCount, querySpectrum, thrsMatchFactor, shiftTol);
//                return Assembly.assembleSeq(rankedSSCLibrary, j, minMatchingSphereCount, querySpectrum, thrsMatchFactor, shiftTol);
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
                .forEach(tempHashMap -> {
                    solutions.putAll(tempHashMap);
                });
        // shut down the executor service
        Utils.stopExecuter(executor, 5);
        
        return solutions;
    } 

    public static boolean isFinalSSC(final SSC ssc, final Spectrum querySpectrum, final double shiftTol, final double thrsMatchFactor) {
        System.out.println("\nno more unsaturated atoms left? -> " + !ssc.hasUnsaturatedAtoms());
        if(ssc.hasUnsaturatedAtoms()){
            return false;
        }
        System.out.println("query spectrum size reached? -> " + ssc.getSubspectrum().getSignalCount() + " == " + querySpectrum.getSignalCount() + " -> " + (ssc.getSubspectrum().getSignalCount() == querySpectrum.getSignalCount()));
        if((ssc.getSubspectrum().getSignalCount() != querySpectrum.getSignalCount())){
            return false;
        }
        System.out.println("isValidSpectrum? -> " + Assembly.isValidSubspectrum(ssc.getSubspectrum(), querySpectrum, shiftTol, thrsMatchFactor));
        if(!Assembly.isValidSubspectrum(ssc.getSubspectrum(), querySpectrum, shiftTol, thrsMatchFactor)){
            return false;
        }

        try {
            Kekulization.kekulize(ssc.getSubstructure());
            System.out.println("kekulization? -> true");
        } catch (CDKException e) {
            System.out.println("kekulization? -> false");
            return false;
        }
        
        return true;
    }
    
    private static HashSet<Integer> getMissingAtomIndices(final SSC ssc, final Set<Integer> mappedAtomIndices){
        final HashSet<Integer> missingAtomIndices = new HashSet<>();
        for (int k = 0; k < ssc.getAtomCount(); k++) {
            if (!mappedAtomIndices.contains(k)) {
                missingAtomIndices.add(k);
            }
        }
        
        return missingAtomIndices;
    }
    
    /**
     *
     * @param ssc1 SSC to extend 
     * @param ssc2
     * @param atomMappings indices mappings from {@code ssc1} (key) to
     * {@code ssc2} (value)
     * @return
     * @throws CDKException
     * @throws CloneNotSupportedException
     */
    public static boolean addMissingBonds(final SSC ssc1, final SSC ssc2, final HashMap<Integer, Integer> atomMappings) throws CDKException, CloneNotSupportedException {
        int mappedAtomIndexSSC1, mappedAtomIndexSSC2;
        IAtom mappedAtomSSC1, mappedAtomSSC2, connectedAtomSSC1;
        IBond bondToAdd;
        boolean addedAnyBond = false;
        final HashMap<Integer, Integer> reversedAtomMappings = Assembly.reverseHashMap(atomMappings);
        
        for (final Entry<Integer, Integer> entry : reversedAtomMappings.entrySet()) {
            mappedAtomIndexSSC2 = entry.getKey();
            mappedAtomIndexSSC1 = entry.getValue();            
            if ((ssc1.isUnsaturated(mappedAtomIndexSSC1) == null) || !ssc1.isUnsaturated(mappedAtomIndexSSC1)) {
                continue;
            }            
//            System.out.println("\nmapping: " + mappedAtomIndexSSC2 + ", " + mappedAtomIndexSSC1
//                    + "\n-> " + ssc1.getHOSECode(mappedAtomIndexSSC1) + "\n-> " + ssc2.getHOSECode(mappedAtomIndexSSC2));
            
            mappedAtomSSC2 = ssc2.getSubstructure().getAtom(mappedAtomIndexSSC2);
            mappedAtomSSC1 = ssc1.getSubstructure().getAtom(mappedAtomIndexSSC1);
            for (final IAtom connectedAtomSSC2 : ssc2.getSubstructure().getConnectedAtomsList(mappedAtomSSC2)) {
                if(!reversedAtomMappings.containsKey(connectedAtomSSC2.getIndex()) || reversedAtomMappings.get(connectedAtomSSC2.getIndex()) == mappedAtomIndexSSC1){
                    continue;
                }
                connectedAtomSSC1 = ssc1.getSubstructure().getAtom(reversedAtomMappings.get(connectedAtomSSC2.getIndex()));
                if(ssc1.getSubstructure().getBond(mappedAtomSSC1, connectedAtomSSC1) == null){
                    bondToAdd = ssc2.getSubstructure().getBond(mappedAtomSSC2, connectedAtomSSC2).clone();
                    bondToAdd.setAtom(mappedAtomSSC1, 0);
                    bondToAdd.setAtom(connectedAtomSSC1, 1);
                    if(Assembly.isValidBondAddition(ssc1.getSubstructure(), mappedAtomIndexSSC1, bondToAdd)
                            && Assembly.isValidBondAddition(ssc1.getSubstructure(), connectedAtomSSC1.getIndex(), bondToAdd)){
//                        if(bondToAdd.isAromatic()){
//                            bondToAdd.setOrder(IBond.Order.SINGLE);
//                        }
                        ssc1.getSubstructure().addBond(bondToAdd);
                        ssc1.update();  
                        addedAnyBond = true;           
                        System.out.println("in addMissingBonds: new bond (" + bondToAdd.getOrder().numeric() + ", " + bondToAdd.isAromatic() + ") between " + mappedAtomIndexSSC1 + ", " + connectedAtomSSC1.getIndex() + " added!!!");
                    }                    
                }
            }
        }
        
        return addedAnyBond;
    }
    
//    public static SSC assemblyCore(final SSC ssc1, final SSC ssc2, final Spectrum querySpectrum, final int minMatchingSphereCount, final double shiftTol, final double thrsMatchFactor) throws Exception {
//        int mappedAtomIndexSSC2;
//        HashMap<Integer, Integer> mappedAtomIndices;
//        HashSet<Integer> missingAtomIndicesSSC2;
//
//        // 1. check for partial structural identity (overlaps)
//        // 1.1 map structural overlaps in all spheres separately
//        final HashMap<Integer, ArrayList<Double[]>> overlapsInSpheres = Assembly.getOverlapsHOSECode(ssc2, ssc1, minMatchingSphereCount, shiftTol);
//        // 1.2 insertion of all valid mappings in different spheres into one map
//        mappedAtomIndices = Assembly.getAtomMappingsHOSECode(ssc2, ssc1, overlapsInSpheres);
//        if (mappedAtomIndices.isEmpty()) {
//            return null;
//        }
//        missingAtomIndicesSSC2 = Assembly.getMissingAtomIndices(ssc2, mappedAtomIndices.keySet());
//        System.out.println(" -> mapped atom indices: \t" + mappedAtomIndices);
//        System.out.println(" -> missing atom indices:\t" + missingAtomIndicesSSC2);
//        // 1.3 for subspectra combination: add signals from SSC2's subspectrum associated with mapping (overlapping) atoms to a new subspectrum to combine with subspectrum of SSC1
//        final LinkedHashSet<Integer> signalIndicesToIgnoreFromSubspectrumSSC2 = new LinkedHashSet<>();
//        for (final Entry<Integer, Integer> entry : mappedAtomIndices.entrySet()) {
//            mappedAtomIndexSSC2 = entry.getKey();
//            if (ssc2.getSubstructure().getAtom(mappedAtomIndexSSC2).getSymbol().equals(ssc1.getSubspectrumAtomType())) {
//                signalIndicesToIgnoreFromSubspectrumSSC2.add(ssc2.getAssignments().getSignalIndex(0, mappedAtomIndexSSC2));
//            }
//        }
//        final Spectrum subspectrumToAddSSC2 = new Spectrum(ssc2.getSubspectrum().getNuclei());
//        for (int s = 0; s < ssc2.getSubspectrum().getSignalCount(); s++) {
//            if (!signalIndicesToIgnoreFromSubspectrumSSC2.contains(s)) {
//                subspectrumToAddSSC2.addSignal(ssc2.getSubspectrum().getSignal(s).getClone(), ssc2.getSubspectrum().getEquivalence(s));
//            }
//        }
//        // 2. combine subspectra of SSC to extend and matched SSC and validate it
//        final Spectrum combinedSpectrum = Matcher.combineSpectra(ssc1.getSubspectrum(), subspectrumToAddSSC2, 0, 0, Start.EQUIV_SIGNAL_THRS);
////            System.out.println("\nspectrum1:\t" + ssc1.getSubspectrum().getShifts(0));
////            System.out.println("spectrum2:\t" + ssc2.getSubspectrum().getShifts(0));
////            System.out.println("spectrum2a:\t" + subspectrumToAddSSC2.getShifts(0));
////            System.out.println("spectrum3:\t" + combinedSpectrum.getShifts(0));
////            System.out.println("spectrum3 equ:\t" + combinedSpectrum.getEquivalences());
//
//        // if no valid spectrum could be built then go to next pairwise SSC comparison
//        if ((combinedSpectrum == null) || !Assembly.isValidSubspectrum(combinedSpectrum, querySpectrum, shiftTol, thrsMatchFactor)) {
////                System.out.println("-> no valid combined spectrum!");
//            return null;
//        }
//        System.out.println("-> !!!valid combined spectrum!!!");
//
//        // 3. substructure extension in SSC1
//        // 3.1 try to add atoms and bonds if some are missing
//        boolean extended = false;
//        if (!missingAtomIndicesSSC2.isEmpty()) {
//            mappedAtomIndices = Assembly.extendSSC(ssc1, ssc2, mappedAtomIndices, shiftTol);
//            if ((mappedAtomIndices == null) || mappedAtomIndices.isEmpty()) {
////                    System.out.println("---> could not extend ssc1 -> set to backup SSC!");
//                return null;
//            }
//            missingAtomIndicesSSC2 = Assembly.getMissingAtomIndices(ssc2, mappedAtomIndices.keySet());
//            System.out.println(" -> mapped atoms indices at the end:\t" + mappedAtomIndices);
//            System.out.println(" -> missing atom indices at the end:\t" + missingAtomIndicesSSC2);
//            if (!missingAtomIndicesSSC2.isEmpty()) {
////                System.out.println("---> could not add all missing atoms!");
//                return null;
//            }
//            extended = true;
//        }
//        // 3.2 check for potential missing bonds (i.e. to close rings)
//        if (!extended && !Assembly.addMissingBonds(ssc1, ssc2, mappedAtomIndices)) {
////                System.out.println("nothing added -> skip");
//            return null;
//        }
//        // the whole (sub)structure is saturated but the query spectrum is not fully covered -> skip
//        if (!ssc1.hasUnsaturatedAtoms() && (ssc1.getSubspectrum().getSignalCount() != querySpectrum.getSignalCount())) {
////                System.out.println("built a saturated but not valid structure -> skip");
//            return null;
//        }
//
//        return ssc1;
//    }
    
    public static HashMap<String, SSC> assembleBFS(final SSCLibrary rankedSSCLibrary, final long startSSCIndex, final int minMatchingSphereCount, 
            final Spectrum querySpectrum, final double thrsMatchFactor, final double shiftTol) throws Exception {
         
        final SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Absolute);
        String structureAsSMILES;
        final HashMap<String, SSC> solutions = new HashMap<>();
        SSC intermediate, newIntermediate, ssc2;
        // notice: clone (!!!) the SSC contents only; don't use the object (reference) itself because of modifications
        intermediate = rankedSSCLibrary.getSSC(startSSCIndex).getClone();
        // check whether the current SSC is already a final SSC
        if (Assembly.isFinalSSC(intermediate, querySpectrum, shiftTol, thrsMatchFactor)) {
            structureAsSMILES = smilesGenerator.create(intermediate.getSubstructure());
            if (!solutions.containsKey(structureAsSMILES)) {
                solutions.put(structureAsSMILES, intermediate);
                System.out.println("--> new solution found!!! -> " + solutions.size() + " -> " + structureAsSMILES);
                System.out.println("-> atom count: " + intermediate.getAtomCount() + ", bond count: " + intermediate.getBondCount());
                System.out.println("-> query spectrum:\t" + querySpectrum.getShifts(0));
                System.out.println("-> equivalences:\t" + querySpectrum.getEquivalences());
                System.out.println("-> pred. spectrum:\t" + intermediate.getSubspectrum().getShifts(0));
                System.out.println("-> equivalences:\t" + intermediate.getSubspectrum().getEquivalences());
            }
            return solutions;
        }
                
        LinkedHashSet<Long> path = new LinkedHashSet<>(), newPath;
        path.add(startSSCIndex);
        // create queue with initial state for BFS
        final Queue<Object[]> intermediates = new LinkedList<>();
        intermediates.add(new Object[]{intermediate, path});
        
        
        while (!intermediates.isEmpty()) {            
            intermediate = ((SSC) intermediates.peek()[0]).getClone();
            path = new LinkedHashSet<>((LinkedHashSet<Long>) intermediates.peek()[1]);
            System.out.println("--> for path: " + path + "\nnext ssc index: " + (Collections.max(path) + 1));
            System.out.println("--> ranked SSC indices: " + rankedSSCLibrary.getSSCIndices());
            for (long i = Collections.max(path) + 1; i < rankedSSCLibrary.getSSCCount(); i++) {

                ssc2 = rankedSSCLibrary.getSSC(i);

                System.out.println("\n\n-------------------------------- " + path + ", " + i + " --------------------------------");                                
                
//            backupSSC = intermediate.getClone();

                newIntermediate = Assembly.assemblyCore(intermediate.getClone(), ssc2, querySpectrum, minMatchingSphereCount, shiftTol, thrsMatchFactor);
                if (newIntermediate == null) {
//                intermediate = backupSSC.getClone();
                    continue;
                }

//                try {
//                    Utils.generatePicture(newIntermediate.getSubstructure(), "results/temp_" + path + "_" + i + ".png");
//                } catch (Exception e) {
//
//                }

                if (Assembly.isFinalSSC(newIntermediate, querySpectrum, shiftTol, thrsMatchFactor)) {
                    structureAsSMILES = smilesGenerator.create(newIntermediate.getSubstructure());
                    if (!solutions.containsKey(structureAsSMILES)) {
                        solutions.put(structureAsSMILES, newIntermediate);
                        System.out.println("--> new solution found!!! -> " + solutions.size() + " -> " + structureAsSMILES);
                        System.out.println("-> atom count: " + newIntermediate.getAtomCount() + ", bond count: " + newIntermediate.getBondCount());
                        System.out.println("-> query spectrum:\t" + querySpectrum.getShifts(0));
                        System.out.println("-> equivalences:\t" + querySpectrum.getEquivalences());
                        System.out.println("-> pred. spectrum:\t" + newIntermediate.getSubspectrum().getShifts(0));
                        System.out.println("-> equivalences:\t" + newIntermediate.getSubspectrum().getEquivalences());
                    }

//                intermediate = startSSC.getClone();//backupSSC1.getClone();
                    continue;
                }
                
                newPath = new LinkedHashSet<>(path);
                newPath.add(i);
                intermediates.add(new Object[]{newIntermediate, newPath});
            }
            
            intermediates.poll();
        }
        
        return solutions;
    }


    public static HashMap<String, SSC> assembleDFS(final SSCLibrary rankedSSCLibrary, final long startSSCIndex, final int minMatchingSphereCount,
                                                   final Spectrum querySpectrum, final double thrsMatchFactor, final double shiftTol, final String pathToOutputsFolder, final long querySpectrumCounter) throws Exception {

        final BufferedWriter bw = new BufferedWriter(new FileWriter(pathToOutputsFolder + "/results_" + querySpectrumCounter + "_temp.smiles"));
        final SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Absolute);
        String structureAsSMILES;
        final HashMap<String, SSC> solutions = new HashMap<>();
        SSC intermediate, newIntermediate, ssc2;
        // notice: clone (!!!) the SSC contents only; don't use the object (reference) itself because of modifications
        intermediate = rankedSSCLibrary.getSSC(startSSCIndex).getClone();
        // check whether the current SSC is already a final SSC
        if (Assembly.isFinalSSC(intermediate, querySpectrum, shiftTol, thrsMatchFactor)) {
            structureAsSMILES = smilesGenerator.create(intermediate.getSubstructure());
            if (!solutions.containsKey(structureAsSMILES)) {
                solutions.put(structureAsSMILES, intermediate);
                System.out.println("--> new solution found!!! -> " + solutions.size() + " -> " + structureAsSMILES);
                System.out.println("-> atom count: " + intermediate.getAtomCount() + ", bond count: " + intermediate.getBondCount());
                System.out.println("-> query spectrum:\t" + querySpectrum.getShifts(0));
                System.out.println("-> equivalences:\t" + querySpectrum.getEquivalences());
                System.out.println("-> pred. spectrum:\t" + intermediate.getSubspectrum().getShifts(0));
                System.out.println("-> equivalences:\t" + intermediate.getSubspectrum().getEquivalences());

                bw.append(structureAsSMILES);
                bw.newLine();
                bw.flush();
            }
            return solutions;
        }

        LinkedHashSet<Long> path = new LinkedHashSet<>(), newPath;
        path.add(startSSCIndex);
        // create queue with initial state for BFS
        final Stack<Object[]> intermediates = new Stack<>();
        intermediates.push(new Object[]{intermediate, path});

        long j = Collections.max(path) + 1;
        while (!intermediates.isEmpty()) {

            while ((j < rankedSSCLibrary.getSSCCount() - 1) && !intermediates.isEmpty()) {

                for (long i = j; i < rankedSSCLibrary.getSSCCount(); i++) {

                    if(intermediates.isEmpty()){
                        break;
                    }

                    intermediate = ((SSC) intermediates.peek()[0]).getClone();
                    path = new LinkedHashSet<>((LinkedHashSet<Long>) intermediates.peek()[1]);
                    System.out.println("--> for path: " + path + "\nnext ssc index: " + i + "/" + (rankedSSCLibrary.getSSCCount() - 1));
//                System.out.println("--> ranked SSC indices: " + rankedSSCLibrary.getSSCIndices());
                    ssc2 = rankedSSCLibrary.getSSC(i);

                    System.out.println("\n\n-------------------------------- " + path + ", " + i + " --------------------------------");

                    newIntermediate = Assembly.assemblyCore(intermediate.getClone(), ssc2, querySpectrum, minMatchingSphereCount, shiftTol, thrsMatchFactor);
                    if (newIntermediate == null) {
                        continue;
                    }

                try {
                    Utils.generatePicture(newIntermediate.getSubstructure(), "results/temp_" + path + "_" + i + ".png");
                } catch (Exception e) {

                }

                    if (Assembly.isFinalSSC(newIntermediate, querySpectrum, shiftTol, thrsMatchFactor)) {
                        structureAsSMILES = smilesGenerator.create(newIntermediate.getSubstructure());
                        if (!solutions.containsKey(structureAsSMILES)) {
                            solutions.put(structureAsSMILES, newIntermediate);
                            System.out.println("--> new solution found!!! -> " + solutions.size() + " -> " + structureAsSMILES);
                            System.out.println("-> atom count: " + newIntermediate.getAtomCount() + ", bond count: " + newIntermediate.getBondCount());
                            System.out.println("-> query spectrum:\t" + querySpectrum.getShifts(0));
                            System.out.println("-> equivalences:\t" + querySpectrum.getEquivalences());
                            System.out.println("-> pred. spectrum:\t" + newIntermediate.getSubspectrum().getShifts(0));
                            System.out.println("-> equivalences:\t" + newIntermediate.getSubspectrum().getEquivalences());

                            bw.append(structureAsSMILES);
                            bw.newLine();
                            bw.flush();

//                            return solutions;
                        }
                        continue;
                    }

//                    if(i < (rankedSSCLibrary.getSSCCount() - 1)){
                        newPath = new LinkedHashSet<>(path);
                        newPath.add(i);
                        intermediates.push(new Object[]{newIntermediate, newPath});
//                    }
                }
                j = Collections.max((LinkedHashSet<Long>) intermediates.peek()[1]) + 1;
                intermediates.pop();
            }

            if(intermediates.isEmpty()){
                break;
            }
            j = Collections.max((LinkedHashSet<Long>) intermediates.peek()[1]) + 1;
            intermediates.pop();
        }

        bw.close();

        return solutions;
    }
    
    
    public static HashMap<String, SSC> assembleSeq(final SSCLibrary rankedSSCLibrary, final long startSSCIndex, final int minMatchingSphereCount,
                                                   final Spectrum querySpectrum, final double thrsMatchFactor, final double shiftTol) throws Exception {

        final SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Absolute);
        String structureAsSMILES;
        final HashMap<String, SSC> solutions = new HashMap<>();
        SSC intermediate, backupSSC, ssc2, startSSC;
        // notice: clone (!!!) the SSC contents only; don't use the object (reference) itself because of modifications
        startSSC = rankedSSCLibrary.getSSC(startSSCIndex).getClone();
        intermediate = startSSC.getClone();
        intermediate.setIndex(startSSCIndex);
        // check whether the current SSC is already a final SSC
        if (Assembly.isFinalSSC(intermediate, querySpectrum, shiftTol, thrsMatchFactor)) {
            structureAsSMILES = smilesGenerator.create(intermediate.getSubstructure());
            solutions.put(structureAsSMILES, intermediate);
            System.out.println("--> new solution found!!! -> " + solutions.size() + " -> " + structureAsSMILES);
            System.out.println("-> atom count: " + intermediate.getAtomCount() + ", bond count: " + intermediate.getBondCount());
            System.out.println("-> query spectrum:\t" + querySpectrum.getShifts(0));
            System.out.println("-> equivalences:\t" + querySpectrum.getEquivalences());
            System.out.println("-> pred. spectrum:\t" + intermediate.getSubspectrum().getShifts(0));
            System.out.println("-> equivalences:\t" + intermediate.getSubspectrum().getEquivalences());

            if(solutions.containsKey(null)){
                System.out.println(" NULL key at 1");
            }
            if(solutions.containsValue(null)){
                System.out.println(" NULL value at 1");
            }
            return solutions;
        }
        
        for (long i = 0; i < rankedSSCLibrary.getSSCCount(); i++) {
            if(i == startSSCIndex){
                continue;
            }            
            
            System.out.println("\n\n-------------------------------- " + startSSCIndex + ", " + i + " --------------------------------");
            backupSSC = intermediate.getClone();
            ssc2 = rankedSSCLibrary.getSSC(i);
            intermediate = Assembly.assemblyCore(intermediate.getClone(), ssc2, querySpectrum, minMatchingSphereCount, shiftTol, thrsMatchFactor);
            if (intermediate == null) {
                intermediate = backupSSC.getClone();
                intermediate.setIndex(startSSCIndex);
                continue;
            }

            try {
                Utils.generatePicture(intermediate.getSubstructure(), "results/temp_" + startSSCIndex + "_" + i + ".png");
            } catch (Exception e) {

            }

            if (Assembly.isFinalSSC(intermediate, querySpectrum, shiftTol, thrsMatchFactor)) {
                structureAsSMILES = smilesGenerator.create(intermediate.getSubstructure());
                if (!solutions.containsKey(structureAsSMILES)) {
                    solutions.put(structureAsSMILES, intermediate);
                    System.out.println("--> new solution found!!! -> " + solutions.size() + " -> " + structureAsSMILES);
                    System.out.println("-> atom count: " + intermediate.getAtomCount() + ", bond count: " + intermediate.getBondCount());
                    System.out.println("-> query spectrum:\t" + querySpectrum.getShifts(0));
                    System.out.println("-> equivalences:\t" + querySpectrum.getEquivalences());
                    System.out.println("-> pred. spectrum:\t" + intermediate.getSubspectrum().getShifts(0));
                    System.out.println("-> equivalences:\t" + intermediate.getSubspectrum().getEquivalences());
                }

                intermediate = startSSC.getClone();//backupSSC1.getClone();
                intermediate.setIndex(startSSCIndex);
//                continue;
            }
            
        }

        if(solutions.containsKey(null)){
            System.out.println(" NULL key at 2");
        }
        if(solutions.containsValue(null)){
            System.out.println(" NULL value at 2");
        }
        return solutions;
    }











    private static HashMap<Integer, Integer> getAtomMappingsMCSS(final SSC ssc1, final SSC ssc2, final List<RMap> subgraph) throws CDKException {
        final HashMap<Integer, Integer> atomMappingsTemp = new HashMap<>(), atomMappings = new HashMap<>();
        final UniversalIsomorphismTester universalIsomorphismTester = new UniversalIsomorphismTester();
        final IAtomContainer overlap = UniversalIsomorphismTester.project(subgraph, ssc1.getSubstructure(),0);

        for (final RMap rMap : universalIsomorphismTester.getSubgraphAtomsMap(ssc1.getSubstructure(), overlap)){
            atomMappingsTemp.put(rMap.getId1(), rMap.getId2());
        }
        for (final RMap rMap : universalIsomorphismTester.getSubgraphAtomsMap(ssc2.getSubstructure(), overlap)){
            for (final int atomIndexSSC1 : atomMappingsTemp.keySet()) {
                if (atomMappingsTemp.get(atomIndexSSC1) == rMap.getId2()) {
                    atomMappings.put(atomIndexSSC1, rMap.getId1());
                    break;
                }
            }
        }

        return atomMappings;
    }

    private static boolean isValidAtomMapping(final SSC ssc1, final SSC ssc2, final int atomIndexSSC1, final int atomIndexSSC2, final double shiftTol){
        final Signal signalSSC1 = ssc1.getSubspectrum().getSignal(ssc1.getAssignments().getIndex(0, atomIndexSSC1));
        final Signal signalSSC2 = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getIndex(0, atomIndexSSC2));
        // check for similar signal properties
        // atoms with assigned signals
        if((signalSSC1 != null) && (signalSSC2 != null)){
            // check multiplicities
            if(!signalSSC1.getMultiplicity().equals(signalSSC2.getMultiplicity())){
                System.out.println("---> atom ssc1: " + atomIndexSSC1 + " (" + signalSSC1.getMultiplicity() + "), atom ssc2: " + atomIndexSSC2 + "(" + signalSSC2.getMultiplicity() + ") have not the same multiplicity!!!!");
                System.out.println("!!!current subgraph (structural overlap) is not valid -> skip!!!");
                return false;
            }
            // @TODO add solvent effect deviation to shiftTol?

            // check shift deviations
            if(Math.abs(signalSSC1.getShift(0) - signalSSC2.getShift(0)) > shiftTol){
                System.out.println("---> atom ssc1: " + atomIndexSSC1 + " (" + signalSSC1.getShift(0) + ", " + ssc1.getSubspectrum().getSolvent() + "), atom ssc2: " + atomIndexSSC2 + "(" + signalSSC2.getShift(0) + ", " + ssc2.getSubspectrum().getSolvent() + ") differ in their shifts too much!!!!");
                System.out.println("!!!current subgraph (structural overlap) is not valid -> skip!!!");
                return false;
            }
        } else if((ssc1.getSubstructure().getAtom(atomIndexSSC1).getImplicitHydrogenCount() != null)
                && (ssc2.getSubstructure().getAtom(atomIndexSSC2).getImplicitHydrogenCount() != null)
                && (ssc1.getSubstructure().getAtom(atomIndexSSC1).getImplicitHydrogenCount() != ssc2.getSubstructure().getAtom(atomIndexSSC2).getImplicitHydrogenCount())){
            // atoms without assigned signals
            System.out.println("---> hetero atom ssc1: " + atomIndexSSC1 + " (" + ssc1.getSubstructure().getAtom(atomIndexSSC1).getImplicitHydrogenCount() + ", " + ssc1.getSubspectrum().getSolvent() + "), atom ssc2: " + atomIndexSSC2 + "(" + ssc2.getSubstructure().getAtom(atomIndexSSC2).getImplicitHydrogenCount() + ", " + ssc2.getSubspectrum().getSolvent() + ") differ in their shifts too much!!!!");
            System.out.println("!!!current subgraph (structural overlap) is not valid -> skip!!!");
            return false;
        }

        return true;
    }

    private static HashMap<Integer, Integer> getInvalidAtomMappings(final SSC ssc1, final SSC ssc2, final HashMap<Integer, Integer> atomMappings, final double shiftTol){
        final HashMap<Integer, Integer> invalidAtomMappings = new HashMap<>();
        // check for same multiplicity for two mapped (carbon) atoms in current subgraph (atom mappings)
        for (final Entry<Integer, Integer> entry : atomMappings.entrySet()) {
            if(!Assembly.isValidAtomMapping(ssc1, ssc2, entry.getKey(), entry.getValue(), shiftTol)){
                invalidAtomMappings.put(entry.getKey(), entry.getValue());
            }
        }

        return invalidAtomMappings;
    }

    private static boolean isValidAtomMappings(final SSC ssc1, final SSC ssc2, final HashMap<Integer, Integer> atomMappings, final double shiftTol){

        return Assembly.getInvalidAtomMappings(ssc1, ssc2, atomMappings, shiftTol).isEmpty();
    }

    private static HashMap<Integer, Integer> reverseHashMap(final HashMap<Integer, Integer> map){
        final HashMap<Integer, Integer> reversedMap = new HashMap<>();
        for (final Entry<Integer, Integer> entry : map.entrySet()){
            reversedMap.put(entry.getValue(), entry.getKey());
        }

        return reversedMap;
    }

    private static HashMap<Integer, Integer> extendSSC(final SSC ssc1, final SSC ssc2, final HashMap<Integer, Integer> atomMappings, final double shiftTol) throws CDKException, CloneNotSupportedException{
        int mappedAtomIndexSSC1, mappedAtomIndexSSC2, equivalentSignalIndex;
        IBond bondToAdd;
        IAtom atomToAdd, parentAtomSSC1;
        Signal signalToAddSSC2;
        ConnectionTree connectionTreeToAddFromSSC2;
        final HashMap<Integer, Integer> reversedAtomMappings = Assembly.reverseHashMap(atomMappings);
        final HashMap<Integer, Integer> reversedAtomMappingsTemp = new HashMap<>(reversedAtomMappings);

        // for each atom pair map
        for (final Entry<Integer, Integer> entry : reversedAtomMappingsTemp.entrySet()) {
            mappedAtomIndexSSC2 = entry.getKey();
            mappedAtomIndexSSC1 = entry.getValue();
            // check whether the current mapped atom in SSC1 (to extend) is an open site (unsaturated);
            // if no then skip because there is no chance to add atoms and bonds
            if ((ssc1.isUnsaturated(mappedAtomIndexSSC1) == null) || !ssc1.isUnsaturated(mappedAtomIndexSSC1)) {
                continue;
            }
            System.out.println("\nmapping: " + mappedAtomIndexSSC1 + ", " + mappedAtomIndexSSC2
                    + "\n-> " + ssc1.getHOSECode(mappedAtomIndexSSC1) + "\n-> " + ssc2.getHOSECode(mappedAtomIndexSSC2));

            // BFS to build connection tree which contains atoms in SSC2 to add to SSC1;
            // all mapped atom indices in SSC2 are used as list of visited atoms and then used
            // during BFS for building the connection tree;
            // that means that connected but unmapped atoms in SSC2 should exist in the resulting
            // connection tree to add to SSC1
            connectionTreeToAddFromSSC2 = HOSECodeBuilder.buildConnectionTree(ssc2.getSubstructure(), mappedAtomIndexSSC2, null, new HashSet<>(reversedAtomMappings.keySet()));
            System.out.println(" -> BFS: " + connectionTreeToAddFromSSC2.toString()
                    + "\n -> maxSphere: " + connectionTreeToAddFromSSC2.getMaxSphere());
            // traverse connection tree and try to add missing (not mapped) atoms from SSC2 to SSC1
            for (int s = 1; s <= connectionTreeToAddFromSSC2.getMaxSphere(); s++) {
                // traverse via spheres of connection tree
                for (final ConnectionTreeNode connectedNodeInSphereToAddSSC2 : connectionTreeToAddFromSSC2.getNodesInSphere(s)) {

                    if (connectedNodeInSphereToAddSSC2.isRingClosureNode()) {
                        continue;
                    }

                    bondToAdd = connectedNodeInSphereToAddSSC2.getBondsToParents().get(0).clone();
                    System.out.println(" -> key: " + connectedNodeInSphereToAddSSC2.getParentNodes().get(0).getKey() + " -> " + reversedAtomMappings.get(connectedNodeInSphereToAddSSC2.getParentNodes().get(0).getKey()));
                    parentAtomSSC1 = ssc1.getSubstructure().getAtom(reversedAtomMappings.get(connectedNodeInSphereToAddSSC2.getParentNodes().get(0).getKey()));
                    System.out.println("in s: " + s + " -> bond: (" + connectedNodeInSphereToAddSSC2.getKey() + ") to " + reversedAtomMappings.get(connectedNodeInSphereToAddSSC2.getParentNodes().get(0).getKey()));

                    if (Assembly.isValidBondAddition(ssc1.getSubstructure(), parentAtomSSC1.getIndex(), bondToAdd)) {

                        atomToAdd = connectedNodeInSphereToAddSSC2.getAtom().clone();
                        ssc1.getSubstructure().addAtom(atomToAdd);

                        bondToAdd.setAtom(atomToAdd, 0);
                        bondToAdd.setAtom(parentAtomSSC1, 1);
                        if(bondToAdd.isAromatic()){
                            bondToAdd.setOrder(IBond.Order.SINGLE);
                        }
                        ssc1.getSubstructure().addBond(bondToAdd);
                        System.out.println("in extendSSC: new bond (" + bondToAdd.getOrder().numeric() + ", " + bondToAdd.isAromatic() + ") between " + parentAtomSSC1.getIndex() + ", " + (ssc1.getAtomCount() - 1) + " added!!!");

                        // add belonging signal from SSC2 to SSC1
                        if (atomToAdd.getSymbol().equals(ssc1.getSubspectrumAtomType())) {
                            signalToAddSSC2 = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getIndex(0, connectedNodeInSphereToAddSSC2.getKey()));
                            if (signalToAddSSC2 == null) {
                                return null;
                            }
                            equivalentSignalIndex = ssc1.getSubspectrum().pickClosestSignal(signalToAddSSC2.getShift(0), 0, shiftTol);
                            ssc1.getSubspectrum().addSignal(signalToAddSSC2, equivalentSignalIndex);
                            ssc1.getAssignments().addAssignment(new int[]{ssc1.getAtomCount() - 1});
                        }

                        reversedAtomMappings.put(connectedNodeInSphereToAddSSC2.getKey(), ssc1.getAtomCount() - 1);
                        System.out.println("-> added atom " + connectedNodeInSphereToAddSSC2.getKey() + " from SSC2 to SSC1");


                    } else {
                        System.out.println("not a valid bond addition!!!");
                        return null;
                    }
                }
            }
            ssc1.update();
        }
//        System.out.println("\n\n -> mapped atom indices at the end: \t" + mappedAtomIndices2);

        return Assembly.reverseHashMap(reversedAtomMappings);
    }

    public static List<RMap> getMaximumCommonSubgraph(final SSC ssc1, final SSC ssc2, final double shiftTol){
        final UniversalIsomorphismTester universalIsomorphismTester = new UniversalIsomorphismTester();
        final List<List<RMap>> subgraphList;
        try {
            subgraphList = universalIsomorphismTester.search(ssc1.getSubstructure(), ssc2.getSubstructure(), new BitSet(), new BitSet(), true, false);
        } catch (CDKException e) {
            return null;
        }

        int counter = 0;
        int maxValidSubgraphIndex = -1;
        HashMap<Integer, Integer> atomMappingsTemp;
        for (final List<RMap> subgraph : subgraphList){
            try {
                atomMappingsTemp = Assembly.getAtomMappingsMCSS(ssc1, ssc2, subgraph);
            } catch (CDKException e) {
                continue;
            }
//            System.out.println("atom mappings for subgraph " + counter + " -> " + atomMappingsTemp);
//            HashMap<Integer, Integer> invalidAtomMappings = Assembly.getInvalidAtomMappings(ssc1, ssc2, atomMappingsTemp, shiftTol);
//            System.out.println("-> invalid mappings                        -> " + invalidAtomMappings);

            /*
             @ToDo what if matched atoms are just in wrong order? (multiplicities/shifts are not correct but structural overlap itself is valid)
             -> use HOSE code information somehow to redirect/correct atom mappings ?
            */

            if((subgraph.size() > maxValidSubgraphIndex) && Assembly.isValidAtomMappings(ssc1, ssc2, atomMappingsTemp, shiftTol)){
                maxValidSubgraphIndex = counter;
//                System.out.println(" -> is valid? -> " + Assembly.isValidAtomMappings(ssc1, ssc2, atomMappingsTemp, shiftTol));
            }
            counter++;
        }
        if(maxValidSubgraphIndex < 0){
            System.out.println("-> !!! no valid subgraph found!!! -> skip");
            return null;
        }
        System.out.println(" --> subgraoh index to use: " + maxValidSubgraphIndex);

        return subgraphList.get(maxValidSubgraphIndex);
    }

    public static SSC assemblyCore(final SSC ssc1, final SSC ssc2, final Spectrum querySpectrum, final int minMatchingSphereCount, final double shiftTol, final double thrsMatchFactor) throws Exception {

//        // potential approach via spectral matching first (or only)
//        final Assignment assignment = Matcher.matchSpectra(ssc2.getSubspectrum(), ssc1.getSubspectrum(), 0, 0, shiftTol);
//        System.out.println(" --> match assignments: " + assignment.isFullyAssigned(0) + " -> " + assignment.getAssignments(0));
////        System.out.println(" ---> atom pairs: ");
//        HashMap<Integer, Integer> atomMappingsTemp = new HashMap<>();
//        for (int i = 0; i < assignment.getAssignmentsCount(); i++) { // iterate over ssc2 signal count via assignments
////            System.out.println(ssc2.getAssignments().getAssignment(0, i) + " - " + ssc1.getAssignments().getAssignment(0, assignment.getAssignment(0, i)));
//            if(ssc1.getAssignments().getAssignment(0, assignment.getAssignment(0, i)) != null){
//                atomMappingsTemp.put(ssc1.getAssignments().getAssignment(0, assignment.getAssignment(0, i)), ssc2.getAssignments().getAssignment(0, i));
//            }
//        }
//        for (final Entry<Integer, Integer> entry : atomMappingsTemp.entrySet()){
//            System.out.println("-> mapped atoms     : " + entry.getKey() + " - " + entry.getValue());
//            if(ssc1.isUnsaturated(entry.getKey())){
//                System.out.println("-> is unsaturated!");
//            }
//        }
//        HashMap<Integer, Integer> atomMappings = Assembly.getAtomMappingsHOSECode(ssc1, ssc2, Assembly.getOverlapsHOSECodeCore(ssc1, ssc2, atomMappingsTemp, minMatchingSphereCount, shiftTol), atomMappingsTemp);



//        // 1. check for partial structural identity (overlaps) via MCSS
//        // atom mapping by maximum common (valid) subgraph
//        HashMap<Integer, Integer> atomMappings = Assembly.getAtomMappingsMCSS(ssc1, ssc2, Assembly.getMaximumCommonSubgraph(ssc1, ssc2, shiftTol));


        HashMap<Integer, ArrayList<Integer[]>> overlapsHOSECodeNew = Assembly.getOverlapsHOSECodeNew(ssc1, ssc2, minMatchingSphereCount, shiftTol);
        if(overlapsHOSECodeNew.isEmpty()){
            return ssc1;
        }

        final ArrayList<SSC> validSSCExtensions = new ArrayList<>();
        ArrayList<Integer[]> overlapsHOSECodeInSphere;
        // for each max. matching sphere (key); starting with highest
        for (int s = Collections.max(overlapsHOSECodeNew.keySet()); s >= minMatchingSphereCount; s--){
            if(!overlapsHOSECodeNew.keySet().contains(s)){
                continue;
            }

            overlapsHOSECodeInSphere = overlapsHOSECodeNew.get(s);
            System.out.println(" --> s = " + s + ": ");
            for (final Integer[] indices : overlapsHOSECodeInSphere){
                System.out.println(" -> " + Arrays.toString(indices));
            }

            // for each overlapping atom pairs in SSC1 and SSC2 in sphere (maybe in a certain order?)
            for (int k = 0; k < overlapsHOSECodeInSphere.size(); k++) {
                // reset ssc1Extended to original SSC1
                SSC ssc1Extended = ssc1.getClone();
                HashMap<Integer, Integer> atomMappingsTemp = new HashMap<>();
                System.out.println("\n -> k: " + k);
                System.out.println(" -> HOSE code maps:" + Arrays.toString(overlapsHOSECodeInSphere.get(k)));
                // indices of matched (root) atoms
                int i = overlapsHOSECodeInSphere.get(k)[0];
                int j = overlapsHOSECodeInSphere.get(k)[1];

                // @TODO check this constraint again
                if(ssc1Extended.isUnsaturated(i) || ssc2.isUnsaturated(j)){
                    System.out.println(" atom " + i + " in SSC1 or atom " + j + " in SSC2 is unsaturated and not allowed as overlap root");
                    continue;
                }

                ConnectionTree completeConnectionTreeSSC2 = HOSECodeBuilder.buildConnectionTree(ssc2.getSubstructure(), j, null);
                ConnectionTree maxSphereConnectionTreeSSC2 = HOSECodeBuilder.buildConnectionTree(ssc2.getSubstructure(), j, ssc2.getMaxSphere());
                ConnectionTree matchingConnectionTreeSSC1 = HOSECodeBuilder.buildConnectionTree(ssc1Extended.getSubstructure(), i, s);
                ConnectionTree matchingConnectionTreeSSC2 = HOSECodeBuilder.buildConnectionTree(ssc2.getSubstructure(), j, s);
                System.out.println(" --> atoms in SSC1: " + matchingConnectionTreeSSC1.getKeys(true));
                System.out.println(" --> atoms in SSC2: " + matchingConnectionTreeSSC2.getKeys(true));
                System.out.println(ssc1Extended.getHOSECode(i) + "\n" + ssc2.getHOSECode(j));
                System.out.println(ssc1Extended.getConnectionTree(i) + " -> " + ssc1Extended.getConnectionTree(i).getNodesCountInSphere(s) + "\n" + ssc2.getConnectionTree(j) + " -> " + ssc2.getConnectionTree(j).getNodesCountInSphere(s));
                ArrayList<Integer> unsaturatedAtomsSSC1 = new ArrayList<>();
                ArrayList<Integer> connectionTreeKeysSSC1 = new ArrayList<>(matchingConnectionTreeSSC1.getKeys(false));
                ArrayList<Integer> connectionTreeKeysSSC2 = new ArrayList<>(matchingConnectionTreeSSC2.getKeys(false));
                if(connectionTreeKeysSSC1.size() != connectionTreeKeysSSC2.size()){
                    System.out.println(" <--- max. matching sphere is not the same anymore!!!! --->");
                    continue;
                }
                for (int l = 0; l < connectionTreeKeysSSC1.size(); l++) {
                    atomMappingsTemp.put(connectionTreeKeysSSC1.get(l), connectionTreeKeysSSC2.get(l));
                }
                for (final int key : connectionTreeKeysSSC1){
                    if(key < 0){
                        continue;
                    }
                    if(ssc1Extended.isUnsaturated(key)){
                        unsaturatedAtomsSSC1.add(key);
                    }
                }
                if(unsaturatedAtomsSSC1.isEmpty()){
                    continue;
                }
                Signal signalToAdd;
                System.out.println("predicted spectrum orig.: " + ssc1Extended.getSubspectrum().getShifts(0));
                System.out.println(" -----> unsaturated atoms SSC1:" + unsaturatedAtomsSSC1);

                // check for each found unsaturated atom in SSC1 whether there is a valid extension possible from SSC2 via its children
                for (final int unsaturatedAtomKeySSC1 : unsaturatedAtomsSSC1){
                    int indexUntilMaxMatchingSphere = connectionTreeKeysSSC1.indexOf(unsaturatedAtomKeySSC1);
                    int nodeKeyInCompleteConnectionTreeSSC1 = connectionTreeKeysSSC1.get(indexUntilMaxMatchingSphere);
                    int nodeKeyInCompleteConnectionTreeSSC2 = connectionTreeKeysSSC2.get(indexUntilMaxMatchingSphere);
                    System.out.println(" ---> unsaturated atoms in max. matching sphere: " + nodeKeyInCompleteConnectionTreeSSC1 + " - " + nodeKeyInCompleteConnectionTreeSSC2);
                    
                    ArrayList<ConnectionTreeNode> childNodesToAppend = maxSphereConnectionTreeSSC2.getNode(nodeKeyInCompleteConnectionTreeSSC2).getChildNodes();
                    for (final ConnectionTreeNode childNodeToAppend : childNodesToAppend){
                        System.out.println(" --> to append from SSC2: " + childNodeToAppend.getKey());

                        // if ring closure node
                        if(childNodeToAppend.isRingClosureNode()){

//                            int index = connectionTreeKeysSSC2.indexOf(childNodeToAppend.getParentNodes().get(0).getKey());
//                            int parentKeyInSSC1 = connectionTreeKeysSSC1.get(index);
//                            int parentKeyInSSC2 = connectionTreeKeysSSC2.get(index);
//
//                            IBond bondToAdd = completeConnectionTreeSSC2.getBond(parentKeyInSSC2, childNodeToAppend.getKey()).clone();
//
//                            if(Assembly.isValidBondAddition(ssc1Extended.getSubstructure(), parentKeyInSSC1, bondToAdd)){
//
////                                System.out.println(" ->  bond addition for ring closure " + parentKeyInSSC1
////                                        + " (" + parentKeyInSSC2 + " - " + childNodeToAppend.getKey() + ") : " + " valid!!!");
//                                // parent node of child node shouldn have more than one parent node;
//                                // we now search for them in complete conn. tree of SSC2 and check the connection
//                                // between the grand parent and parent to complete the ring closure by adding the bond
//
////                                for (final ConnectionTreeNode grandParentNodeSSC2 : completeConnectionTreeSSC2.getNode(parentKeyInSSC2).getParentNodes()){
////                                    if(atomMappingsTemp.containsValue(parentKeyInSSC2)
////                                            && atomMappingsTemp.containsValue(grandParentNodeSSC2.getKey())){
////                                        int parentAtomIndexSSC1 = -1, grandParentAtomIndexSSC1 = -1;
////                                        for (final Entry<Integer, Integer> entry : atomMappingsTemp.entrySet()){
////                                            if(entry.getValue() == parentKeyInSSC2){
////                                                parentAtomIndexSSC1 = entry.getKey();
////                                            } else if(entry.getValue() == grandParentNodeSSC2.getKey()) {
////                                                grandParentAtomIndexSSC1 = entry.getKey();
////                                            }
////                                            if((parentAtomIndexSSC1 != -1)
////                                                    && (grandParentAtomIndexSSC1 != -1)){
////                                                break;
////                                            }
////                                        }
////                                        if(((parentAtomIndexSSC1 != -1)
////                                                && (grandParentAtomIndexSSC1 != -1)
////                                                && sscToExtend.getSubstructure().getBond(sscToExtend.getSubstructure().getAtom(parentAtomIndexSSC1), sscToExtend.getSubstructure().getAtom(grandParentAtomIndexSSC1)) == null)){
////
////                                            if(Assembly.isValidBondAddition(sscToExtend.getSubstructure(), grandParentAtomIndexSSC1, bondToAdd)){
////                                                bondToAdd.setAtom(sscToExtend.getSubstructure().getAtom(parentAtomIndexSSC1), 0);
////                                                bondToAdd.setAtom(sscToExtend.getSubstructure().getAtom(grandParentAtomIndexSSC1), 1);
////                                                sscToExtend.getSubstructure().addBond(bondToAdd);
////                                            }
////                                        }
////                                    }
////                                }
//                            }


                        } else {
                            // 1. add signal from SSC2 to SSC1 and check for predicted (sub)spectrum validation
                            // 2. check valid bond addition; add atom and bond to parent
                            // 3. take connection tree of nodeKeyInCompleteConnectionTreeSSC2 as root and add all non-visited child nodes
                            //  -> add atom and bond of child nodes as well as ring closures




//                            if(completeConnectionTreeSSC2.getBond(nodeKeyInCompleteConnectionTreeSSC2, childNodeToAppend.getKey()) == null){
//                                System.out.println("completeConnectionTreeSSC2: " + completeConnectionTreeSSC2);
//                                System.out.println("!!! bond between " + nodeKeyInCompleteConnectionTreeSSC2 + " and " + childNodeToAppend.getKey() + " does not exist -> !!!! SKIP !!!!");
//                                continue;
//                            }
                            IBond bondToAdd = completeConnectionTreeSSC2.getBond(nodeKeyInCompleteConnectionTreeSSC2, childNodeToAppend.getKey()).clone();
                            if(Assembly.isValidBondAddition(ssc1Extended.getSubstructure(), nodeKeyInCompleteConnectionTreeSSC1, bondToAdd)){

                                ConnectionTree substructureTreeToAddSSC2 = hose.Utils.getSubtree(completeConnectionTreeSSC2, childNodeToAppend.getKey());
                                System.out.println(" JOOO: " + (substructureTreeToAddSSC2 == null));
                                System.out.println(" JOOO 2: " + (substructureTreeToAddSSC2.getNodesCount()));
                                System.out.println("-> subtree would be: " + substructureTreeToAddSSC2 + "\n" + substructureTreeToAddSSC2.getKeys(true));
                                int counter = 1;
                                for (final int substructureTreeNodeKeySSC2 : substructureTreeToAddSSC2.getKeys(true)) {
                                    atomMappingsTemp.put((ssc1Extended.getAtomCount() - 1) + counter, substructureTreeNodeKeySSC2);
                                    counter++;
                                }
                                System.out.println(" -> atom mappings temp: " + atomMappingsTemp);

                                IAtomContainer substructureToAddSSC2 = HOSECodeBuilder.buildAtomContainer(substructureTreeToAddSSC2);
                                System.out.println("subtree size: " + substructureToAddSSC2.getAtomCount());
                                ssc1Extended.getSubstructure().add(substructureToAddSSC2);

                                bondToAdd.setAtom(ssc1Extended.getSubstructure().getAtom((ssc1Extended.getAtomCount() - 1) - (substructureToAddSSC2.getAtomCount() - 1)), 0);
                                bondToAdd.setAtom(ssc1Extended.getSubstructure().getAtom(nodeKeyInCompleteConnectionTreeSSC1), 1);
                                System.out.println("is valid bond addition? -> true");
                                ssc1Extended.getSubstructure().addBond(bondToAdd);

                                // add signals and assignments from added subtructure tree in SSC2 to SSC1
                                for (final int nodeKeyInSubtreeSSC2 : substructureTreeToAddSSC2.getKeys(true)){
                                    signalToAdd = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getIndex(0, nodeKeyInSubtreeSSC2));
                                    if(signalToAdd != null){
                                        ssc1Extended.getSubspectrum().addSignal(signalToAdd);
                                        int index = new ArrayList<>(substructureTreeToAddSSC2.getKeys(true)).indexOf(nodeKeyInSubtreeSSC2);
                                        index += ((ssc1Extended.getAtomCount() - 1) - (substructureToAddSSC2.getAtomCount() - 1));
                                        ssc1Extended.getAssignments().addAssignment(new int[]{index});
                                    }
                                }

                                ConnectionTree completeConnectionTreeSSC1 = HOSECodeBuilder.buildConnectionTree(ssc1Extended.getSubstructure(), i, null);
                                // try to close rings directly from added subtree of SSC2
                                ConnectionTreeNode nodeSSC1, parentNodeSSC1;
                                // for each node in added subtree
                                for (final int nodeKeyInSubtreeSSC2 : substructureTreeToAddSSC2.getKeys(true)){
                                    // if node is ring closure point
                                    if(completeConnectionTreeSSC2.getNode(nodeKeyInSubtreeSSC2).getParentNodes().size() > 1){
                                        // then check for each (additional) parent node whether there is a bond missing
                                        for (final ConnectionTreeNode parentNodeSSC2 : completeConnectionTreeSSC2.getNode(nodeKeyInSubtreeSSC2).getParentNodes()){
//                                            // skip current parent node (ancestor of root of subtree)
//                                            if(parentNodeSSC2.getKey() == childNodeToAppend.getKey()){
//                                                continue;
//                                            }
                                            System.out.println(" FOR " + nodeKeyInSubtreeSSC2 + " -> PARENT NODE: " + parentNodeSSC2.getKey());
                                            nodeSSC1 = null;
                                            parentNodeSSC1 = null;
                                            if(atomMappingsTemp.containsValue(parentNodeSSC2.getKey())){
                                                System.out.println(" HUHU");
                                                for (final Entry<Integer, Integer> entry : atomMappingsTemp.entrySet()){
                                                    if(entry.getValue() == parentNodeSSC2.getKey()){
                                                        System.out.println(" HUHU 1");
                                                        parentNodeSSC1 = completeConnectionTreeSSC1.getNode(entry.getKey());
                                                    }
                                                    if(entry.getValue() == nodeKeyInSubtreeSSC2){
                                                        System.out.println(" HUHU 2");
                                                        nodeSSC1 = completeConnectionTreeSSC1.getNode(entry.getKey());
                                                    }
                                                    if((nodeSSC1 != null) && (parentNodeSSC1 != null)){
                                                        System.out.println(" HUHU 3");
                                                        System.out.println(" --> EQU: " + nodeSSC1.getKey() + " - " + parentNodeSSC1.getKey());
                                                        break;
                                                    }
                                                }
                                            }
                                            if((nodeSSC1 != null) && (parentNodeSSC1 != null)
                                                    && (ssc1Extended.getSubstructure().getBond(
                                                    ssc1Extended.getSubstructure().getAtom(nodeSSC1.getKey()),
                                                    ssc1Extended.getSubstructure().getAtom(parentNodeSSC1.getKey())) == null)){
                                                bondToAdd = ssc2.getSubstructure().getBond(ssc2.getSubstructure().getAtom(nodeKeyInSubtreeSSC2), ssc2.getSubstructure().getAtom(parentNodeSSC2.getKey())).clone();
                                                if(Assembly.isValidBondAddition(ssc1Extended.getSubstructure(), nodeSSC1.getKey(), bondToAdd)
                                                        && Assembly.isValidBondAddition(ssc1Extended.getSubstructure(), parentNodeSSC1.getKey(), bondToAdd)){
                                                    bondToAdd.setAtom(ssc1Extended.getSubstructure().getAtom(nodeSSC1.getKey()), 0);
                                                    bondToAdd.setAtom(ssc1Extended.getSubstructure().getAtom(parentNodeSSC1.getKey()), 1);
                                                    ssc1Extended.getSubstructure().addBond(bondToAdd);

                                                    System.out.println(" --> new BOND: " + nodeSSC1.getKey() + " - " + parentNodeSSC1.getKey());
                                                }
                                            }
                                        }
                                    }
                                }


                            } else {
                                System.out.println("is valid bond addition? -> false");
                            }

                        }
                    }
                }
                System.out.println("predicted spectrum modi.: " + ssc1Extended.getSubspectrum().getShifts(0));

                // for current HOSE code matching atom pair
                // if a valid substructure and also subspectrum could be assembled
                if(Assembly.isValidSubspectrum(ssc1Extended.getSubspectrum(), querySpectrum, shiftTol, thrsMatchFactor)){

                    ssc1Extended.update();
                    if((ssc1Extended.getAtomCount() == ssc1.getAtomCount())
                            && (ssc1Extended.getBondCount() == ssc1.getBondCount())){
                        continue;
                    }

                    System.out.println("\nvalid SSC built: ");
                    System.out.println("substructure: " + ssc1Extended.getSubstructure().getAtomCount());
                    System.out.println("subspectrum : " + ssc1Extended.getSubspectrum().getShifts(0));
                    System.out.println("assignments : " + ssc1Extended.getAssignments().getAssignments(0) + "\n");
                    System.out.println("unsaturated atoms: " + ssc1Extended.getUnsaturatedAtomIndices());

                    validSSCExtensions.add(ssc1Extended);

//                    Utils.generatePicture(ssc1Extended.getSubstructure(), "results/temp_" + 0 + "_" + ssc2.getIndex() + " -> " + i + "_" + j + ".png");
                }
            }
            System.out.println("\n");
        }

        if(validSSCExtensions.isEmpty()){
            return null;
        }

        validSSCExtensions.sort((validExtendedSSC1, validExtendedSSC2) -> {
            final int atomCountComp = -1 * Integer.compare(validExtendedSSC1.getAtomCount(), validExtendedSSC2.getAtomCount());
            if(atomCountComp != 0){
                return atomCountComp;
            }
            final int bondCountComp = -1 * Integer.compare(validExtendedSSC1.getBondCount(), validExtendedSSC2.getBondCount());
            if(bondCountComp != 0){
                return bondCountComp;
            }

            return  Double.compare(Matcher.calculateAverageDeviation(validExtendedSSC1.getSubspectrum(), querySpectrum, 0, 0, shiftTol), Matcher.calculateAverageDeviation(validExtendedSSC2.getSubspectrum(), querySpectrum, 0, 0, shiftTol));
        });

        for (int i = 0; i < validSSCExtensions.size(); i++){
            System.out.println(" -------> for valid extension " + i + " with size " + validSSCExtensions.get(i).getAtomCount() + " and " + validSSCExtensions.get(i).getBondCount() + " and " + Matcher.calculateAverageDeviation(validSSCExtensions.get(i).getSubspectrum(), querySpectrum, 0, 0, shiftTol));
        }

        return validSSCExtensions.get(0);//return ssc1;





//        // 1. check for partial structural identity (overlaps) via HOSE codes
//        // map structural overlaps in all spheres separately and
//        // insertion of all valid mappings in different spheres into one map
//        final HashMap<Integer, Integer> preAtomMappings = new HashMap<>();
//        HashMap<Integer, Integer> atomMappings = Assembly.getAtomMappingsHOSECode(ssc1, ssc2, Assembly.getOverlapsHOSECode(ssc1, ssc2, minMatchingSphereCount, shiftTol), preAtomMappings);
//        System.out.println("-> mapped atoms     : " + atomMappings);
//
//        // @TODO use parameter "shiftTol" from solvent effect deviations?
//        if(atomMappings.isEmpty() || !Assembly.isValidAtomMappings(ssc1, ssc2, atomMappings, shiftTol)){
//            System.out.println(" no (valid) atom mappings found!!! -> skip");
//            return null;
//        }
//        HashSet<Integer> missingAtomsSSC2 = Assembly.getMissingAtomIndices(ssc2, new HashSet<>(atomMappings.values()));
//        System.out.println("missing atoms: " + missingAtomsSSC2);
//
//        int mappedAtomIndexSSC2;
//        // spectra combination and validation
//        // 2.1 for subspectra combination: add signals from SSC2's subspectrum associated with mapping (overlapping) atoms to a new subspectrum to combine with subspectrum of SSC1
//        final LinkedHashSet<Integer> signalIndicesToIgnoreFromSubspectrumSSC2 = new LinkedHashSet<>();
//        for (final Entry<Integer, Integer> entry : atomMappings.entrySet()) {
//            mappedAtomIndexSSC2 = entry.getValue();
//            if (ssc2.getSubstructure().getAtom(mappedAtomIndexSSC2).getSymbol().equals(ssc1.getSubspectrumAtomType())) {
//                signalIndicesToIgnoreFromSubspectrumSSC2.add(ssc2.getAssignments().getIndex(0, mappedAtomIndexSSC2));
//            }
//        }
//        final Spectrum subspectrumToAddSSC2 = new Spectrum(ssc2.getSubspectrum().getNuclei());
//        for (int s = 0; s < ssc2.getSubspectrum().getSignalCount(); s++) {
//            if (!signalIndicesToIgnoreFromSubspectrumSSC2.contains(s)) {
//                subspectrumToAddSSC2.addSignal(ssc2.getSubspectrum().getSignal(s).getClone(), ssc2.getSubspectrum().getEquivalence(s));
//            }
//        }
//        // 2.2 combine subspectra of SSC to extend and matched SSC and validate it
//        final Spectrum combinedSpectrum = Matcher.combineSpectra(ssc1.getSubspectrum(), subspectrumToAddSSC2, 0, 0, Start.EQUIV_SIGNAL_THRS);
////            System.out.println("\nspectrum1:\t" + ssc1.getSubspectrum().getShifts(0));
////            System.out.println("spectrum2:\t" + ssc2.getSubspectrum().getShifts(0));
////            System.out.println("spectrum2a:\t" + subspectrumToAddSSC2.getShifts(0));
////            System.out.println("spectrum3:\t" + combinedSpectrum.getShifts(0));
////            System.out.println("spectrum3 equ:\t" + combinedSpectrum.getEquivalences());
//
//        // if no valid spectrum could be built then go to next pairwise SSC comparison
//        if ((combinedSpectrum == null) || !Assembly.isValidSubspectrum(combinedSpectrum, querySpectrum, shiftTol, thrsMatchFactor)) {
//                System.out.println("-> no valid combined spectrum!");
//            return null;
//        }
//        System.out.println("-> !!!valid combined spectrum!!!");
//
//        // 3. substructure extension in SSC1
//        // 3.1 try to add atoms and bonds if some are missing
//        boolean extended = false;
//        if (!missingAtomsSSC2.isEmpty()) {
//            atomMappings = Assembly.extendSSC(ssc1, ssc2, atomMappings, shiftTol);
//            if ((atomMappings == null) || atomMappings.isEmpty()) {
//                    System.out.println("---> could not extend ssc1 -> set to backup SSC!");
//                return null;
//            }
//            missingAtomsSSC2 = Assembly.getMissingAtomIndices(ssc2, new HashSet<>(atomMappings.values()));
//            System.out.println(" -> mapped atoms indices at the end:\t" + atomMappings);
//            System.out.println(" -> missing atom indices at the end:\t" + missingAtomsSSC2);
//            if (!missingAtomsSSC2.isEmpty()) {
//                System.out.println("---> could not add all missing atoms!");
//                return null;
//            }
//            extended = true;
//        }
//        // 3.2 check for potential missing bonds (i.e. to close rings)
//        if (!extended && !Assembly.addMissingBonds(ssc1, ssc2, atomMappings)) {
//                System.out.println("nothing added -> skip");
//            return null;
//        }
//        // the whole (sub)structure is saturated but the query spectrum is not fully covered -> skip
//        if (!ssc1.hasUnsaturatedAtoms() && (ssc1.getSubspectrum().getSignalCount() != querySpectrum.getSignalCount())) {
//                System.out.println("built a saturated but not valid structure -> skip");
//            return null;
//        }
//
//        return ssc1;
    }



    //    private static HashMap<Integer, Integer> extendSSC(final SSC ssc1, final SSC ssc2, final HashMap<Integer, Integer> mappedAtomIndices, final double shiftTol) throws CDKException, CloneNotSupportedException{
//        int mappedAtomIndexSSC1, mappedAtomIndexSSC2, equivalentSignalIndex;
//        IBond bondToAdd;
//        IAtom atomToAdd, parentAtomSSC1, atomInRingSSC1_1, atomInRingSSC1_2, atomSSC2;
//        Signal signalToAddSSC2;
//        ConnectionTree connectionTreeToAddSSC2;
//        final HashMap<Integer, Integer> mappedAtomIndices2 = new HashMap<>(mappedAtomIndices);
//
//        // for each atom pair map
//        for (final Entry<Integer, Integer> entry : mappedAtomIndices.entrySet()) {
//            mappedAtomIndexSSC2 = entry.getKey();
//            mappedAtomIndexSSC1 = entry.getValue();
//            // check whether the current mapped atom in SSC1 (to extend) is an open site (unsaturated);
//            // if no then skip because there is no chance to add atoms and bonds
//            if (/*(ssc1.isUnsaturated(mappedAtomIndexSSC1) == null) || */!ssc1.isUnsaturated(mappedAtomIndexSSC1)) {
//                continue;
//            }
//            System.out.println("\nmapping: " + mappedAtomIndexSSC2 + ", " + mappedAtomIndexSSC1
//                    + "\n-> " + ssc1.getHOSECode(mappedAtomIndexSSC1) + "\n-> " + ssc2.getHOSECode(mappedAtomIndexSSC2));
//
//            // BFS to build connection tree which contains atoms in SSC2 to add to SSC1;
//            // all mapped atom indices in SSC2 are used as list of visited atoms and then used
//            // during BFS for building the connection tree;
//            // that means that connected but unmapped atoms in SSC2 should exist in the resulting
//            // connection tree to add to SSC1
//            connectionTreeToAddSSC2 = HOSECodeBuilder.buildConnectionTree(ssc2.getSubstructure(), mappedAtomIndexSSC2, null, new HashSet<>(mappedAtomIndices2.keySet()));
//            System.out.println(" -> BFS: " + connectionTreeToAddSSC2.toString()
//                    + "\n -> maxSphere: " + connectionTreeToAddSSC2.getMaxSphere());
//            // traverse connection tree and try to add missing (not mapped) atoms from SSC2 to SSC1
//            for (int s = 1; s <= connectionTreeToAddSSC2.getMaxSphere(); s++) {
//                // traverse via spheres of connection tree
//                for (final ConnectionTreeNode connectedNodeInSphereToAddSSC2 : connectionTreeToAddSSC2.getNodesInSphere(s)) {
//                    // in case of ring closure node try to close a ring
//                    if (connectedNodeInSphereToAddSSC2.isRingClosureNode()) {
//                        System.out.println("ringClosureNode at: " + connectedNodeInSphereToAddSSC2.getKey());
//                        ConnectionTreeNode parentNodeSSC2 = connectedNodeInSphereToAddSSC2.getParentNodes().get(0);
//                        ConnectionTreeNode ringClosureParentNodeSSC2 = connectedNodeInSphereToAddSSC2.getParentNodes().get(0).getParentNodes().get(1);
//                        if (!mappedAtomIndices2.containsKey(parentNodeSSC2.getKey())
//                                || !mappedAtomIndices2.containsKey(ringClosureParentNodeSSC2.getKey())) {
//                            continue;
//                        }
//                        // try to add the missing bond in ring
//                        bondToAdd = parentNodeSSC2.getBondsToParents().get(1).clone();
//                        if (Assembly.isValidBondAddition(ssc1.getSubstructure(), mappedAtomIndices2.get(parentNodeSSC2.getKey()), bondToAdd)
//                                && Assembly.isValidBondAddition(ssc1.getSubstructure(), mappedAtomIndices2.get(ringClosureParentNodeSSC2.getKey()), bondToAdd)) {
//                            bondToAdd.setAtom(ssc1.getSubstructure().getAtom(mappedAtomIndices2.get(parentNodeSSC2.getKey())), 0);
//                            bondToAdd.setAtom(ssc1.getSubstructure().getAtom(mappedAtomIndices2.get(ringClosureParentNodeSSC2.getKey())), 1);
//                            if(bondToAdd.isAromatic()){
//                                bondToAdd.setOrder(IBond.Order.SINGLE);
//                            }
//                            ssc1.getSubstructure().addBond(bondToAdd);
//                        }
//
//                        continue;
//                    }
//                    bondToAdd = connectedNodeInSphereToAddSSC2.getBondsToParents().get(0).clone();
//                    System.out.println(" -> key: " + connectedNodeInSphereToAddSSC2.getParentNodes().get(0).getKey() + " -> " + mappedAtomIndices2.get(connectedNodeInSphereToAddSSC2.getParentNodes().get(0).getKey()));
//                    parentAtomSSC1 = ssc1.getSubstructure().getAtom(mappedAtomIndices2.get(connectedNodeInSphereToAddSSC2.getParentNodes().get(0).getKey()));
//                    System.out.println("in s: " + s + " -> bond: (" + connectedNodeInSphereToAddSSC2.getKey() + ") to " + mappedAtomIndices2.get(connectedNodeInSphereToAddSSC2.getParentNodes().get(0).getKey()));
//
//                    if (Assembly.isValidBondAddition(ssc1.getSubstructure(), parentAtomSSC1.getIndex(), bondToAdd)) {
//
//                        atomToAdd = connectedNodeInSphereToAddSSC2.getAtom().clone();
//                        ssc1.getSubstructure().addAtom(atomToAdd);
//
//                        bondToAdd.setAtom(atomToAdd, 0);
//                        bondToAdd.setAtom(parentAtomSSC1, 1);
//                        if(bondToAdd.isAromatic()){
//                            bondToAdd.setOrder(IBond.Order.SINGLE);
//                        }
//                        ssc1.getSubstructure().addBond(bondToAdd);
//                        System.out.println("in extendSSC: new bond (" + bondToAdd.getOrder().numeric() + ", " + bondToAdd.isAromatic() + ") between " + parentAtomSSC1.getIndex() + ", " + (ssc1.getAtomCount() - 1) + " added!!!");
//
//                        // add belonging signal from SSC2 to SSC1
//                        if (atomToAdd.getSymbol().equals(ssc1.getSubspectrumAtomType())) {
//                            signalToAddSSC2 = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getSignalIndex(0, connectedNodeInSphereToAddSSC2.getKey()));
//                            if (signalToAddSSC2 == null) {
//                                return null;
//                            }
//                            equivalentSignalIndex = ssc1.getSubspectrum().pickClosestSignal(signalToAddSSC2.getShift(0), 0, shiftTol);
//                            ssc1.getSubspectrum().addSignal(signalToAddSSC2, equivalentSignalIndex);
//                            ssc1.getAssignments().addAssignment(new int[]{ssc1.getAtomCount() - 1});
//                        }
//
//                        mappedAtomIndices2.put(connectedNodeInSphereToAddSSC2.getKey(), ssc1.getAtomCount() - 1);
//                        System.out.println("-> added atom " + connectedNodeInSphereToAddSSC2.getKey() + " from SSC2 to SSC1");
//
//                        // ring closure directly from connection tree information
//                        atomSSC2 = ssc2.getSubstructure().getAtom(connectedNodeInSphereToAddSSC2.getKey());
//                        for (final IAtom connectedAtomSSC2 : ssc2.getSubstructure().getConnectedAtomsList(atomSSC2)) {
//                            if ((connectedAtomSSC2.getIndex() != connectedNodeInSphereToAddSSC2.getParentNodes().get(0).getKey())
//                                    && mappedAtomIndices2.containsKey(connectedAtomSSC2.getIndex())
//                                    && !connectionTreeToAddSSC2.containsKey(connectedAtomSSC2.getIndex())) {
//                                System.out.println("connected but mapped atom found -> ring closure!!!"
//                                        + "\n -> " + connectedNodeInSphereToAddSSC2.getKey() + ", " + connectedAtomSSC2.getIndex());
//                                bondToAdd = ssc2.getSubstructure().getBond(atomSSC2, connectedAtomSSC2).clone();
//                                if (Assembly.isValidBondAddition(ssc1.getSubstructure(), ssc1.getAtomCount() - 1, bondToAdd)
//                                        && Assembly.isValidBondAddition(ssc1.getSubstructure(), mappedAtomIndices2.get(connectedAtomSSC2.getIndex()), bondToAdd)) {
//                                    System.out.println("!!! ring closure !!!");
//                                    atomInRingSSC1_1 = ssc1.getSubstructure().getAtom(ssc1.getAtomCount() - 1);
//                                    atomInRingSSC1_2 = ssc1.getSubstructure().getAtom(mappedAtomIndices2.get(connectedAtomSSC2.getIndex()));
//
//                                    bondToAdd.setAtom(atomInRingSSC1_1, 0);
//                                    bondToAdd.setAtom(atomInRingSSC1_2, 1);
//                                    if(bondToAdd.isAromatic()){
//                                        bondToAdd.setOrder(IBond.Order.SINGLE);
//                                    }
//                                    ssc1.getSubstructure().addBond(bondToAdd);
//                                }
//                            }
//                        }
//                    } else {
//                        System.out.println("not a valid bond addition!!!");
//                        return null;
//                    }
//                }
//            }
//            ssc1.update();
//        }
////        System.out.println("\n\n -> mapped atom indices at the end: \t" + mappedAtomIndices2);
//
//        return mappedAtomIndices2;
//    }
}
