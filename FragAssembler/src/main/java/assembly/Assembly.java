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
import match.Match;
import model.ConnectionTree;
import model.ConnectionTreeNode;
import model.SSC;
import model.SSCLibrary;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.Bond;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Assembly {  
    
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
    public static HashMap<Integer, HashMap<Integer, Integer>> getStructuralOverlapsHOSECode(final SSC ssc1, final SSC ssc2, final int minSphereMatchCount, final double shiftTol) throws CDKException {
        final HashMap<Integer, HashMap<Integer, Integer>> mappedAtomIndicesInSpheres = new HashMap<>();

        String HOSECodeSSC1, HOSECodeSSC2;
        Signal signalRootAtomSSC1, signalRootAtomSSC2;
        int maxMatchingSphere, mappedAtomIndexSSC1, mappedAtomIndexSSC2;
        ConnectionTree connectionTreeSSC1;
        ConnectionTree connectionTreeSSC2;
        Signal signalSSC1, signalSSC2;
        ArrayList<Integer> nodeKeysInSphereToCheckSSC1, nodeKeysInSphereToCheckSSC2;
        HashMap<Integer, Integer> mappedEqualNodesInSphere;
        boolean stop, notValid, containsUnsaturatedAtomsSSC2;        
        for (int i = 0; i < ssc1.getAtomCount(); i++) {
            for (int j = 0; j < ssc2.getAtomCount(); j++) {
                // check for same atom types
                if(!ssc1.getSubstructure().getAtom(i).getSymbol().equals(ssc2.getSubstructure().getAtom(j).getSymbol())){
                    continue;
                }
                // init max. matching sphere as default to -1 (no match)
                maxMatchingSphere = -1;
                // get signals of both atoms to compare (if available) and check for same signal properties
                // to set a starting point with higher probability to be correct
                signalRootAtomSSC1 = ssc1.getSubspectrum().getSignal(ssc1.getAssignments().getSignalIndex(0, i));
                signalRootAtomSSC2 = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getSignalIndex(0, j));                      
                if ((signalRootAtomSSC1 != null) && (signalRootAtomSSC2 != null) 
                        && (!signalRootAtomSSC1.getMultiplicity().equals(signalRootAtomSSC2.getMultiplicity()) || Math.abs(signalRootAtomSSC1.getShift(0) - signalRootAtomSSC2.getShift(0)) > shiftTol)) {
                    continue;                    
                }
                System.out.println("\n-> i: " + i + ", j: " + j);
                // check for pure structural identity via HOSE code
                for (int s = 0; s <= Integer.min(ssc1.getMaxSphere(), ssc2.getMaxSphere()); s++) {                    
                    HOSECodeSSC1 = HOSECodeBuilder.buildHOSECode(ssc1.getSubstructure(), i, s, false);
                    HOSECodeSSC2 = HOSECodeBuilder.buildHOSECode(ssc2.getSubstructure(), j, s, false);

                    if(!HOSECodeSSC1.equals(HOSECodeSSC2)){
                        break;                        
                    }
                    maxMatchingSphere = s;
                    System.out.println(" --> in s: " + s + " -> " + HOSECodeSSC1 + " vs. " + HOSECodeSSC2);                                        
                }    
                System.out.println("-> maxMatchingSphere: " + maxMatchingSphere);
                // skip non-matching atom pairs
                if (maxMatchingSphere == -1) {
                    continue;
                }
                // skip atom pairs which do not match at least in given minimum number of spheres
                if((maxMatchingSphere < minSphereMatchCount) 
                        // but atoms at open sites are allowed because those sites could enable 
                        // an extension if no other atom mappings exist
                        && (!ssc1.isUnsaturated(i) && !ssc2.isUnsaturated(j))
                        ){
                        continue;
                }    
                // if a atom pair matched succesfully via HOSE codes (structurally) then check now for spectral identity/similarity
                connectionTreeSSC1 = HOSECodeBuilder.buildConnectionTree(ssc1.getSubstructure(), i, maxMatchingSphere);
                connectionTreeSSC2 = HOSECodeBuilder.buildConnectionTree(ssc2.getSubstructure(), j, maxMatchingSphere);
                notValid = false;
                // check all spheres until max. matching sphere for spectral similarity
                for (int s = 0; s <= maxMatchingSphere; s++) {
                    nodeKeysInSphereToCheckSSC1 = connectionTreeSSC1.getNodeKeysInSphere(s);
                    nodeKeysInSphereToCheckSSC2 = connectionTreeSSC2.getNodeKeysInSphere(s);
                    System.out.println("in s: " + s + " -> node keys SSC1: " + nodeKeysInSphereToCheckSSC1);
                    System.out.println("in s: " + s + " -> node keys SSC2: " + nodeKeysInSphereToCheckSSC2);
                    stop = false;
                     // check for same/similar signal properties in same order as the belonging HOSE codes until max. matching sphere
                    for (int k = 0; k < nodeKeysInSphereToCheckSSC1.size(); k++) {
                        signalSSC1 = ssc1.getSubspectrum().getSignal(ssc1.getAssignments().getSignalIndex(0, nodeKeysInSphereToCheckSSC1.get(k)));
                        signalSSC2 = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getSignalIndex(0, nodeKeysInSphereToCheckSSC2.get(k)));
                        if ((signalSSC1 != null) && (signalSSC2 != null)
                                && (!signalSSC1.getMultiplicity().equals(signalSSC2.getMultiplicity()) || Math.abs(signalSSC1.getShift(0) - signalSSC2.getShift(0)) > shiftTol)) {
                            stop = true;
                            break;
                        }                                               
                    }
                    // it could happen that a structural identity until a certain sphere is there but the actual order of atoms 
                    // within a HOSE code is different than in the second HOSE codes;
                    // then one could search for a possible atom mapping via spectral properties 
                    // within a sphere in different order than HOSE code
                    if (stop) {
                        mappedEqualNodesInSphere = Match.mapEqualNodesInSphere(ssc1, ssc2, i, j, s, shiftTol);
                        System.out.println("in s: " + s + " --> maps in sphere: " + mappedEqualNodesInSphere);
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
                    System.out.println("-> maxMatchingSphere is NOT valid!!!");
                    continue;
                }
                System.out.println("-> maxMatchingSphere " + maxMatchingSphere + " is valid!!!");
                // create new key for the found max. matching sphere if it's not existing
                if (!mappedAtomIndicesInSpheres.containsKey(maxMatchingSphere)) {
                    mappedAtomIndicesInSpheres.put(maxMatchingSphere, new HashMap<>());
                }

//                connectionTreeSSC1 = HOSECodeBuilder.buildConnectionTree(ssc1.getSubstructure(), i, maxMatchingSphere);
//                connectionTreeSSC2 = HOSECodeBuilder.buildConnectionTree(ssc2.getSubstructure(), j, maxMatchingSphere);
//                ArrayList<Integer> nodeKeysSSC1 = connectionTreeSSC1.getKeysInOrder(true);
//                ArrayList<Integer> nodeKeysSSC2 = connectionTreeSSC2.getKeysInOrder(true);                
//                for (int k = 0; k < nodeKeysSSC1.size(); k++) {
//                    mappedAtomIndexSSC1 = nodeKeysSSC1.get(k);
//                    mappedAtomIndexSSC2 = nodeKeysSSC2.get(k);
//                    if (!mappedAtomIndicesInSpheres.get(maxMatchingSphere).containsKey(mappedAtomIndexSSC1)
//                            && !mappedAtomIndicesInSpheres.get(maxMatchingSphere).containsValue(mappedAtomIndexSSC2)) {                        
//                        
//                        
////                        // check whether at least one mapped atom in SSC2 is unsaturated
////                        containsUnsaturatedAtomsSSC2 = false;
////                        for (final int nodeKey : ssc2.getConnectionTree(mappedAtomIndexSSC2).getKeysInOrder(true)) {
////                            if (ssc2.isUnsaturated(nodeKey)) {
////                                containsUnsaturatedAtomsSSC2 = true;
////                                break;
////                            }
////                        }
////                        if (!containsUnsaturatedAtomsSSC2) {
////                            System.out.println("\nmap index ssc1: " + mappedAtomIndexSSC1
////                                    + ", map index ssc2: " + mappedAtomIndexSSC2
////                                    + " -> nothing unsaturated -> skip this overlap");
////                            continue;
////                        }
//                        
//                        
//                        mappedAtomIndicesInSpheres.get(maxMatchingSphere).put(mappedAtomIndexSSC1, mappedAtomIndexSSC2);
//                    }
//                }
                
                // insert matching atom pair into mappings hashmap with found max. matching sphere
                mappedAtomIndicesInSpheres.get(maxMatchingSphere).put(i, j);
            }
        }
                        
        System.out.println("all maps: " + mappedAtomIndicesInSpheres);
        
//        // insertion of all valid mappings in different spheres into one map, 
//        // starting with highest matching sphere (higher chance of correctness)
//        if(mappedAtomIndicesInSpheres.isEmpty()){
//            return new HashMap<>();
//        }
//        if(mappedAtomIndicesInSpheres.size() == 1){
//            return mappedAtomIndicesInSpheres.get(new ArrayList<>(mappedAtomIndicesInSpheres.keySet()).get(0));            
//        } 
//        
//        final HashMap<Integer, Integer> mappedAtomIndices = new HashMap<>();
//        ArrayList<Integer> nodeKeysSSC1, nodeKeysSSC2;
//        int rootMappedAtomIndexSSC1, rootMappedAtomIndexSSC2;
//        // for all found spheres try to insert the mappings
//        for (int m = Collections.max(mappedAtomIndicesInSpheres.keySet()); m >= Collections.min(mappedAtomIndicesInSpheres.keySet()); m--) {
//            if (!mappedAtomIndicesInSpheres.containsKey(m)) {
//                continue;
//            }
//            // for each atom pair as valid starting point insert its belonging mapped atom pairs until its max. matching sphere
//            for (final Map.Entry<Integer, Integer> entry : mappedAtomIndicesInSpheres.get(m).entrySet()) {
//                rootMappedAtomIndexSSC1 = entry.getKey();
//                rootMappedAtomIndexSSC2 = entry.getValue();
//                connectionTreeSSC1 = HOSECodeBuilder.buildConnectionTree(ssc1.getSubstructure(), rootMappedAtomIndexSSC1, m);
//                connectionTreeSSC2 = HOSECodeBuilder.buildConnectionTree(ssc2.getSubstructure(), rootMappedAtomIndexSSC2, m);
//                nodeKeysSSC1 = connectionTreeSSC1.getKeysInOrder(true);
//                nodeKeysSSC2 = connectionTreeSSC2.getKeysInOrder(true);
//                for (int k = 0; k < nodeKeysSSC1.size(); k++) {
//                    mappedAtomIndexSSC1 = nodeKeysSSC1.get(k);
//                    mappedAtomIndexSSC2 = nodeKeysSSC2.get(k);
//                    if (!mappedAtomIndices.containsKey(mappedAtomIndexSSC1)
//                            && !mappedAtomIndices.containsValue(mappedAtomIndexSSC2)) {
//
//                        // check whether at least one mapped atom in SSC2 (SSC to extend) is unsaturated
//                        containsUnsaturatedAtomsSSC2 = false;
//                        for (final int nodeKey : ssc2.getConnectionTree(mappedAtomIndexSSC2).getKeysInOrder(true)) {
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
//                        
//                        mappedAtomIndices.put(mappedAtomIndexSSC1, mappedAtomIndexSSC2);
//                    }
//                }                                        
//            }                        
//        }
//        
//        return mappedAtomIndices;

        return mappedAtomIndicesInSpheres;
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
     * @param mappedAtomIndicesInSpheres HashMap containing the maximum 
     * matching spheres as keys and the belonging mapped atom pairs as values
     * @return 
     * @throws CDKException
     * 
     * @see #getStructuralOverlapsHOSECode(model.SSC, model.SSC, int, double) 
     */
    public static HashMap<Integer, Integer> mapStructuralOverlapsHOSECode(final SSC ssc1, final SSC ssc2, final HashMap<Integer, HashMap<Integer, Integer>> mappedAtomIndicesInSpheres) throws CDKException{
        // insertion of all valid mappings in different spheres into one map, 
        // starting with highest matching sphere (higher chance of correctness)
        if (mappedAtomIndicesInSpheres.isEmpty()) {
            return new HashMap<>();
        }
        if (mappedAtomIndicesInSpheres.size() == 1) {
            return mappedAtomIndicesInSpheres.get(new ArrayList<>(mappedAtomIndicesInSpheres.keySet()).get(0));
        }

        final HashMap<Integer, Integer> mappedAtomIndices = new HashMap<>();
        ArrayList<Integer> nodeKeysSSC1, nodeKeysSSC2;
        int rootMappedAtomIndexSSC1, rootMappedAtomIndexSSC2, mappedAtomIndexSSC1, mappedAtomIndexSSC2;
        ConnectionTree connectionTreeSSC1, connectionTreeSSC2;
        boolean containsUnsaturatedAtomsSSC2;
        // for all found spheres try to insert the mappings
        for (int m = Collections.max(mappedAtomIndicesInSpheres.keySet()); m >= Collections.min(mappedAtomIndicesInSpheres.keySet()); m--) {
            if (!mappedAtomIndicesInSpheres.containsKey(m)) {
                continue;
            }
            // for each atom pair as valid starting point insert its belonging mapped atom pairs until its max. matching sphere
            for (final Map.Entry<Integer, Integer> entry : mappedAtomIndicesInSpheres.get(m).entrySet()) {
                rootMappedAtomIndexSSC1 = entry.getKey();
                rootMappedAtomIndexSSC2 = entry.getValue();
                connectionTreeSSC1 = HOSECodeBuilder.buildConnectionTree(ssc1.getSubstructure(), rootMappedAtomIndexSSC1, m);
                connectionTreeSSC2 = HOSECodeBuilder.buildConnectionTree(ssc2.getSubstructure(), rootMappedAtomIndexSSC2, m);
                nodeKeysSSC1 = connectionTreeSSC1.getKeysInOrder(true);
                nodeKeysSSC2 = connectionTreeSSC2.getKeysInOrder(true);
                for (int k = 0; k < nodeKeysSSC1.size(); k++) {
                    mappedAtomIndexSSC1 = nodeKeysSSC1.get(k);
                    mappedAtomIndexSSC2 = nodeKeysSSC2.get(k);
                    if (!mappedAtomIndices.containsKey(mappedAtomIndexSSC1)
                            && !mappedAtomIndices.containsValue(mappedAtomIndexSSC2)) {

                        // check whether at least one mapped atom in SSC2 (SSC to extend) is unsaturated
                        containsUnsaturatedAtomsSSC2 = false;
                        for (final int nodeKey : ssc2.getConnectionTree(mappedAtomIndexSSC2).getKeysInOrder(true)) {
                            if (ssc2.isUnsaturated(nodeKey)) {
                                containsUnsaturatedAtomsSSC2 = true;
                                break;
                            }
                        }
                        if (!containsUnsaturatedAtomsSSC2) {
                            System.out.println("\nmap index ssc1: " + mappedAtomIndexSSC1
                                    + ", map index ssc2: " + mappedAtomIndexSSC2
                                    + " -> nothing unsaturated -> skip this overlap");
                            continue;
                        }

                        mappedAtomIndices.put(mappedAtomIndexSSC1, mappedAtomIndexSSC2);
                    }
                }
            }
        }

        return mappedAtomIndices;
    }
    
    
    /**
     *
     * @param HOSECode1
     * @param HOSECode2
     * @return
     * 
     * @deprecated 
     */
    public static int getMaximumMatchingSphereHOSECode(final String HOSECode1, final String HOSECode2) {
        final ArrayList<String> HOSECodeIntoSpheresSSC1 = HOSECodeBuilder.splitHOSECodeIntoSpheres(HOSECode1);
        final ArrayList<String> HOSECodeIntoSpheresSSC2 = HOSECodeBuilder.splitHOSECodeIntoSpheres(HOSECode2);
        int maxMatchingSphere = -1;
        for (int s = 0; s < Integer.min(HOSECodeIntoSpheresSSC1.size(), HOSECodeIntoSpheresSSC2.size()); s++) {
            if (!HOSECodeIntoSpheresSSC1.get(s).equals(HOSECodeIntoSpheresSSC2.get(s))) {
                break;
            }
            maxMatchingSphere = s;
        }
        return maxMatchingSphere;
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
    
    public static boolean isValidSubspectrum(final Spectrum subspectrum, final Spectrum querySpectrum, final double pickPrecisiom, final double thrsMatchFactor){
        if(subspectrum == null){
            return false;
        }
        final Assignment matchAssignments = Match.matchSpectra(subspectrum, querySpectrum, pickPrecisiom);    
        // filter for unset assignments
        if (!matchAssignments.isFullyAssigned(0)){//matchAssignments.getSetAssignmentsCount(0) < matchAssignments.getAssignmentsCount()) {
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
        if (Match.getMatchFactor(subspectrum, querySpectrum, pickPrecisiom) > thrsMatchFactor) {
//            System.out.println("-> match factor not allowed!!!");
            return false;
        }
        
        return true;
    }    
        
    
    public static HashMap<String, SSC> assemble(final long nStartSSCs, final int nThreads, final SSCLibrary rankedSSCLibrary, final int minMatchingSphereCount, 
            final Spectrum querySpectrum, final double shiftTol, final double thrsMatchFactor, final double pickPrecision) throws InterruptedException{
        final HashMap<String, SSC> solutions = new HashMap<>();
        // initialize an executor for parallelization
        final ExecutorService executor = Utils.initExecuter(nThreads);
        final ArrayList<Callable<HashMap<String, SSC>>> callables = new ArrayList<>();
        // add all task to do
        for (int i = 0; i < nStartSSCs; i++) {
            final int j = i;
            callables.add((Callable<HashMap<String, SSC>>) () -> {                                
                return Assembly.assembleHOSECode(rankedSSCLibrary, j, minMatchingSphereCount, querySpectrum, shiftTol, thrsMatchFactor, pickPrecision);
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
    
    /**
     *
     * @param rankedSSCLibrary
     * @param startSSCIndex
     * @param minMatchingSphereCount
     * @param querySpectrum
     * @param shiftTol
     * @param thrsMatchFactor
     * @param pickPrecision
     * @return
     * @throws CloneNotSupportedException
     * @throws CDKException
     * @throws IOException
     * 
     * @deprecated 
     */
    public static HashMap<String, SSC> assemble(final SSCLibrary rankedSSCLibrary, final int startSSCIndex, final int minMatchingSphereCount, 
            final Spectrum querySpectrum, final double shiftTol, final double thrsMatchFactor, final double pickPrecision) throws CloneNotSupportedException, CDKException, IOException {
        
        final SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Absolute);
//        final SmilesParser smilesParser = new SmilesParser(new SilentChemObjectBuilder());
        String structureAsSMILES;
        final HashMap<String, SSC> solutions = new HashMap<>();
        SSC ssc2, ssc1Backup;
        // clone (!!!) the SSC1 contents only; don't use the object (reference) itself because of modifications
        SSC ssc1 = rankedSSCLibrary.getSSC(startSSCIndex).getClone();
        
        int maxMatchingSphere, mappedAtomIndexSSC1, mappedAtomIndexSSC2;
        Signal signalToRemoveSSC2;
        Spectrum subspectrumSSC2Clone, combinedSpectrum;
        HashMap<Integer, Integer> mappedAtomIndices;
        HashSet<Integer> missingAtomIndicesSSC2;
        HashMap<Integer, ArrayList<Integer[]>> overlaps;
        
        // check whether the current SSC is already a final SSC
        if (Assembly.isFinalSSC(ssc1, querySpectrum)) {
            structureAsSMILES = smilesGenerator.create(ssc1.getSubstructure());
            solutions.put(structureAsSMILES, ssc1);
            System.out.println("--> new solution found!!! -> " + solutions.size() + " -> " + structureAsSMILES);
            System.out.println("-> atom count: " + ssc1.getAtomCount() + ", bond count: " + ssc1.getBondCount());
            System.out.println("-> query spectrum:\t" + querySpectrum.getShifts(0));
            System.out.println("-> equivalences:\t" + querySpectrum.getEquivalences());
            System.out.println("-> pred. spectrum:\t" + ssc1.getSubspectrum().getShifts(0));
            System.out.println("-> equivalences:\t" + ssc1.getSubspectrum().getEquivalences());
            return solutions;
        }
        
        
        for (final long ssc2Index : rankedSSCLibrary.getSSCIndices()) {
            if (ssc2Index == startSSCIndex) {
                continue;
            }
            ssc2 = rankedSSCLibrary.getSSC(ssc2Index);

            System.out.println("\n\n-------------------------------- " + startSSCIndex + ", " + ssc2Index + " --------------------------------");

            
            
            
            
            
//            System.out.println("\n\nMCSS:\n");
//            final UniversalIsomorphismTester universalIsomorphismTester = new UniversalIsomorphismTester();
//            List<IAtomContainer> structuralOverlaps = universalIsomorphismTester.getOverlaps(ssc1.getSubstructure(), ssc2.getSubstructure());
//            structuralOverlaps.sort(new Comparator<IAtomContainer>() {
//                @Override
//                public int compare(final IAtomContainer ac1, final IAtomContainer ac2) {
//                    return -1 * Integer.compare(ac1.getAtomCount(), ac2.getAtomCount());
//                }
//            });
//            System.out.println("-> " + startSSCIndex + ", " + ssc2Index + " -> " + structuralOverlaps.size());
//    
//            mappedAtomIndices = new HashMap<>();
//            for (int j = 0; j < structuralOverlaps.size(); j++){
//                System.out.println(" -> j: " + j);
//                Utils.generatePicture(ssc1.getSubstructure(), "/Users/mwenk/Downloads/outputs/overlap_" + startSSCIndex + "_" + ssc2Index + "_" + j + "_p1.png");
//                Utils.generatePicture(ssc2.getSubstructure(), "/Users/mwenk/Downloads/outputs/overlap_" + startSSCIndex + "_" + ssc2Index + "_" + j + "_p2.png");
//                Utils.generatePicture(structuralOverlaps.get(j), "/Users/mwenk/Downloads/outputs/overlap_" + startSSCIndex + "_" + ssc2Index + "_" + j + ".png");
//                List<RMap> atomMapsSSC1 = universalIsomorphismTester.getSubgraphAtomsMap(ssc1.getSubstructure(), structuralOverlaps.get(j));
//                List<RMap> atomMapsSSC2 = universalIsomorphismTester.getSubgraphAtomsMap(ssc2.getSubstructure(), structuralOverlaps.get(j));
//                for (final RMap atomMapSSC1 : atomMapsSSC1) {
//                    for (final RMap atomMapSSC2 : atomMapsSSC2) {
//                        if (atomMapSSC1.getId2() == atomMapSSC2.getId2()
//                                && !mappedAtomIndices.containsKey(atomMapSSC2.getId1())
//                                && !mappedAtomIndices.containsValue(atomMapSSC1.getId1())) {
//                            mappedAtomIndices.put(atomMapSSC2.getId1(), atomMapSSC1.getId1());
//                        }
//                    }
//
//                }
//                missingAtomIndicesSSC2 = Assembly.getMissingAtomIndices(ssc2, mappedAtomIndices.keySet());
//                System.out.println("----> mapped atom indices: " + mappedAtomIndices);
//                System.out.println("----> missing atom indices SSC2: " + missingAtomIndicesSSC2);
//            }
//            System.out.println("\nMCSS ends:\n\n");
//            if (true) {
//                continue;
//            }
            
            
            
            
            
            
            // 1. check for partial structural identity (overlaps); via HOSE code or connection tree comparison                                        
            overlaps = Assembly.getStructuralOverlaps(ssc2, ssc1, shiftTol, minMatchingSphereCount, ssc1.getSubspectrumAtomType());
            // if there is no structural identity then skip that SSC pair comparison            
            if (overlaps.isEmpty()) {
                continue;
            }            
            ssc1Backup = ssc1.getClone();
            subspectrumSSC2Clone = ssc2.getSubspectrum().getClone();                                      
            
            
            maxMatchingSphere = Collections.max(overlaps.keySet());
            System.out.println("-> maxMatchingSphere: " + maxMatchingSphere);
            
            // 1.1 map structural overlaps
            mappedAtomIndices = mapStructuralOverlaps(ssc1, ssc2, overlaps, shiftTol);
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
            combinedSpectrum = Match.combineSpectra(ssc1.getSubspectrum(), subspectrumSSC2Clone, pickPrecision);
//            System.out.println("\nspectrum1:\t" + ssc1.getSubspectrum().getShifts(0));
//            System.out.println("spectrum2:\t" + ssc2.getSubspectrum().getShifts(0));
//            System.out.println("spectrum2a:\t" + subspectrumSSC2Clone.getShifts(0));
//            System.out.println("spectrum3:\t" + combinedSpectrum.getShifts(0));
//            System.out.println("spectrum3 equ:\t" + combinedSpectrum.getEquivalences());

            // if no valid spectrum could be built then go to next pairwise SSC comparison
            if ((combinedSpectrum == null) || !Assembly.isValidSubspectrum(combinedSpectrum, querySpectrum, shiftTol, thrsMatchFactor)) {
//                System.out.println("-> no valid combined spectrum!");
                continue;
            }
            System.out.println("-> !!!valid combined spectrum!!!");

//            System.out.println("entering substructure extension!!!");
            // 3. substructure extension in SSC1
            // 3. a) try to add atoms and bonds because atoms are missing
            if (!missingAtomIndicesSSC2.isEmpty()) {
                mappedAtomIndices = Assembly.extendSSC(ssc1, ssc2, mappedAtomIndices, pickPrecision);
                if (mappedAtomIndices == null) {
                    System.out.println("---> could not extend ssc1 -> set to backup SSC!");
                    ssc1 = ssc1Backup.getClone();
                    continue;
                }
                if (mappedAtomIndices.isEmpty()) {
                    ssc1 = ssc1Backup.getClone();
                    continue;
                }

                missingAtomIndicesSSC2 = Assembly.getMissingAtomIndices(ssc2, mappedAtomIndices.keySet());
                System.out.println(" -> mapped atoms indices  at the end:\t" + mappedAtomIndices);
                System.out.println(" -> missing atom indices  at the end:\t" + missingAtomIndicesSSC2);
                if (!missingAtomIndicesSSC2.isEmpty()) {
//                System.out.println("---> could not add all missing atoms!");
                    ssc1 = ssc1Backup.getClone();
                    continue;
                }
            }
            // 3. b) check for potential bond additions
            Assembly.checkBonds(ssc1, ssc2, mappedAtomIndices);
            
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
            
            
//            Utils.generatePicture(smilesParser.parseSmiles(smilesGenerator.create(ssc1.getSubstructure())), "outputs/temp_" + startSSCIndex + "_" + ssc2Index + ".png");
//            Utils.generatePicture(ssc1.getSubstructure(), "results/temp_" + startSSCIndex + "_" + ssc2Index + ".png");
            
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
                }
                ssc1 = ssc1Backup.getClone();
                continue;
            }
            
//            // in case there are two hetero atoms the same (added two times) and unsaturated, then one could remove the first atom and link it to the neighbor(s) atom(s) of the second one
//            // this step should not be necessary anymore if one use the molecula formula of the unknown compound to elucidate
//            if((ssc1.getSubspectrum().getSignalCount() == querySpectrum.getSignalCount())){
//                System.out.println("query spectrum covered but unsaturated atoms are left: " + ssc1.getUnsaturatedAtomIndices());
//                System.out.println(" -> unsaturated atoms: " + ssc1.getUnsaturatedAtomIndices());
//                if(ssc1.getUnsaturatedAtomIndices().size() == 2){
//                    int unsaturatedAtomIndex1 = ssc1.getUnsaturatedAtomIndices().get(0);
//                    int unsaturatedAtomIndex2 = ssc1.getUnsaturatedAtomIndices().get(1);                    
//                    ConnectionTreeNode node1 = ssc1.getConnectionTree(unsaturatedAtomIndex1).getRootNode();
//                    ConnectionTreeNode node2 = ssc1.getConnectionTree(unsaturatedAtomIndex2).getRootNode();
//                    Signal signal1 = ssc1.getSubspectrum().getSignal(ssc1.getAssignments().getSignalIndex(0, unsaturatedAtomIndex1));
//                    Signal signal2 = ssc1.getSubspectrum().getSignal(ssc1.getAssignments().getSignalIndex(0, unsaturatedAtomIndex2));
//                    if(Match.isEqualNode(node1, node2, signal1, signal2, shiftTol)){
//                        System.out.println(" -> link!!!");
//                        IAtom unsaturatedAtom1 = ssc1.getSubstructure().getAtom(unsaturatedAtomIndex1);
//                        IAtom unsaturatedAtom2 = ssc1.getSubstructure().getAtom(unsaturatedAtomIndex2);
//                        for (final IBond bond2 : ssc1.getSubstructure().getConnectedBondsList(unsaturatedAtom2)) {
//                            if(isValidBondAddition(ssc1.getSubstructure(), unsaturatedAtomIndex1, bond2)){
//                                bond2.setAtom(bond2.getOther(unsaturatedAtom2), 0);   
//                                bond2.setAtom(unsaturatedAtom1, 1);
//                            }
//                        }
//                        ssc1.getSubstructure().removeAtom(unsaturatedAtom2);
//                        ssc1.update();
//                        
//                        if (Assembly.isFinalSSC(ssc1, querySpectrum)) {
//                            structureAsSMILES = smilesGenerator.create(ssc1.getSubstructure());
//                            if (!solutions.containsKey(structureAsSMILES)) {
//                                solutions.put(structureAsSMILES, ssc1);
//                                System.out.println("--> new solution found!!! -> " + solutions.size() + " -> " + structureAsSMILES);
//                                System.out.println("-> atom count: " + ssc1.getAtomCount() + ", bond count: " + ssc1.getBondCount());
//                                System.out.println("-> query spectrum:\t" + querySpectrum.getShifts(0));
//                                System.out.println("-> equivalences:\t" + querySpectrum.getEquivalences());
//                                System.out.println("-> pred. spectrum:\t" + ssc1.getSubspectrum().getShifts(0));
//                                System.out.println("-> equivalences:\t" + ssc1.getSubspectrum().getEquivalences());
//                            }
//                            break;
//                        }
//                    }
//                }
//            }
        }
        
        return solutions;
    }
    
    public static boolean isFinalSSC(final SSC ssc, final Spectrum querySpectrum){
//        System.out.println("\n\nno more unsaturated atoms left? -> " + !ssc.hasUnsaturatedAtoms());
//        System.out.println("query spectrum size reached? -> " + ssc.getSubspectrum().getSignalCount() + " == " + querySpectrum.getSignalCount() + " -> " + (ssc.getSubspectrum().getSignalCount() == querySpectrum.getSignalCount()));
        return !ssc.hasUnsaturatedAtoms() && (ssc.getSubspectrum().getSignalCount() == querySpectrum.getSignalCount());
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
    
    
    private static HashMap<Integer, Integer> extendSSC(final SSC ssc1, final SSC ssc2, final HashMap<Integer, Integer> mappedAtomIndices, final double pickPrecision) throws CDKException, CloneNotSupportedException{
        int mappedAtomIndexSSC1, mappedAtomIndexSSC2, equivalentSignalIndex;
        IBond bondToAdd;
        IAtom atomToAdd, parentAtomSSC1;
        Signal signalToAddSSC2;
        ConnectionTree connectionTreeToAddSSC2;
        final HashMap<Integer, Integer> mappedAtomIndices2 = new HashMap<>(mappedAtomIndices);
        
        // for each atom pair map
        for (final Map.Entry<Integer, Integer> entry : mappedAtomIndices.entrySet()) {
            mappedAtomIndexSSC2 = entry.getKey();
            mappedAtomIndexSSC1 = entry.getValue();
            // check whether the current mapped atom in SSC1 (to extend) is an open site (unsaturated);
            // if no then skip because there is no chance to add atoms and bonds
            if ((ssc1.isUnsaturated(mappedAtomIndexSSC1) == null) || !ssc1.isUnsaturated(mappedAtomIndexSSC1)) {
                continue;
            }
            System.out.println("\nmapping: " + mappedAtomIndexSSC2 + ", " + mappedAtomIndexSSC1
                    + "\n-> " + ssc1.getHOSECode(mappedAtomIndexSSC1) + "\n-> " + ssc2.getHOSECode(mappedAtomIndexSSC2));

            // BFS to build connection tree which contains atoms in SSC2 to add to SSC1;
            // all mapped atom indices in SSC2 are used as list of visited atoms and then used 
            // during BFS for building the connection tree; 
            // that means that connected but unmapped atoms in SSC2 should exist in the resulting 
            // connection tree to add to SSC1
            connectionTreeToAddSSC2 = HOSECodeBuilder.buildConnectionTree(ssc2.getSubstructure(), mappedAtomIndexSSC2, null, new HashSet<>(mappedAtomIndices2.keySet()));
            System.out.println(" -> BFS: " + connectionTreeToAddSSC2.toString()
                    + "\n -> maxSphere: " + connectionTreeToAddSSC2.getMaxSphere());
            // traverse connection tree and try to add missing (not mapped) atoms from SSC2 to SSC1
            for (int s = 1; s <= connectionTreeToAddSSC2.getMaxSphere(); s++) {
                // traverse via spheres of connection tree
                for (final ConnectionTreeNode connectedNodeInSphereToAddSSC2 : connectionTreeToAddSSC2.getNodesInSphere(s)) {
                    // in case of ring closure node try to close a ring
                    if (connectedNodeInSphereToAddSSC2.isRingClosureNode()) {
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
//                            ssc1.getSubstructure().addBond(bondToAdd);
//                        }

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
        
        return mappedAtomIndices2;        
    }   
    
    static void checkBonds(final SSC ssc1, final SSC ssc2, final HashMap<Integer, Integer> mappedAtomIndices) throws CDKException, CloneNotSupportedException {
        int mappedAtomIndexSSC1, mappedAtomIndexSSC2;
        IAtom mappedAtomSSC1, mappedAtomSSC2, connectedAtomSSC1;
        IBond bondToAdd;

        for (final Map.Entry<Integer, Integer> entry : mappedAtomIndices.entrySet()) {
            mappedAtomIndexSSC2 = entry.getKey();
            mappedAtomIndexSSC1 = entry.getValue();            
            if ((ssc1.isUnsaturated(mappedAtomIndexSSC1) == null) || !ssc1.isUnsaturated(mappedAtomIndexSSC1)) {
                continue;
            }            
            System.out.println("\nmapping: " + mappedAtomIndexSSC2 + ", " + mappedAtomIndexSSC1
                    + "\n-> " + ssc1.getHOSECode(mappedAtomIndexSSC1) + "\n-> " + ssc2.getHOSECode(mappedAtomIndexSSC2));
            
            mappedAtomSSC2 = ssc2.getSubstructure().getAtom(mappedAtomIndexSSC2);
            mappedAtomSSC1 = ssc1.getSubstructure().getAtom(mappedAtomIndexSSC1);
            for (final IAtom connectedAtomSSC2 : ssc2.getSubstructure().getConnectedAtomsList(mappedAtomSSC2)) {
                if(mappedAtomIndices.get(connectedAtomSSC2.getIndex()) == mappedAtomIndexSSC1){
                    continue;
                }
                connectedAtomSSC1 = ssc1.getSubstructure().getAtom(mappedAtomIndices.get(connectedAtomSSC2.getIndex()));
                if(ssc1.getSubstructure().getBond(mappedAtomSSC1, connectedAtomSSC1) == null){
                    bondToAdd = new Bond(mappedAtomSSC1, connectedAtomSSC1, ssc2.getSubstructure().getBond(mappedAtomSSC2, connectedAtomSSC2).getOrder());
                    if(Assembly.isValidBondAddition(ssc1.getSubstructure(), mappedAtomIndexSSC1, bondToAdd)
                            && Assembly.isValidBondAddition(ssc1.getSubstructure(), connectedAtomSSC1.getIndex(), bondToAdd)){
                        ssc1.getSubstructure().addBond(bondToAdd);
                        ssc1.update();                        
                    }                    
                }
            }
        }
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
     * 
     * @deprecated 
     */
    public static HashMap<Integer, ArrayList<Integer[]>> getStructuralOverlaps(final SSC ssc1, final SSC ssc2, final double shiftTol, final int minSphereMatchCount, final String atomType) throws CloneNotSupportedException {
        final HashMap<Integer, ArrayList<Integer[]>> overlaps = new HashMap<>();
        int maxMatchingSphere;
        for (int i = 0; i < ssc1.getAtomCount(); i++) {
            if (!ssc1.getSubstructure().getAtom(i).getSymbol().equals(atomType)) {
                continue;
            }
            for (int j = 0; j < ssc2.getAtomCount(); j++) {
                if (!ssc2.getSubstructure().getAtom(j).getSymbol().equals(atomType)) {
                    continue;
                }
                maxMatchingSphere = -1;
                if (Match.isEqualNode(ssc1.getConnectionTree(i).getRootNode(), ssc2.getConnectionTree(j).getRootNode(), 
                        ssc1.getSubspectrum().getSignal(ssc1.getAssignments().getSignalIndex(0, i)), 
                        ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getSignalIndex(0, j)), shiftTol)) {
                    maxMatchingSphere = Match.getMaximumMatchingSphere(ssc1, ssc2, i, j, shiftTol);
                }
                if (maxMatchingSphere >= minSphereMatchCount) {
                    if (!overlaps.containsKey(maxMatchingSphere)) {
                        overlaps.put(maxMatchingSphere, new ArrayList<>());
                    }
                    overlaps.get(maxMatchingSphere).add(new Integer[]{i, j});
                }
            }
        }
        return overlaps;
    }

    /**
     *
     * @param ssc1
     * @param ssc2
     * @param overlaps
     * @param shiftTol
     * @return
     * @throws CDKException
     * @throws CloneNotSupportedException
     * 
     * @deprecated 
     */
    public static HashMap<Integer, Integer> mapStructuralOverlaps(final SSC ssc1, final SSC ssc2, final HashMap<Integer, ArrayList<Integer[]>> overlaps, final double shiftTol) throws CDKException, CloneNotSupportedException {
        final HashMap<Integer, Integer> mappedAtomIndices = new HashMap<>();
        HashMap<Integer, Integer> mappedEqualNodesKeys;
        int overlapAtomIndexSSC1;
        int overlapAtomIndexSSC2;
        int mappedAtomIndexSSC1;
        int mappedAtomIndexSSC2;
        boolean containsUnsaturatedAtomsSSC1;
        // for all found max. spheres
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
                    if (ssc1.isUnsaturated(nodeKey)) {
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
                        if (ssc1.getConnectionTree(overlapAtomIndexSSC1).getNode(mappedAtomIndexSSC1).isRingClosureNode()
                                || (mappedAtomIndexSSC1 == -1)
                                || ssc2.getConnectionTree(overlapAtomIndexSSC2).getNode(mappedAtomIndexSSC2).isRingClosureNode()) {
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
                    mappedEqualNodesKeys = Match.mapEqualNodesInSphere(ssc2, ssc1, overlapAtomIndexSSC2, overlapAtomIndexSSC1, s, shiftTol);
                    //                    System.out.println(" -> in s: " + s + " -> mapped equal node:  " + mappedEqualNodesKeys);
                    for (final Map.Entry<Integer, Integer> entry : mappedEqualNodesKeys.entrySet()) {
                        mappedAtomIndexSSC2 = entry.getKey();
                        mappedAtomIndexSSC1 = entry.getValue();
                        if (ssc2.getConnectionTree(overlapAtomIndexSSC2).getNode(mappedAtomIndexSSC2).isRingClosureNode()
                                || (mappedAtomIndexSSC1 == -1)
                                || ssc1.getConnectionTree(overlapAtomIndexSSC1).getNode(mappedAtomIndexSSC1).isRingClosureNode()) {
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
    
    
     public static HashMap<String, SSC> assembleHOSECode(final SSCLibrary rankedSSCLibrary, final int startSSCIndex, final int minMatchingSphereCount, 
            final Spectrum querySpectrum, final double shiftTol, final double thrsMatchFactor, final double pickPrecision) throws CloneNotSupportedException, CDKException, IOException {
        
//         int counter = 0;         
//         for (final long sscIndex : rankedSSCLibrary.getSSCIndices()) {
//             if (counter <= 50) {
//                 try {
//                    Utils.generatePicture(rankedSSCLibrary.getSSC(sscIndex).getSubstructure(),  "/Users/mwenk/Downloads/results/out_" + counter + ".png");
//                 } catch (Exception e) {
//                     System.out.println("could not depict for ranked SSC index " + sscIndex + ": " + e.getMessage());
//                 }
//                 counter++;
//             } else {
//                 break;
//             }
//         }
         
        final SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Absolute);
        String structureAsSMILES;
        final HashMap<String, SSC> solutions = new HashMap<>();
        SSC ssc2, ssc1Backup;
        // clone (!!!) the SSC1 contents only; don't use the object (reference) itself because of modifications
        SSC ssc1 = rankedSSCLibrary.getSSC(startSSCIndex).getClone();
        
        int mappedAtomIndexSSC2;
        Signal signalToRemoveSSC2;
        Spectrum subspectrumSSC2Clone, combinedSpectrum;
        HashMap<Integer, Integer> mappedAtomIndices;
        HashMap<Integer, HashMap<Integer, Integer>> mappedAtomIndicesInSpheres;
        HashSet<Integer> missingAtomIndicesSSC2;
        HashMap<Integer, ArrayList<Integer[]>> overlaps;
        
        // check whether the current SSC is already a final SSC
        if (Assembly.isFinalSSC(ssc1, querySpectrum)) {
            structureAsSMILES = smilesGenerator.create(ssc1.getSubstructure());
            solutions.put(structureAsSMILES, ssc1);
            System.out.println("--> new solution found!!! -> " + solutions.size() + " -> " + structureAsSMILES);
            System.out.println("-> atom count: " + ssc1.getAtomCount() + ", bond count: " + ssc1.getBondCount());
            System.out.println("-> query spectrum:\t" + querySpectrum.getShifts(0));
            System.out.println("-> equivalences:\t" + querySpectrum.getEquivalences());
            System.out.println("-> pred. spectrum:\t" + ssc1.getSubspectrum().getShifts(0));
            System.out.println("-> equivalences:\t" + ssc1.getSubspectrum().getEquivalences());
            return solutions;
        }
        
        
        for (final long ssc2Index : rankedSSCLibrary.getSSCIndices()) {
            if (ssc2Index == startSSCIndex) {
                continue;
            }
            ssc2 = rankedSSCLibrary.getSSC(ssc2Index);

            System.out.println("\n\n-------------------------------- " + startSSCIndex + ", " + ssc2Index + " --------------------------------");                       
                   
            ssc1Backup = ssc1.getClone();
            subspectrumSSC2Clone = ssc2.getSubspectrum().getClone();                                      
            // 1. check for partial structural identity (overlaps)
            // 1.1 map structural overlaps in all spheres separately
            mappedAtomIndicesInSpheres = Assembly.getStructuralOverlapsHOSECode(ssc2, ssc1, minMatchingSphereCount, pickPrecision);
            // 1.2 insertion of all valid mappings in different spheres into one map
            mappedAtomIndices = Assembly.mapStructuralOverlapsHOSECode(ssc2, ssc1, mappedAtomIndicesInSpheres);
            if (mappedAtomIndices.isEmpty()) {
                continue;
            }
            missingAtomIndicesSSC2 = Assembly.getMissingAtomIndices(ssc2, mappedAtomIndices.keySet());
            System.out.println(" -> mapped atom indices: \t" + mappedAtomIndices);
            System.out.println(" -> missing atom indices:\t" + missingAtomIndicesSSC2);                      
            // 1.3 remove signals from SSC2's subspectrum associated with mapping (overlapping) atoms; for subspectra combination
            for (final Map.Entry<Integer, Integer> entry : mappedAtomIndices.entrySet()) {
                mappedAtomIndexSSC2 = entry.getKey();
                if (ssc2.getSubstructure().getAtom(mappedAtomIndexSSC2).getSymbol().equals(ssc1.getSubspectrumAtomType())) {
                    signalToRemoveSSC2 = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getSignalIndex(0, mappedAtomIndexSSC2));
                    subspectrumSSC2Clone.removeSignal(signalToRemoveSSC2);
                }
            }

            // 2. combine subspectra of SSC to extend and matched SSC and validate it
            combinedSpectrum = Match.combineSpectra(ssc1.getSubspectrum(), subspectrumSSC2Clone, pickPrecision);
//            System.out.println("\nspectrum1:\t" + ssc1.getSubspectrum().getShifts(0));
//            System.out.println("spectrum2:\t" + ssc2.getSubspectrum().getShifts(0));
//            System.out.println("spectrum2a:\t" + subspectrumSSC2Clone.getShifts(0));
//            System.out.println("spectrum3:\t" + combinedSpectrum.getShifts(0));
//            System.out.println("spectrum3 equ:\t" + combinedSpectrum.getEquivalences());

            // if no valid spectrum could be built then go to next pairwise SSC comparison
            if ((combinedSpectrum == null) || !Assembly.isValidSubspectrum(combinedSpectrum, querySpectrum, shiftTol, thrsMatchFactor)) {
                System.out.println("-> no valid combined spectrum!");
                continue;
            }
            System.out.println("-> !!!valid combined spectrum!!!");

            // 3. substructure extension in SSC1
            // 3.1 try to add atoms and bonds if some are missing
            if(!missingAtomIndicesSSC2.isEmpty()){
                mappedAtomIndices = Assembly.extendSSC(ssc1, ssc2, mappedAtomIndices, pickPrecision);
                if (mappedAtomIndices == null) {
                    System.out.println("---> could not extend ssc1 -> set to backup SSC!");
                    ssc1 = ssc1Backup.getClone();
                    continue;
                }
                if (mappedAtomIndices.isEmpty()) {
                    ssc1 = ssc1Backup.getClone();
                    continue;
                }
                
                missingAtomIndicesSSC2 = Assembly.getMissingAtomIndices(ssc2, mappedAtomIndices.keySet());
                System.out.println(" -> mapped atoms indices  at the end:\t" + mappedAtomIndices);
                System.out.println(" -> missing atom indices  at the end:\t" + missingAtomIndicesSSC2);
                if (!missingAtomIndicesSSC2.isEmpty()) {
//                System.out.println("---> could not add all missing atoms!");
                    ssc1 = ssc1Backup.getClone();
                    continue;
                }
            } 
            // 3.2 check for potential missing bonds (i.e. to close rings)
            Assembly.checkBonds(ssc1, ssc2, mappedAtomIndices);
            
//            Utils.generatePicture(smilesParser.parseSmiles(smilesGenerator.create(ssc1.getSubstructure())), "outputs/temp_" + startSSCIndex + "_" + ssc2Index + ".png");
//            try {
//                Utils.generatePicture(ssc1.getSubstructure(), "results/temp_" + startSSCIndex + "_" + ssc2Index + ".png");
//            } catch (Exception e) {
//                
//            }
            
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
                }
//                break;
                ssc1 = ssc1Backup.getClone();
                continue;
            }
            
//            if(!ssc1.hasUnsaturatedAtoms()){
//                ssc1 = ssc1Backup.getClone();
//                continue;
//            }
             
        }
        
        return solutions;
    }
}
