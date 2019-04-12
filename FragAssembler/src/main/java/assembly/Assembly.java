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

import hose.HOSECodeBuilder;
import casekit.NMR.Utils;
import casekit.NMR.model.Assignment;
import casekit.NMR.model.Signal;
import casekit.NMR.model.Spectrum;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Queue;
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
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import start.Start;

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
    public static HashMap<Integer, ArrayList<Double[]>> getOverlapsHOSECode(final SSC ssc1, final SSC ssc2, final int minSphereMatchCount, final double shiftTol) throws CDKException {
        final HashMap<Integer, ArrayList<Double[]>> overlapsInSpheres = new HashMap<>();

        IAtom rootAtomSSC1, rootAtomSSC2;
        Signal signalRootAtomSSC1, signalRootAtomSSC2;        
        int maxMatchingSphere, overlappingAtomsCount;
        ConnectionTree connectionTreeSSC1;
        ConnectionTree connectionTreeSSC2;
        Signal signalSSC1, signalSSC2;
        ArrayList<Integer> nodeKeysInSphereToCheckSSC1, nodeKeysInSphereToCheckSSC2;
        HashMap<Integer, Integer> mappedEqualNodesInSphere;
        boolean stop, notValid;        
        for (int i = 0; i < ssc1.getAtomCount(); i++) {
            rootAtomSSC1 = ssc1.getSubstructure().getAtom(i);
            signalRootAtomSSC1 = ssc1.getSubspectrum().getSignal(ssc1.getAssignments().getSignalIndex(0, i));
            for (int j = 0; j < ssc2.getAtomCount(); j++) {                
                rootAtomSSC2 = ssc2.getSubstructure().getAtom(j);
                // check for same atom types
                if(!rootAtomSSC1.getSymbol().equals(rootAtomSSC2.getSymbol())){
                    continue;
                }                
                // get signals of both atoms to compare (if available) and check for same signal properties
                // to set a starting point with higher probability to be correct                
                signalRootAtomSSC2 = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getSignalIndex(0, j));                                  
                if ((signalRootAtomSSC1 != null) && (signalRootAtomSSC2 != null) 
                        &&  (!signalRootAtomSSC1.getMultiplicity().equals(signalRootAtomSSC2.getMultiplicity()) 
                            || Math.abs(signalRootAtomSSC1.getShift(0) - signalRootAtomSSC2.getShift(0)) > shiftTol
                            )
                        ) {
                    continue;                    
                }
                System.out.println("\n-> i: " + i + ", j: " + j);
                // check for pure structural identity via HOSE code for each further sphere
                maxMatchingSphere = Match.getMaximumMatchingSphereHOSECode(ssc1, ssc2, i, j);
                System.out.println("-> maxMatchingSphere: " + maxMatchingSphere);
                // skip non-matching atom pairs
                if (maxMatchingSphere == -1) {
                    continue;
                }
                // skip atom pairs which do not match at least in given minimum number of spheres
                if(maxMatchingSphere < minSphereMatchCount) {
//                        // but atoms at open sites are allowed because those sites could enable 
//                        // an extension if no other atom mappings exist
                        if(!ssc1.isUnsaturated(i) && !ssc2.isUnsaturated(j)){
                            continue;
                        }
                }                    
                connectionTreeSSC1 = HOSECodeBuilder.buildConnectionTree(ssc1.getSubstructure(), i, maxMatchingSphere);
                connectionTreeSSC2 = HOSECodeBuilder.buildConnectionTree(ssc2.getSubstructure(), j, maxMatchingSphere);                
                // if a atom pair matched succesfully via HOSE codes (structurally) then check now for spectral identity/similarity
                notValid = false;
                // check all spheres until max. matching sphere for spectral similarity/identity
                for (int s = 0; s <= maxMatchingSphere; s++) {
                    nodeKeysInSphereToCheckSSC1 = connectionTreeSSC1.getNodeKeysInSphere(s);
                    nodeKeysInSphereToCheckSSC2 = connectionTreeSSC2.getNodeKeysInSphere(s);
//                    System.out.println("in s: " + s + " -> node keys SSC1: " + nodeKeysInSphereToCheckSSC1);
//                    System.out.println("in s: " + s + " -> node keys SSC2: " + nodeKeysInSphereToCheckSSC2);                                    
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
                    // within a HOSE code is different than in the second HOSE code;
                    // then one could search for a possible atom mapping via spectral properties 
                    // within a sphere in different order than HOSE code order
                    if (stop) {
                        mappedEqualNodesInSphere = Match.mapEqualNodesInSphere(ssc1, ssc2, i, j, s, shiftTol);
                        System.out.println("---> in s: " + s + " (spectral matching) --> maps in sphere: " + mappedEqualNodesInSphere);
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
                    continue;
                }
                System.out.println("-> maxMatchingSphere " + maxMatchingSphere + " is valid!!!");
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
                        signalSSC1 = ssc1.getSubspectrum().getSignal(ssc1.getAssignments().getSignalIndex(0, nodesInSphereSSC1.get(k).getKey()));
                        signalSSC2 = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getSignalIndex(0, nodesInSphereSSC2.get(k).getKey()));
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
     * matching spheres as keys and the belonging mapped atom pairs and 
     * overlapping atoms count as values
     * @return 
     * @throws CDKException
     * 
     * @see #getOverlapsHOSECode(model.SSC, model.SSC, int, double) 
     */
    public static HashMap<Integer, Integer> mapOverlapsHOSECode(final SSC ssc1, final SSC ssc2, final HashMap<Integer, ArrayList<Double[]>> overlapsInSpheres) throws CDKException{
        // insertion of all valid mappings in different spheres into one map, 
        // starting with highest matching sphere (higher chance of correctness)
        if (overlapsInSpheres.isEmpty()) {
            return new HashMap<>();
        }        
        
        final HashMap<Integer, Integer> mappedAtomIndices = new HashMap<>();
        ArrayList<Integer> nodeKeysSSC1, nodeKeysSSC2;
        int rootMappedAtomIndexSSC1, rootMappedAtomIndexSSC2, mappedAtomIndexSSC1, mappedAtomIndexSSC2;
        ConnectionTree connectionTreeSSC1, connectionTreeSSC2;
        boolean containsUnsaturatedAtomsSSC2;
        ArrayList<Double[]> overlapsInSphere;
        // for all found spheres try to insert the mappings
        for (int m = Collections.max(overlapsInSpheres.keySet()); m >= Collections.min(overlapsInSpheres.keySet()); m--) {
            if (!overlapsInSpheres.containsKey(m)) {
                continue;
            }
            // rank overlaps in current max. sphere by (1) mean of deviations and (2) total number of overlapping atoms
            overlapsInSphere = overlapsInSpheres.get(m);
            overlapsInSphere.sort(new Comparator<Double[]>() {
                @Override
                public int compare(final Double[] entry1, final Double[] entry2) {
                    int deviationsComp = ((entry1[3] == null) || (entry2[3] == null)) ? 0 : Double.compare(entry1[3], entry2[3]);
                    if(deviationsComp != 0){
                        return deviationsComp;
                    }
                    return -1 * Integer.compare(entry1[2].intValue(), entry2[2].intValue());
                }
            });
            // for each atom pair as valid starting point insert its belonging mapped atom pairs until its max. matching sphere            
            for (final Double[] entry : overlapsInSphere) {
                rootMappedAtomIndexSSC1 = entry[0].intValue();
                rootMappedAtomIndexSSC2 = entry[1].intValue();
                connectionTreeSSC1 = HOSECodeBuilder.buildConnectionTree(ssc1.getSubstructure(), rootMappedAtomIndexSSC1, m);
                connectionTreeSSC2 = HOSECodeBuilder.buildConnectionTree(ssc2.getSubstructure(), rootMappedAtomIndexSSC2, m);
                nodeKeysSSC1 = new ArrayList<>(connectionTreeSSC1.getKeys(true));
                nodeKeysSSC2 = new ArrayList<>(connectionTreeSSC2.getKeys(true));
                for (int k = 0; k < nodeKeysSSC1.size(); k++) {
                    mappedAtomIndexSSC1 = nodeKeysSSC1.get(k);
                    mappedAtomIndexSSC2 = nodeKeysSSC2.get(k);
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
                    }
                }
            }
        }

        return mappedAtomIndices;
    }
    
    
    public static boolean isValidBondAddition(final IAtomContainer ac, final int atomIndex, final IBond bondToAdd){        

        System.out.println(atomIndex + " --> " + Utils.getBondOrderSum(ac, atomIndex, true) + " + " + Utils.getBondOrderAsNumeric(bondToAdd) + " = " + (Utils.getBondOrderSum(ac, atomIndex, true) + Utils.getBondOrderAsNumeric(bondToAdd)) + " <= " + ac.getAtom(atomIndex).getValency() + " ? -> " + ((Utils.getBondOrderSum(ac, atomIndex, true) + Utils.getBondOrderAsNumeric(bondToAdd)) <= ac.getAtom(atomIndex).getValency()));        
        
        return (Utils.getBondOrderSum(ac, atomIndex, true) + Utils.getBondOrderAsNumeric(bondToAdd)) <= ac.getAtom(atomIndex).getValency();
    }
    
    /**
     *
     * @param HOSECodeLookupTable
     * @param HOSECode
     * @return
     * 
     * @see casekit.NMR.Utils#getRMS(java.util.ArrayList) 
     * 
     * @deprecated 
     */
    public static Double predictShift(final HashMap<String, ArrayList<Double>> HOSECodeLookupTable, final String HOSECode){
        if(HOSECodeLookupTable.containsKey(HOSECode)){
            return Utils.getRMS(HOSECodeLookupTable.get(HOSECode));
        }
        
        return null;
    }
    
    /**
     * Specified for carbons (13C) only. Not generic at the moment because of 
     * usage of {@link casekit.NMR.Utils#getMultiplicityFromHydrogenCount(int)} 
     * with {@code hCount}.
     *
     * @param HOSECodeLookupTable
     * @param ac
     * @param atomIndex
     * @param maxSphere
     * @param nucleus
     * @param hCount
     * @return
     * @throws CDKException
     * 
     * @see #predictShift(java.util.HashMap, java.lang.String) 
     * 
     * @deprecated 
     */
    public static Signal predictSignal(final HashMap<String, ArrayList<Double>> HOSECodeLookupTable, final IAtomContainer ac, final int atomIndex, final Integer maxSphere, final String nucleus, final Integer hCount) throws CDKException{        
        if(!Utils.checkIndexInAtomContainer(ac, atomIndex) || (hCount == null)){
            return null;
        }        
        final Double predictedShift = Assembly.predictShift(HOSECodeLookupTable, HOSECodeBuilder.buildHOSECode(ac, atomIndex, maxSphere, false));
        if (predictedShift == null) {
            return null;
        }
        return new Signal(
                new String[]{nucleus},
                new Double[]{predictedShift},
                Utils.getMultiplicityFromHydrogenCount(hCount),
                null                
        );
    }
    
    /**
     * Specified for carbons (13C) only. Not generic at the moment because of 
     * {@link casekit.NMR.Utils#getMultiplicityFromHydrogenCount(int)}.
     *
     * @param HOSECodeLookupTable
     * @param ac
     * @param maxSphere
     * @param atomType
     * @param nucleus
     * @return
     * @throws org.openscience.cdk.exception.CDKException
     * 
     * @see #predictSignal(java.util.HashMap, org.openscience.cdk.interfaces.IAtomContainer, int, java.lang.Integer, java.lang.String, java.lang.Integer) 
     * 
     * @deprecated 
     */
    public static Spectrum predictSpectrum(final HashMap<String, ArrayList<Double>> HOSECodeLookupTable, final IAtomContainer ac, final Integer maxSphere, final String atomType, final String nucleus) throws CDKException{
        final Spectrum predictedSpectrum = new Spectrum(new String[]{nucleus});
        for (final IAtom atom: ac.atoms()) {    
            if(atom.getSymbol().equals(atomType)){
                predictedSpectrum.addSignal(Assembly.predictSignal(HOSECodeLookupTable, ac, atom.getIndex(), maxSphere, nucleus, atom.getImplicitHydrogenCount()));
            }
        }
        
        return predictedSpectrum;
    } 
    
    public static boolean isValidSubspectrum(final Spectrum subspectrum, final Spectrum querySpectrum, final double shiftTol, final double thrsMatchFactor){
        if((subspectrum == null) || (subspectrum.getSignalCount() > querySpectrum.getSignalCount())){
            return false;
        }
        final Assignment matchAssignments = Match.matchSpectra(subspectrum, querySpectrum, shiftTol);    
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
        if (Utils.roundDouble(Match.calculateMatchFactor(subspectrum, querySpectrum, shiftTol), Start.DECIMAL_PLACES) > thrsMatchFactor) {
//            System.out.println("-> match factor not allowed!!!");
            return false;
        }
        
        return true;
    }    
        
    
    public static HashMap<String, SSC> assemble(final long nStarts, final int nThreads, final SSCLibrary rankedSSCLibrary, final int minMatchingSphereCount, 
            final Spectrum querySpectrum, final double thrsMatchFactor, final double shiftTol) throws InterruptedException {
        
//        long counter = 0;
//        for (final long sscIndex : rankedSSCLibrary.getSSCIndices()) {
//            if (counter <= 200) {
//                try {
//                    Utils.generatePicture(rankedSSCLibrary.getSSC(sscIndex).getSubstructure(), "results/out_" + counter + ".png");
//                } catch (Exception e) {
//                    System.out.println("could not depict for ranked SSC index " + sscIndex + ": " + e.getMessage());
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
            callables.add((Callable<HashMap<String, SSC>>) () -> {                                
//                return Assembly.assembleBFS(rankedSSCLibrary, j, minMatchingSphereCount, querySpectrum, thrsMatchFactor, shiftTol);
                return Assembly.assemble(rankedSSCLibrary, j, minMatchingSphereCount, querySpectrum, thrsMatchFactor, shiftTol);
//              
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
    
//    /**
//     *
//     * @param rankedSSCLibrary
//     * @param startSSCIndex
//     * @param minMatchingSphereCount
//     * @param querySpectrum
//     * @param shiftTol
//     * @param thrsMatchFactor
//     * @return
//     * @throws CloneNotSupportedException
//     * @throws CDKException
//     * @throws IOException
//     * 
//     * @deprecated 
//     */
//    public static HashMap<String, SSC> assemble(final SSCLibrary rankedSSCLibrary, final int startSSCIndex, final int minMatchingSphereCount, 
//            final Spectrum querySpectrum, final double shiftTol, final double thrsMatchFactor) throws CloneNotSupportedException, CDKException, IOException {
//        
//        final SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Absolute);
////        final SmilesParser smilesParser = new SmilesParser(new SilentChemObjectBuilder());
//        String structureAsSMILES;
//        final HashMap<String, SSC> solutions = new HashMap<>();
//        SSC ssc2, ssc1Backup;
//        // clone (!!!) the SSC1 contents only; don't use the object (reference) itself because of modifications
//        SSC ssc1 = rankedSSCLibrary.getSSC(startSSCIndex).getClone();
//        
//        int maxMatchingSphere, mappedAtomIndexSSC1, mappedAtomIndexSSC2;
//        Signal signalToRemoveSSC2;
//        Spectrum subspectrumSSC2Clone, combinedSpectrum;
//        HashMap<Integer, Integer> mappedAtomIndices;
//        HashSet<Integer> missingAtomIndicesSSC2;
//        HashMap<Integer, ArrayList<Integer[]>> overlaps;
//        
//        // check whether the current SSC is already a final SSC
//        if (Assembly.isFinalSSC(ssc1, querySpectrum, shiftTol, thrsMatchFactor)) {
//            structureAsSMILES = smilesGenerator.create(ssc1.getSubstructure());
//            solutions.put(structureAsSMILES, ssc1);
//            System.out.println("--> new solution found!!! -> " + solutions.size() + " -> " + structureAsSMILES);
//            System.out.println("-> atom count: " + ssc1.getAtomCount() + ", bond count: " + ssc1.getBondCount());
//            System.out.println("-> query spectrum:\t" + querySpectrum.getShifts(0));
//            System.out.println("-> equivalences:\t" + querySpectrum.getEquivalences());
//            System.out.println("-> pred. spectrum:\t" + ssc1.getSubspectrum().getShifts(0));
//            System.out.println("-> equivalences:\t" + ssc1.getSubspectrum().getEquivalences());
//            return solutions;
//        }
//        
//        
//        for (final long ssc2Index : rankedSSCLibrary.getSSCIndices()) {
//            if (ssc2Index == startSSCIndex) {
//                continue;
//            }
//            ssc2 = rankedSSCLibrary.getSSC(ssc2Index);
//
//            System.out.println("\n\n-------------------------------- " + startSSCIndex + ", " + ssc2Index + " --------------------------------");
//
//            
//            
//            
//            
//            
////            System.out.println("\n\nMCSS:\n");
////            final UniversalIsomorphismTester universalIsomorphismTester = new UniversalIsomorphismTester();
////            List<IAtomContainer> structuralOverlaps = universalIsomorphismTester.getOverlaps(ssc1.getSubstructure(), ssc2.getSubstructure());
////            structuralOverlaps.sort(new Comparator<IAtomContainer>() {
////                @Override
////                public int compare(final IAtomContainer ac1, final IAtomContainer ac2) {
////                    return -1 * Integer.compare(ac1.getAtomCount(), ac2.getAtomCount());
////                }
////            });
////            System.out.println("-> " + startSSCIndex + ", " + ssc2Index + " -> " + structuralOverlaps.size());
////    
////            mappedAtomIndices = new HashMap<>();
////            for (int j = 0; j < structuralOverlaps.size(); j++){
////                System.out.println(" -> j: " + j);
////                Utils.generatePicture(ssc1.getSubstructure(), "/Users/mwenk/Downloads/outputs/overlap_" + startSSCIndex + "_" + ssc2Index + "_" + j + "_p1.png");
////                Utils.generatePicture(ssc2.getSubstructure(), "/Users/mwenk/Downloads/outputs/overlap_" + startSSCIndex + "_" + ssc2Index + "_" + j + "_p2.png");
////                Utils.generatePicture(structuralOverlaps.get(j), "/Users/mwenk/Downloads/outputs/overlap_" + startSSCIndex + "_" + ssc2Index + "_" + j + ".png");
////                List<RMap> atomMapsSSC1 = universalIsomorphismTester.getSubgraphAtomsMap(ssc1.getSubstructure(), structuralOverlaps.get(j));
////                List<RMap> atomMapsSSC2 = universalIsomorphismTester.getSubgraphAtomsMap(ssc2.getSubstructure(), structuralOverlaps.get(j));
////                for (final RMap atomMapSSC1 : atomMapsSSC1) {
////                    for (final RMap atomMapSSC2 : atomMapsSSC2) {
////                        if (atomMapSSC1.getId2() == atomMapSSC2.getId2()
////                                && !mappedAtomIndices.containsKey(atomMapSSC2.getId1())
////                                && !mappedAtomIndices.containsValue(atomMapSSC1.getId1())) {
////                            mappedAtomIndices.put(atomMapSSC2.getId1(), atomMapSSC1.getId1());
////                        }
////                    }
////
////                }
////                missingAtomIndicesSSC2 = Assembly.getMissingAtomIndices(ssc2, mappedAtomIndices.keySet());
////                System.out.println("----> mapped atom indices: " + mappedAtomIndices);
////                System.out.println("----> missing atom indices SSC2: " + missingAtomIndicesSSC2);
////            }
////            System.out.println("\nMCSS ends:\n\n");
////            if (true) {
////                continue;
////            }
//            
//            
//            
//            
//            
//            
//            // 1. check for partial structural identity (overlaps); via HOSE code or connection tree comparison                                        
//            overlaps = Assembly.getStructuralOverlaps(ssc2, ssc1, shiftTol, minMatchingSphereCount, ssc1.getSubspectrumAtomType());
//            // if there is no structural identity then skip that SSC pair comparison            
//            if (overlaps.isEmpty()) {
//                continue;
//            }            
//            ssc1Backup = ssc1.getClone();
//            subspectrumSSC2Clone = ssc2.getSubspectrum().getClone();                                      
//            
//            
//            maxMatchingSphere = Collections.max(overlaps.keySet());
//            System.out.println("-> maxMatchingSphere: " + maxMatchingSphere);
//            
//            // 1.1 map structural overlaps
//            mappedAtomIndices = mapStructuralOverlaps(ssc1, ssc2, overlaps, shiftTol);
//            if (mappedAtomIndices.isEmpty()) {
//                continue;
//            }
//
//            missingAtomIndicesSSC2 = Assembly.getMissingAtomIndices(ssc2, mappedAtomIndices.keySet());
//            System.out.println(" -> mapped atom indices: \t" + mappedAtomIndices);
//            System.out.println(" -> missing atom indices:\t" + missingAtomIndicesSSC2);
//            if (missingAtomIndicesSSC2.isEmpty()) {
//                continue;
//            }
//
//            // 1.2 remove signals from SSC2's subspectrum associated with mapping (overlapping) atoms; for subspectra combination
//            for (final Map.Entry<Integer, Integer> entry : mappedAtomIndices.entrySet()) {
//                mappedAtomIndexSSC2 = entry.getKey();
//                if (ssc2.getSubstructure().getAtom(mappedAtomIndexSSC2).getSymbol().equals(ssc1.getSubspectrumAtomType())) {
//                    signalToRemoveSSC2 = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getSignalIndex(0, mappedAtomIndexSSC2));
//                    subspectrumSSC2Clone.removeSignal(signalToRemoveSSC2);
//                }
//            }
//
//            // 2. combine subspectra of SSC to extend and matched SSC and validate it
//            combinedSpectrum = Match.combineSpectra(ssc1.getSubspectrum(), subspectrumSSC2Clone, shiftTol);
////            System.out.println("\nspectrum1:\t" + ssc1.getSubspectrum().getShifts(0));
////            System.out.println("spectrum2:\t" + ssc2.getSubspectrum().getShifts(0));
////            System.out.println("spectrum2a:\t" + subspectrumSSC2Clone.getShifts(0));
////            System.out.println("spectrum3:\t" + combinedSpectrum.getShifts(0));
////            System.out.println("spectrum3 equ:\t" + combinedSpectrum.getEquivalences());
//
//            // if no valid spectrum could be built then go to next pairwise SSC comparison
//            if ((combinedSpectrum == null) || !Assembly.isValidSubspectrum(combinedSpectrum, querySpectrum, shiftTol, thrsMatchFactor)) {
////                System.out.println("-> no valid combined spectrum!");
//                continue;
//            }
//            System.out.println("-> !!!valid combined spectrum!!!");
//
////            System.out.println("entering substructure extension!!!");
//            // 3. substructure extension in SSC1
//            // 3. a) try to add atoms and bonds because atoms are missing
//            if (!missingAtomIndicesSSC2.isEmpty()) {
//                mappedAtomIndices = Assembly.extendSSC(ssc1, ssc2, mappedAtomIndices, shiftTol);
//                if (mappedAtomIndices == null) {
//                    System.out.println("---> could not extend ssc1 -> set to backup SSC!");
//                    ssc1 = ssc1Backup.getClone();
//                    continue;
//                }
//                if (mappedAtomIndices.isEmpty()) {
//                    ssc1 = ssc1Backup.getClone();
//                    continue;
//                }
//
//                missingAtomIndicesSSC2 = Assembly.getMissingAtomIndices(ssc2, mappedAtomIndices.keySet());
//                System.out.println(" -> mapped atoms indices  at the end:\t" + mappedAtomIndices);
//                System.out.println(" -> missing atom indices  at the end:\t" + missingAtomIndicesSSC2);
//                if (!missingAtomIndicesSSC2.isEmpty()) {
////                System.out.println("---> could not add all missing atoms!");
//                    ssc1 = ssc1Backup.getClone();
//                    continue;
//                }
//            }
//            // 3. b) check for potential bond additions
//            Assembly.addMissingBonds(ssc1, ssc2, mappedAtomIndices);
//            
////            // update overlaps
////            overlaps = Assembly.getStructuralOverlaps(ssc2, ssc1, shiftTol, minMatchingSphereCount, maxMatchingSphereCount);
////            // map structural overlaps again
////            mappedAtomIndices = Assembly.mapStructuralOverlaps(ssc1, ssc2, overlaps, shiftTol);
//                        
//            missingAtomIndicesSSC2 = Assembly.getMissingAtomIndices(ssc2, mappedAtomIndices.keySet());
//            System.out.println(" -> mapped atoms indices  at the end:\t" + mappedAtomIndices);
//            System.out.println(" -> missing atom indices  at the end:\t" + missingAtomIndicesSSC2);
//            if (!missingAtomIndicesSSC2.isEmpty()) {
////                System.out.println("---> could not add all missing atoms!");
//                ssc1 = ssc1Backup.getClone();
//                continue;
//            }
//            
//            
////            Utils.generatePicture(smilesParser.parseSmiles(smilesGenerator.create(ssc1.getSubstructure())), "outputs/temp_" + startSSCIndex + "_" + ssc2Index + ".png");
////            Utils.generatePicture(ssc1.getSubstructure(), "results/temp_" + startSSCIndex + "_" + ssc2Index + ".png");
//            
//            if(Assembly.isFinalSSC(ssc1, querySpectrum, shiftTol, thrsMatchFactor)){
//                structureAsSMILES = smilesGenerator.create(ssc1.getSubstructure());
//                if (!solutions.containsKey(structureAsSMILES)) {
//                    solutions.put(structureAsSMILES, ssc1);
//                    System.out.println("--> new solution found!!! -> " + solutions.size() + " -> " + structureAsSMILES);
//                    System.out.println("-> atom count: " + ssc1.getAtomCount() + ", bond count: " + ssc1.getBondCount());
//                    System.out.println("-> query spectrum:\t" + querySpectrum.getShifts(0));
//                    System.out.println("-> equivalences:\t" + querySpectrum.getEquivalences());
//                    System.out.println("-> pred. spectrum:\t" + ssc1.getSubspectrum().getShifts(0));
//                    System.out.println("-> equivalences:\t" + ssc1.getSubspectrum().getEquivalences());
//                }
//                ssc1 = ssc1Backup.getClone();
//                continue;
//            }
//            
////            // in case there are two hetero atoms the same (added two times) and unsaturated, then one could remove the first atom and link it to the neighbor(s) atom(s) of the second one
////            // this step should not be necessary anymore if one use the molecula formula of the unknown compound to elucidate
////            if((ssc1.getSubspectrum().getSignalCount() == querySpectrum.getSignalCount())){
////                System.out.println("query spectrum covered but unsaturated atoms are left: " + ssc1.getUnsaturatedAtomIndices());
////                System.out.println(" -> unsaturated atoms: " + ssc1.getUnsaturatedAtomIndices());
////                if(ssc1.getUnsaturatedAtomIndices().size() == 2){
////                    int unsaturatedAtomIndex1 = ssc1.getUnsaturatedAtomIndices().get(0);
////                    int unsaturatedAtomIndex2 = ssc1.getUnsaturatedAtomIndices().get(1);                    
////                    ConnectionTreeNode node1 = ssc1.getConnectionTree(unsaturatedAtomIndex1).getRootNode();
////                    ConnectionTreeNode node2 = ssc1.getConnectionTree(unsaturatedAtomIndex2).getRootNode();
////                    Signal signal1 = ssc1.getSubspectrum().getSignal(ssc1.getAssignments().getSignalIndex(0, unsaturatedAtomIndex1));
////                    Signal signal2 = ssc1.getSubspectrum().getSignal(ssc1.getAssignments().getSignalIndex(0, unsaturatedAtomIndex2));
////                    if(Match.isEqualNode(node1, node2, signal1, signal2, shiftTol)){
////                        System.out.println(" -> link!!!");
////                        IAtom unsaturatedAtom1 = ssc1.getSubstructure().getAtom(unsaturatedAtomIndex1);
////                        IAtom unsaturatedAtom2 = ssc1.getSubstructure().getAtom(unsaturatedAtomIndex2);
////                        for (final IBond bond2 : ssc1.getSubstructure().getConnectedBondsList(unsaturatedAtom2)) {
////                            if(isValidBondAddition(ssc1.getSubstructure(), unsaturatedAtomIndex1, bond2)){
////                                bond2.setAtom(bond2.getOther(unsaturatedAtom2), 0);   
////                                bond2.setAtom(unsaturatedAtom1, 1);
////                            }
////                        }
////                        ssc1.getSubstructure().removeAtom(unsaturatedAtom2);
////                        ssc1.update();
////                        
////                        if (Assembly.isFinalSSC(ssc1, querySpectrum)) {
////                            structureAsSMILES = smilesGenerator.create(ssc1.getSubstructure());
////                            if (!solutions.containsKey(structureAsSMILES)) {
////                                solutions.put(structureAsSMILES, ssc1);
////                                System.out.println("--> new solution found!!! -> " + solutions.size() + " -> " + structureAsSMILES);
////                                System.out.println("-> atom count: " + ssc1.getAtomCount() + ", bond count: " + ssc1.getBondCount());
////                                System.out.println("-> query spectrum:\t" + querySpectrum.getShifts(0));
////                                System.out.println("-> equivalences:\t" + querySpectrum.getEquivalences());
////                                System.out.println("-> pred. spectrum:\t" + ssc1.getSubspectrum().getShifts(0));
////                                System.out.println("-> equivalences:\t" + ssc1.getSubspectrum().getEquivalences());
////                            }
////                            break;
////                        }
////                    }
////                }
////            }
//        }
//        
//        return solutions;
//    }
    
    public static boolean isFinalSSC(final SSC ssc, final Spectrum querySpectrum, final double shiftTol, final double thrsMatchFactor) throws CDKException{
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
        Assembly.correctAromaticBondOrders(ssc);
        System.out.println("could all added aromatic bonds be adjusted? -> " + Assembly.isValidAromaticBondOrdersSet(ssc) + "\n");
        if (!Assembly.isValidAromaticBondOrdersSet(ssc)) {
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
    
    
    private static HashMap<Integer, Integer> extendSSC(final SSC ssc1, final SSC ssc2, final HashMap<Integer, Integer> mappedAtomIndices, final double shiftTol) throws CDKException, CloneNotSupportedException{
        int mappedAtomIndexSSC1, mappedAtomIndexSSC2, equivalentSignalIndex;
        IBond bondToAdd;
        IAtom atomToAdd, parentAtomSSC1, atomInRingSSC1_1, atomInRingSSC1_2, atomSSC2;
        Signal signalToAddSSC2;
        ConnectionTree connectionTreeToAddSSC2;
        final HashMap<Integer, Integer> mappedAtomIndices2 = new HashMap<>(mappedAtomIndices);
        
        // for each atom pair map
        for (final Map.Entry<Integer, Integer> entry : mappedAtomIndices.entrySet()) {
            mappedAtomIndexSSC2 = entry.getKey();
            mappedAtomIndexSSC1 = entry.getValue();
            // check whether the current mapped atom in SSC1 (to extend) is an open site (unsaturated);
            // if no then skip because there is no chance to add atoms and bonds
            if (/*(ssc1.isUnsaturated(mappedAtomIndexSSC1) == null) || */!ssc1.isUnsaturated(mappedAtomIndexSSC1)) {
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
                        System.out.println("ringClosureNode at: " + connectedNodeInSphereToAddSSC2.getKey());
                        ConnectionTreeNode parentNodeSSC2 = connectedNodeInSphereToAddSSC2.getParentNodes().get(0);
                        ConnectionTreeNode ringClosureParentNodeSSC2 = connectedNodeInSphereToAddSSC2.getParentNodes().get(0).getParentNodes().get(1);
                        if (!mappedAtomIndices2.containsKey(parentNodeSSC2.getKey())
                                || !mappedAtomIndices2.containsKey(ringClosureParentNodeSSC2.getKey())) {
                            continue;
                        }
                        // try to add the missing bond in ring
                        bondToAdd = parentNodeSSC2.getBondsToParents().get(1).clone();
                        if (Assembly.isValidBondAddition(ssc1.getSubstructure(), mappedAtomIndices2.get(parentNodeSSC2.getKey()), bondToAdd)
                                && Assembly.isValidBondAddition(ssc1.getSubstructure(), mappedAtomIndices2.get(ringClosureParentNodeSSC2.getKey()), bondToAdd)) {
                            bondToAdd.setAtom(ssc1.getSubstructure().getAtom(mappedAtomIndices2.get(parentNodeSSC2.getKey())), 0);
                            bondToAdd.setAtom(ssc1.getSubstructure().getAtom(mappedAtomIndices2.get(ringClosureParentNodeSSC2.getKey())), 1);
                            if(bondToAdd.isAromatic()){
                                bondToAdd.setOrder(IBond.Order.SINGLE);
                            }
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
                        if(bondToAdd.isAromatic()){
                            bondToAdd.setOrder(IBond.Order.SINGLE);
                        }
                        ssc1.getSubstructure().addBond(bondToAdd); 
                        System.out.println("in extendSSC: new bond (" + bondToAdd.getOrder().numeric() + ", " + bondToAdd.isAromatic() + ") between " + parentAtomSSC1.getIndex() + ", " + (ssc1.getAtomCount() - 1) + " added!!!");
                        
                        // add belonging signal from SSC2 to SSC1 
                        if (atomToAdd.getSymbol().equals(ssc1.getSubspectrumAtomType())) {
                            signalToAddSSC2 = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getSignalIndex(0, connectedNodeInSphereToAddSSC2.getKey()));
                            if (signalToAddSSC2 == null) {
                                return null;
                            }
                            equivalentSignalIndex = ssc1.getSubspectrum().pickClosestSignal(signalToAddSSC2.getShift(0), 0, shiftTol);
                            ssc1.getSubspectrum().addSignal(signalToAddSSC2, equivalentSignalIndex);
                            ssc1.getAssignments().addAssignment(new int[]{ssc1.getAtomCount() - 1});
                        }

                        mappedAtomIndices2.put(connectedNodeInSphereToAddSSC2.getKey(), ssc1.getAtomCount() - 1);
                        System.out.println("-> added atom " + connectedNodeInSphereToAddSSC2.getKey() + " from SSC2 to SSC1");

                        // ring closure directly from connection tree information
                        atomSSC2 = ssc2.getSubstructure().getAtom(connectedNodeInSphereToAddSSC2.getKey());
                        for (final IAtom connectedAtomSSC2 : ssc2.getSubstructure().getConnectedAtomsList(atomSSC2)) {
                            if ((connectedAtomSSC2.getIndex() != connectedNodeInSphereToAddSSC2.getParentNodes().get(0).getKey())
                                    && mappedAtomIndices2.containsKey(connectedAtomSSC2.getIndex())
                                    && !connectionTreeToAddSSC2.containsKey(connectedAtomSSC2.getIndex())) {
                                System.out.println("connected but mapped atom found -> ring closure!!!"
                                        + "\n -> " + connectedNodeInSphereToAddSSC2.getKey() + ", " + connectedAtomSSC2.getIndex());
                                bondToAdd = ssc2.getSubstructure().getBond(atomSSC2, connectedAtomSSC2).clone();                                
                                if (Assembly.isValidBondAddition(ssc1.getSubstructure(), ssc1.getAtomCount() - 1, bondToAdd)
                                        && Assembly.isValidBondAddition(ssc1.getSubstructure(), mappedAtomIndices2.get(connectedAtomSSC2.getIndex()), bondToAdd)) {
                                    System.out.println("!!! ring closure !!!");
                                    atomInRingSSC1_1 = ssc1.getSubstructure().getAtom(ssc1.getAtomCount() - 1);
                                    atomInRingSSC1_2 = ssc1.getSubstructure().getAtom(mappedAtomIndices2.get(connectedAtomSSC2.getIndex()));
                                    
                                    bondToAdd.setAtom(atomInRingSSC1_1, 0);
                                    bondToAdd.setAtom(atomInRingSSC1_2, 1);
                                    if(bondToAdd.isAromatic()){                                       
                                        bondToAdd.setOrder(IBond.Order.SINGLE);
                                    }
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
    
    /**
     *
     * @param ssc1 SSC to extend 
     * @param ssc2
     * @param mappedAtomIndices indices mappings from {@code ssc2} (key) to 
     * {@code ssc1} (value)
     * @return
     * @throws CDKException
     * @throws CloneNotSupportedException
     */
    public static boolean addMissingBonds(final SSC ssc1, final SSC ssc2, final HashMap<Integer, Integer> mappedAtomIndices) throws CDKException, CloneNotSupportedException {
        int mappedAtomIndexSSC1, mappedAtomIndexSSC2;
        IAtom mappedAtomSSC1, mappedAtomSSC2, connectedAtomSSC1;
        IBond bondToAdd;
        boolean addedAnyBond = false;
        
        for (final Map.Entry<Integer, Integer> entry : mappedAtomIndices.entrySet()) {
            mappedAtomIndexSSC2 = entry.getKey();
            mappedAtomIndexSSC1 = entry.getValue();            
            if (/*(ssc1.isUnsaturated(mappedAtomIndexSSC1) == null) || */!ssc1.isUnsaturated(mappedAtomIndexSSC1)) {
                continue;
            }            
//            System.out.println("\nmapping: " + mappedAtomIndexSSC2 + ", " + mappedAtomIndexSSC1
//                    + "\n-> " + ssc1.getHOSECode(mappedAtomIndexSSC1) + "\n-> " + ssc2.getHOSECode(mappedAtomIndexSSC2));
            
            mappedAtomSSC2 = ssc2.getSubstructure().getAtom(mappedAtomIndexSSC2);
            mappedAtomSSC1 = ssc1.getSubstructure().getAtom(mappedAtomIndexSSC1);
            for (final IAtom connectedAtomSSC2 : ssc2.getSubstructure().getConnectedAtomsList(mappedAtomSSC2)) {
                if(!mappedAtomIndices.containsKey(connectedAtomSSC2.getIndex()) || mappedAtomIndices.get(connectedAtomSSC2.getIndex()) == mappedAtomIndexSSC1){
                    continue;
                }
                connectedAtomSSC1 = ssc1.getSubstructure().getAtom(mappedAtomIndices.get(connectedAtomSSC2.getIndex()));                
                if(ssc1.getSubstructure().getBond(mappedAtomSSC1, connectedAtomSSC1) == null){
                    bondToAdd = ssc2.getSubstructure().getBond(mappedAtomSSC2, connectedAtomSSC2).clone();
                    bondToAdd.setAtom(mappedAtomSSC1, 0);
                    bondToAdd.setAtom(connectedAtomSSC1, 1);
                    if(Assembly.isValidBondAddition(ssc1.getSubstructure(), mappedAtomIndexSSC1, bondToAdd)
                            && Assembly.isValidBondAddition(ssc1.getSubstructure(), connectedAtomSSC1.getIndex(), bondToAdd)){
                        if(bondToAdd.isAromatic()){
                            bondToAdd.setOrder(IBond.Order.SINGLE);
                        }
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
    
    public static void correctAromaticBondOrders(final SSC ssc) throws CDKException{
        IAtom atom1, atom2;
        int bondOrderSumAtom1, bondOrderSumAtom2;        
        for (final IBond bond : ssc.getSubstructure().bonds()) {
            if(!bond.isAromatic() || (bond.getOrder() != IBond.Order.SINGLE)){
                continue;
            }
            atom1 = bond.getAtom(0);            
            bondOrderSumAtom1 = 0;            
            for (final IBond bondAtom1 : ssc.getSubstructure().getConnectedBondsList(atom1)) {
                bondOrderSumAtom1 += bondAtom1.getOrder().numeric();
            }
            if (atom1.getImplicitHydrogenCount() != null) {
                bondOrderSumAtom1 += atom1.getImplicitHydrogenCount();
            }
            // atom is already saturated, so nothing to adjust and skip that atom
            if (bondOrderSumAtom1 >= atom1.getValency()) {
                continue;
            }
            atom2 = bond.getAtom(1);
            bondOrderSumAtom2 = 0;
            for (final IBond bondAtom2 : ssc.getSubstructure().getConnectedBondsList(atom2)) {
                bondOrderSumAtom2 += bondAtom2.getOrder().numeric();
            }
            if (atom2.getImplicitHydrogenCount() != null) {
                bondOrderSumAtom2 += atom2.getImplicitHydrogenCount();
            }
            // atom is already saturated, so nothing to adjust and skip that atom
            if (bondOrderSumAtom2 >= atom2.getValency()) {
                continue;
            }
            // adjust aromatic bond from single to double bond order if possible
            if(((bondOrderSumAtom1 - IBond.Order.SINGLE.numeric()) + IBond.Order.DOUBLE.numeric() <= atom1.getValency())
                    && ((bondOrderSumAtom2 - IBond.Order.SINGLE.numeric()) + IBond.Order.DOUBLE.numeric() <= atom2.getValency())
                    ){
                bond.setOrder(IBond.Order.DOUBLE);
            }
        }
        ssc.update();                
    }
    
    public static boolean isValidAromaticBondOrdersSet(final SSC ssc){
        // check whether each atom connected to an aromatic bond is saturated now, without considering aromatic bond orders by 1.5
        IAtom atom;
        int bondOrderSum;
        boolean isConnectedToAnyAromaticBond;
        for (int i = 0; i < ssc.getAtomCount(); i++) {
            bondOrderSum = 0;
            isConnectedToAnyAromaticBond = false;
            atom = ssc.getSubstructure().getAtom(i);
            for (final IBond bond : ssc.getSubstructure().getConnectedBondsList(atom)) {
                if (bond.isAromatic()) {
                    isConnectedToAnyAromaticBond = true;
                }
                bondOrderSum += bond.getOrder().numeric();
            }
            // atom is not connected to any aromatic bond, so it is not to check and skip that atom
            if (!isConnectedToAnyAromaticBond) {
                continue;
            }
            if (atom.getImplicitHydrogenCount() != null) {
                bondOrderSum += atom.getImplicitHydrogenCount();
            }
            // atom is still not saturated (without considering aromatic bond orders by 1.5), so return false
            if (bondOrderSum < atom.getValency()) {
                return false;
            }
        }

        return true;
    }

//    /**
//     * Returns pairwise atom indices of two structural overlapping atoms in both
//     * substructures, including the maximum matching sphere counts as keys.
//     *
//     * @param ssc1 first SSC
//     * @param ssc2 second SSC
//     * @param shiftTol
//     * @param minSphereMatchCount number of minimum matching spheres
//     * @param atomType
//     * @return pairs of atom indices; first in SSC1, second in SSC2
//     * @throws java.lang.CloneNotSupportedException
//     * 
//     * @deprecated 
//     */
//    public static HashMap<Integer, ArrayList<Integer[]>> getStructuralOverlaps(final SSC ssc1, final SSC ssc2, final double shiftTol, final int minSphereMatchCount, final String atomType) throws CloneNotSupportedException {
//        final HashMap<Integer, ArrayList<Integer[]>> overlaps = new HashMap<>();
//        int maxMatchingSphere;
//        for (int i = 0; i < ssc1.getAtomCount(); i++) {
//            if (!ssc1.getSubstructure().getAtom(i).getSymbol().equals(atomType)) {
//                continue;
//            }
//            for (int j = 0; j < ssc2.getAtomCount(); j++) {
//                if (!ssc2.getSubstructure().getAtom(j).getSymbol().equals(atomType)) {
//                    continue;
//                }
//                maxMatchingSphere = -1;
//                if (Match.isEqualNode(ssc1.getConnectionTree(i).getRootNode(), ssc2.getConnectionTree(j).getRootNode(), 
//                        ssc1.getSubspectrum().getSignal(ssc1.getAssignments().getSignalIndex(0, i)), 
//                        ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getSignalIndex(0, j)), shiftTol)) {
//                    maxMatchingSphere = Match.getMaximumMatchingSphere(ssc1, ssc2, i, j, shiftTol);
//                }
//                if (maxMatchingSphere >= minSphereMatchCount) {
//                    if (!overlaps.containsKey(maxMatchingSphere)) {
//                        overlaps.put(maxMatchingSphere, new ArrayList<>());
//                    }
//                    overlaps.get(maxMatchingSphere).add(new Integer[]{i, j});
//                }
//            }
//        }
//        return overlaps;
//    }

//    /**
//     *
//     * @param ssc1
//     * @param ssc2
//     * @param overlaps
//     * @param shiftTol
//     * @return
//     * @throws CDKException
//     * @throws CloneNotSupportedException
//     * 
//     * @deprecated 
//     */
//    public static HashMap<Integer, Integer> mapStructuralOverlaps(final SSC ssc1, final SSC ssc2, final HashMap<Integer, ArrayList<Integer[]>> overlaps, final double shiftTol) throws CDKException, CloneNotSupportedException {
//        final HashMap<Integer, Integer> mappedAtomIndices = new HashMap<>();
//        HashMap<Integer, Integer> mappedEqualNodesKeys;
//        int overlapAtomIndexSSC1;
//        int overlapAtomIndexSSC2;
//        int mappedAtomIndexSSC1;
//        int mappedAtomIndexSSC2;
//        boolean containsUnsaturatedAtomsSSC1;
//        // for all found max. spheres
//        for (int m = Collections.max(overlaps.keySet()); m >= Collections.min(overlaps.keySet()); m--) {
//            if (!overlaps.containsKey(m)) {
//                continue;
//            }
//            System.out.println("\n -> m: " + m);
//            for (final Integer[] overlap : overlaps.get(m)) {
//                overlapAtomIndexSSC2 = overlap[0];
//                overlapAtomIndexSSC1 = overlap[1];
//                System.out.println("\n--> overlap: " + overlapAtomIndexSSC2 + ", " + overlapAtomIndexSSC1);
//                System.out.println("HOSE code ssc1:\t" + HOSECodeBuilder.buildHOSECode(ssc1.getConnectionTree(overlapAtomIndexSSC1), false));
//                System.out.println("HOSE code ssc2:\t" + HOSECodeBuilder.buildHOSECode(ssc2.getConnectionTree(overlapAtomIndexSSC2), false));
//                containsUnsaturatedAtomsSSC1 = false;
//                for (final int nodeKey : ssc1.getConnectionTree(overlapAtomIndexSSC1).getKeysInOrder(true)) {
//                    if (ssc1.isUnsaturated(nodeKey)) {
//                        containsUnsaturatedAtomsSSC1 = true;
//                        break;
//                    }
//                }
//                if (!containsUnsaturatedAtomsSSC1) {
//                    System.out.println("\nnothing unsaturated -> skip this overlap");
//                    continue;
//                }
//                // add mappings of all atoms between SSC1 and SSC2 which have an identical structural overlap
//                for (int s = 0; s <= m; s++) {
//                    for (int i = 0; i < ssc2.getConnectionTree(overlapAtomIndexSSC2).getNodesCountInSphere(s); i++) {
//                        mappedAtomIndexSSC2 = ssc2.getConnectionTree(overlapAtomIndexSSC2).getNodeKeysInSphere(s).get(i);
//                        mappedAtomIndexSSC1 = ssc1.getConnectionTree(overlapAtomIndexSSC1).getNodeKeysInSphere(s).get(i);
//                        if (ssc1.getConnectionTree(overlapAtomIndexSSC1).getNode(mappedAtomIndexSSC1).isRingClosureNode()
//                                || (mappedAtomIndexSSC1 == -1)
//                                || ssc2.getConnectionTree(overlapAtomIndexSSC2).getNode(mappedAtomIndexSSC2).isRingClosureNode()) {
//                            continue;
//                        }
//                        if (mappedAtomIndices.containsKey(mappedAtomIndexSSC2) && (mappedAtomIndices.get(mappedAtomIndexSSC2) != -1)) {
//                            //                            if (mappedAtomIndices.get(mappedAtomIndexSSC2) != mappedAtomIndexSSC1) {
//                            //                                System.out.println(" in s: " + s + " -> !!!tried to set mappedAtomIndexSSC1 more than one time!!! -> "
//                            //                                        + mappedAtomIndexSSC2 + " : " + mappedAtomIndices.get(mappedAtomIndexSSC2) + " vs. " + mappedAtomIndexSSC1);
//                            //                            }
//                            continue;
//                        }
//                        if (mappedAtomIndices.containsValue(mappedAtomIndexSSC1)) {
//                            //                            System.out.println(" in s: " + s + " -> !!!tried to set mappedAtomIndexSSC2 more than one time!!! -> "
//                            //                                    + mappedAtomIndexSSC2 + " : " + mappedAtomIndexSSC1);
//                            continue;
//                        }
//                        //                        System.out.println(" in s: " + s + " -> new map: " + mappedAtomIndexSSC2 + ", " + mappedAtomIndexSSC1);
//                        mappedAtomIndices.put(mappedAtomIndexSSC2, mappedAtomIndexSSC1);
//                    }
//                }
//                // add mappings of all atoms between SSC1 and SSC2 which have to be assigned to each other (not structurally identical)
//                for (int s = m + 1; s <= Integer.min(ssc2.getConnectionTree(overlapAtomIndexSSC2).getMaxSphere(), ssc1.getConnectionTree(overlapAtomIndexSSC1).getMaxSphere()); s++) {
//                    mappedEqualNodesKeys = Match.mapEqualNodesInSphere(ssc2, ssc1, overlapAtomIndexSSC2, overlapAtomIndexSSC1, s, shiftTol);
//                    //                    System.out.println(" -> in s: " + s + " -> mapped equal node:  " + mappedEqualNodesKeys);
//                    for (final Map.Entry<Integer, Integer> entry : mappedEqualNodesKeys.entrySet()) {
//                        mappedAtomIndexSSC2 = entry.getKey();
//                        mappedAtomIndexSSC1 = entry.getValue();
//                        if (ssc2.getConnectionTree(overlapAtomIndexSSC2).getNode(mappedAtomIndexSSC2).isRingClosureNode()
//                                || (mappedAtomIndexSSC1 == -1)
//                                || ssc1.getConnectionTree(overlapAtomIndexSSC1).getNode(mappedAtomIndexSSC1).isRingClosureNode()) {
//                            continue;
//                        }                                               
//                        if (mappedAtomIndices.containsKey(mappedAtomIndexSSC2) && (mappedAtomIndices.get(mappedAtomIndexSSC2) != -1)) {
//                            //                            if (mappedAtomIndices.get(mappedAtomIndexSSC2) != mappedAtomIndexSSC1) {
//                            //                                System.out.println(" in s: " + s + " -> !!!tried to set mappedAtomIndexSSC1 more than one time!!! -> "
//                            //                                        + mappedAtomIndexSSC2 + " : " + mappedAtomIndices.get(mappedAtomIndexSSC2) + " vs. " + mappedAtomIndexSSC1);
//                            //                            }
//                            continue;
//                        }
//                        if (mappedAtomIndices.containsValue(mappedAtomIndexSSC1)) {
//                            //                            System.out.println(" in s: " + s + " -> !!!tried to set mappedAtomIndexSSC2 more than one time!!! -> "
//                            //                                    + mappedAtomIndexSSC2 + " : " + mappedAtomIndexSSC1);
//                            continue;
//                        }
//                        //                        System.out.println(" in s: " + s + " -> new map: " + mappedAtomIndexSSC2 + ", " + mappedAtomIndexSSC1);
//                        mappedAtomIndices.put(mappedAtomIndexSSC2, mappedAtomIndexSSC1);
//                    }
//                }
//            }
//            System.out.println(" -> mapped atom indices after m: " + m + ":\t" + mappedAtomIndices);
//            if (!mappedAtomIndices.isEmpty()) {
//                System.out.println(" -> not empty");
//                break;
//            }
//        }
//        return mappedAtomIndices;
//    }
    
    public static SSC assemblyCore(final SSC ssc1, final SSC ssc2, final Spectrum querySpectrum, final int minMatchingSphereCount, final double shiftTol, final double thrsMatchFactor) throws CDKException, CloneNotSupportedException{        
        int mappedAtomIndexSSC2;
        HashMap<Integer, Integer> mappedAtomIndices;
        HashSet<Integer> missingAtomIndicesSSC2;
        
        // 1. check for partial structural identity (overlaps)
        // 1.1 map structural overlaps in all spheres separately
        final HashMap<Integer, ArrayList<Double[]>> overlapsInSpheres = Assembly.getOverlapsHOSECode(ssc2, ssc1, minMatchingSphereCount, shiftTol);
        // 1.2 insertion of all valid mappings in different spheres into one map
        mappedAtomIndices = Assembly.mapOverlapsHOSECode(ssc2, ssc1, overlapsInSpheres);
        if (mappedAtomIndices.isEmpty()) {
//            continue;
            return null;
        }
        missingAtomIndicesSSC2 = Assembly.getMissingAtomIndices(ssc2, mappedAtomIndices.keySet());
        System.out.println(" -> mapped atom indices: \t" + mappedAtomIndices);
        System.out.println(" -> missing atom indices:\t" + missingAtomIndicesSSC2);
        // 1.3 for subspectra combination: add signals from SSC2's subspectrum associated with mapping (overlapping) atoms to a new subspectrum to combine with subspectrum of SSC1
        final LinkedHashSet<Integer> signalIndicesToIgnoreFromSubspectrumSSC2 = new LinkedHashSet<>();
        for (final Map.Entry<Integer, Integer> entry : mappedAtomIndices.entrySet()) {
            mappedAtomIndexSSC2 = entry.getKey();
            if (ssc2.getSubstructure().getAtom(mappedAtomIndexSSC2).getSymbol().equals(ssc1.getSubspectrumAtomType())) {
                signalIndicesToIgnoreFromSubspectrumSSC2.add(ssc2.getAssignments().getSignalIndex(0, mappedAtomIndexSSC2));
            }
        }
        final Spectrum subspectrumToAddSSC2 = new Spectrum(ssc2.getSubspectrum().getNuclei());
        for (int s = 0; s < ssc2.getSubspectrum().getSignalCount(); s++) {
            if (!signalIndicesToIgnoreFromSubspectrumSSC2.contains(s)) {
                subspectrumToAddSSC2.addSignal(ssc2.getSubspectrum().getSignal(s).getClone(), ssc2.getSubspectrum().getEquivalence(s));
            }
        }
        // 2. combine subspectra of SSC to extend and matched SSC and validate it
        final Spectrum combinedSpectrum = Match.combineSpectra(ssc1.getSubspectrum(), subspectrumToAddSSC2, Start.EQUIV_SIGNAL_THRS);
//            System.out.println("\nspectrum1:\t" + ssc1.getSubspectrum().getShifts(0));
//            System.out.println("spectrum2:\t" + ssc2.getSubspectrum().getShifts(0));
//            System.out.println("spectrum2a:\t" + subspectrumToAddSSC2.getShifts(0));
//            System.out.println("spectrum3:\t" + combinedSpectrum.getShifts(0));
//            System.out.println("spectrum3 equ:\t" + combinedSpectrum.getEquivalences());

        // if no valid spectrum could be built then go to next pairwise SSC comparison
        if ((combinedSpectrum == null) || !Assembly.isValidSubspectrum(combinedSpectrum, querySpectrum, shiftTol, thrsMatchFactor)) {
//                System.out.println("-> no valid combined spectrum!");
//            continue;
            return null;
        }
        System.out.println("-> !!!valid combined spectrum!!!");

        // 3. substructure extension in SSC1
        // 3.1 try to add atoms and bonds if some are missing
        boolean extended = false;
        if (!missingAtomIndicesSSC2.isEmpty()) {
            mappedAtomIndices = Assembly.extendSSC(ssc1, ssc2, mappedAtomIndices, shiftTol);
            if ((mappedAtomIndices == null) || mappedAtomIndices.isEmpty()) {
//                    System.out.println("---> could not extend ssc1 -> set to backup SSC!");
//                ssc1 = backupSSC1.getClone();
//                continue;
                return null;
            }
            missingAtomIndicesSSC2 = Assembly.getMissingAtomIndices(ssc2, mappedAtomIndices.keySet());
            System.out.println(" -> mapped atoms indices at the end:\t" + mappedAtomIndices);
            System.out.println(" -> missing atom indices at the end:\t" + missingAtomIndicesSSC2);
            if (!missingAtomIndicesSSC2.isEmpty()) {
//                System.out.println("---> could not add all missing atoms!");
//                ssc1 = backupSSC1.getClone();
//                continue;
                return null;
            }
            extended = true;
        }
        // 3.2 check for potential missing bonds (i.e. to close rings)
        if (!extended && !Assembly.addMissingBonds(ssc1, ssc2, mappedAtomIndices)) {
//                System.out.println("nothing added -> skip");
//            continue;
            return null;
        }
        // the whole (sub)structure is saturated but the query spectrum is not fully covered -> skip
        if (!ssc1.hasUnsaturatedAtoms() && (ssc1.getSubspectrum().getSignalCount() != querySpectrum.getSignalCount())) {
//                System.out.println("built a saturated but not valid structure -> skip");
//            ssc1 = backupSSC1.getClone();
//            continue;
            return null;
        }
        
        return ssc1;
    }
    
    public static HashMap<String, SSC> assembleBFS(final SSCLibrary rankedSSCLibrary, final long startSSCIndex, final int minMatchingSphereCount, 
            final Spectrum querySpectrum, final double thrsMatchFactor, final double shiftTol) throws CloneNotSupportedException, CDKException, IOException {                
         
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
    
    
    
    public static HashMap<String, SSC> assemble(final SSCLibrary rankedSSCLibrary, final long startSSCIndex, final int minMatchingSphereCount,
            final Spectrum querySpectrum, final double thrsMatchFactor, final double shiftTol) throws CloneNotSupportedException, CDKException, IOException {

        final SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Absolute);
        String structureAsSMILES;
        final HashMap<String, SSC> solutions = new HashMap<>();
        SSC intermediate, backupSSC, ssc2, startSSC;
        // notice: clone (!!!) the SSC contents only; don't use the object (reference) itself because of modifications
        startSSC = rankedSSCLibrary.getSSC(startSSCIndex).getClone();
        intermediate = startSSC.getClone();
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
//                continue;
            }
            
        }
        
        return solutions;
    }
}
