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
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import parallel.ParallelTasks;
import start.Start;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.*;
import java.util.Map.Entry;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Assembly {

    public static boolean isValidSubspectrum(final Spectrum subspectrum, final Spectrum querySpectrum, final double shiftTol, final double thrsMatchFactor){
        if((subspectrum == null) || (subspectrum.getSignalCount() > querySpectrum.getSignalCount())){
//            System.out.println("-> subspectrum == null or signal count subspectrum > signal count query spectrum!!!");
            return false;
        }
        final Assignment matchAssignments = Matcher.matchSpectra(subspectrum, querySpectrum, 0, 0, shiftTol);
        // filter for unset assignments
        if (!matchAssignments.isFullyAssigned(0)){
//            System.out.println("-> set assigments not allowed!!! -> " + matchAssignments.getAssignments(0));
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
        // @TODO add/enable filter for intensities
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

    public static HashMap<Integer, ArrayList<Integer[]>> getOverlaps(final SSC ssc1, final SSC ssc2, final int minMatchingSphereCount, final double shiftTol){
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
    
    public static ConcurrentHashMap<String, SSC> assemble(final long nStarts, final int nThreads, final SSCLibrary rankedSSCLibrary, final int minMatchingSphereCount,
            final Spectrum querySpectrum, final double thrsMatchFactor, final double shiftTol, final String pathToOutputsFolder, final long querySpectrumCounter) throws InterruptedException {

//        long counter = 0;
//        for (final long sscindex : rankedSSCLibrary.getSSCIndices()) {
//            if (counter <= 200) {
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


        final ConcurrentHashMap<String, SSC> solutions = new ConcurrentHashMap<>();
        final ArrayList<Callable<HashMap<String, SSC>>> callables = new ArrayList<>();
        // add all task to do
        for (int i = 0; i < nStarts; i++) {
            final int j = i;
            callables.add(() -> {
                return Assembly.assembleDFS(rankedSSCLibrary, j, minMatchingSphereCount, querySpectrum, thrsMatchFactor, shiftTol, pathToOutputsFolder, querySpectrumCounter);
//                return Assembly.assembleBFS(rankedSSCLibrary, j, minMatchingSphereCount, querySpectrum, thrsMatchFactor, shiftTol, pathToOutputsFolder, querySpectrumCounter);
//                return Assembly.assembleSeq(rankedSSCLibrary, j, minMatchingSphereCount, querySpectrum, thrsMatchFactor, shiftTol);
            });
        }

        ParallelTasks.processTasks(callables, tempHashMap -> solutions.putAll(tempHashMap), nThreads);
        
        return solutions;
    }

    public static HashMap<String, SSC> assembleBFS(final SSCLibrary rankedSSCLibrary, final long startSSCIndex, final int minMatchingSphereCount, 
            final Spectrum querySpectrum, final double thrsMatchFactor, final double shiftTol, final String pathToOutputsFolder, final long querySpectrumCounter) throws Exception {

        final BufferedWriter bw = new BufferedWriter(new FileWriter(pathToOutputsFolder + "/results_" + querySpectrumCounter + "_temp_" + startSSCIndex + ".smiles"));
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
                bw.close();
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

                newIntermediate = Assembly.assemblyCore(intermediate.getClone(), ssc2, querySpectrum, minMatchingSphereCount, shiftTol, thrsMatchFactor);
                if (newIntermediate == null) {
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

                        bw.append(structureAsSMILES);
                        bw.newLine();
                        bw.flush();

//                    return solutions;
                    }
                    continue;
                }
                
                newPath = new LinkedHashSet<>(path);
                newPath.add(i);
                intermediates.add(new Object[]{newIntermediate, newPath});
            }
            
            intermediates.poll();
        }
        bw.close();
        
        return solutions;
    }


    public static HashMap<String, SSC> assembleDFS(final SSCLibrary rankedSSCLibrary, final long startSSCIndex, final int minMatchingSphereCount,
                                                   final Spectrum querySpectrum, final double thrsMatchFactor, final double shiftTol, final String pathToOutputsFolder, final long querySpectrumCounter) throws Exception {

        final BufferedWriter bw = new BufferedWriter(new FileWriter(pathToOutputsFolder + "/results_" + querySpectrumCounter + "_temp_" + startSSCIndex + ".smiles"));
        final SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Absolute);
        String structureAsSMILES;
        final HashMap<String, SSC> solutions = new HashMap<>();
        SSC intermediate, newIntermediate, ssc2;

        // swap the start SSC with the first SSC
        final ArrayList<SSC> rankedSSCs = new ArrayList<>(rankedSSCLibrary.getSSCs());
        rankedSSCs.add(0, rankedSSCs.get((int) startSSCIndex));
        rankedSSCs.remove((int) startSSCIndex + 1);

        // notice: clone (!!!) the SSC contents only; don't use the object (reference) itself because of modifications
        intermediate = rankedSSCs.get(0).getClone();
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
                bw.close();
            }
            return solutions;
        }

        LinkedHashSet<Long> path = new LinkedHashSet<>(), newPath;
        path.add((long) 0);
        // create stack with initial state for DFS
        final Stack<Object[]> intermediates = new Stack<>();
        intermediates.push(new Object[]{intermediate, path});

        long j = 1;
        while (!intermediates.isEmpty()) {
            for (long i = j; i < rankedSSCs.size(); i++) {
                intermediate = ((SSC) intermediates.peek()[0]).getClone();
                path = new LinkedHashSet<>((LinkedHashSet<Long>) intermediates.peek()[1]);
                System.out.println("--> for path: " + path + "\nnext ssc index: " + i + "/" + (rankedSSCs.size() - 1));

                ssc2 = rankedSSCs.get((int) i);
                System.out.println("\n\n-------------------------------- " + startSSCIndex + "_" + path + ", " + i + " --------------------------------");

                newIntermediate = Assembly.assemblyCore(intermediate.getClone(), ssc2, querySpectrum, minMatchingSphereCount, shiftTol, thrsMatchFactor);
                if ((newIntermediate == null)){
                    continue;
                }

//                try {
//                    Utils.generatePicture(newIntermediate.getSubstructure(), "results/temp_" + startSSCIndex + "_" + path + "_" + i + ".png");
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

                        bw.append(structureAsSMILES);
                        bw.newLine();
                        bw.flush();

//                        return solutions;
                    }
                    continue;
                }
                // valid and extended SSC built but not final, so add it to stack
                newPath = new LinkedHashSet<>(path);
                newPath.add(i);
                intermediates.push(new Object[]{newIntermediate, newPath});
            }
            if(!intermediates.isEmpty()){
                j = Collections.max((LinkedHashSet<Long>) intermediates.peek()[1]) + 1;
                intermediates.pop();
            }
        }
        bw.close();

        return solutions;
    }


    public static SSC assemblyCore(final SSC ssc1, final SSC ssc2, final Spectrum querySpectrum, final int minMatchingSphereCount, final double shiftTol, final double thrsMatchFactor) throws Exception {

//        // 1. check for partial structural identity (overlaps) via MCSS
//        // atom mapping by maximum common (valid) subgraph
//        HashMap<Integer, Integer> atomMappings = Assembly.getAtomMappingsMCSS(ssc1, ssc2, Assembly.getMaximumCommonSubgraph(ssc1, ssc2, shiftTol));


        // @TODO usage of shift tolerance window from min-max information and solvent effect information while searching for overlaps?
        // @TODO use combinations of valid overlaps to achieve closures of rings or cyclic structure
        HashMap<Integer, ArrayList<Integer[]>> overlapsHOSECodeNew = Assembly.getOverlaps(ssc1, ssc2, minMatchingSphereCount, shiftTol);
        if(overlapsHOSECodeNew.isEmpty()){
            return null;
        }

        final ArrayList<Object[]> validSSCExtensions = new ArrayList<>();
        final ArrayList<Object[]> invalidSSCExtensions = new ArrayList<>();
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

                // @TODO check this constraint again but should make sense
                if(ssc1Extended.isUnsaturated(i) || ssc2.isUnsaturated(j)){
                    System.out.println(" atom " + i + " in SSC1 or atom " + j + " in SSC2 is unsaturated and not allowed as overlap root");
                    continue;
                }

                ConnectionTree completeConnectionTreeSSC2 = HOSECodeBuilder.buildConnectionTree(ssc2.getSubstructure(), j, null);
//                ConnectionTree maxSphereConnectionTreeSSC2 = HOSECodeBuilder.buildConnectionTree(ssc2.getSubstructure(), j, ssc2.getMaxSphere());
                ConnectionTree matchingConnectionTreeSSC1 = HOSECodeBuilder.buildConnectionTree(ssc1Extended.getSubstructure(), i, s);
                ConnectionTree matchingConnectionTreeSSC2 = HOSECodeBuilder.buildConnectionTree(ssc2.getSubstructure(), j, s);

                if(!Match.rearrangeInvalidNodePairsInSphere(ssc1Extended, ssc2, matchingConnectionTreeSSC1, matchingConnectionTreeSSC2, s, shiftTol)){
                    continue;
                }

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
                System.out.println("predicted spectrum orig.: " + ssc1Extended.getSubspectrum().getShifts(0));
                System.out.println(" -----> unsaturated atoms SSC1:" + unsaturatedAtomsSSC1);

                // check for each found unsaturated atom in SSC1 whether there is a valid extension possible from SSC2 via its children
                for (final int unsaturatedAtomKeySSC1 : unsaturatedAtomsSSC1){
                    int indexUntilMaxMatchingSphere = connectionTreeKeysSSC1.indexOf(unsaturatedAtomKeySSC1);
                    int nodeKeyInCompleteConnectionTreeSSC1 = connectionTreeKeysSSC1.get(indexUntilMaxMatchingSphere);
                    int nodeKeyInCompleteConnectionTreeSSC2 = connectionTreeKeysSSC2.get(indexUntilMaxMatchingSphere);
                    System.out.println(" ---> for unsaturated atom in SSC1: " + nodeKeyInCompleteConnectionTreeSSC1 + " -> " + nodeKeyInCompleteConnectionTreeSSC2);
                    
                    ArrayList<ConnectionTreeNode> childNodesToAppend = completeConnectionTreeSSC2.getNode(nodeKeyInCompleteConnectionTreeSSC2).getChildNodes();
                    for (final ConnectionTreeNode childNodeToAppend : childNodesToAppend){
                        SSC ssc1ExtendedBackup = ssc1Extended.getClone();
                        System.out.println(" --> to append from SSC2: " + childNodeToAppend.getKey());

                        // if ring closure node
                        if(childNodeToAppend.isRingClosureNode()){
                            continue;
                        }

                        IBond bondToAdd = completeConnectionTreeSSC2.getBond(nodeKeyInCompleteConnectionTreeSSC2, childNodeToAppend.getKey()).clone();
                        if(!Utils.isValidBondAddition(ssc1Extended.getSubstructure(), nodeKeyInCompleteConnectionTreeSSC1, bondToAdd)){
                          continue;
                        }

                        ConnectionTree connectionTreeToAddSSC2 = hose.Utils.getSubtree(completeConnectionTreeSSC2, childNodeToAppend.getKey());
                        System.out.println("-> subtree would be: " + connectionTreeToAddSSC2 + "\n" + connectionTreeToAddSSC2.getKeys(true));
                        int counter = 1;
                        for (final int substructureTreeNodeKeySSC2 : connectionTreeToAddSSC2.getKeys(true)) {
                            atomMappingsTemp.put((ssc1Extended.getAtomCount() - 1) + counter, substructureTreeNodeKeySSC2);
                            counter++;
                        }
                        System.out.println(" -> atom mappings temp: " + atomMappingsTemp);
                        System.out.println("is valid bond addition? -> true");

                        ssc1Extended = Assembly.extendSSC(ssc1Extended, ssc2, connectionTreeToAddSSC2, nodeKeyInCompleteConnectionTreeSSC1, bondToAdd);
                        System.out.println(" -> EXTENSION method done!!!");

                        // if one child node extends already existing or just invalid atoms then skip this child node and set to backup SSC
                        if(!Assembly.isValidSubspectrum(ssc1Extended.getSubspectrum(), querySpectrum, shiftTol, thrsMatchFactor)){
                            ssc1Extended = ssc1ExtendedBackup.getClone();
                            continue;
                        }
                        System.out.println("predicted spectrum temp.: " + ssc1Extended.getSubspectrum().getShifts(0));
                        System.out.println(" -> after EXTENSION -> valid subspectrum!!!");


                        // try to close rings directly from added (and mapped) subtree of SSC2 to mapped atoms in SSC1
                        ConnectionTree completeConnectionTreeSSC1 = HOSECodeBuilder.buildConnectionTree(ssc1Extended.getSubstructure(), i, null);
                        ConnectionTreeNode nodeSSC1, parentNodeSSC1;
                        // for each node in added subtree
                        for (final int nodeKeyInSubtreeSSC2 : connectionTreeToAddSSC2.getKeys(true)){
                            // if node is ring closure point
                            if(ConnectionTree.isAtRingClosure(completeConnectionTreeSSC2.getNode(nodeKeyInSubtreeSSC2))){
                                // then check for each (additional) parent node whether there is a bond missing
                                for (final ConnectionTreeNode parentNodeSSC2 : completeConnectionTreeSSC2.getNode(nodeKeyInSubtreeSSC2).getParentNodes()){
//                                            // skip current parent node (ancestor of root of subtree)
//                                            if(parentNodeSSC2.getKey() == childNodeToAppend.getKey()){
//                                                continue;
//                                            }
                                    nodeSSC1 = null;
                                    parentNodeSSC1 = null;
                                    if(atomMappingsTemp.containsValue(parentNodeSSC2.getKey())){
                                        for (final Entry<Integer, Integer> entry : atomMappingsTemp.entrySet()){
                                            if(entry.getValue() == parentNodeSSC2.getKey()){
                                                parentNodeSSC1 = completeConnectionTreeSSC1.getNode(entry.getKey());
                                            }
                                            if(entry.getValue() == nodeKeyInSubtreeSSC2){
                                                nodeSSC1 = completeConnectionTreeSSC1.getNode(entry.getKey());
                                            }
                                            if((nodeSSC1 != null) && (parentNodeSSC1 != null)){
                                                break;
                                            }
                                        }
                                    }
                                    if((nodeSSC1 != null) && (parentNodeSSC1 != null)
                                            && (ssc1Extended.getSubstructure().getBond(
                                            ssc1Extended.getSubstructure().getAtom(nodeSSC1.getKey()),
                                            ssc1Extended.getSubstructure().getAtom(parentNodeSSC1.getKey())) == null)){
                                        bondToAdd = ssc2.getSubstructure().getBond(ssc2.getSubstructure().getAtom(nodeKeyInSubtreeSSC2), ssc2.getSubstructure().getAtom(parentNodeSSC2.getKey())).clone();
                                        if(Utils.isValidBondAddition(ssc1Extended.getSubstructure(), nodeSSC1.getKey(), bondToAdd)
                                                && Utils.isValidBondAddition(ssc1Extended.getSubstructure(), parentNodeSSC1.getKey(), bondToAdd)){
                                            bondToAdd.setAtom(ssc1Extended.getSubstructure().getAtom(nodeSSC1.getKey()), 0);
                                            bondToAdd.setAtom(ssc1Extended.getSubstructure().getAtom(parentNodeSSC1.getKey()), 1);
                                            ssc1Extended.getSubstructure().addBond(bondToAdd);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                System.out.println("predicted spectrum modi.: " + ssc1Extended.getSubspectrum().getShifts(0));

                // for current HOSE code matching atom pair
                // if a valid substructure and also subspectrum could be assembled
                if(Assembly.isValidSubspectrum(ssc1Extended.getSubspectrum(), querySpectrum, shiftTol, thrsMatchFactor)){

                    System.out.println("\nvalid SSC built: ");
                    ssc1Extended.update();
                    if((ssc1Extended.getAtomCount() <= ssc1.getAtomCount())
                            && (ssc1Extended.getBondCount() <= ssc1.getBondCount())){
                        System.out.println(" --> but no extension happened");
                        continue;
                    }

                    System.out.println("substructure: " + ssc1Extended.getSubstructure().getAtomCount());
                    System.out.println("subspectrum : " + ssc1Extended.getSubspectrum().getShifts(0));
                    System.out.println("assignments : " + ssc1Extended.getAssignments().getAssignments(0) + "\n");
                    System.out.println("unsaturated atoms: " + ssc1Extended.getUnsaturatedAtomIndices());

                    validSSCExtensions.add(new Object[]{ssc1Extended.getClone(), new HashMap<>(atomMappingsTemp)});
                } else {
                    invalidSSCExtensions.add(new Object[]{ssc1Extended.getClone(), new HashMap<>(atomMappingsTemp)});
                }
            }
            System.out.println("\n");
        }


        final ArrayList<Integer[]> atomPairsToLinkSSC1 = new ArrayList<>();
        final ArrayList<Object[]> sscExtensions = new ArrayList<>();
        sscExtensions.addAll(validSSCExtensions);
        sscExtensions.addAll(invalidSSCExtensions);
        HashMap<Integer, Integer> atomMappingsExtendedSSC1, atomMappingsExtendedSSC2;
        for (int i = 0; i < sscExtensions.size(); i++) {
            atomMappingsExtendedSSC1 = (HashMap<Integer, Integer>) sscExtensions.get(i)[1];
            for (int j = i + 1; j < sscExtensions.size(); j++) {
                atomMappingsExtendedSSC2 = (HashMap<Integer, Integer>) sscExtensions.get(j)[1];
                for (final Entry<Integer, Integer> entryExtendedSSC1 : atomMappingsExtendedSSC1.entrySet()){
                    for (final Entry<Integer, Integer> entryExtendedSSC2 : atomMappingsExtendedSSC2.entrySet()){
                        if(entryExtendedSSC1.getValue().equals(entryExtendedSSC2.getValue())){
//                            System.out.println("\n\n at " + i + ", " + j + " -> " + entryExtendedSSC1.getValue());
                            if((Utils.checkIndexInAtomContainer(ssc1.getSubstructure(), entryExtendedSSC1.getKey())) && ssc1.isUnsaturated(entryExtendedSSC1.getKey())){
//                                System.out.println("  --> to link: " + entryExtendedSSC1.getKey());
                                atomPairsToLinkSSC1.add(new Integer[]{entryExtendedSSC1.getKey(), entryExtendedSSC1.getValue()});
                            }
                        }
                    }
                }
            }
        }
        System.out.println("\n\n ----> atoms link: ");
        for (final Integer[] atomPair : atomPairsToLinkSSC1){
            System.out.println(" -> " + Arrays.toString(atomPair));
        }
        System.out.println("\n\n");


        if(validSSCExtensions.isEmpty()){
            return null;
        }
        if (validSSCExtensions.size() == 1){
            return (SSC) validSSCExtensions.get(0)[0];
        }

        validSSCExtensions.sort((validExtendedSSC1Array, validExtendedSSC2Array) -> {
            final SSC validExtendedSSC1 = (SSC) validExtendedSSC1Array[0];
            final SSC validExtendedSSC2 = (SSC) validExtendedSSC2Array[0];
            final int atomCountComp = -1 * Integer.compare(validExtendedSSC1.getAtomCount(), validExtendedSSC2.getAtomCount());
            if(atomCountComp != 0){
                return atomCountComp;
            }
            final int bondCountComp = -1 * Integer.compare(validExtendedSSC1.getBondCount(), validExtendedSSC2.getBondCount());
            if(bondCountComp != 0){
                return bondCountComp;
            }

            return Double.compare(Matcher.calculateAverageDeviation(validExtendedSSC1.getSubspectrum(), querySpectrum, 0, 0, shiftTol), Matcher.calculateAverageDeviation(validExtendedSSC2.getSubspectrum(), querySpectrum, 0, 0, shiftTol));
        });

        for (int i = 0; i < validSSCExtensions.size(); i++){
            System.out.println(" -------> for valid extension " + i + " with size " + ((SSC) validSSCExtensions.get(i)[0]).getAtomCount() + " and " + ((SSC) validSSCExtensions.get(i)[0]).getBondCount() + " and " + Matcher.calculateAverageDeviation(((SSC) validSSCExtensions.get(i)[0]).getSubspectrum(), querySpectrum, 0, 0, shiftTol));
        }

        return (SSC) validSSCExtensions.get(0)[0];
    }

    public static SSC extendSSC(final SSC ssc1, final SSC ssc2, final ConnectionTree connectionTreeToAddSSC, final Integer parentAtomIndexToLinkSSC1, final IBond bondToLinkSSC2){

        final SSC ssc1Backup;
        try {
            ssc1Backup = ssc1.getClone();
        } catch (Exception e) {
            return ssc1;
        }
        ConnectionTreeNode nodeInSphere;
        // add root atom of connection tree to add to SSC1 and link it via a given bond to parent atom in SSC1
        if(!ssc1.addAtom(connectionTreeToAddSSC.getRootNode().getAtom(), bondToLinkSSC2, parentAtomIndexToLinkSSC1)){
            return ssc1Backup;
        }
        Signal signalToAddSSC2 = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getIndex(0, connectionTreeToAddSSC.getRootNode().getKey()));
        if(signalToAddSSC2 != null){
            if(!ssc1.addSignal(signalToAddSSC2, ssc1.getAtomCount() - 1)){
                return ssc1Backup;
            }
        }

        // for each sphere: add the atom which is stored as node to atom container and set bonds between parent nodes
        for (int s = 1; s <= connectionTreeToAddSSC.getMaxSphere(); s++) {
            // first add all atoms and its parents (previous sphere only) to structure
            for (int i = 0; i < connectionTreeToAddSSC.getNodesInSphere(s).size(); i++) {
                nodeInSphere = connectionTreeToAddSSC.getNodesInSphere(s).get(i);
                if(nodeInSphere.isRingClosureNode()){
                   continue;
                }
                for (final ConnectionTreeNode parentNode : nodeInSphere.getParentNodes()){
                    if(parentNode.getSphere() < nodeInSphere.getSphere()){
                        if(!ssc1.addAtom(nodeInSphere.getAtom(), connectionTreeToAddSSC.getBond(parentNode.getKey(), nodeInSphere.getKey()), ssc1.getSubstructure().indexOf(parentNode.getAtom()))){
                            return ssc1Backup;
                        }
                        signalToAddSSC2 = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getIndex(0, nodeInSphere.getKey()));
                        if(signalToAddSSC2 != null){
                            if(!ssc1.addSignal(signalToAddSSC2, ssc1.getAtomCount() - 1)){
                                return ssc1Backup;
                            }
                        }
                    }
                }
            }
        }
        ConnectionTreeNode parentNode;
        IBond bondToParent;
        for (int s = 1; s <= connectionTreeToAddSSC.getMaxSphere(); s++) {
            // and as second add the remaining bonds (ring closures) to structure
            for (int i = 0; i < connectionTreeToAddSSC.getNodesInSphere(s).size(); i++) {
                nodeInSphere = connectionTreeToAddSSC.getNodesInSphere(s).get(i);
                if(nodeInSphere.isRingClosureNode()){
                    continue;
                }
                for (int j = 0; j < nodeInSphere.getParentNodes().size(); j++) {
                    parentNode = nodeInSphere.getParentNodes().get(j);
                    if(ssc1.getSubstructure().getBond(ssc1.getSubstructure().getAtom(ssc1.getSubstructure().indexOf(nodeInSphere.getAtom())),
                            ssc1.getSubstructure().getAtom(ssc1.getSubstructure().indexOf(parentNode.getAtom()))) == null){
                        bondToParent = nodeInSphere.getBondsToParents().get(j);
                        ssc1.addBond(bondToParent, ssc1.getSubstructure().indexOf(nodeInSphere.getAtom()), ssc1.getSubstructure().indexOf(parentNode.getAtom()));
                    }
                }
            }
        }

        return ssc1;
    }

    /**
     * @param rankedSSCLibrary
     * @param startSSCIndex
     * @param minMatchingSphereCount
     * @param querySpectrum
     * @param thrsMatchFactor
     * @param shiftTol
     * @return
     * @throws Exception
     *
     * @deprecated
     */
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
}
