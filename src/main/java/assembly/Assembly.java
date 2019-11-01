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
import model.SSC;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import start.Start;
import utils.Compare;
import utils.ParallelTasks;

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

    public static boolean isValidExtension(final Spectrum subspectrum, final Spectrum querySpectrum, final double shiftTol, final double thrsMatchFactor, final IAtomContainer substructure, final IMolecularFormula molecularFormula, final IAtomContainer prevSubstructure){

        return ((substructure.getAtomCount() > prevSubstructure.getAtomCount())
                || (substructure.getBondCount() > prevSubstructure.getBondCount()))
                && Assembly.isValidSSC(subspectrum, querySpectrum, shiftTol, thrsMatchFactor, substructure, molecularFormula);
    }

    public static boolean isValidSSC(final Spectrum subspectrum, final Spectrum querySpectrum, final double shiftTol, final double thrsMatchFactor, final IAtomContainer substructure,final IMolecularFormula molecularFormula){

        return isValidSubspectrum(subspectrum, querySpectrum, shiftTol, thrsMatchFactor) && isValidSubstructure(substructure, molecularFormula);
    }

    public static boolean isFinalSSC(final SSC ssc, final Spectrum querySpectrum, final double shiftTol, final double thrsMatchFactor, final IMolecularFormula molecularFormula) {
//        System.out.println("\nno more unsaturated atoms left? -> " + !ssc.hasUnsaturatedAtoms() + " -> " + ssc.getUnsaturatedAtomIndices());
        if(ssc.hasUnsaturatedAtoms()){
            return false;
        }
//        System.out.println("query spectrum size reached? -> " + ssc.getSubspectrum().getSignalCount() + " == " + querySpectrum.getSignalCount() + " -> " + (ssc.getSubspectrum().getSignalCount() == querySpectrum.getSignalCount()));
        if((ssc.getSubspectrum().getSignalCount() != querySpectrum.getSignalCount())){
            return false;
        }
//        System.out.println("isValidSpectrum? -> " + Assembly.isValidSubspectrum(ssc.getSubspectrum(), querySpectrum, shiftTol, thrsMatchFactor));
//        if(!Assembly.isValidSubspectrum(ssc.getSubspectrum(), querySpectrum, shiftTol, thrsMatchFactor)){
        if(!Assembly.isValidSSC(ssc.getSubspectrum(), querySpectrum, shiftTol, thrsMatchFactor, ssc.getSubstructure(), molecularFormula)){
            return false;
        }
        if((molecularFormula != null) && (MolecularFormulaManipulator.getAtomCount(molecularFormula) != MolecularFormulaManipulator.getAtomCount(MolecularFormulaManipulator.getMolecularFormula(ssc.getSubstructure())))){
            return false;
        }
        try {
            Kekulization.kekulize(ssc.getSubstructure());
//            System.out.println("kekulization? -> true");
        } catch (CDKException e) {
//            System.out.println("kekulization? -> false");
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
                    maxMatchingSphere = Compare.getMaximumMatchingSphereHOSECode(ssc1, ssc2, i, j, shiftTol);
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
    
    public static ConcurrentHashMap<String, SSC> assemble(final long nStarts, final int nThreads, final ArrayList<SSC> rankedSSCList, final int minMatchingSphereCount,
            final Spectrum querySpectrum, final double thrsMatchFactor, final double shiftTol, final IMolecularFormula molecularFormula, final String pathToOutputsFolder, final long querySpectrumCounter) throws InterruptedException {

//        long counter = 0;
//        for (final SSC rankedSSC : rankedSSCList) {
//            if (counter <= 1) {
//                int rootAtomIndex = 0;
//                final ExtendedConnectionMatrixHOSECode extendedConnectionMatrixHOSECode = new ExtendedConnectionMatrixHOSECode(rankedSSC.getSubstructure(), rootAtomIndex, rankedSSC.getMaxSphere());
//                System.out.println("\nsize -> " + extendedConnectionMatrixHOSECode.getAtomCount());
//                System.out.println(extendedConnectionMatrixHOSECode.toString());
//                int atomIndex1 = 7, atomIndex2 = 9;
//                System.out.println(" -----> in " + rootAtomIndex + ": " + atomIndex1 + " and " + atomIndex2 + " form a ring closure: " + extendedConnectionMatrixHOSECode.formRingClosure(rootAtomIndex, atomIndex1, atomIndex2));
//                atomIndex1 = 6; atomIndex2 = 12;
//                System.out.println(" -----> in " + rootAtomIndex + ": " + atomIndex1 + " and " + atomIndex2 + " form a ring closure: " + extendedConnectionMatrixHOSECode.formRingClosure(rootAtomIndex, atomIndex1, atomIndex2));
//                atomIndex1 = 8; atomIndex2 = 9;
//                System.out.println(" -----> in " + rootAtomIndex + ": " + atomIndex1 + " and " + atomIndex2 + " form a ring closure: " + extendedConnectionMatrixHOSECode.formRingClosure(rootAtomIndex, atomIndex1, atomIndex2));
//
//                extendedConnectionMatrixHOSECode.addAtom("Cl", 0,2,0, false, false, IAtomType.Hybridization.SP3);
//                boolean addedBond = extendedConnectionMatrixHOSECode.addBond(8, extendedConnectionMatrixHOSECode.getAtomCount() - 1, 1.0, false, false);
//                if(addedBond){
//                    System.out.println("!!!bond addition successful!!!");
//                } else {
//                    System.out.println("!!!bond addition NOT successful!!!");
//                }
//                addedBond = extendedConnectionMatrixHOSECode.addBond(8, extendedConnectionMatrixHOSECode.getAtomCount() - 1, 1.0, false, false);
//                if(addedBond){
//                    System.out.println("!!!bond addition successful!!!");
//                } else {
//                    System.out.println("!!!bond addition NOT successful!!!");
//                }
//                System.out.println("\nsize -> " + extendedConnectionMatrixHOSECode.getAtomCount());
//                System.out.println(extendedConnectionMatrixHOSECode.toString());
//            }
//            counter++;
//        }

//        int counter = 0;
//        for (final SSC rankedSSC : rankedSSCList) {
//            if (counter <= 200) {
//                try {
//                    Utils.generatePicture(rankedSSC.getSubstructure(), "results/out_" + counter + ".png");
//                } catch (Exception e) {
//                    System.out.println("could not depict for ranked ssc index " + rankedSSC.getIndex() + ": " + e.getMessage());
//                }
//                counter++;
//            } else {
//                break;
//            }
//        }


//        final SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
//        IAtomContainer ac;
//        String HOSECode, expHOSECode;
//        final HOSECodeGenerator hcg = new HOSECodeGenerator();
//        System.out.println("\n");
//        try {
//            // Bremser paper first example with two rings
//            ac = parser.parseSmiles("CC1=C(Cl)C=C(C=C1)C(=O)C1CCCCC1N");
//            Utils.setAromaticity(ac);
//            HOSECode = HOSECodeBuilder.buildHOSECode(ac, 5, 4, true);
//            expHOSECode = "C-3;*C*CC(*C,*C,=OC/*CX,*&,,CC/*&C,,CN,C)";
//
//            System.out.println("expected: " + expHOSECode);
//            System.out.println("build   : " + HOSECode);
//            System.out.println("build 2 : " + hcg.getHOSECode(ac, ac.getAtom(5),4));
//
//            try {
//                Utils.generatePicture(HOSECodeBuilder.buildAtomContainer(HOSECode, true), "results/out_HOSECodeBackwards" + 1 + ".png");
//            } catch (Exception e) {
//            }
//        } catch (CDKException e) {
//            e.printStackTrace();
//        }
//        System.out.println("");
//        try {
//            // Bremser paper second example with sulfur
//            // example where CDK HOSECodeGenerator fails
//            ac  = SilentChemObjectBuilder.getInstance().newAtomContainer();
//            ac.addAtom(new Atom(16));
//            ac.addAtom(new Atom(6, 1));
//            ac.addBond(0, 1, IBond.Order.SINGLE);
//            ac.addAtom(new Atom(6, 1));
//            ac.addBond(0, 2, IBond.Order.SINGLE);
//            ac.addAtom(new Atom(6, 0));
//            ac.addAtom(new Atom(6, 1));
//            ac.addBond(1, 3, IBond.Order.DOUBLE);
//            ac.addBond(2, 4, IBond.Order.DOUBLE);
//            ac.addBond(3, 4, IBond.Order.SINGLE);
//            ac.addAtom(new Atom(9));
//            ac.addBond(5, 3, IBond.Order.SINGLE);
//            HOSECode = HOSECodeBuilder.buildHOSECode(ac, 0, 4, true);
//            expHOSECode = "S-2;CC(=C,=C/F&,&/)";
//
//            System.out.println("expected: " + expHOSECode);
//            System.out.println("build   : " + HOSECode);
//            System.out.println("build 2 : " + hcg.getHOSECode(ac, ac.getAtom(0),4));
//
//            try {
//                Utils.generatePicture(HOSECodeBuilder.buildAtomContainer(HOSECode, true), "results/out_HOSECodeBackwards" + 2 + ".png");
//            } catch (Exception e) {
//            }
//        } catch (CDKException e) {
//            e.printStackTrace();
//        }
//        System.out.println("");
//        try {
//            // SpecInfo example 2:
//            ac  = parser.parseSmiles("C\\C(\\C=C/C=C)=C/C=C");
//            HOSECode = HOSECodeBuilder.buildHOSECode(ac, 1, 4, true);
//            expHOSECode = "C-3;=CCC(C,=C,/=C,C/,=C)"; // with error in article; was written as "C-3;=CCC(C,=C,/=C,C/C,=C)"
//
//            System.out.println("expected: " + expHOSECode);
//            System.out.println("build   : " + HOSECode);
//            System.out.println("build 2 : " + hcg.getHOSECode(ac, ac.getAtom(1),4));
//
//            try {
//                Utils.generatePicture(HOSECodeBuilder.buildAtomContainer(HOSECode, true), "results/out_HOSECodeBackwards" + 3 + ".png");
//            } catch (Exception e) {
//            }
//        } catch (CDKException e) {
//            e.printStackTrace();
//        }
//        System.out.println("");
//        try {
//            // example where CDK HOSECodeGenerator fails
//            ac  = parser.parseSmiles("CC12C[N]3(CON3)(C1)C2");
//            Utils.setAromaticity(ac);
//            HOSECode = HOSECodeBuilder.buildHOSECode(ac, 0, 6, true);
//            expHOSECode = "C-4;C(CCC/N,&,&/CN&&/O,&/&)";
//
//            System.out.println("expected: " + expHOSECode);
//            System.out.println("build   : " + HOSECode);
//            System.out.println("build 2 : " + hcg.getHOSECode(ac, ac.getAtom(0),6));
//
//            try {
//                Utils.generatePicture(HOSECodeBuilder.buildAtomContainer(HOSECode, true), "results/out_HOSECodeBackwards" + 4 + ".png");
//            } catch (Exception e) {
//            }
//        } catch (CDKException e) {
//            e.printStackTrace();
//        }
//        System.out.println("");
//        try {
//            // example where CDK HOSECodeGenerator fails
//            ac  = parser.parseSmiles("C1CCC2(CC1)CCCCC2");
//            HOSECode = HOSECodeBuilder.buildHOSECode(ac, 3, 4, false);
//            expHOSECode = "C-4;CCCC(C,C,C,C/C,C,&,&/&,&)";
//
//            System.out.println("expected: " + expHOSECode);
//            System.out.println("build   : " + HOSECode);
//            System.out.println("build 2 : " + hcg.getHOSECode(ac, ac.getAtom(3),4));
//
//            try {
//                Utils.generatePicture(HOSECodeBuilder.buildAtomContainer(HOSECode, false), "results/out_HOSECodeBackwards" + 5 + ".png");
//            } catch (Exception e) {
//            }
//        } catch (CDKException e) {
//            e.printStackTrace();
//        }
//        System.out.println("");
//        try {
//            // example where CDK HOSECodeGenerator fails
//            ac  = parser.parseSmiles("C1CCC2(C1)CCCCC2");
//            HOSECode = HOSECodeBuilder.buildHOSECode(ac, 3, 4, false);
//            expHOSECode = "C-4;CCCC(C,C,C,C/C,&,&,&/&)";
//
//            System.out.println("expected: " + expHOSECode);
//            System.out.println("build   : " + HOSECode);
//            System.out.println("build 2 : " + hcg.getHOSECode(ac, ac.getAtom(3),4));
//
//            try {
//                Utils.generatePicture(HOSECodeBuilder.buildAtomContainer(HOSECode, false), "results/out_HOSECodeBackwards" + 6 + ".png");
//            } catch (Exception e) {
//            }
//        } catch (CDKException e) {
//            e.printStackTrace();
//        }
//        System.out.println("\n");


        final ConcurrentHashMap<String, SSC> solutions = new ConcurrentHashMap<>();
        final ConcurrentHashMap<String, HashSet<Long>> invalidExtensions = new ConcurrentHashMap<>();
        final ArrayList<Callable<HashMap<String, SSC>>> callables = new ArrayList<>();
        // add all task to do
        for (int i = 0; i < nStarts; i++) {
            final int j = i;
            callables.add(() -> {
                return Assembly.assembleDFS(rankedSSCList, j, minMatchingSphereCount, querySpectrum, thrsMatchFactor, shiftTol, molecularFormula, invalidExtensions, pathToOutputsFolder, querySpectrumCounter);
//                return Assembly.assembleBFS(rankedSSCLibrary, j, minMatchingSphereCount, querySpectrum, thrsMatchFactor, shiftTol, pathToOutputsFolder, querySpectrumCounter);
//                return Assembly.assembleSeq(rankedSSCLibrary, j, minMatchingSphereCount, querySpectrum, thrsMatchFactor, shiftTol);
            });
        }
        ParallelTasks.processTasks(callables, tempHashMap -> solutions.putAll(tempHashMap), nThreads);
        
        return solutions;
    }

    /**
     * @param rankedSSCList
     * @param startSSCIndex
     * @param minMatchingSphereCount
     * @param querySpectrum
     * @param thrsMatchFactor
     * @param shiftTol
     * @param molecularFormula
     * @param pathToOutputsFolder
     * @param querySpectrumCounter
     * @return
     * @throws Exception
     *
     * @deprecated
     */
    public static HashMap<String, SSC> assembleBFS(final ArrayList<SSC> rankedSSCList, final long startSSCIndex, final int minMatchingSphereCount,
            final Spectrum querySpectrum, final double thrsMatchFactor, final double shiftTol, final IMolecularFormula molecularFormula, final String pathToOutputsFolder, final long querySpectrumCounter) throws Exception {

        final BufferedWriter bw = new BufferedWriter(new FileWriter(pathToOutputsFolder + "/results_" + querySpectrumCounter + "_temp_" + startSSCIndex + ".smiles"));
        final SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Absolute);
        String structureAsSMILES;
        final HashMap<String, SSC> solutions = new HashMap<>();
        SSC intermediate, newIntermediate, ssc2;
        // notice: clone (!!!) the SSC contents only; don't use the object (reference) itself because of modifications
        intermediate = rankedSSCList.get((int) startSSCIndex).getClone(false);
        // check whether the current SSC is already a final SSC
        if (Assembly.isFinalSSC(intermediate, querySpectrum, shiftTol, thrsMatchFactor, molecularFormula)) {
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
            intermediate = ((SSC) intermediates.peek()[0]).getClone(false);
            path = new LinkedHashSet<>((LinkedHashSet<Long>) intermediates.peek()[1]);
            System.out.println("--> for path: " + path + "\nnext ssc index: " + (Collections.max(path) + 1));
//            System.out.println("--> ranked SSC indices: " + rankedSSCLibrary.getSSCIndices());
            for (long i = Collections.max(path) + 1; i < rankedSSCList.size(); i++) {

                ssc2 = rankedSSCList.get((int) i);

                System.out.println("\n\n-------------------------------- " + path + ", " + i + " --------------------------------");                                

                newIntermediate = null;//Assembly.assemblyCore(intermediate.getClone(false), ssc2, querySpectrum, minMatchingSphereCount, shiftTol, thrsMatchFactor);
                if (newIntermediate == null) {
                    continue;
                }

//                try {
//                    Utils.generatePicture(newIntermediate.getSubstructure(), "results/temp_" + path + "_" + i + ".png");
//                } catch (Exception e) {
//
//                }

                if (Assembly.isFinalSSC(newIntermediate, querySpectrum, shiftTol, thrsMatchFactor, molecularFormula)) {
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


    public static HashMap<String, SSC> assembleDFS(final ArrayList<SSC> rankedSSCList, final long startSSCIndex, final int minMatchingSphereCount,
                                                   final Spectrum querySpectrum, final double thrsMatchFactor, final double shiftTol, final IMolecularFormula molecularFormula, final ConcurrentHashMap<String, HashSet<Long>> invalidExtensions, final String pathToOutputsFolder, final long querySpectrumCounter) throws Exception {

        final BufferedWriter bw = new BufferedWriter(new FileWriter(pathToOutputsFolder + "/results_" + querySpectrumCounter + "_temp_" + startSSCIndex + ".smiles"));
        final SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Absolute);
        String structureAsSMILES;
        final HashMap<String, SSC> solutions = new HashMap<>();
        SSC intermediate, newIntermediate;
        ArrayList<SSC> newIntermediates;

        // notice: clone (!!!) the SSC contents only; don't use the object (reference) itself because of modifications
        intermediate = rankedSSCList.get((int) startSSCIndex).getClone(false);
        // check whether the current SSC is already a final SSC
        if (Assembly.isFinalSSC(intermediate, querySpectrum, shiftTol, thrsMatchFactor, molecularFormula)) {
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
        // create stack with initial state for DFS
        final Stack<Object[]> intermediates = new Stack<>();
        LinkedHashSet<Float> path = new LinkedHashSet<>(), newPath;
        path.add((float) 0);
        intermediates.push(new Object[]{path, intermediate});//, addedAtomsCount, addedBondsCount});
        // HashMap to store invalid extensions
//        final HashMap<String, HashSet<Long>> invalidExtensions = new HashMap<>();

        String extendedHOSECode;
        long j = 1, i;
        while (!intermediates.isEmpty()) {
            i = j;
//            for (long i = j; i < rankedSSCs.size(); i++) {
            while (i < rankedSSCList.size()){
                // do not compare intermediate with its start structure
                if(i == startSSCIndex){
                    i++;
                    continue;
                }
                path = (LinkedHashSet<Float>) intermediates.peek()[0];
                intermediate = (SSC) intermediates.peek()[1];
                extendedHOSECode = Compare.getExtendedHOSECode(intermediate);

                System.out.println("\n\n--> for path: " + path + "\n -> " + extendedHOSECode + "\n -> " + intermediate.getConnectionTree(intermediate.getRootAtomIndex()) + "\nnext ssc index: " + i + "/" + (rankedSSCList.size() - 1));
                System.out.println("-------------------------------- " + startSSCIndex + "_" + path + "_" + i + " --------------------------------");

                synchronized (invalidExtensions){
                    if(invalidExtensions.containsKey(extendedHOSECode) && invalidExtensions.get(extendedHOSECode).contains(i)){
                        System.out.println(" -> extension by " + i + " already failed before");
                        for (long k = i + 1; k < rankedSSCList.size(); k++) {
                            i = k;
                            if(!invalidExtensions.get(extendedHOSECode).contains(k)){
                                break;
                            }
                        }
                        if(i < rankedSSCList.size()){
                            if(i == (rankedSSCList.size() - 1)){
                                System.out.println(" -> last SSC reached -> break and go to next intermediate in stack");
                                break;
                            }
                            System.out.println(" -> skip and continue with next unchecked SSC: " + i);
                            continue;
                        }
                        System.out.println(" -> all next SSC extensions are not valid -> break and go to next intermediate in stack");
                        break;
                    }
                }


                newIntermediates = Assembly.assemblyCore(intermediate.getClone(false), rankedSSCList.get((int) i), querySpectrum, minMatchingSphereCount, shiftTol, thrsMatchFactor, molecularFormula);
                if ((newIntermediates.isEmpty())){
                    synchronized (invalidExtensions){
                        if(!invalidExtensions.containsKey(extendedHOSECode)){
                            invalidExtensions.put(extendedHOSECode, new HashSet<>());
                        }
                        invalidExtensions.get(extendedHOSECode).add(i);
                    }

                    i++;
                    continue;
                }

                // traverse over valid new intermediates in reversed order (pushes to stack) because output from AssemblyCore is sorted by similarity to query
                for (int k = newIntermediates.size() - 1; k >= 0; k--) {
                    newIntermediate = newIntermediates.get(k);

//                    try {
//                        Utils.generatePicture(newIntermediate.getSubstructure(), "results/temp_" + startSSCIndex + "_" + path + "-" + i + ".png");
//                    } catch (Exception e) {
//
//                    }

                    if (Assembly.isFinalSSC(newIntermediate, querySpectrum, shiftTol, thrsMatchFactor, molecularFormula)) {
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
//                        i++;
                        continue;
                    }
                    // valid and extended SSC built but not final, so add it to stack
                    newPath = new LinkedHashSet<>(path);
                    newPath.add(Float.parseFloat(i + "." + k));
                    intermediates.push(new Object[]{newPath, newIntermediates.get(k)});// addedAtomsCount, addedBondsCount});
                }

                i++;
            }
            if(!intermediates.isEmpty()){
                j = (long) (Collections.max((LinkedHashSet<Float>) intermediates.peek()[0]) + 1);
                intermediates.pop();

//                System.gc();
            }
        }
        bw.close();

        return solutions;
    }


    public static ArrayList<SSC> assemblyCore(final SSC ssc1, final SSC ssc2, final Spectrum querySpectrum, final int minMatchingSphereCount, final double shiftTol, final double thrsMatchFactor, final IMolecularFormula molecularFormula) throws Exception {

        // @TODO usage of shift tolerance window from min-max information and solvent effect information while searching for overlaps?
        // @TODO use combinations of valid overlaps to achieve closures of rings or cyclic structure
        HashMap<Integer, ArrayList<Integer[]>> overlapsHOSECodeNew = Assembly.getOverlaps(ssc1, ssc2, minMatchingSphereCount, shiftTol);
        if(overlapsHOSECodeNew.isEmpty()){
            return new ArrayList<>();
        }

        final ArrayList<SSC> validSSCExtensions = new ArrayList<>();
        ArrayList<Integer[]> overlapsHOSECodeInSphere;
        SSC ssc1Extended, ssc1ExtendedBackup;
        HashMap<Integer, Integer> atomMappings;
        ConnectionTree matchingConnectionTreeSSC1, matchingConnectionTreeSSC2, completeConnectionTreeSSC1, completeConnectionTreeSSC2, connectionTreeToAddSSC2;
        ArrayList<Integer> unsaturatedAtomsSSC1, connectionTreeKeysSSC1, connectionTreeKeysSSC2;
        ArrayList<ConnectionTreeNode> childNodesToAppend;
        int i, j, counter, indexUntilMaxMatchingSphere, nodeKeyInCompleteConnectionTreeSSC1, nodeKeyInCompleteConnectionTreeSSC2;
        IBond bondToAdd;

        // for each max. matching sphere (key); starting with highest
        for (int s = Collections.max(overlapsHOSECodeNew.keySet()); s >= minMatchingSphereCount; s--){
            if(!overlapsHOSECodeNew.containsKey(s)){
                continue;
            }

            overlapsHOSECodeInSphere = overlapsHOSECodeNew.get(s);
//            System.out.println("\n\n --> s = " + s + ": ");
//            for (final Integer[] indices : overlapsHOSECodeInSphere){
//                System.out.println(" -> " + Arrays.toString(indices));
//            }

            // for each overlapping atom pairs in SSC1 and SSC2 in sphere (maybe in a certain order?)
            for (int k = 0; k < overlapsHOSECodeInSphere.size(); k++) {
                // reset ssc1Extended to original SSC1
                ssc1Extended = ssc1.getClone(false);
                atomMappings = new HashMap<>();
//                System.out.println("\n -> k: " + k);
//                System.out.println(" -> HOSE code maps:" + Arrays.toString(overlapsHOSECodeInSphere.get(k)));
                // indices of matched (root) atoms
                i = overlapsHOSECodeInSphere.get(k)[0];
                j = overlapsHOSECodeInSphere.get(k)[1];

                if(ssc1Extended.isUnsaturated(i) || ssc2.isUnsaturated(j)){
//                    System.out.println(" atom " + i + " in SSC1 or atom " + j + " in SSC2 is unsaturated and not allowed as overlap root");
                    continue;
                }

                completeConnectionTreeSSC2 = HOSECodeBuilder.buildConnectionTree(ssc2.getSubstructure(), j, null);
//                ConnectionTree maxSphereConnectionTreeSSC2 = HOSECodeBuilder.buildConnectionTree(ssc2.getSubstructure(), j, ssc2.getMaxSphere());
                matchingConnectionTreeSSC1 = HOSECodeBuilder.buildConnectionTree(ssc1Extended.getSubstructure(), i, s);
                matchingConnectionTreeSSC2 = HOSECodeBuilder.buildConnectionTree(ssc2.getSubstructure(), j, s);

                if(!Compare.rearrangeInvalidNodePairsInSphere(ssc1Extended, ssc2, matchingConnectionTreeSSC1, matchingConnectionTreeSSC2, s, shiftTol)){
                    continue;
                }

//                System.out.println(" --> atoms in SSC1: " + matchingConnectionTreeSSC1.getKeys(true));
//                System.out.println(" --> atoms in SSC2: " + matchingConnectionTreeSSC2.getKeys(true));
//                System.out.println(ssc1Extended.getHOSECode(i) + "\n" + ssc2.getHOSECode(j));
//                System.out.println(ssc1Extended.getConnectionTree(i) + " -> " + ssc1Extended.getConnectionTree(i).getNodesCountInSphere(s) + "\n" + ssc2.getConnectionTree(j) + " -> " + ssc2.getConnectionTree(j).getNodesCountInSphere(s));
                unsaturatedAtomsSSC1 = new ArrayList<>();
                connectionTreeKeysSSC1 = new ArrayList<>(matchingConnectionTreeSSC1.getKeys());
                connectionTreeKeysSSC2 = new ArrayList<>(matchingConnectionTreeSSC2.getKeys());
                if(connectionTreeKeysSSC1.size() != connectionTreeKeysSSC2.size()){
//                    System.out.println(" <--- max. matching sphere is not the same anymore!!!! --->");
                    continue;
                }
                for (int l = 0; l < connectionTreeKeysSSC1.size(); l++) {
                    atomMappings.put(connectionTreeKeysSSC1.get(l), connectionTreeKeysSSC2.get(l));
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
//                System.out.println("predicted spectrum orig.: " + ssc1Extended.getSubspectrum().getShifts(0));
//                System.out.println(" -----> unsaturated atoms SSC1:" + unsaturatedAtomsSSC1);

                // check for each found unsaturated atom in SSC1 whether there is a valid extension possible from SSC2 via its children
                for (final int unsaturatedAtomKeySSC1 : unsaturatedAtomsSSC1){
                    indexUntilMaxMatchingSphere = connectionTreeKeysSSC1.indexOf(unsaturatedAtomKeySSC1);
                    nodeKeyInCompleteConnectionTreeSSC1 = connectionTreeKeysSSC1.get(indexUntilMaxMatchingSphere);
                    nodeKeyInCompleteConnectionTreeSSC2 = connectionTreeKeysSSC2.get(indexUntilMaxMatchingSphere);
//                    System.out.println(" ---> for unsaturated atom in SSC1: " + nodeKeyInCompleteConnectionTreeSSC1 + " -> " + nodeKeyInCompleteConnectionTreeSSC2);
                    
                    childNodesToAppend = completeConnectionTreeSSC2.getNode(nodeKeyInCompleteConnectionTreeSSC2).getChildNodes();
                    for (final ConnectionTreeNode childNodeToAppend : childNodesToAppend){
                        ssc1ExtendedBackup = ssc1Extended.getClone(false);
//                        System.out.println(" --> to append from SSC2: " + childNodeToAppend.getKey());

                        // if ring closure node
                        if(childNodeToAppend.isRingClosureNode()){

//                            System.out.println("ring closure node: " + "-" + childNodeToAppend.getParent().getKey() + " -> ring closure parent: " + childNodeToAppend.getRingClosureParent().getKey() + " -> mapped?: " + atomMappings.containsValue(childNodeToAppend.getRingClosureParent().getKey()));
                            if(atomMappings.containsValue(childNodeToAppend.getRingClosureParent().getKey())){
                                Integer parentkeySSC1 = null, ringClosureParentKeySSC1 = null;
                                for (final Entry<Integer, Integer> entry : atomMappings.entrySet()){
                                    if(entry.getValue().equals(childNodeToAppend.getParent().getKey())){
                                        parentkeySSC1 = entry.getKey();
                                    } else if(entry.getValue().equals(childNodeToAppend.getRingClosureParent().getKey())){
                                        ringClosureParentKeySSC1 = entry.getKey();
                                    }
                                    if((parentkeySSC1 != null) && (ringClosureParentKeySSC1 != null)){
                                        ssc1Extended.addBond(ssc2.getBond(childNodeToAppend.getParent().getKey(), childNodeToAppend.getRingClosureParent().getKey()), parentkeySSC1, ringClosureParentKeySSC1);
                                        break;
                                    }
                                }
                            }

                            continue;
                        }
                        bondToAdd = completeConnectionTreeSSC2.getBond(nodeKeyInCompleteConnectionTreeSSC2, childNodeToAppend.getKey());
                        if(!Utils.isValidBondAddition(ssc1Extended.getSubstructure(), nodeKeyInCompleteConnectionTreeSSC1, bondToAdd)){
                          continue;
                        }

                        connectionTreeToAddSSC2 = ConnectionTree.buildSubtree(completeConnectionTreeSSC2, childNodeToAppend.getKey());
//                        System.out.println("-> subtree would be: " + connectionTreeToAddSSC2 + "\n" + connectionTreeToAddSSC2.getKeys());
                        counter = 1;
                        for (final int substructureTreeNodeKeySSC2 : connectionTreeToAddSSC2.getKeys()) {
                            atomMappings.put((ssc1Extended.getAtomCount() - 1) + counter, substructureTreeNodeKeySSC2);
                            counter++;
                        }

                        ssc1Extended = Assembly.extendSSC(ssc1Extended, ssc2, connectionTreeToAddSSC2, nodeKeyInCompleteConnectionTreeSSC1, bondToAdd);
//                        System.out.println(" -> EXTENSION method done!!!");

                        // if one child node extends already existing or just invalid atoms then skip this child node and set to backup SSC
                        // @TODO what if valid the "right" extension comes after a "wrong" one? (i.e. at a ring node where to extend to the "outside")
                        if(!Assembly.isValidSSC(ssc1Extended.getSubspectrum(), querySpectrum, shiftTol, thrsMatchFactor, ssc1Extended.getSubstructure(), molecularFormula)){
                            ssc1Extended = ssc1ExtendedBackup.getClone(false);
                            continue;
                        }
//                        System.out.println("predicted spectrum temp.: " + ssc1Extended.getSubspectrum().getShifts(0));
//                        System.out.println(" -> after EXTENSION -> valid SSC!!!");

//                        System.out.println("atom mappings: " + atomMappings);

                        // try to close rings directly from added (and mapped) subtree of SSC2 to mapped atoms in SSC1
                        completeConnectionTreeSSC1 = HOSECodeBuilder.buildConnectionTree(ssc1Extended.getSubstructure(), i, null);
                        // for each node in added subtree
                        for (final int nodeKeyInSubtreeSSC2 : connectionTreeToAddSSC2.getKeys()){
                            ConnectionTreeNode nodeSSC2 = completeConnectionTreeSSC2.getNode(nodeKeyInSubtreeSSC2);
                            // if node is ring closure point
                            if(ConnectionTree.isAtRingClosure(nodeSSC2)){
                                // then check for each (additional) parent node whether there is a bond missing
                                for (final ConnectionTreeNode childNodeSSC2 : nodeSSC2.getChildNodes()){
                                    if(!childNodeSSC2.isRingClosureNode()){
                                        continue;
                                    }

                                    ConnectionTreeNode nodeSSC1 = null;
                                    ConnectionTreeNode ringClosureParentNodeSSC1 = null;
                                    ConnectionTreeNode ringClosureParentNodeSSC2 = childNodeSSC2.getRingClosureParent();
                                    if (atomMappings.containsValue(nodeSSC2.getKey())) {
                                        for (final Entry<Integer, Integer> entry : atomMappings.entrySet()) {
                                            if (entry.getValue() == nodeSSC2.getKey()) {
                                                nodeSSC1 = completeConnectionTreeSSC1.getNode(entry.getKey());
                                            }
                                            if (entry.getValue() == ringClosureParentNodeSSC2.getKey()) {
                                                ringClosureParentNodeSSC1 = completeConnectionTreeSSC1.getNode(entry.getKey());
                                            }
                                            if ((nodeSSC1 != null) && (ringClosureParentNodeSSC1 != null)) {
                                                break;
                                            }
                                        }
                                    }
                                    if ((nodeSSC1 != null) && (ringClosureParentNodeSSC1 != null)) {
                                        bondToAdd = ssc2.getBond(nodeSSC2.getKey(), ringClosureParentNodeSSC2.getKey());
                                        ssc1Extended.addBond(bondToAdd, nodeSSC1.getKey(), ringClosureParentNodeSSC1.getKey());
                                    }
                                }
                            }
                        }

                    }
                }
//                System.out.println("predicted spectrum modi.: " + ssc1Extended.getSubspectrum().getShifts(0));

                // for current HOSE code matching atom pair
                // if a valid substructure and also subspectrum could be assembled
                if(Assembly.isValidExtension(ssc1Extended.getSubspectrum(), querySpectrum, shiftTol, thrsMatchFactor, ssc1Extended.getSubstructure(), molecularFormula, ssc1.getSubstructure())){
//                    System.out.println("\nvalid SSC extension  built: ");
//                    System.out.println("HOSE code   : " + ssc1Extended.toHOSECode());
//                    System.out.println("substructure: " + ssc1Extended.getSubstructure().getAtomCount());
//                    System.out.println("subspectrum : " + ssc1Extended.getSubspectrum().getShifts(0));
//                    System.out.println("assignments : " + ssc1Extended.getAssignments().getAssignments(0) + "\n");
//                    System.out.println("unsat. atoms: " + ssc1Extended.getUnsaturatedAtomIndices());

                    if(validSSCExtensions.isEmpty()){
                        validSSCExtensions.add(ssc1Extended.getClone(false));
                    } else {
                        // check whether a same SSC was already built here
                        for (int l = 0; l < validSSCExtensions.size(); l++) {
                            if(Compare.compareSSC(validSSCExtensions.get(l), ssc1Extended)){
                                continue;
                            }
                            if(l == (validSSCExtensions.size() - 1)){
                                validSSCExtensions.add(ssc1Extended.getClone(false));
                            }
                        }
                    }

                }
            }
//            System.out.println("\n");
        }

        if (validSSCExtensions.size() > 1){
            validSSCExtensions.sort((validExtendedSSC1, validExtendedSSC2) -> {
//                // ranking by number of overlapping signals
//                final int setAssignmentsCountComp = -1 * Integer.compare(
//                        Matcher.matchSpectra(querySpectrum, validExtendedSSC1.getSubspectrum(), 0,0, shiftTol).getSetAssignmentsCount(0),
//                        Matcher.matchSpectra(querySpectrum, validExtendedSSC2.getSubspectrum(), 0,0, shiftTol).getSetAssignmentsCount(0));
//                if (setAssignmentsCountComp != 0) {
//                    return setAssignmentsCountComp;
//                }
                // ranking by atom counts
                final int atomCountComp = -1 * Integer.compare(validExtendedSSC1.getAtomCount(), validExtendedSSC2.getAtomCount());
                if(atomCountComp != 0){
                    return atomCountComp;
                }
                // ranking by bond counts
                final int bondCountComp = -1 * Integer.compare(validExtendedSSC1.getBondCount(), validExtendedSSC2.getBondCount());
                if(bondCountComp != 0){
                    return bondCountComp;
                }
                // ranking by average shift deviations (match factor)
                return Double.compare(
                        Matcher.calculateAverageDeviation(validExtendedSSC1.getSubspectrum(), querySpectrum, 0, 0, shiftTol),
                        Matcher.calculateAverageDeviation(validExtendedSSC2.getSubspectrum(), querySpectrum, 0, 0, shiftTol));
            });
        }

//        for (int k = 0; k < validSSCExtensions.size(); k++){
//            System.out.println("\n -------> for valid extension " + k + " with size " + validSSCExtensions.get(k).getAtomCount() + " and " + validSSCExtensions.get(k).getBondCount() + " and " + Matcher.calculateAverageDeviation(validSSCExtensions.get(k).getSubspectrum(), querySpectrum, 0, 0, shiftTol));
//            System.out.println("HOSE code   : " + validSSCExtensions.get(k).toHOSECode());
//            System.out.println("subspectrum : " + validSSCExtensions.get(k).getSubspectrum().getShifts(0));
//            System.out.println("multipl.    : " + validSSCExtensions.get(k).getSubspectrum().getMultiplicities());
//            System.out.println("unsat. atoms: " + validSSCExtensions.get(k).getUnsaturatedAtomIndices());
//        }

//        if(!validSSCExtensions.isEmpty()){
//            final ArrayList<SSC> validSSCExtensionsTemp = new ArrayList<>();
//            validSSCExtensionsTemp.add(validSSCExtensions.get(0));
//
//            return validSSCExtensionsTemp;
//        }
        return validSSCExtensions;
    }

    public static SSC extendSSC(final SSC ssc1, final SSC ssc2, final ConnectionTree connectionTreeToAddSSC, final Integer parentAtomIndexToLinkSSC1, final IBond bondToLinkSSC2) {

        final SSC ssc1Backup;
        try {
            ssc1Backup = ssc1.getClone(false);
        } catch (Exception e) {
            return ssc1;
        }

        final HashMap<Integer, Integer> insertedAtomMappings = new HashMap<>();
        // add root atom of connection tree to add to SSC1 and link it via a given bond to parent atom in SSC1
        Signal signalToAddSSC2 = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getIndex(0, connectionTreeToAddSSC.getRootNode().getKey()));
        if(!ssc1.addAtom(connectionTreeToAddSSC.getRootNode().getAtom(), bondToLinkSSC2, parentAtomIndexToLinkSSC1, signalToAddSSC2)){
            return ssc1Backup;
        }
        insertedAtomMappings.put(connectionTreeToAddSSC.getRootNode().getKey(), ssc1.getAtomCount() - 1);
        ArrayList<ConnectionTreeNode> nodesInSphere;
        ConnectionTreeNode nodeInSphere;
        IBond bondToAdd;
        // for each sphere: add the atom which is stored as node to atom container and set bonds between parent nodes
        for (int s = 1; s <= connectionTreeToAddSSC.getMaxSphere(); s++) {
            // first add all atoms and its parents (previous sphere only) to structure
            nodesInSphere = connectionTreeToAddSSC.getNodesInSphere(s, false);
            for (int i = 0; i < nodesInSphere.size(); i++) {
                nodeInSphere = nodesInSphere.get(i);
                bondToAdd = connectionTreeToAddSSC.getBond(nodeInSphere.getParent().getKey(), nodeInSphere.getKey());
                signalToAddSSC2 = ssc2.getSubspectrum().getSignal(ssc2.getAssignments().getIndex(0, nodeInSphere.getKey()));
                if(!ssc1.addAtom(nodeInSphere.getAtom(), bondToAdd, insertedAtomMappings.get(nodeInSphere.getParent().getKey()), signalToAddSSC2)){
                    return ssc1Backup;
                }
                insertedAtomMappings.put(nodeInSphere.getKey(), ssc1.getAtomCount() - 1);
            }
        }
        for (int s = 1; s <= connectionTreeToAddSSC.getMaxSphere(); s++) {
            // and as second add the remaining bonds (ring closures) to structure
            nodesInSphere = connectionTreeToAddSSC.getNodesInSphere(s, true);
            for (int i = 0; i < nodesInSphere.size(); i++) {
                nodeInSphere = nodesInSphere.get(i);
                if(!nodeInSphere.isRingClosureNode()){
                    continue;
                }
                ssc1.addBond(nodeInSphere.getBondToParent(), insertedAtomMappings.get(nodeInSphere.getRingClosureParent().getKey()), insertedAtomMappings.get(nodeInSphere.getParent().getKey()));
            }
        }
//        System.out.println("inserted atom mappings: " + insertedAtomMappings);

        return ssc1;
    }

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
//        System.out.println("\n" + subspectrum.getShifts(0) + "\n" + querySpectrum.getShifts(0));
        if(Utils.roundDouble(Matcher.calculateAverageDeviation(subspectrum, querySpectrum, 0, 0, shiftTol),  Start.DECIMAL_PLACES) > thrsMatchFactor) {
//            System.out.println("-> match factor not allowed!!!");
            return false;
        }

        return true;
    }

    public static boolean isValidSubstructure(final IAtomContainer substructure, final IMolecularFormula molecularFormula) {
        return Compare.compareWithMolecularFormula(substructure, molecularFormula);
    }

//    /**
//     * @param rankedSSCList
//     * @param startSSCIndex
//     * @param minMatchingSphereCount
//     * @param querySpectrum
//     * @param thrsMatchFactor
//     * @param shiftTol
//     * @return
//     * @throws Exception
//     *
//     * @deprecated
//     */
//    public static HashMap<String, SSC> assembleSeq(final ArrayList<SSC> rankedSSCList, final long startSSCIndex, final int minMatchingSphereCount,
//                                                   final Spectrum querySpectrum, final double thrsMatchFactor, final double shiftTol) throws Exception {
//
//        final SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Absolute);
//        String structureAsSMILES;
//        final HashMap<String, SSC> solutions = new HashMap<>();
//        SSC intermediate, backupSSC, ssc2, startSSC;
//        // notice: clone (!!!) the SSC contents only; don't use the object (reference) itself because of modifications
//        startSSC = rankedSSCList.get((int) startSSCIndex).getClone(false);
//        intermediate = startSSC.getClone(false);
//        intermediate.setIndex(startSSCIndex);
//        // check whether the current SSC is already a final SSC
//        if (Assembly.isFinalSSC(intermediate, querySpectrum, shiftTol, thrsMatchFactor)) {
//            structureAsSMILES = smilesGenerator.create(intermediate.getSubstructure());
//            solutions.put(structureAsSMILES, intermediate);
//            System.out.println("--> new solution found!!! -> " + solutions.size() + " -> " + structureAsSMILES);
//            System.out.println("-> atom count: " + intermediate.getAtomCount() + ", bond count: " + intermediate.getBondCount());
//            System.out.println("-> query spectrum:\t" + querySpectrum.getShifts(0));
//            System.out.println("-> equivalences:\t" + querySpectrum.getEquivalences());
//            System.out.println("-> pred. spectrum:\t" + intermediate.getSubspectrum().getShifts(0));
//            System.out.println("-> equivalences:\t" + intermediate.getSubspectrum().getEquivalences());
//
//            if(solutions.containsKey(null)){
//                System.out.println(" NULL key at 1");
//            }
//            if(solutions.containsValue(null)){
//                System.out.println(" NULL value at 1");
//            }
//            return solutions;
//        }
//
//        for (long i = 0; i < rankedSSCList.size(); i++) {
//            if(i == startSSCIndex){
//                continue;
//            }
//
//            System.out.println("\n\n-------------------------------- " + startSSCIndex + ", " + i + " --------------------------------");
//            backupSSC = intermediate.getClone(false);
//            ssc2 = rankedSSCList.get((int) i);
//            intermediate = null;//Assembly.assemblyCore(intermediate.getClone(false), ssc2, querySpectrum, minMatchingSphereCount, shiftTol, thrsMatchFactor);
//            if (intermediate == null) {
//                intermediate = backupSSC.getClone(false);
//                intermediate.setIndex(startSSCIndex);
//                continue;
//            }
//
//            try {
//                Utils.generatePicture(intermediate.getSubstructure(), "results/temp_" + startSSCIndex + "_" + i + ".png");
//            } catch (Exception e) {
//
//            }
//
//            if (Assembly.isFinalSSC(intermediate, querySpectrum, shiftTol, thrsMatchFactor)) {
//                structureAsSMILES = smilesGenerator.create(intermediate.getSubstructure());
//                if (!solutions.containsKey(structureAsSMILES)) {
//                    solutions.put(structureAsSMILES, intermediate);
//                    System.out.println("--> new solution found!!! -> " + solutions.size() + " -> " + structureAsSMILES);
//                    System.out.println("-> atom count: " + intermediate.getAtomCount() + ", bond count: " + intermediate.getBondCount());
//                    System.out.println("-> query spectrum:\t" + querySpectrum.getShifts(0));
//                    System.out.println("-> equivalences:\t" + querySpectrum.getEquivalences());
//                    System.out.println("-> pred. spectrum:\t" + intermediate.getSubspectrum().getShifts(0));
//                    System.out.println("-> equivalences:\t" + intermediate.getSubspectrum().getEquivalences());
//                }
//
//                intermediate = startSSC.getClone(false);//backupSSC1.getClone();
//                intermediate.setIndex(startSSCIndex);
////                continue;
//            }
//
//        }
//
//        if(solutions.containsKey(null)){
//            System.out.println(" NULL key at 2");
//        }
//        if(solutions.containsValue(null)){
//            System.out.println(" NULL value at 2");
//        }
//        return solutions;
//    }
}
