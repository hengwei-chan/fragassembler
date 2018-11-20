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

import casekit.NMR.Utils;
import casekit.NMR.model.Assignment;
import casekit.NMR.model.Signal;
import casekit.NMR.model.Spectrum;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import model.SSC;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.Bond;
import search.SSCLibraryHolder;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Assembly {  
    
    /**
     * Assumes that the root atom index (start of HOSE code build) of each SSC is 0.
     * 
     * @param ssc1
     * @param ssc2
     * @return
     */
    public static ArrayList<Integer[]> getStructuralIdentities(final SSC ssc1, final SSC ssc2) {
        final ArrayList<Integer[]> identities = new ArrayList<>();
        for (int i = 0; i < ssc1.getAtomCount(); i++) {
            for (int j = 0; j < ssc2.getAtomCount(); j++) {
                if (ssc2.getHOSECode(j).startsWith(ssc1.getHOSECode(i).substring(0, ssc1.getHOSECode(i).length() - 2))) {
                    identities.add(new Integer[]{i, j});
                }
            }
        }
        
        return identities;                                 
    }
    
    /**
     * Assumes that the root atom index (start of HOSE code build) of each SSC is 0.
     *
     * @param ssc1
     * @param ssc2
     * @return
     */
    public static boolean hasStructuralIdentity(final SSC ssc1, final SSC ssc2) {
        boolean hasIdentity = false;
        for (int i = 0; i < ssc1.getAtomCount(); i++) {
            for (int j = 0; j < ssc2.getAtomCount(); j++) {
                if (ssc2.getHOSECode(j).startsWith(ssc1.getHOSECode(i).substring(0, ssc1.getHOSECode(i).length() - 2))) {
                    hasIdentity = true;
                    break;
                }
            }
        }
        
        return hasIdentity;                                 
    }
    
    
    public static SSC extendSSC(final SSC ssc1, final SSC ssc2, final int atomIndexInSSC1, final int atomIndexInSSC2, 
            final SSCLibraryHolder sscLibraryHandler, final Spectrum querySpectrum, 
            final double tol, final int minOverlap, final double thrsMatchFactor) throws CloneNotSupportedException, CDKException {        
        
        SSC SSCToExtend = ssc1.clone();
        // SSC1 is still used for SSCToExtend in the following code lines 
        int parentAtomIndexInSphereSSC2, indexOfParentInSSC1, indexOfAddedAtomInSSC1;
        IAtom atomToAddInSSC1;
        Integer signalIndexToAddFromSSC2;
        IBond.Order bondOrderToAddInSSC1;
        IBond bondToAddInSSC1;
        if (ssc2.getHOSECode(atomIndexInSSC2).startsWith(SSCToExtend.getHOSECode(atomIndexInSSC1).substring(0, SSCToExtend.getHOSECode(atomIndexInSSC1).length() - 2))) {
                System.out.println("\n--> i, j: " + atomIndexInSSC1 + ", " + atomIndexInSSC2 + " -> SSC1: " + SSCToExtend.getHOSECode(atomIndexInSSC1).substring(0, SSCToExtend.getHOSECode(atomIndexInSSC1).length() - 2) + " -> ssc2: " + ssc2.getHOSECode(atomIndexInSSC2));
            for (int s = 0; s < SSCToExtend.getMaxSphere(); s++) {
                ArrayList<Integer> atomIndicesInHOSECodeSphereSSC1 = SSCToExtend.getAtomIndicesInHOSECodeSpheres(atomIndexInSSC1, s+1);
                ArrayList<Integer> parentAtomIndicesInHOSECodeSphereSSC1 = SSCToExtend.getParentAtomIndicesInHOSECodeSpheres(atomIndexInSSC1, s+1);
                ArrayList<Integer> atomIndicesInHOSECodeSphereSSC2 = ssc2.getAtomIndicesInHOSECodeSpheres(atomIndexInSSC2, s+1);
                ArrayList<Integer> parentAtomIndicesInHOSECodeSphereSSC2 = ssc2.getParentAtomIndicesInHOSECodeSpheres(atomIndexInSSC2, s+1);
                   
                System.out.println("sphere " + (s+1) + ": ssc1 indices current:\t\t" + atomIndicesInHOSECodeSphereSSC1);
                System.out.println("sphere " + (s+1) + ": ssc1 indices parents:\t\t" + parentAtomIndicesInHOSECodeSphereSSC1);
                System.out.println("sphere " + (s+1) + ": ssc2 indices current:\t\t" + atomIndicesInHOSECodeSphereSSC2);
                System.out.println("sphere " + (s+1) + ": ssc2 indices parents:\t\t" + parentAtomIndicesInHOSECodeSphereSSC2 + "\n");

                for (int k = 0; k < atomIndicesInHOSECodeSphereSSC1.size(); k++) {
                    if ((atomIndicesInHOSECodeSphereSSC1.get(k) != null) 
                            && (atomIndicesInHOSECodeSphereSSC1.get(k) == -1) // a known but missing atom in sphere in SSC1 was found 
                            ) {                            
                        if (    (k < atomIndicesInHOSECodeSphereSSC2.size())
                                && (atomIndicesInHOSECodeSphereSSC2.get(k) != null)
                                && (atomIndicesInHOSECodeSphereSSC2.get(k) >= 0)
                                && (k < parentAtomIndicesInHOSECodeSphereSSC2.size())
                                && (parentAtomIndicesInHOSECodeSphereSSC2.get(k) != null)
                                && ((parentAtomIndicesInHOSECodeSphereSSC2.get(k) >= 0) || (parentAtomIndicesInHOSECodeSphereSSC2.get(k) == -3)) // second as special case for dummy parent index of root node
                                ) {
                            System.out.println("\n-> query:");
                            System.out.println("--> query signal shifts:\t" + querySpectrum.getShifts(0));
                            System.out.println("--> query signal mults:\t" + querySpectrum.getMultiplicities());
                            System.out.println("\n-> atom count before:\t" + SSCToExtend.getAtomCount() + ", bond count: " + SSCToExtend.getBondCount());
                            System.out.println("-> signal shifts before:\t" + SSCToExtend.getSubspectrum().getShifts(0));
                            System.out.println("-> signal mults before:\t" + SSCToExtend.getSubspectrum().getMultiplicities());
                            System.out.println("-> assignments before:\t" + Arrays.toString(SSCToExtend.getAssignments().getAtomIndices(0)));
                            System.out.println("-> unsaturated atoms before:\t" + SSCToExtend.getUnsaturatedAtomIndices());
                            parentAtomIndexInSphereSSC2 = parentAtomIndicesInHOSECodeSphereSSC2.get(k);
                            if((parentAtomIndexInSphereSSC2 == -3) 
                                    || !SSCToExtend.getSubstructure().contains(ssc2.getSubstructure().getAtom(parentAtomIndexInSphereSSC2))
                                    && (s+1 >= 2)){
                                // in case of a root atom in SSC2 to add in SSC1, there is no parent atom index which one can use for setting a new atom-bond-connection in SSC1 (just a dummy value of -3)
                                // so we are looking here in the first sphere of root atom in SSC2 for same atom in previous sphere (s) of SSC1 then the actually current used sphere (s+1)
                                boolean stop = false;
                                for (final IAtom atomInFirstHOSECodeSphereSSC2 : ssc2.getAtomsInHOSECodeSpheres(atomIndicesInHOSECodeSphereSSC2.get(k), 1)) {
                                    if ((atomInFirstHOSECodeSphereSSC2 != null) && (SSCToExtend.getAtomsInHOSECodeSpheres(atomIndexInSSC1, s) != null)) {
                                        for (final IAtom atomInSthHOSECodeSphereSSC1 : SSCToExtend.getAtomsInHOSECodeSpheres(atomIndexInSSC1, s)) {
                                            if ((atomInSthHOSECodeSphereSSC1 != null) && atomInFirstHOSECodeSphereSSC2.equals(atomInSthHOSECodeSphereSSC1)) {
                                                parentAtomIndexInSphereSSC2 = ssc2.getAtomIndicesInHOSECodeSpheres(atomIndicesInHOSECodeSphereSSC2.get(k), 1).
                                                        get(ssc2.getAtomsInHOSECodeSpheres(atomIndicesInHOSECodeSphereSSC2.get(k), 1).indexOf(atomInFirstHOSECodeSphereSSC2));
                                                stop = true;
                                                break;
                                            }
                                        }
                                    }
                                    if (stop) {
                                        break;
                                    }
                                }
                                // if a parent could not be found then skip
                                if(parentAtomIndexInSphereSSC2 < 0){
                                    System.out.println("parentAtomIndexInSphereSSC2 still under 0!!!");
                                    continue;
                                }
                            }
                            // if the parent for the new atom, which is to add in SSC1, could not be found then skip
                            indexOfParentInSSC1 = SSCToExtend.getSubstructure().indexOf(ssc2.getSubstructure().getAtom(parentAtomIndexInSphereSSC2));
                            if(indexOfParentInSSC1 == -1){
                                System.out.println("indexOfParentInSSC1 == -1 !!!");
                                continue;
                            }
                            atomToAddInSSC1 = ssc2.getSubstructure().getAtom(atomIndicesInHOSECodeSphereSSC2.get(k));
                            // add the known but missing atom to substructure of SSC1 from SSC2
                            if (!SSCToExtend.getSubstructure().contains(atomToAddInSSC1)) {
                                SSCToExtend.getSubstructure().addAtom(atomToAddInSSC1);
                                indexOfAddedAtomInSSC1 = SSCToExtend.getSubstructure().indexOf(atomToAddInSSC1);
                                // set HOSE code to new added atom
                                SSCToExtend.setHOSECode(indexOfAddedAtomInSSC1, ssc2.getHOSECode(atomIndicesInHOSECodeSphereSSC2.get(k)));
                                for (int s2 = 0; s2 < SSCToExtend.getMaxSphere(); s2++) {
                                    // set atoms (of origin structure in SSC2) in HOSE code spheres of new added atom
                                    SSCToExtend.setAtomsInHOSECodeSpheres(indexOfAddedAtomInSSC1, ssc2.getAtomsInHOSECodeSpheres(atomIndicesInHOSECodeSphereSSC2.get(k), s2+1), s2+1);
                                    // update atom indices and parents atom indices of atoms in sphere
                                    SSCToExtend.updateAtomIndicesInHOSECodeSpheres(atomIndexInSSC1, s2+1);
                                }      
                                // add signal to subspectrum and assignment for that to new attached atom
                                signalIndexToAddFromSSC2 = ssc2.getAssignments().getSignalIndex(0, atomIndicesInHOSECodeSphereSSC2.get(k));
                                if (signalIndexToAddFromSSC2 != null) { // if atom to add to SSC1 is not from subspectrum atom type then skip
                                    SSCToExtend.getSubspectrum().addSignal(ssc2.getSubspectrum().getSignal(signalIndexToAddFromSSC2));
                                    SSCToExtend.getAssignments().addAssignment(new int[]{indexOfAddedAtomInSSC1});
                                }
                                // update the atom type indices in SSC
                                SSCToExtend.updateAtomTypeIndices();
                                // update list of unsaturated atoms in SSC1
                                SSCToExtend.updateUnsaturatedAtomIndices();
                                // update present multiplicities and shifts inm SSC1
                                SSCToExtend.updatePresenceMultiplicities();
                                System.out.println("--> atom added: " + indexOfAddedAtomInSSC1 + "(" + atomToAddInSSC1.getSymbol() + ")");

                                // does predicted spectrum matches with query spectrum? -> move to Assembly.extendSSC() function in s2 section
                                if (!Assembly.isValidSSC(sscLibraryHandler, SSCToExtend, querySpectrum, tol, minOverlap, thrsMatchFactor)) {
//                                    return ssc1;
                                    System.out.println("-> !!!extension is not valid!!!");
                                    System.out.println("-> atom count after:\t" + SSCToExtend.getAtomCount() + ", bond count: " + SSCToExtend.getBondCount());
                                    System.out.println("-> signal shifts after:\t" + SSCToExtend.getSubspectrum().getShifts(0));
                                    System.out.println("-> signal mults after:\t" + SSCToExtend.getSubspectrum().getMultiplicities());
                                    System.out.println("-> assignments after:\t" + Arrays.toString(SSCToExtend.getAssignments().getAtomIndices(0)));
                                    System.out.println("-> unsaturated atoms after:\t" + SSCToExtend.getUnsaturatedAtomIndices() + "\n");
                                    continue;
                                }   
                            }                                 
                            // add new bond if such bond not exists and if both belonging atoms are not already saturated
                            if ((SSCToExtend.getSubstructure().getBond(atomToAddInSSC1, SSCToExtend.getSubstructure().getAtom(indexOfParentInSSC1)) == null)
                                    && (SSCToExtend.getUnsaturatedAtomIndices().contains(SSCToExtend.getSubstructure().indexOf(atomToAddInSSC1)))
                                    && (SSCToExtend.getUnsaturatedAtomIndices().contains(SSCToExtend.getSubstructure().indexOf(SSCToExtend.getSubstructure().getAtom(indexOfParentInSSC1))))
                                    && (ssc2.getSubstructure().getBond(atomToAddInSSC1, ssc2.getSubstructure().getAtom(parentAtomIndexInSphereSSC2)) != null)
                                    ){                                
                                bondOrderToAddInSSC1 = ssc2.getSubstructure().getBond(atomToAddInSSC1,
                                        ssc2.getSubstructure().getAtom(parentAtomIndexInSphereSSC2)).getOrder();                                    
                                bondToAddInSSC1 = new Bond(atomToAddInSSC1, SSCToExtend.getSubstructure().getAtom(indexOfParentInSSC1), 
                                        bondOrderToAddInSSC1);
                                SSCToExtend.getSubstructure().addBond(bondToAddInSSC1);
                                // update list of unsaturated atoms in SSC1
                                SSCToExtend.updateUnsaturatedAtomIndices();
                                    System.out.println("--> bond added: " + bondOrderToAddInSSC1 + " " + atomToAddInSSC1.getSymbol());
                            }  
                            System.out.println("\n-> !!!extension is valid!!!");
                            System.out.println("-> atom count after:\t" + SSCToExtend.getAtomCount() + ", bond count: " + SSCToExtend.getBondCount());
                            System.out.println("-> signal shifts after:\t" + SSCToExtend.getSubspectrum().getShifts(0));
                            System.out.println("-> signal mults after:\t" + SSCToExtend.getSubspectrum().getMultiplicities());
                            System.out.println("-> assignments after:\t" + Arrays.toString(SSCToExtend.getAssignments().getAtomIndices(0)));
                            System.out.println("-> unsaturated atoms after:\t" + SSCToExtend.getUnsaturatedAtomIndices());
                        }
                    }
                }
            }
        }
                    
        return SSCToExtend;
    }        
       
    public static Spectrum predictSpectrum(final SSC ssc, final SSCLibraryHolder sscLibraryHandler){
        final Spectrum predictedSpectrum = new Spectrum(ssc.getSubspectrum().getNuclei());
        int atomIndexCounter = 0;
        for (final int i: ssc.getAtomTypeIndices().get(ssc.getAtomType())) {
            if(sscLibraryHandler.getHOSECodeLookupTable().containsKey(ssc.getHOSECode(i))){
                predictedSpectrum.addSignal(new Signal(
                        ssc.getSubspectrum().getNuclei(), 
                        new Double[]{sscLibraryHandler.getHOSECodeLookupTableRMS().get(ssc.getHOSECode(i))},
                        ssc.getSubspectrum().getIntensity(atomIndexCounter),
                        ssc.getSubspectrum().getMultiplicity(atomIndexCounter)
                ));
            } else {
                return null;
            }
            atomIndexCounter++;
        }
        
        return predictedSpectrum;
    } 
    
    public static boolean isValidSSC(final SSCLibraryHolder sscLibraryHandler, final SSC ssc, final Spectrum querySpectrum, final double tol, final int minOverlap, final double thrsMatchFactor){
        final Spectrum predictedSpectrum = Assembly.predictSpectrum(ssc, sscLibraryHandler); // usage of predicted shifts by means of HOSE codes    
        final Assignment matchAssignments = Assembly.findMatches(predictedSpectrum, querySpectrum, tol);
        // filter for unset assignments
        if(matchAssignments.getSetAssignmentsCount(0) < matchAssignments.getAssignmentsCount()){
//            System.out.println("-> set assigments!!!");
            return false;
        }
        // filter for multiple assignments
        // there might be multiple assignments to same signals, so check for possible symmetry
        for (final int matchedSignalIndexInQuerySpectrum : matchAssignments.getAtomIndices(0)) {
            if (Collections.frequency(Utils.ArrayToArrayList(matchAssignments.getAtomIndices(0)), matchedSignalIndexInQuerySpectrum) 
                    > querySpectrum.getEquivalentSignals(matchedSignalIndexInQuerySpectrum).size() + 1) {                
                    return false;
//                }                
            }
        }   
        // filter for multiplicities and intensities, actually not necessary because of the HOSE code guideness
//        for (int i = 0; i < matchAssignments.getAssignmentsCount(); i++) {
//            if((predictedSpectrum.getMultiplicity(i) == null) 
//                    || (querySpectrum.getMultiplicity(matchAssignments.getAtomIndex(0, i)) == null) 
//                    || !predictedSpectrum.getMultiplicity(i).equals(querySpectrum.getMultiplicity(matchAssignments.getAtomIndex(0, i)))
//                    ){
//                return false;
//            }
////            if(Math.abs(predictedSpectrum.getIntensity(i) - this.querySpectrum.getIntensity(matchAssignments.getAtomIndex(0, i))) > this.toleranceIntensityMatching){
////                return false;
////            }
//        }   
        // filter for match factor
        if(Assembly.getMatchFactor(predictedSpectrum, querySpectrum, tol, minOverlap) > thrsMatchFactor){
            return false;
        }
        
        return true;
    }
        
    
    /**
     * Returns the closest shift matches between two spectra as an 
     * Assignment object.
     * Despite intensities are given, they are still not considered here.
     *
     * @param spectrum1
     * @param spectrum2
     * @param tol Tolerance value [ppm] used during shift matching
     * @return Assignments with signal indices of spectrum1 and matched indices 
     * in spectrum2; assignments are unset (-1) if spectrum1 has more signals 
     * then spectrum2
     */
    public static Assignment findMatches(final Spectrum spectrum1, final Spectrum spectrum2, final double tol) {
        final Assignment matchAssignments = new Assignment(spectrum1);
        // first nuclei in both spectra are not the same [or size of spectrum1 is bigger then in spectrum2]
        if (!spectrum1.getNuclei()[0].equals(spectrum2.getNuclei()[0]) 
//                || spectrum1.getSignalCount() > spectrum2.getSignalCount()
                ) {
            return matchAssignments;
        }
        final HashSet<Integer> pickedSignalIndices = new HashSet<>();
        int pickedSignalIndexSpectrum2;
        for (int i = 0; i < spectrum1.getSignalCount(); i++) {
            if(spectrum1.getShift(i, 0) == null){
                pickedSignalIndexSpectrum2 = -1;
            } else {           
                pickedSignalIndexSpectrum2 = spectrum2.pickClosestSignal(spectrum1.getShift(i, 0), 0, tol);
                // if matched signal is already assigned, then consider symmetries (equiv. signals)
                if (pickedSignalIndices.contains(pickedSignalIndexSpectrum2)) {
                    // symmetry exists
                    if (spectrum2.hasEquivalences(pickedSignalIndexSpectrum2)) {
                        // assign the next signal in equivalence list
                        for (final int equivalentSignalIndexSpectrum2 : spectrum2.getEquivalentSignals(pickedSignalIndexSpectrum2)) {
                            if(!pickedSignalIndices.contains(equivalentSignalIndexSpectrum2)){
                                pickedSignalIndexSpectrum2 = equivalentSignalIndexSpectrum2;
                                break;
                            }
                        }
                    } else {
                        // not symmetric signals but same HOSE code, 
                        // so the same (predicted) or very similar shifts and multiple assignments to catch
                        // -> still open
                        pickedSignalIndexSpectrum2 = -1;
                    }
                }
                // check multiplicity 
                if((spectrum1.getMultiplicity(i) == null) 
                        || (spectrum2.getMultiplicity(pickedSignalIndexSpectrum2) == null)
                        || !spectrum1.getMultiplicity(i).equals(spectrum2.getMultiplicity(pickedSignalIndexSpectrum2))) {
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
        
//        System.out.println("--> assignments before:\t" + Utils.ArrayToArrayList(matchAssignments.getAtomIndices(0)));
        // try to assign the still unassigned shifts in spectrum1 tp shifts in spectrum2
        for (int i = 0; i < matchAssignments.getAssignmentsCount(); i++) {
            if((matchAssignments.getAtomIndex(0, i) == -1) && (spectrum1.getShift(i, 0) != null)){
                for (final int pickedSignalIndexInSpectrum2 : spectrum2.pickSignals(spectrum1.getShift(i, 0), 0, tol)) {
                    if (!pickedSignalIndices.contains(pickedSignalIndexInSpectrum2)
                            && (spectrum1.getMultiplicity(i) != null)
                            && (spectrum2.getMultiplicity(pickedSignalIndexInSpectrum2) != null)
                            && spectrum1.getMultiplicity(i).equals(spectrum2.getMultiplicity(pickedSignalIndexInSpectrum2))
                            ){                            
                        matchAssignments.setAssignment(0, i, pickedSignalIndexInSpectrum2);
                        break;
                    }
                }
            }
        }
//        System.out.println("--> assignments after:\t" + Utils.ArrayToArrayList(matchAssignments.getAtomIndices(0)));

        return matchAssignments;
    }    

    /**
     * Returns deviatons between matched shifts in SSC and query query spectrum.
     * The matching procedure is already included here.
     *
     * @param spectrum1 
     * @param spectrum2
     * @param tol
     * @return
     *
     * @see #findMatches(casekit.NMR.model.Spectrum, casekit.NMR.model.Spectrum, double) 
     */
    public static Double[] getDeviations(final Spectrum spectrum1, final Spectrum spectrum2, final double tol) {
        final Double[] deviations = new Double[spectrum1.getSignalCount()];
        final Assignment matchAssignments = Assembly.findMatches(spectrum1, spectrum2, tol);

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
     * Returns the average of all deviations of matched shifts between two 
     * spectra.
     * The calculation of deviations is already included here.
     *
     * @param spectrum1 
     * @param spectrum2
     * @param tol Tolerance value [ppm] used during shift matching
     * @param minOverlap Minimum overlap threshold
     * @return
     *
     * @see #getDeviations(casekit.NMR.model.Spectrum, casekit.NMR.model.Spectrum, double) 
     * @see casekit.NMR.Utils#getMatchFactor(java.lang.Double[], int) 
     */
    public static Double getMatchFactor(final Spectrum spectrum1, final Spectrum spectrum2, final double tol, final int minOverlap) {
        return Utils.getMatchFactor(Assembly.getDeviations(spectrum1, spectrum2, tol), minOverlap);
    }

}
