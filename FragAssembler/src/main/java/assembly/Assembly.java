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

import java.util.ArrayList;
import model.SSC;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.Bond;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Assembly {    
    
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
    
    
    public static SSC extendSSC(final SSC ssc1, final SSC ssc2) throws CloneNotSupportedException, CDKException{                
        final SSC newSSC = ssc1.clone();
        // SSC1 is still used for newSSC in the following code lines 
        int i, j;
        int parentAtomIndexInSphereSSC2, indexOfParentInSSC1, indexOfAddedAtomInSSC1;
        IAtom atomToAddInSSC1;
        Integer signalIndexToAddFromSSC2;
        IBond.Order bondOrderToAddInSSC1;
        IBond bondToAddInSSC1;
        for (final Integer[] identity : Assembly.getStructuralIdentities(newSSC, ssc2)) {
            i = identity[0]; j = identity[1];
            if (ssc2.getHOSECode(j).startsWith(newSSC.getHOSECode(i).substring(0, newSSC.getHOSECode(i).length() - 2))) {
//                System.out.println("\n--> i, j: " + i + ", " + j + " -> SSC1: " + newSSC.getHOSECode(i).substring(0, newSSC.getHOSECode(i).length() - 2) + " -> ssc2: " + ssc2.getHOSECode(j));
                for (int s = 0; s < newSSC.getMaxSphere(); s++) {
                    ArrayList<Integer> atomIndicesInHOSECodeSphereSSC1 = newSSC.getAtomIndicesInHOSECodeSpheres(i, s+1);
//                    ArrayList<Integer> parentAtomIndicesInHOSECodeSphereSSC1 = newSSC.getParentAtomIndicesInHOSECodeSpheres(i, s+1);
                    ArrayList<Integer> atomIndicesInHOSECodeSphereSSC2 = ssc2.getAtomIndicesInHOSECodeSpheres(j, s+1);
                    ArrayList<Integer> parentAtomIndicesInHOSECodeSphereSSC2 = ssc2.getParentAtomIndicesInHOSECodeSpheres(j, s+1);
//                    System.out.println("\nsphere " + (s+1) + ": ssc1 atoms current:\t\t" + newSSC.getAtomsInHOSECodeSpheres(i, s+1).size());
//                    System.out.println("sphere " + (s+1) + ": ssc1 indices current:\t\t" + atomIndicesInHOSECodeSphereSSC1);
//                    System.out.println("sphere " + (s+1) + ": ssc1 indices parents:\t\t" + parentAtomIndicesInHOSECodeSphereSSC1);
//                    System.out.println("sphere " + (s+1) + ": ssc2 atoms current:\t\t" + ssc2.getAtomsInHOSECodeSpheres(j, s+1).size());
//                    System.out.println("sphere " + (s+1) + ": ssc2 indices current:\t\t" + atomIndicesInHOSECodeSphereSSC2);
//                    System.out.println("sphere " + (s+1) + ": ssc2 indices parents:\t\t" + parentAtomIndicesInHOSECodeSphereSSC2 + "\n");

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
//                                System.out.println("\n-> atom count before:\t" + newSSC.getAtomCount() + ", bond count: " + newSSC.getBondCount());
//                                System.out.println("-> signals before:\t" + newSSC.getSubspectrum().getShifts(0));
//                                System.out.println("-> assignments before:\t" + Arrays.toString(newSSC.getAssignments().getAtomIndices(0)));
//                                System.out.println("-> unsaturated atoms before:\t" + newSSC.getUnsaturatedAtomIndices());
                                parentAtomIndexInSphereSSC2 = parentAtomIndicesInHOSECodeSphereSSC2.get(k);
                                if((parentAtomIndexInSphereSSC2 == -3) || !newSSC.getSubstructure().contains(ssc2.getSubstructure().getAtom(parentAtomIndexInSphereSSC2))
                                        && (s+1 >= 2)){
                                    // in case of a root atom in SSC2 to add in SSC1, there is no parent atom index which one can use for setting a new atom-bond-connection in SSC1 (just a dummy value of -3)
                                    // so we are looking here in the first sphere of root atom in SSC2 for same atom in previous sphere (s) of SSC1 then the actually current used sphere (s+1)
                                    for (final IAtom atomInFirstHOSECodeSphereSSC2 : ssc2.getAtomsInHOSECodeSpheres(atomIndicesInHOSECodeSphereSSC2.get(k), 1)) {
                                        if ((atomInFirstHOSECodeSphereSSC2 != null) && (newSSC.getAtomsInHOSECodeSpheres(i, s) != null)) {
                                            for (final IAtom atomInSthHOSECodeSphereSSC1 : newSSC.getAtomsInHOSECodeSpheres(i, s)) {
                                                if ((atomInSthHOSECodeSphereSSC1 != null) && atomInFirstHOSECodeSphereSSC2.equals(atomInSthHOSECodeSphereSSC1)) {
                                                    parentAtomIndexInSphereSSC2 = ssc2.getAtomIndicesInHOSECodeSpheres(atomIndicesInHOSECodeSphereSSC2.get(k), 1).get(ssc2.getAtomsInHOSECodeSpheres(atomIndicesInHOSECodeSphereSSC2.get(k), 1).indexOf(atomInFirstHOSECodeSphereSSC2));
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                    // if a parent could not be found then skip
                                    if(parentAtomIndexInSphereSSC2 < 0){
                                        continue;
                                    }
                                }
                                // if the parent for the new atom, which is to add in SSC1, could not be found then skip
                                indexOfParentInSSC1 = newSSC.getSubstructure().indexOf(ssc2.getSubstructure().getAtom(parentAtomIndexInSphereSSC2));
                                if(indexOfParentInSSC1 == -1){
                                    continue;
                                }
                                atomToAddInSSC1 = ssc2.getSubstructure().getAtom(atomIndicesInHOSECodeSphereSSC2.get(k));
                                // add the known but missing atom to substructure of SSC1 from SSC2
                                if (!newSSC.getSubstructure().contains(atomToAddInSSC1)) {
                                    newSSC.getSubstructure().addAtom(atomToAddInSSC1);
                                    indexOfAddedAtomInSSC1 = newSSC.getSubstructure().indexOf(atomToAddInSSC1);
                                    // set HOSE code to new added atom
                                    newSSC.setHOSECode(indexOfAddedAtomInSSC1, ssc2.getHOSECode(atomIndicesInHOSECodeSphereSSC2.get(k)));
                                    for (int s2 = 0; s2 < newSSC.getMaxSphere(); s2++) {
                                        // set atoms (of origin structure in SSC2) in HOSE code spheres of new added atom
                                        newSSC.setAtomsInHOSECodeSpheres(indexOfAddedAtomInSSC1, ssc2.getAtomsInHOSECodeSpheres(atomIndicesInHOSECodeSphereSSC2.get(k), s2+1), s2+1);
                                        // update atom indices and parents atom indices of atoms in sphere
                                        newSSC.updateAtomIndicesInHOSECodeSpheres(i, s2+1);
                                    }                                    
                                    // add signal to subspectrum and assignment for that to new attached atom
                                    signalIndexToAddFromSSC2 = ssc2.getAssignments().getSignalIndex(0, atomIndicesInHOSECodeSphereSSC2.get(k));
                                    if (signalIndexToAddFromSSC2 != null) { // if atom to add to SSC1 is not from subspectrum atom type then skip
                                        newSSC.getSubspectrum().addSignal(ssc2.getSubspectrum().getSignal(signalIndexToAddFromSSC2));
                                        newSSC.getAssignments().addAssignment(new int[]{indexOfAddedAtomInSSC1});
                                    }    
                                    // update the atom type indices in SSC
                                    newSSC.updateAtomTypeIndices();
                                    // update list of unsaturated atoms in SSC1
                                    newSSC.updateUnsaturatedAtomIndices();
                                    // update present multiplicities and shifts inm SSC1
                                    newSSC.updatePresenceMultiplicities();
//                                    System.out.println("-> atom added: " + indexOfAddedAtomInSSC1 + "(" + atomToAddInSSC1.getSymbol() + ")");
                                }                                 
                                
                                // add new bond if such bond not exists and if both belonging atoms are not already saturated
                                if ((newSSC.getSubstructure().getBond(atomToAddInSSC1, newSSC.getSubstructure().getAtom(indexOfParentInSSC1)) == null)
                                        && (newSSC.getUnsaturatedAtomIndices().contains(newSSC.getSubstructure().indexOf(atomToAddInSSC1)))
                                        && (newSSC.getUnsaturatedAtomIndices().contains(newSSC.getSubstructure().indexOf(newSSC.getSubstructure().getAtom(indexOfParentInSSC1))))
                                        ){
                                    bondOrderToAddInSSC1 = ssc2.getSubstructure().getBond(atomToAddInSSC1,
                                            ssc2.getSubstructure().getAtom(parentAtomIndexInSphereSSC2)).getOrder();                                    
                                    bondToAddInSSC1 = new Bond(atomToAddInSSC1, newSSC.getSubstructure().getAtom(indexOfParentInSSC1), 
                                            bondOrderToAddInSSC1);
                                    newSSC.getSubstructure().addBond(bondToAddInSSC1);
                                    // update list of unsaturated atoms in SSC1
                                    newSSC.updateUnsaturatedAtomIndices();
//                                    System.out.println("-> bond added: " + bondOrderToAddInSSC1 + " " + atomToAddInSSC1.getSymbol());
                                }  
//                                System.out.println("\n-> atom count after:\t" + newSSC.getAtomCount() + ", bond count: " + newSSC.getBondCount());
//                                System.out.println("-> signals after:\t" + newSSC.getSubspectrum().getShifts(0));
//                                System.out.println("-> assignments after:\t" + Arrays.toString(newSSC.getAssignments().getAtomIndices(0)));
//                                System.out.println("-> unsaturated atoms after:\t" + newSSC.getUnsaturatedAtomIndices());
                            }
                        }
                    }
                }
            }
        }
        
        return newSSC;
    }        
   
}
