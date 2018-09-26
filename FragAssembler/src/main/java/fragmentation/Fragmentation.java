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
package fragmentation;

import model.SSC;
import casekit.NMR.Utils;
import casekit.NMR.model.Assignment;
import casekit.NMR.model.Signal;
import casekit.NMR.model.Spectrum;
import java.util.ArrayList;
import java.util.HashMap;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.HOSECodeGenerator;


public class Fragmentation {
    
    /**
     * Builds a set of substructure-subspectrum-correlations (SSC objects) from a 
     * structure for all its atoms by using a breadth first search 
     * with spherical limit. 
     *
     * @param ac IAtomContainer as structure
     * @param maxNoOfSpheres Spherical limit for building a substructure into 
     * all directions
     * @param atomType Atom type (element) used subspectrum creation. It has 
     * to be the same type as in used spectrum property.
     * @param NMRShiftDBSpectrumProp Spectrum property name/ID belonging 
     * to given atom container.
     * @return
     * @throws org.openscience.cdk.exception.CDKException
     * @see Fragmentation#buildSSC(org.openscience.cdk.interfaces.IAtomContainer, int, int, java.lang.String, java.lang.String) 
     */
    public static ArrayList<SSC> buildSSCs(final IAtomContainer ac, final int maxNoOfSpheres, final String atomType, final String NMRShiftDBSpectrumProp) throws CDKException{
        
        final ArrayList<SSC> SSCs = new ArrayList<>();
        for (int i = 0; i < ac.getAtomCount(); i++) {
            SSCs.add(Fragmentation.buildSSC(ac, i, maxNoOfSpheres, atomType, NMRShiftDBSpectrumProp));
        }
        
        return SSCs;
    }
    
    /**
     * Builds a substructure-subspectrum-correlation (SSC) object from a 
     * structure by using a breadth first search 
     * with spherical limit and starting point. 
     *
     * @param ac IAtomContainer as structure
     * @param rootAtomIndex Index of start atom
     * @param maxNoOfSpheres Spherical limit for building a substructure into 
     * all directions
     * @param atomType Atom type (element) used subspectrum creation. It has 
     * to be the same type as in used spectrum property.
     * @param NMRShiftDBSpectrumProp Spectrum property name/ID belonging 
     * to given atom container.
     * @return
     * @throws org.openscience.cdk.exception.CDKException
     * @see Fragmentation#BFS(org.openscience.cdk.interfaces.IAtomContainer, int, int, int, org.openscience.cdk.interfaces.IAtomContainer, int) 
     */
    public static SSC buildSSC(final IAtomContainer ac, final int rootAtomIndex, final int maxNoOfSpheres, final String atomType, final String NMRShiftDBSpectrumProp) throws CDKException{
        
        // spectrum property string has to be changed for general cases
        if(!Utils.setNMRShiftDBShiftsToAtomContainer(ac, NMRShiftDBSpectrumProp)){
            return null;
        }        
        final IAtomContainer substructure = Fragmentation.buildSubstructure(ac, rootAtomIndex, maxNoOfSpheres);
        
        final Spectrum subspectrum = Fragmentation.createSubspectrum(substructure, atomType);
        final Assignment assignment = Fragmentation.createAssignments(subspectrum, substructure, atomType);
        
        final HOSECodeGenerator hcg = new HOSECodeGenerator();
        final HashMap<String, ArrayList<Double>> hoseShifts = new HashMap<>();
        final HashMap<String, ArrayList<Integer>> hoseIndices = new HashMap<>();
        String hose;
        for (int i = 0; i < substructure.getAtomCount(); i++) {
            hose = hcg.getHOSECode(substructure, substructure.getAtom(i), maxNoOfSpheres);
            if (!hoseShifts.containsKey(hose)){
                hoseShifts.put(hose, new ArrayList<>());
                hoseIndices.put(hose, new ArrayList<>());
            }
            if (substructure.getAtom(i).getSymbol().equals(atomType) && (substructure.getAtom(i).getProperty(Utils.getNMRShiftConstant(atomType)) != null)) {
                hoseShifts.get(hose).add(substructure.getAtom(i).getProperty(Utils.getNMRShiftConstant(atomType)));
                hoseIndices.get(hose).add(i);
            }
        }
        
        return new SSC(subspectrum, assignment, substructure, hoseShifts, hoseIndices);
    }
    
    /**
     * Builds a substructure from a structure using a breadth first search 
     * with spherical limit and starting point. 
     *
     * @param ac IAtomContainer as structure
     * @param rootAtomIndex Index of start atom
     * @param maxNoOfSpheres Spherical limit for building a substructure into 
     * all directions
     * @return
     * @see Fragmentation#BFS(org.openscience.cdk.interfaces.IAtomContainer, int, int, int, org.openscience.cdk.interfaces.IAtomContainer, int) 
     */
    public static IAtomContainer buildSubstructure(final IAtomContainer ac, final int rootAtomIndex, final int maxNoOfSpheres) {

        return Fragmentation.BFS(ac, rootAtomIndex, -1, maxNoOfSpheres + 1, SilentChemObjectBuilder.getInstance().newAtomContainer(), 0);
    }

    /**
     * Builds a substructure from a structure using a breadth first search 
     * with spherical limit and starting point and its predecessor as well as 
     * search depth. 
     * This function is used for substructure building in 
     * {@link Fragmentation#buildSubstructure(org.openscience.cdk.interfaces.IAtomContainer, int, int)}.
     *
     * @param ac IAtomContainer as structure
     * @param currentAtomIndex Index of start atom
     * @param predecessorAtomIndex Index of predecessor atom of start atom
     * @param maxNoOfSpheres Spherical limit for building a substructure into 
     * all directions
     * @param fragmentToBuild until now build substructure to extend
     * @param depth In which depth this search is
     * @return
     */
    public static IAtomContainer BFS(final IAtomContainer ac, final int currentAtomIndex, final int predecessorAtomIndex, final int maxNoOfSpheres, IAtomContainer fragmentToBuild, final int depth) {

        if ((currentAtomIndex < 0 || currentAtomIndex >= ac.getAtomCount()) 
                || (depth >= maxNoOfSpheres)
                || ac.getAtom(currentAtomIndex).getSymbol().equals("H")) {
            return fragmentToBuild;
        }
        
        if (fragmentToBuild.indexOf(ac.getAtom(currentAtomIndex)) == -1){            
            fragmentToBuild.addAtom(ac.getAtom(currentAtomIndex));
        }        
        if ((predecessorAtomIndex != -1) && (fragmentToBuild.getBond(fragmentToBuild.getAtom(fragmentToBuild.indexOf(ac.getAtom(predecessorAtomIndex))), fragmentToBuild.getAtom(fragmentToBuild.indexOf(ac.getAtom(currentAtomIndex)))) == null)){
            fragmentToBuild.addBond(fragmentToBuild.indexOf(ac.getAtom(predecessorAtomIndex)), fragmentToBuild.indexOf(ac.getAtom(currentAtomIndex)), ac.getBond(ac.getAtom(predecessorAtomIndex), ac.getAtom(currentAtomIndex)).getOrder());
        }
        
        for (final IAtom connectedAtom : ac.getConnectedAtomsList(ac.getAtom(currentAtomIndex))) {
            fragmentToBuild = BFS(ac, ac.indexOf(connectedAtom), currentAtomIndex, maxNoOfSpheres, fragmentToBuild, depth + 1);
        }
        
        return fragmentToBuild;
    }
    
    /**
     * Creates a subspectrum for an atom container from its atoms, including set
     * NMR shift values in {@link IAtomContainer#getProperty(java.lang.Object)}.
     *
     * @param substructure IAtomContainer as substructure
     * @param atomType Atom type for which elements in atom container should be considered
     * @return
     */
    public static Spectrum createSubspectrum(final IAtomContainer substructure, final String atomType){
        final Spectrum subspectrum = new Spectrum(new String[]{Utils.getIsotopeIdentifier(atomType)});  
        for (int i = 0; i < substructure.getAtomCount(); i++) {
            if(substructure.getAtom(i).getSymbol().equals(atomType) && (substructure.getAtom(i).getProperty(Utils.getNMRShiftConstant(atomType)) != null)){                
                subspectrum.addSignal(new Signal(new String[]{Utils.getIsotopeIdentifier(atomType)}, new Double[]{substructure.getAtom(i).getProperty(Utils.getNMRShiftConstant(atomType))}));
            }
        }
        
        return subspectrum;
    }
    
    /**
     * Creates assignments between an substructure and its subspectrum.
     *
     * @param subspectrum  Spectrum as subspectrum
     * @param substructure IAtomContainer as substructure
     * @param atomType Atom type for which elements in atom container should be considered
     * @return
     */
    public static Assignment createAssignments(final Spectrum subspectrum, final IAtomContainer substructure, final String atomType){
        int counter = 0;
        final Assignment assignment = new Assignment(subspectrum);
        for (int i = 0; i < substructure.getAtomCount(); i++) {
            if (substructure.getAtom(i).getSymbol().equals(atomType) && (substructure.getAtom(i).getProperty(Utils.getNMRShiftConstant(atomType)) != null)) {
                assignment.setAssignment(0, counter, i);
                counter++;
            }
        }
        
        return assignment;
    }
}
