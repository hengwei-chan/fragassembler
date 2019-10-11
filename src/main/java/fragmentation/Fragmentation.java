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
package fragmentation;

import casekit.NMR.Utils;
import casekit.NMR.model.Assignment;
import casekit.NMR.model.Spectrum;
import hose.HOSECodeBuilder;
import hose.model.ConnectionTree;
import model.SSC;
import model.SSCLibrary;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import parallel.ParallelTasks;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentLinkedQueue;


public class Fragmentation {
    
    /**
     * Builds a set of substructure-subspectrum-correlations (SSC objects) from
     * an atom container set for all its molecules and atoms by using a 
     * breadth first search with spherical limit.
     *
     * @param SSCComponentsSet
     * @param maxSphere Spherical limit for building a substructure into 
     * all directions
     * @param nThreads Number of threads to use for parallelization
     *
     * @return
     *
     * @throws java.lang.InterruptedException
     *
     * @see Fragmentation#buildSSCs(Object[], int)
     */
    public static SSCLibrary buildSSCLibrary(final HashMap<Integer, Object[]> SSCComponentsSet, final int maxSphere, final int nThreads) throws InterruptedException {

        final ConcurrentLinkedQueue<SSC> buildSSCs = new ConcurrentLinkedQueue<>();
        final ArrayList<Callable<SSCLibrary>> callables = new ArrayList<>();
        // add all task to do
        for (final int index: SSCComponentsSet.keySet()) {
            callables.add(() -> Fragmentation.buildSSCs(SSCComponentsSet.get(index), maxSphere));
        }
        ParallelTasks.processTasks(callables, sscLibraryTemp -> buildSSCs.addAll(sscLibraryTemp.getSSCs()), nThreads);

        final SSCLibrary sscLibrary = new SSCLibrary(nThreads);
        sscLibrary.extend(buildSSCs);

        return sscLibrary;
    }
    
    /**
     * Builds a set of substructure-subspectrum-correlations (SSC objects) from one
     * structure for all its atoms by using a breadth first search 
     * with spherical limit. 
     *
     * @param SSCComponentsSet
     * @param maxSphere Spherical limit for building a substructure into 
     * all directions
     * to be the same type as in used spectrum property.
     *
     * @return
     * @see Fragmentation#buildSSC(IAtomContainer, Spectrum, Assignment, int, int)
     */
    private static SSCLibrary buildSSCs(final Object[] SSCComponentsSet, final int maxSphere) throws CDKException {
        final IAtomContainer structure = (IAtomContainer) SSCComponentsSet[0];
        final Spectrum spectrum = (Spectrum) SSCComponentsSet[1];
        final Assignment assignment = (Assignment) SSCComponentsSet[2];
        // if the structure contains explicit hydrogens atoms then, after 
        // removing them, the assignments have to be corrected 
        if(Utils.containsExplicitHydrogens(structure)){
            final HashMap<IAtom, Integer> prevAtomIndices = Utils.convertExplicitToImplicitHydrogens(structure);  
            for (final IAtom atom : prevAtomIndices.keySet()) {
                // consider only atoms with same atom type as first spectrum nucleus type
                if(!atom.getSymbol().equals(Utils.getAtomTypeFromSpectrum(spectrum, 0))){
                    continue;
                }
                for (int i = 0; i < structure.getAtomCount(); i++) {
                    if((structure.getAtom(i) == atom) && (i != prevAtomIndices.get(atom))){
                        assignment.setAssignment(0, assignment.getIndex(0, prevAtomIndices.get(atom)), i);
                        break;
                    } 
                }
            }
        }
        final SSCLibrary sscLibrary = new SSCLibrary();        
        SSC ssc;
        for (int i = 0; i < structure.getAtomCount(); i++) {      
            ssc = Fragmentation.buildSSC(structure, spectrum, assignment, i, maxSphere);
            // if one part of the structure could not be built then skip the whole structure
            // and return an empty SSC library
            if (ssc == null) {
                return new SSCLibrary();
            }
            sscLibrary.insert(ssc);
        }
       
        return sscLibrary;
    }        
    
    /**
     * Builds a substructure-subspectrum-correlation ({@link model.SSC}) object 
     * from a structure, a spectrum and signal to atom assignments.
     * The structure fragmentation is done by using breadth first 
     * search with a spherical limit and each atom as starting point. 
     *
     * @param structure structure to fragment into substructures
     * @param spectrum spectrum to split into subspectra
     * @param assignment signal to atom assignments
     * @param rootAtomIndex Index of start atom
     * @param maxSphere Spherical limit for building a substructure into 
     * all directions
     * 
     * @return
     * @throws org.openscience.cdk.exception.CDKException
     * @see Fragmentation#buildSubstructure(org.openscience.cdk.interfaces.IAtomContainer, int, int)
     */
    public static SSC buildSSC(final IAtomContainer structure, final Spectrum spectrum, final Assignment assignment, final int rootAtomIndex, final int maxSphere) throws CDKException {
        Utils.setAromaticityAndKekulize(structure);
        final ArrayList<Integer> substructureAtomIndices = new ArrayList<>(Fragmentation.buildSubstructureAtomIndicesSet(structure, rootAtomIndex, maxSphere));
        final IAtomContainer substructure = Fragmentation.buildSubstructure(structure, rootAtomIndex, maxSphere);
        final Spectrum subspectrum = new Spectrum(spectrum.getNuclei());
        final Assignment subassignment = new Assignment(subspectrum);
        IAtom atomInStructure;
        for (int j = 0; j < substructure.getAtomCount(); j++) {
            atomInStructure = structure.getAtom(substructureAtomIndices.get(j));
            if(atomInStructure.getSymbol().equals(Utils.getAtomTypeFromSpectrum(subspectrum, 0))){
                if((assignment.getIndex(0, substructureAtomIndices.get(j)) == null) || (spectrum.getSignal(assignment.getIndex(0, substructureAtomIndices.get(j))) == null)){
                    return null;
                }
                subspectrum.addSignal(spectrum.getSignal(assignment.getIndex(0, substructureAtomIndices.get(j))));
                subassignment.addAssignment(new int[]{j});                
            }
        }
        subspectrum.setSolvent(spectrum.getSolvent());
        subspectrum.setSpectrometerFrequency(spectrum.getSpectrometerFrequency());
        subspectrum.detectEquivalences();
        // tries to return a valid SSC with all complete information
        // if something is missing/incomplete then null will be returned 
        try {
            return new SSC(subspectrum, subassignment, substructure, 0, maxSphere);
        } catch (Exception e) {
            return null;
        }         
    }
    
    /**
     * Builds a substructure from a structure using a breadth first search 
     * with spherical limit, starting point as well as HOSE code priority order
     * of next neighbor atoms. 
     *
     * @param structure IAtomContainer as structure
     * @param rootAtomIndex Index of start atom
     * @param maxSphere Spherical limit for building a substructure into 
     * all directions
     * @return
     * @throws org.openscience.cdk.exception.CDKException
     * @see HOSECodeBuilder#buildConnectionTree(IAtomContainer, int, Integer)
     * @see HOSECodeBuilder#buildAtomContainer(ConnectionTree)
     */
    public static IAtomContainer buildSubstructure(final IAtomContainer structure, final int rootAtomIndex, final int maxSphere) throws CDKException {
        return HOSECodeBuilder.buildAtomContainer(HOSECodeBuilder.buildConnectionTree(structure, rootAtomIndex, maxSphere));
    }
    
    /**
     * Builds a set of substructure atom indices from a structure using a 
     * breadth first search with spherical limit, starting point as well as 
     * HOSE code priority order of next neighbor atoms. 
     *
     * @param structure IAtomContainer as structure
     * @param rootAtomIndex Index of start atom
     * @param maxSphere Spherical limit for building a substructure into 
     * all directions
     * @return
     * @throws org.openscience.cdk.exception.CDKException
     * @see HOSECodeBuilder#buildConnectionTree(IAtomContainer, int, Integer)
     */
    public static LinkedHashSet<Integer> buildSubstructureAtomIndicesSet(final IAtomContainer structure, final int rootAtomIndex, final int maxSphere) throws CDKException {
        return HOSECodeBuilder.buildConnectionTree(structure, rootAtomIndex, maxSphere).getKeys(true);
    }
}
