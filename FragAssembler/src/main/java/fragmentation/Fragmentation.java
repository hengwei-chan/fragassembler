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
import casekit.NMR.model.Spectrum;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import model.SSCLibrary;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;


public class Fragmentation {
    
    /**
     * Builds a set of substructure-subspectrum-correlations (SSC objects) from  
     * an atom container set for all its molecules and atoms by using a 
     * breadth first search with spherical limit.
     *
     * @param SSCComponentsSet
     * @param maxNoOfSpheres Spherical limit for building a substructure into 
     * all directions
     * @param nThreads Number of threads to use for parallelization
     * @return
     * @throws java.lang.InterruptedException
     * @throws org.openscience.cdk.exception.CDKException
     * @throws java.lang.CloneNotSupportedException
     * @see Fragmentation#buildSSCs(org.openscience.cdk.interfaces.IAtomContainer, int, java.lang.String, java.lang.String) 
     */
    public static SSCLibrary buildSSCLibrary(final HashMap<Integer, Object[]> SSCComponentsSet, final int maxNoOfSpheres, final int nThreads) throws InterruptedException, CDKException, CloneNotSupportedException {                
        
        return Fragmentation.buildSSCLibrary(SSCComponentsSet, maxNoOfSpheres, nThreads, 0);
    }
    
    /**
     * Builds a set of substructure-subspectrum-correlations (SSC objects) from  
     * an atom container set for all its molecules and atoms by using a 
     * breadth first search with spherical limit.
     *
     * @param SSCComponentsSet
     * @param maxNoOfSpheres Spherical limit for building a substructure into 
     * all directions
     * @param nThreads Number of threads to use for parallelization
     * @param offset offset value as starting point for indexing the SSCs
     * @return
     * @throws java.lang.InterruptedException
     * @throws org.openscience.cdk.exception.CDKException
     * @throws java.lang.CloneNotSupportedException
     * @see Fragmentation#buildSSCs(org.openscience.cdk.interfaces.IAtomContainer, int, java.lang.String, java.lang.String) 
     */
    public static SSCLibrary buildSSCLibrary(final HashMap<Integer, Object[]> SSCComponentsSet, final int maxNoOfSpheres, final int nThreads, final long offset) throws InterruptedException, CDKException, CloneNotSupportedException {
        // initialize an executor
        final ExecutorService executor = Utils.initExecuter(nThreads);
        final SSCLibrary sscLibrary = new SSCLibrary();
        final ArrayList<Callable<SSCLibrary>> callables = new ArrayList<>();
        // add all task to do        
        long offsetSSCIndex = offset;
        for (final int index: SSCComponentsSet.keySet()) {
            final long offsetSSCIndexFinalCopy = offsetSSCIndex;           
            callables.add((Callable<SSCLibrary>) () -> Fragmentation.buildSSCs(SSCComponentsSet.get(index), maxNoOfSpheres, offsetSSCIndexFinalCopy));
            offsetSSCIndex += ((IAtomContainer) SSCComponentsSet.get(index)[0]).getAtomCount();
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
                .forEach((sscs) -> {
                    sscLibrary.extend(sscs);                    
                });
        // shut down the executor service
        Utils.stopExecuter(executor, 5);                        
        
        
        return sscLibrary;
    }
    
    /**
     * Builds a set of substructure-subspectrum-correlations (SSC objects) from one 
     * structure for all its atoms by using a breadth first search 
     * with spherical limit. 
     *
     * @param SSCComponentsSet
     * @param maxNoOfSpheres Spherical limit for building a substructure into 
     * all directions
     * to be the same type as in used spectrum property.
     * @param offsetSSCIndex  offset value for indexing each new built SSC
     * @return
     * @throws java.lang.CloneNotSupportedException
     * @see Fragmentation#buildSSC(org.openscience.cdk.interfaces.IAtomContainer, int, int, java.lang.String, java.lang.String) 
     */
    private static SSCLibrary buildSSCs(final Object[] SSCComponentsSet, final int maxNoOfSpheres, final long offsetSSCIndex) throws CloneNotSupportedException, CDKException {
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
                        assignment.setAssignment(0, assignment.getSignalIndex(0, prevAtomIndices.get(atom)), i);
                        break;
                    } 
                }
            }
        }
        final SSCLibrary sscLibrary = new SSCLibrary();        
        SSC ssc;
        for (int i = 0; i < structure.getAtomCount(); i++) {      
            ssc = Fragmentation.buildSSC(structure, spectrum, assignment, i, maxNoOfSpheres);
            if (ssc == null) {
                return new SSCLibrary();                
            }
            ssc.setIndex(offsetSSCIndex + i);
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
     * @param maxNoOfSpheres Spherical limit for building a substructure into 
     * all directions
     * 
     * @return
     * @throws org.openscience.cdk.exception.CDKException
     * @throws java.lang.CloneNotSupportedException
     * @see Fragmentation#buildSubstructure(org.openscience.cdk.interfaces.IAtomContainer, int, int) 
     */
    public static SSC buildSSC(final IAtomContainer structure, final Spectrum spectrum, final Assignment assignment, final int rootAtomIndex, final int maxNoOfSpheres) throws CDKException, CloneNotSupportedException{
        final ArrayList<Integer> substructureAtomIndices = Fragmentation.buildSubstructureAtomIndicesSet(structure, rootAtomIndex, maxNoOfSpheres);        
        final IAtomContainer substructure = Fragmentation.buildSubstructure(structure, rootAtomIndex, maxNoOfSpheres);
        final Spectrum subspectrum = new Spectrum(spectrum.getNuclei());
        final Assignment subassignment = new Assignment(subspectrum);
        for (int j = 0; j < substructure.getAtomCount(); j++) {
            if(structure.getAtom(substructureAtomIndices.get(j)).getSymbol().equals(Utils.getAtomTypeFromSpectrum(subspectrum, 0))){                
                if((assignment.getSignalIndex(0, substructureAtomIndices.get(j)) == null) || (spectrum.getSignal(assignment.getSignalIndex(0, substructureAtomIndices.get(j))) == null)){                    
                    return null;
                }                
                subspectrum.addSignal(spectrum.getSignal(assignment.getSignalIndex(0, substructureAtomIndices.get(j))));
                subassignment.addAssignment(new int[]{j});
            }            
        }
        subspectrum.setSolvent(spectrum.getSolvent());
        subspectrum.setSpectrometerFrequency(spectrum.getSpectrometerFrequency());
        Utils.setSpectrumEquivalences(subspectrum);

        // tries to return a valid SSC with all complete information
        // if something is missing/incomplete then null will be returned 
        try {
            return new SSC(subspectrum, subassignment, substructure, 0, maxNoOfSpheres);
        } catch (CloneNotSupportedException | CDKException e) {
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
     * @param maxNoOfSpheres Spherical limit for building a substructure into 
     * all directions
     * @return
     * @throws org.openscience.cdk.exception.CDKException
     * @see HOSECodeBuilder#buildConnectionTree(org.openscience.cdk.interfaces.IAtomContainer, int, int) 
     * @see HOSECodeBuilder#buildAtomContainer(model.ConnectionTree) 
     */
    public static IAtomContainer buildSubstructure(final IAtomContainer structure, final int rootAtomIndex, final int maxNoOfSpheres) throws CDKException {
        return HOSECodeBuilder.buildAtomContainer(HOSECodeBuilder.buildConnectionTree(structure, rootAtomIndex, maxNoOfSpheres));
    }
    
    /**
     * Builds a set of substructure atom indices from a structure using a 
     * breadth first search with spherical limit, starting point as well as 
     * HOSE code priority order of next neighbor atoms. 
     *
     * @param structure IAtomContainer as structure
     * @param rootAtomIndex Index of start atom
     * @param maxNoOfSpheres Spherical limit for building a substructure into 
     * all directions
     * @return
     * @throws org.openscience.cdk.exception.CDKException
     * @see HOSECodeBuilder#buildConnectionTree(org.openscience.cdk.interfaces.IAtomContainer, int, int) 
     */
    public static ArrayList<Integer> buildSubstructureAtomIndicesSet(final IAtomContainer structure, final int rootAtomIndex, final int maxNoOfSpheres) throws CDKException {
        return HOSECodeBuilder.buildConnectionTree(structure, rootAtomIndex, maxNoOfSpheres).getKeysInOrder(true);
    }
    
}
