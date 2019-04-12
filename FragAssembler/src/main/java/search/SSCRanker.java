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
package search;

import casekit.NMR.Utils;
import casekit.NMR.model.Assignment;
import casekit.NMR.model.Spectrum;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import match.Match;
import model.SSC;
import model.SSCLibrary;
import org.openscience.cdk.exception.CDKException;
import start.Start;

/**
 * Class to search for matches between a SSC library and query spectra.
 * 
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public final class SSCRanker {

    private final SSCLibrary sscLibrary;
    private int nThreads;
    private final HashMap<Long, Assignment> matchAssignments;
    private final HashMap<Long, Double> matchFactors;
    private final HashMap<Long, Float> tanimotoCoefficients;    
    private final ArrayList<Long> rankedSSCIndices;
    private final SSCLibrary rankedSSCLibrary;
    
    /**
     * Instanciates a new object of this class.
     * The number of threads to use for parallelization is set to 1 
     * by default.
     * 
     *
     * @param sscLibrary HashMap object consisting of SSC and their indices.
     */
    public SSCRanker(final SSCLibrary sscLibrary){
        this(sscLibrary, 1);
    }
    
    /**
     * Instanciates a new object of this class.     
     *
     * @param sscLibrary HashMap object consisting of SSC and their indices.
     * @param nThreads number of threads to use for matching and ranking of 
     * SSCs
     */
    public SSCRanker(final SSCLibrary sscLibrary, final int nThreads){
        this.sscLibrary = sscLibrary;
        this.setNThreads(nThreads);
        this.matchAssignments = new HashMap<>();
        this.matchFactors = new HashMap<>();
        this.tanimotoCoefficients = new HashMap<>();
        this.rankedSSCIndices = new ArrayList<>();
        this.rankedSSCLibrary = new SSCLibrary(this.nThreads);
    }     
    
    /**
     * Returns the used input SSC library which is used for finding hits 
     * for query spectra.
     *
     * @return
     */
    public SSCLibrary getSSCLibrary(){
        return this.sscLibrary;
    }
    
    /**
     * Sets the number of threads to use.
     *
     * @param nThreads number of threads
     * @return false if {@code nThreads < 1}
     */
    public boolean setNThreads(final int nThreads){
        if(nThreads < 1){
            return false;
        }
        this.nThreads = nThreads;
        
        return true;
    }
    
    /**
     * Returns the number of threads to use.
     *
     * @return
     */
    public int getNThreads(){
        return this.nThreads;
    }
        
    /**
     * Returns the match assignments between a query spectrum and SSCs in 
     * SSC library.
     * To create/update these match assignments use the rank function.
     *
     * @return 
     *  
     * @see #rank(casekit.NMR.model.Spectrum, double)  
     */
    public HashMap<Long, Assignment> getMatchAssignments(){
        return this.matchAssignments;
    }
    
    /**
     * Returns the match factors of SSCs in SSC library.
     * To create/update these match factors use the rank function.
     *
     * @return
     *  
     * @see #rank(casekit.NMR.model.Spectrum, double) 
     */
    public HashMap<Long, Double> getMatchFactors(){
        return this.matchFactors;
    }
    
    /**
     * Returns the tanimoto coefficients of SSCs in SSC library.
     * To create/update these coefficients use the rank function.
     *
     * @return
     *  
     * @see #rank(casekit.NMR.model.Spectrum, double) 
     */
    public HashMap<Long, Float> getTanimotoCoefficients(){
        return this.tanimotoCoefficients;
    }    
    
    /**
     * Returns a SSC library containing clones of the current matching SSCs 
     * for a query spectrum in ascending ranked order.
     * To create/update such ranked SSC library use the rank function.
     *
     * @return
     * 
     * @see #rank(casekit.NMR.model.Spectrum, double) 
     */
    public SSCLibrary getRankedSSCLibrary() {
        return this.rankedSSCLibrary;
    }
    
    public long getRankedSSCsCount(){
        return this.getRankedSSCLibrary().getSSCCount();
    }
                    
    /**
     * Ranks the SSCs in SSC library according to their match factor regarding 
     * the query spectrum. 
     * The results are then available in further class functions, see {@code @see}.
     *
     * @param querySpectrum Query spectrum
     * @param shiftTol Tolerance value [ppm] for shift matching
     * @throws java.lang.InterruptedException
     * @throws org.openscience.cdk.exception.CDKException
     * @throws java.lang.CloneNotSupportedException
     * 
     * @see #getMatchFactors()  
     * @see #getMatchAssignments() 
     * @see #getRankedMatchFactors() 
     * @see #getRankedTanimotoCoefficients() 
     * @see #getRankedSSCIndices()  
     * @see #getRankedSSCLibrary() 
     * 
     */
    public void rank(final Spectrum querySpectrum, final double shiftTol) throws InterruptedException, CDKException, CloneNotSupportedException{                               
        this.calculateMatchAssignments(querySpectrum, shiftTol);
        this.calculateMatchFactors(querySpectrum, shiftTol);
        this.calculateTanimotoCoefficients(querySpectrum);
 
        this.rankSSCIndices();
        this.buildRankedSSCLibrary();
    }
    
    private void buildRankedSSCLibrary() throws CDKException, CloneNotSupportedException {
        this.rankedSSCLibrary.removeAll();
        SSC rankedSSC;
        for (final long rankedSSCIndex : this.rankedSSCIndices) {
            rankedSSC = this.getSSCLibrary().getSSC(rankedSSCIndex).getClone();
            rankedSSC.setIndex(this.rankedSSCLibrary.getSSCCount());
            this.rankedSSCLibrary.insert(rankedSSC);
        }
    }
    
    private void calculateMatchAssignments(final Spectrum querySpectrum, final double shiftTol) throws InterruptedException{
        this.matchAssignments.clear();
        // initialize an executor for parallelization
        final ExecutorService executor = Utils.initExecuter(this.nThreads);
        final ArrayList<Callable<HashMap<Long, Assignment>>> callables = new ArrayList<>();
        // add all task to do
        for (final long sscIndex : this.sscLibrary.getSSCIndices()) {
            callables.add((Callable<HashMap<Long, Assignment>>) () -> {
                final HashMap<Long, Assignment> tempHashMap = new HashMap<>();
                tempHashMap.put(sscIndex, Match.matchSpectra(this.sscLibrary.getSSC(sscIndex).getSubspectrum(), querySpectrum, shiftTol));
                return tempHashMap;
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
                    if (!tempHashMap.values().contains(null)) {
                        this.matchAssignments.putAll(tempHashMap);
                    } 
                });
        // shut down the executor service
        Utils.stopExecuter(executor, 5);
    }
    
    private void calculateMatchFactors(final Spectrum querySpectrum, final double shiftTol) throws InterruptedException{
        this.matchFactors.clear();
        // initialize an executor for parallelization
        final ExecutorService executor = Utils.initExecuter(this.nThreads);
        final ArrayList<Callable<HashMap<Long, Double>>> callables = new ArrayList<>();
        // add all task to do
        for (final long sscIndex : this.sscLibrary.getSSCIndices()) {
            callables.add((Callable<HashMap<Long, Double>>) () -> {
                final HashMap<Long, Double> tempHashMap = new HashMap<>();
                final Assignment matchAssignment = matchAssignments.get(sscIndex);
                if (matchAssignment.isFullyAssigned(0)) {
                    tempHashMap.put(sscIndex, Utils.roundDouble(Match.calculateMatchFactor(this.sscLibrary.getSSC(sscIndex).getSubspectrum(), querySpectrum, shiftTol), Start.DECIMAL_PLACES));
                }                
                return tempHashMap;
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
                    if (!tempHashMap.values().contains(null)) {
                        this.matchFactors.putAll(tempHashMap);
                    }
                });
        // shut down the executor service
        Utils.stopExecuter(executor, 5);
    }
    
    private void calculateTanimotoCoefficients(final Spectrum querySpectrum) throws InterruptedException {
        this.tanimotoCoefficients.clear();
        // initialize an executor for parallelization
        final ExecutorService executor = Utils.initExecuter(this.nThreads);
        final ArrayList<Callable<HashMap<Long, Float>>> callables = new ArrayList<>();
        // add all task to do
        for (final long sscIndex : this.sscLibrary.getSSCIndices()) {
            callables.add((Callable<HashMap<Long, Float>>) () -> {
                final HashMap<Long, Float> tempHashMap = new HashMap<>();
                final Assignment matchAssignment = matchAssignments.get(sscIndex);
                if(matchAssignment.isFullyAssigned(0)){
                    final Spectrum matchedQuerySubspectrum = new Spectrum(querySpectrum.getNuclei());
                    for (final int signalIndexInQuerySpectrum : matchAssignment.getAtomIndices(0)) {
                        matchedQuerySubspectrum.addSignal(querySpectrum.getSignal(signalIndexInQuerySpectrum));
                    }
                    tempHashMap.put(sscIndex, Match.calculateTanimotoCoefficient(matchedQuerySubspectrum, this.sscLibrary.getSSC(sscIndex).getSubspectrum(), 0));
                }
                
                return tempHashMap;
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
                    if (!tempHashMap.values().contains(null)) {
                        this.tanimotoCoefficients.putAll(tempHashMap);
                    }
                });
        // shut down the executor service
        Utils.stopExecuter(executor, 5);
    }
    
    /**
     * Ranks SSC indices according their substructure size and match factor. 
     *
     */
    private void rankSSCIndices(){
        this.rankedSSCIndices.clear();
        // use indices of SSCs of the input SSC library which belonging match factor is valid (not null)
        this.rankedSSCIndices.addAll(this.matchFactors.keySet()); 
//        // use indices of SSCs of the input SSC library which belonging tanimoto coefficient was not null
//        this.rankedSSCIndices.addAll(this.tanimotoCoefficients.keySet()); 
        
        Collections.sort(this.rankedSSCIndices, new Comparator<Long>() {
            @Override
            public int compare(final Long indexSSC1, final Long indexSSC2) {     
                // ranking by number of overlapping signals
                final int setAssignmentsCountComp = -1 * Integer.compare(
                        matchAssignments.get(indexSSC1).getSetAssignmentsCount(0),
                        matchAssignments.get(indexSSC2).getSetAssignmentsCount(0));
                if (setAssignmentsCountComp != 0) {
                    return setAssignmentsCountComp;
                }
//                // ranking by tanimoto coefficient 
//                final int tanimotoCoefficientComp = -1 * Float.compare(tanimotoCoefficients.get(indexSSC1), tanimotoCoefficients.get(indexSSC2));
//                if (tanimotoCoefficientComp != 0) {
//                    return tanimotoCoefficientComp;
//                }
                // ranking by match factor 
                final int matchFactorComp = Double.compare(matchFactors.get(indexSSC1), matchFactors.get(indexSSC2));
                if(matchFactorComp != 0){
                    return matchFactorComp;
                }  
                // ranking by total subtructure size
                final int substructureSizeComp = -1 * Integer.compare(
                        sscLibrary.getSSC(indexSSC1).getAtomCount(),
                        sscLibrary.getSSC(indexSSC2).getAtomCount());
//                if (substructureSizeComp != 0) {
                    return substructureSizeComp;
//                }
//                // ranking by unsaturated atoms count (increasing)
//                final int unsaturatedAtomsCountComp = Integer.compare(
//                        sscLibrary.getSSC(indexSSC1).getUnsaturatedAtomIndices().size(),
//                        sscLibrary.getSSC(indexSSC2).getUnsaturatedAtomIndices().size());                
//                return unsaturatedAtomsCountComp;
            }
        });
    }        
}
