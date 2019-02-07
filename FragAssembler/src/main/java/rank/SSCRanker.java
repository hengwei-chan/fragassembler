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
package rank;

import assembly.Assembly;
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
import model.SSC;
import model.SSCLibrary;

/**
 * Class to search for matches between a SSC library and query spectra.
 * 
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public final class SSCRanker {

    private final SSCLibrary sscLibrary;
    private int nThreads;
    private final HashMap<Integer, Assignment> matchAssignments;
    private final HashMap<Integer, Double> matchFactors;    
    private final ArrayList<Integer> rankedSSCIndices;
    
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
        this.rankedSSCIndices = new ArrayList<>();
    }     
    
    /**
     * Returns the used SSC library object.
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
     * To create/update these ranked SSC indices use the rank function.
     *
     * @return 
     *  
     * @see #rank(casekit.NMR.model.Spectrum, double)  
     */
    public HashMap<Integer, Assignment> getMatchAssignments(){
        return this.matchAssignments;
    }
    
    /**
     * Returns the match factors of SSCs in SSC library.
     * To create/update these ranked SSC indices use the rank function.
     *
     * @return
     *  
     * @see #rank(casekit.NMR.model.Spectrum, double) 
     */
    public HashMap<Integer, Double> getMatchFactors(){
        return this.matchFactors;
    }
    
    /**
     * Returns the ascending sorted match factors of SSCs in SSC library.
     * To create/update these ranked SSC indices use the rank function.
     *
     * @return
     * 
     * @see #rank(casekit.NMR.model.Spectrum, double) 
     */
    public ArrayList<Double> getRankedMatchFactors(){
        final ArrayList<Double> rankedMatchFactors = new ArrayList<>();
        for (final int SSCIndex : this.rankedSSCIndices) {
            rankedMatchFactors.add(this.matchFactors.get(SSCIndex));
        }
        
        return rankedMatchFactors;
    }
    
    /**
     * Returns the ascending sorted indices of SSCs in SSC library.
     * To create/update these ranked SSC indices use the rank function.
     *
     * @return
     * 
     * @see #rank(casekit.NMR.model.Spectrum, double) 
     */
    public ArrayList<Integer> getRankedSSCIndices(){
        return this.rankedSSCIndices;
    }
                    
    /**
     * Ranks the SSCs in SSC library according to their match factor regarding 
     * the query spectrum.     
     * The results are then available in further class functions, see {@code @see}.
     * The number of threads to use for this procedure can be set 
     * beforehand with {@see #setNThreads(int)}.
     *
     * @param querySpectrum Query spectrum
     * @param pickPrecision Tolerance value [ppm] for shift matching
     * @throws java.lang.InterruptedException
     * 
     * @see #getMatchFactors()  
     * @see #getMatchAssignments() 
     * @see #getRankedMatchFactors() 
     * @see #getRankedSSCIndices()  
     * 
     */
    public void rank(final Spectrum querySpectrum, final double pickPrecision) throws InterruptedException{                               
        this.calculateMatchAssignments(querySpectrum, pickPrecision);
        this.calculateMatchFactors(querySpectrum, pickPrecision);
 
        this.rankSSCIndices();
    }
    
    private void calculateMatchAssignments(final Spectrum querySpectrum, final double pickPrecision) throws InterruptedException{
        this.matchAssignments.clear();
        // initialize an executor for parallelization
        final ExecutorService executor = Utils.initExecuter(this.nThreads);
        final ArrayList<Callable<HashMap<Integer, Assignment>>> callables = new ArrayList<>();
        // add all task to do
        for (final SSC ssc : this.sscLibrary.getSSCs()) {
            callables.add((Callable<HashMap<Integer, Assignment>>) () -> {
                final HashMap<Integer, Assignment> tempHashMap = new HashMap<>();
                tempHashMap.put(ssc.getIndex(), Assembly.matchSpectra(ssc.getSubspectrum(), querySpectrum, pickPrecision));
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
    
    private void calculateMatchFactors(final Spectrum querySpectrum, final double pickPrecision) throws InterruptedException{
        this.matchFactors.clear();
        // initialize an executor for parallelization
        final ExecutorService executor = Utils.initExecuter(this.nThreads);
        final ArrayList<Callable<HashMap<Integer, Double>>> callables = new ArrayList<>();
        // add all task to do
        for (final SSC ssc : this.sscLibrary.getSSCs()) {
            callables.add((Callable<HashMap<Integer, Double>>) () -> {
                final HashMap<Integer, Double> tempHashMap = new HashMap<>();                
                tempHashMap.put(ssc.getIndex(), Assembly.getMatchFactor(ssc.getSubspectrum(), querySpectrum, pickPrecision));
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
    
    /**
     * Ranks SSC indices according their substructure size and match factor. 
     *
     */
    private void rankSSCIndices(){
        this.rankedSSCIndices.clear();
        this.rankedSSCIndices.addAll(this.matchFactors.keySet()); // use SSC indices of SSC library
        
        Collections.sort(this.rankedSSCIndices, new Comparator<Integer>() {
            @Override
            public int compare(final Integer indexSSC1, final Integer indexSSC2) {          
                // ranking by match factor 
                final int matchFactorComp = Double.compare(matchFactors.get(indexSSC1), matchFactors.get(indexSSC2));
                if(matchFactorComp != 0){
                    return matchFactorComp;
                }                 
                // ranking by number of overlapping signals
                final int setAssignmentsCountComp = -1 * Integer.compare(
                        matchAssignments.get(indexSSC1).getSetAssignmentsCount(0),
                        matchAssignments.get(indexSSC2).getSetAssignmentsCount(0));
                if(setAssignmentsCountComp != 0){
                    return setAssignmentsCountComp; 
                }               
                // ranking by total subtructure size
                final int substructureSizeComp = -1 * Integer.compare(
                        sscLibrary.getSSC(indexSSC1).getAtomCount(), 
                        sscLibrary.getSSC(indexSSC2).getAtomCount());
                if(substructureSizeComp != 0){
                    return substructureSizeComp; 
                }
                // ranking by unsaturated atoms (increasing)
                final int unsaturatedAtomsCountComp = Integer.compare(
                        sscLibrary.getSSC(indexSSC1).getUnsaturatedAtomIndices().size(),
                        sscLibrary.getSSC(indexSSC2).getUnsaturatedAtomIndices().size());                
                return unsaturatedAtomsCountComp;
            }
        });
    }        
    
}
