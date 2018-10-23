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
import casekit.NMR.model.Spectrum;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import model.SSC;

/**
 * Class to search for matches between a SSC library and query spectra.
 * 
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class SSCSearch {

    private final HashMap<Integer, SSC> SSCLibrary;
    private final HashMap<String, ArrayList<Double>> HOSECodeLookupTable;
    private int nThreads, minOverlap;
    private final HashMap<Integer, Double> matchFactors;
    private final ArrayList<Integer> rankedSSCIndices;
    
    /**
     * Instanciates a new object of this class.
     * The number of threads to use for parallelization and minimum number of 
     * overlaps between a SSC subspectrum and query spectrum are set to 1 
     * by default.
     * 
     *
     * @param SSCLibrary HashMap object consisting of SSC and their indices.
     */
    public SSCSearch(final HashMap<Integer, SSC> SSCLibrary){
        this.SSCLibrary = SSCLibrary;
        this.minOverlap = 1;
        this.nThreads = 1;
        this.HOSECodeLookupTable = new HashMap<>();
        this.matchFactors = new HashMap<>();
        this.rankedSSCIndices = new ArrayList<>();
    }
    
    /**
     * Creates a lookup table for generated HOSE codes regarding the given 
     * SSC library. Each HOSE code entry contains a list of belonging shift 
     * values. 
     *
     * @throws java.lang.InterruptedException
     */
    public void createHOSELookupTable() throws InterruptedException{  
        this.HOSECodeLookupTable.clear();
        for (final SSC ssc : this.SSCLibrary.values()) {
            if(ssc != null){
                Utils.combineHashMaps(this.HOSECodeLookupTable, ssc.getHOSELookupShifts());
            }            
        }
    }
    
    /**
     * Returns a lookup table for generated HOSE codes regarding the given 
     * SSC library. Each HOSE code entry contains a list of belonging shift 
     * values. This lookup table has to be created beforehand.
     *
     * @return
     * 
     * @see #createHOSELookupTable()
     */
    public HashMap<String, ArrayList<Double>> getHOSECodeLookupTable(){
        return this.HOSECodeLookupTable;
    }
    
    /**
     * Calculates the root mean square (RMS) values for shift lists in HOSE 
     * code lookup table. This lookup table has to be created beforehand.
     *
     * @return
     * 
     * @see #createHOSELookupTable()   
     */
    public HashMap<String, Double> getHOSECodeLookupTableRMS(){
        return Utils.getRMS(this.getHOSECodeLookupTable());
    }
    
    /**
     * Sets the minimum number of signal matches between SSCs and query spectra.
     *
     * @param minOverlap
     */
    public void setMinOverlap(final int minOverlap){
        this.minOverlap = minOverlap;
    }
    
    /**
     * Returns the minimum number of signal matches between SSCs and query 
     * spectra.
     * 
     * @return
     */
    public int getMinOverlap(){
        return this.minOverlap;
    }
    
    /**
     * Sets the number of threads to use.
     *
     * @param nThreads
     */
    public void setThreadNumber(final int nThreads){
        this.nThreads = nThreads;
    }
    
    /**
     * Returns the number of threads to use.
     *
     * @return
     */
    public int getThreadNumber(){
        return this.nThreads;
    }
        
    /**
     * Returns the match factors of SSCs in SSC library which were not null 
     * itself and its match factor against a query spectrum was not null.
     * To create/update these match factors the usage of the match 
     * function is needed.
     *
     * @return
     *  
     * @see #match(casekit.NMR.model.Spectrum, double) 
     */
    public HashMap<Integer, Double> getMatchFactors(){
        return this.matchFactors;
    }
    
    /**
     * Returns the ascending sorted match factors of SSCs in SSC library which 
     * were not null itself and its match factor against a query spectrum was 
     * not null.
     * To create/update these match factors the usage of the match 
     * function is needed.
     *
     * @return
     * 
     * @see #match(casekit.NMR.model.Spectrum, double) 
     */
    public ArrayList<Double> getRankedMatchFactors(){
        final ArrayList<Double> rankedMatchFactors = new ArrayList<>();
        for (final int SSCIndex : this.rankedSSCIndices) {
            rankedMatchFactors.add(this.matchFactors.get(SSCIndex));
        }
        
        return rankedMatchFactors;
    }
    
    /**
     * Returns the ascending sorted indices of SSCs in SSC library which were 
     * not null itself and its match factor against a query spectrum was not 
     * null.
     * To create/update these ranked SSC indices the usage of the match 
     * function is needed.
     *
     * @return
     * 
     * @see #match(casekit.NMR.model.Spectrum, double) 
     */
    public ArrayList<Integer> getRankedSSCIndices(){
        return this.rankedSSCIndices;
    }
                    
    /**
     * Searches for matches between the SSC library and the query spectrum.
     * Besides that, {@link model.SSC#getMatchFactor(casekit.NMR.model.Spectrum, double)}
     * is used to calculate match factors which will be stored. 
     * The results are then available in further class functions.
     * The number of minimum overlaps between a query spectrum and SSC library 
     * and the number of threads to use for parallelization can be set manually 
     * beforehand.
     *
     * @param querySpectrum Query spectrum
     * @param tol Tolerance value [ppm] for matching
     * @throws java.lang.InterruptedException
     * 
     * @see #getMatchFactors()  
     * @see #getRankedMatchFactors() 
     * @see #getRankedSSCIndices()  
     * @see #setMinOverlap(int) 
     * @see #setThreadNumber(int) 
     */
    public void match(final Spectrum querySpectrum, final double tol) throws InterruptedException{
        this.matchFactors.clear();        
        
        // initialize an executor for parallelization
        final ExecutorService executor = Utils.initExecuter(this.nThreads);  
        final ArrayList<Callable<HashMap<Integer, Double>>> callables = new ArrayList<>();
        // add all task to do
        for (final SSC ssc : this.SSCLibrary.values()) {
            callables.add((Callable<HashMap<Integer, Double>>) () -> {
                final HashMap<Integer, Double> tempHashMap = new HashMap<>(); 
                tempHashMap.put(ssc.getIndex(), ssc.getMatchFactor(querySpectrum, tol, this.minOverlap));
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
                    if(!tempHashMap.values().contains(null)){
                        this.matchFactors.putAll(tempHashMap);
                    }                    
                });
        // shut down the executor service
        Utils.stopExecuter(executor, 3);  
 
        this.rankSSCIndices();
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
            public int compare(final Integer o1, final Integer o2) {    
                if(SSCLibrary.get(o1).getAtomCount() >= SSCLibrary.get(o2).getAtomCount()){
                    if(matchFactors.get(o1) <= matchFactors.get(o2)){
                        return -1;
                    }                    
                } else {//if(SSCLibrary.get(o2).getAtomCount() > SSCLibrary.get(o1).getAtomCount()) {                
                    if(matchFactors.get(o2) < matchFactors.get(o1)){
                        return 1;
                    }
                }
                
                return Double.compare(matchFactors.get(o1), matchFactors.get(o2));
            }
        });
    }
}
