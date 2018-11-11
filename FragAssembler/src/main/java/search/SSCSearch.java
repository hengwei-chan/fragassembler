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
import casekit.NMR.model.Signal;
import casekit.NMR.model.Spectrum;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
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
    private final HashMap<Integer, Assignment> matchAssignments;
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
        this.matchAssignments = new HashMap<>();
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
                Utils.combineHashMaps(this.HOSECodeLookupTable, ssc.getHOSECodeLookupShifts());
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
     * Returns the match assignments between a query spectrum and SSCs in 
     * SSC library which were not null itself.
     * To create/update these match assignments the usage of the match 
     * function is needed.
     *
     * @return 
     *  
     * @see #match(casekit.NMR.model.Spectrum, double) 
     */
    public HashMap<Integer, Assignment> getMatchAssignments(){
        return this.matchAssignments;
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
        this.calculateMatchAssignments(querySpectrum, tol);
        this.calculateMatchFactors(querySpectrum, tol);
 
        this.rankSSCIndices();
    }
    
    private void calculateMatchAssignments(final Spectrum querySpectrum, final double tol) throws InterruptedException{
        this.matchAssignments.clear();
        // initialize an executor for parallelization
        final ExecutorService executor = Utils.initExecuter(this.nThreads);
        final ArrayList<Callable<HashMap<Integer, Assignment>>> callables = new ArrayList<>();
        // add all task to do
        for (final SSC ssc : this.SSCLibrary.values()) {
            callables.add((Callable<HashMap<Integer, Assignment>>) () -> {
                final HashMap<Integer, Assignment> tempHashMap = new HashMap<>();
                tempHashMap.put(ssc.getIndex(), findMatches(ssc, querySpectrum, tol));
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
    
    private void calculateMatchFactors(final Spectrum querySpectrum, final double tol) throws InterruptedException{
        this.matchFactors.clear();
        // initialize an executor for parallelization
        final ExecutorService executor = Utils.initExecuter(this.nThreads);
        final ArrayList<Callable<HashMap<Integer, Double>>> callables = new ArrayList<>();
        // add all task to do
        for (final SSC ssc : this.SSCLibrary.values()) {
            callables.add((Callable<HashMap<Integer, Double>>) () -> {
                final HashMap<Integer, Double> tempHashMap = new HashMap<>();
                tempHashMap.put(ssc.getIndex(), getMatchFactor(getDeviations(ssc, querySpectrum, ssc.getIndex()), this.minOverlap));
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
            public int compare(final Integer o1, final Integer o2) {          
                final int setAssignmentsCountComp = -1 * Integer.compare(
                        matchAssignments.get(o1).getSetAssignmentsCount(0), 
                        matchAssignments.get(o2).getSetAssignmentsCount(0));
                if(setAssignmentsCountComp != 0){
                    return setAssignmentsCountComp;
                }
                return Double.compare(matchFactors.get(o1), matchFactors.get(o2));
            }
        });
    }
    
    /**
     * Returns the full list of matched shifts as atom indices between a SSC
     * and a query spectrum.
     * Despite intensities are given, they are still not considered here.
     *
     * @param ssc
     * @param multiplicity Multiplicity as multiplet string, e.g. "S" (Singlet)
     * @param shift Shift value [ppm]
     * @param intensity Intensity value
     * @param tol Tolerance value [ppm] used during shift matching
     * @return
     */
    public ArrayList<Integer> findMatches(final SSC ssc, final String multiplicity, final double shift, final double intensity, final double tol) {
        if (ssc.getPresenceMultiplicities().get(multiplicity) == null) {
            return new ArrayList<>();
        }
        final HashSet<Integer> hs = new HashSet<>();
        for (int i = (int) (shift - tol); i <= (int) (shift + tol); i++) {
            if (Utils.checkMinMaxValue(ssc.getMinShift(), ssc.getMaxShift(), i)) {
                hs.addAll(ssc.getPresenceMultiplicities().get(multiplicity).get(i));
            }
        }

        return new ArrayList<>(hs);
    }

    /**
     * Returns the closest shift matches between a SSC and a
     * query spectrum as an Assignment object.
     * Despite intensities are given, they are still not considered here.
     *
     * @param ssc
     * @param querySpectrum Query spectrum
     * @param tol Tolerance value [ppm] used during shift matching
     * @return
     */
    public Assignment findMatches(final SSC ssc, final Spectrum querySpectrum, final double tol) {
        final Assignment matchAssignmentsSSC = new Assignment(querySpectrum);
        // wrong nucleus in query spectrum or this SSC contains no atoms of requested atom type
        if (!Utils.getElementIdentifier(querySpectrum.getNuclei()[0]).equals(ssc.getAtomType())
                || (ssc.getAtomTypeIndices().get(ssc.getAtomType()) == null)) {
            return matchAssignmentsSSC;
        }
        Signal signal;
        for (int i = 0; i < querySpectrum.getSignalCount(); i++) {
            signal = querySpectrum.getSignal(i);
            matchAssignmentsSSC.setAssignment(0, i, this.getClosestMatch(ssc, this.findMatches(ssc, signal.getMultiplicity(), signal.getShift(0), signal.getIntensity(), tol), signal.getShift(0), tol));
        }

        // too many atoms of that atom type in SSC, but check for possible symmetry
        if (ssc.getAtomTypeIndices().get(ssc.getAtomType()).size() > querySpectrum.getSignalCount()
                || matchAssignmentsSSC.getSetAssignmentsCount(0) > ssc.getAtomTypeIndices().get(ssc.getAtomType()).size()) {
//            System.out.println("Check for multiple assignments/symmetry!!!");
            return new Assignment(querySpectrum);
        }

        return matchAssignmentsSSC;
    }

    private int getClosestMatch(final SSC ssc, final ArrayList<Integer> matchAtomIndices, final double queryShift, final double tol) {
        int closestMatchIndex = -1;
        double diff = tol, shiftMatchIndex;
        for (final int matchIndex : matchAtomIndices) {
            shiftMatchIndex = ssc.getSubstructure().getAtom(matchIndex).getProperty(Utils.getNMRShiftConstant(ssc.getAtomType()));
            if (Math.abs(shiftMatchIndex - queryShift) < diff) {
                diff = Math.abs(shiftMatchIndex - queryShift);
                closestMatchIndex = matchIndex;
            }
        }

        return closestMatchIndex;
    }

    /**
     * Returns deviatons between matched shifts in SSC and query query spectrum.
     * The matching procedure is already included here.
     *
     * @param ssc
     * @param querySpectrum
     * @param sscIndex
     * @return
     *
     * @see #findMatches(model.SSC, casekit.NMR.model.Spectrum, double) 
     */
    public Double[] getDeviations(final SSC ssc, final Spectrum querySpectrum, final int sscIndex) {
        final Double[] deviations = new Double[this.matchAssignments.get(sscIndex).getAssignmentsCount()];

        double shiftMatchIndex;
        for (int i = 0; i < this.matchAssignments.get(sscIndex).getAssignmentsCount(); i++) {
            if (this.matchAssignments.get(sscIndex).getAtomIndex(0, i) == -1) {
                deviations[i] = null;
            } else {
                shiftMatchIndex = ssc.getSubstructure().getAtom(this.matchAssignments.get(sscIndex).getAtomIndex(0, i)).getProperty(Utils.getNMRShiftConstant(ssc.getAtomType()));
                deviations[i] = Math.abs(shiftMatchIndex - querySpectrum.getShift(i, 0));
            }
        }

        return deviations;
    }

    /**
     * Returns deviatons between matched shifts in SSC and query query spectrum.
     * The matching procedure is already included here.
     *
     * @param ssc
     * @param querySpectrum
     * @param tol
     * @return
     *
     * @see #findMatches(model.SSC, casekit.NMR.model.Spectrum, double) 
     */
    public Double[] getDeviations(final SSC ssc, final Spectrum querySpectrum, final double tol) {
        final Double[] deviations = new Double[querySpectrum.getSignalCount()];
        final Assignment matchAssignmentsSSC = this.findMatches(ssc, querySpectrum, tol);

        double shiftMatchIndex;
        for (int i = 0; i < querySpectrum.getSignalCount(); i++) {
            if (matchAssignmentsSSC.getAtomIndex(0, i) == -1) {
                deviations[i] = null;
            } else {
                shiftMatchIndex = ssc.getSubstructure().getAtom(matchAssignmentsSSC.getAtomIndex(0, i)).getProperty(Utils.getNMRShiftConstant(ssc.getAtomType()));
                deviations[i] = Math.abs(shiftMatchIndex - querySpectrum.getShift(i, 0));
            }
        }

        return deviations;
    }

    /**
     * Returns the average of all deviations of matched shifts between a SSC
     * and a query spectrum.
     * The calculation of deviations is already included here.
     *
     * @param ssc
     * @param querySpectrum Query spectrum
     * @param tol Tolerance value [ppm] used during shift matching
     * @param minOverlap Minimum overlap threshold
     * @return
     *
     * @see #getDeviations(model.SSC, casekit.NMR.model.Spectrum, double) 
     */
    public Double getMatchFactor(final SSC ssc, final Spectrum querySpectrum, final double tol, final int minOverlap) {
        return this.getMatchFactor(this.getDeviations(ssc, querySpectrum, tol), minOverlap);
    }

    /**
     * Returns the average of all deviations of matched shifts between a SSC
     * and a query spectrum.
     * If the minimum overlap threshold is not reached, a null value will be
     * returned.
     *
     * @param deviations Deviations between a SSC subspectrum and query
     * spectrum.
     * @param minOverlap Minimum overlap threshold
     * @return
     *
     * @see #getDeviations(model.SSC, casekit.NMR.model.Spectrum, double) 
     */
    public Double getMatchFactor(final Double[] deviations, final int minOverlap) {
        int hitCounter = 0;
        for (final Double deviation : deviations) {
            if (deviation != null) {
                hitCounter++;
            }
        }

        return (hitCounter >= minOverlap) ? Utils.getMean(deviations) : null;
    }
    
    public Spectrum predictSpectrum(final SSC ssc){
        final Spectrum predictedSpectrum = new Spectrum(ssc.getSubspectrum().getNuclei());
        for (final int i: ssc.getAtomTypeIndices().get(ssc.getAtomType())) {
            if(this.getHOSECodeLookupTable().containsKey(ssc.getHOSECode(i))){
                predictedSpectrum.addSignal(new Signal(ssc.getSubspectrum().getNuclei(), new Double[]{this.getHOSECodeLookupTableRMS().get(ssc.getHOSECode(i))}));
            }
        }
        
        return predictedSpectrum;
    }
}
