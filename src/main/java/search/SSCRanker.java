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
package search;

import analysis.MultiplicitySectionsBuilder;
import casekit.NMR.Utils;
import casekit.NMR.match.Matcher;
import casekit.NMR.model.Assignment;
import casekit.NMR.model.Spectrum;
import model.SSC;
import model.SSCLibrary;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IMolecularFormula;
import start.Start;
import utils.Compare;
import utils.ParallelTasks;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Class to search for matches between a SSC library and query spectra.
 * 
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public final class SSCRanker {

    private final SSCLibrary sscLibrary;
    private int nThreads;
    private final ConcurrentHashMap<Long, Object[]> hits;
    private final ArrayList<Long> rankedSSCIndices;
    private final ArrayList<SSC> rankedSSCList;
    private final static MultiplicitySectionsBuilder MULTIPLICITY_SECTIONS_BUILDER = new MultiplicitySectionsBuilder();
    private final static int NO_OF_CALCULATIONS = 3, ASSIGNMENT_IDX = 0, MATCHFACTOR_IDX = 1, TANIMOTO_COEFFICIENT_IDX = 2;
    
    /**
     * Instantiates a new object of this class.
     * The number of threads to use for parallelization is set to the number in {@code sscLibrary} ({@link SSCLibrary#getNThreads()}).
     * 
     *
     * @param sscLibrary HashMap object consisting of SSC and their indices.
     */
    public SSCRanker(final SSCLibrary sscLibrary) {
        this(sscLibrary, sscLibrary.getNThreads());
    }
    
    /**
     * Instantiates a new object of this class.
     *
     * @param sscLibrary HashMap object consisting of SSC and their indices.
     * @param nThreads number of threads to use for matching and ranking of 
     * SSC
     */
    public SSCRanker(final SSCLibrary sscLibrary, final int nThreads){
        this.sscLibrary = sscLibrary;
        this.setNThreads(nThreads);
        this.hits = new ConcurrentHashMap<>();
        this.rankedSSCIndices = new ArrayList<>();
        this.rankedSSCList = new ArrayList<>();
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
     * Returns the matched SSC indices in this SSC library in ranked order.
     * To create/update these coefficients use the findHits function.
     *
     * @return
     *
     * @see #findHits(Spectrum, double, IMolecularFormula)
     */
    public ArrayList<Long> getRankedSSCIndices() {
        return this.rankedSSCIndices;
    }

    /**
     * Returns the match factors of SSC in this SSC library in ranked order.
     * To create/update these match factors use the findHits function.
     *
     * @return
     *
     * @see #findHits(Spectrum, double, IMolecularFormula) 
     */
    public LinkedHashMap<Long, Double> getRankedMatchFactors(){
        final LinkedHashMap<Long, Double> matchFactors = new LinkedHashMap<>();
        for (final Long sscIndex : this.rankedSSCIndices){
            matchFactors.put(sscIndex, (Double) this.hits.get(sscIndex)[SSCRanker.MATCHFACTOR_IDX]);
        }

        return matchFactors;
    }

    /**
     * Returns the tanimoto coefficients of SSC in this SSC library in ranked order.
     * To create/update these coefficients use the findHits function.
     *
     * @return
     *
     * @see #findHits(Spectrum, double, IMolecularFormula) 
     */
    public LinkedHashMap<Long, Float> getRankedTanimotoCoefficients() {
        final LinkedHashMap<Long, Float> tanimotoCoefficients = new LinkedHashMap<>();
        for (final Long sscIndex :this.rankedSSCIndices){
            tanimotoCoefficients.put(sscIndex, (Float) this.hits.get(sscIndex)[SSCRanker.TANIMOTO_COEFFICIENT_IDX]);
        }

        return tanimotoCoefficients;
    }
    
    /**
     * Returns a SSC list containing clones of the current matching SSC
     * for a query spectrum in ranked order.
     * The SSC indices are set in new order starting at 0, 1, 2 etc. .
     * To create/update such ranked SSC library use the findHits function.
     *
     * @return
     * 
     * @see #findHits(Spectrum, double, IMolecularFormula) 
     */
    public ArrayList<SSC> getHits() {
        return this.rankedSSCList;
    }
    
    public long getHitsCount(){
        return this.getHits().size();
    }
                    
    /**
     * Searches for hits between the SSC in the class SSC library and
     * the query spectrum within a certain shift tolerance window. <br>
     * 1. the number of set assignments (matched signals) (highest) <br>
     * 2. the Tanimoto coefficient regarding the query spectrum (highest) <br>
     * 3. the match factor regarding the query spectrum (lowest) <br>
     * 4. the total substructure size (highest) <br> <br>
     * The ranked results are then available in further class functions, see {@code @see}.
     *
     * @param querySpectrum Query spectrum
     * @param shiftTol Tolerance value [ppm] for shift matching
     *
     * @throws CDKException
     * @throws InterruptedException
     *
     * @see #getHits()
     * @see #getRankedMatchFactors()
     * @see #getRankedTanimotoCoefficients()
     *
     */
    public void findHits(final Spectrum querySpectrum, final double shiftTol, final IMolecularFormula molecularFormula) throws Exception {
        this.filter(querySpectrum, shiftTol, molecularFormula);
        this.rank();
        this.buildRankedSSCList();
    }

    private void filter(final Spectrum querySpectrum, final double shiftTol, final IMolecularFormula molecularFormula) throws InterruptedException, CDKException {
        this.hits.clear();

        final HashMap<String, ArrayList<Integer>> multiplicitySectionsQuerySpectrum = SSCRanker.MULTIPLICITY_SECTIONS_BUILDER.buildMultiplicitySections(querySpectrum);
        final ArrayList<Callable<HashMap<Long, Object[]>>> callables = new ArrayList<>();
        // add all task to do
        for (final long sscIndex : this.sscLibrary.getSSCIndices()) {
            final SSC ssc = this.sscLibrary.getSSC(sscIndex);
            callables.add(() -> {
                // 1. pre-search if a molecular formula is given
                if(molecularFormula != null){
                    if(!Compare.compareWithMolecularFormula(ssc.getSubstructure(), molecularFormula)){
                        return null;
                    }
                }
                // 2. pre-search: quick check for too much signals in at least one of the existing multiplicities;
                // could save time in comparison to Matcher.matchSpectra() (below) which has to search for matches (combinations)
                // @TODO might be to extend via comparisons in single sections of each multiplicity
                // @TODO include that pre-search into casekit: Matcher.matchSpectra (?)
                final HashMap<String, ArrayList<Integer>> multiplicitySectionsSSCSubspectrum = ssc.getMultiplicitySections();
                if((multiplicitySectionsSSCSubspectrum.get("Q").size() > multiplicitySectionsQuerySpectrum.get("Q").size())
                        || (multiplicitySectionsSSCSubspectrum.get("T").size() > multiplicitySectionsQuerySpectrum.get("T").size())
                        || (multiplicitySectionsSSCSubspectrum.get("D").size() > multiplicitySectionsQuerySpectrum.get("D").size())
                        || (multiplicitySectionsSSCSubspectrum.get("S").size() > multiplicitySectionsQuerySpectrum.get("S").size())){
                    return null;
                }
                // calculation of matches and check whether each signal in SSC subspectrum has exactly one match in query spectrum
                final Assignment matchAssignment = Matcher.matchSpectra(ssc.getSubspectrum(), querySpectrum, 0, 0, shiftTol);
                if (!matchAssignment.isFullyAssigned(0)) {
                    return null;
                }
                // collect matched signals from query spectrum for Tanimoto coefficient calculation
                final Spectrum matchedQuerySubspectrum = new Spectrum(querySpectrum.getNuclei());
                for (final int signalIndexInQuerySpectrum : matchAssignment.getAssignments(0)) {
                    matchedQuerySubspectrum.addSignal(querySpectrum.getSignal(signalIndexInQuerySpectrum));
                }
                // further calculations: match factor (average deviation) and Tanimoto coefficient
                final HashMap<Long, Object[]> tempHashMap = new HashMap<>();
                final Object[] calculations = new Object[SSCRanker.NO_OF_CALCULATIONS];
                calculations[SSCRanker.ASSIGNMENT_IDX] = matchAssignment;
                calculations[SSCRanker.MATCHFACTOR_IDX] = Utils.roundDouble(Matcher.calculateAverageDeviation(ssc.getSubspectrum(), querySpectrum, 0, 0, shiftTol), Start.DECIMAL_PLACES);
                calculations[SSCRanker.TANIMOTO_COEFFICIENT_IDX] = Matcher.calculateTanimotoCoefficient(matchedQuerySubspectrum, ssc.getSubspectrum(), 0, 0);
                tempHashMap.put(sscIndex, calculations);

                return tempHashMap;
            });
        }
        ParallelTasks.processTasks(callables, tempHashMap -> {
            if ((tempHashMap != null) && !tempHashMap.containsValue(null)) {
                this.hits.putAll(tempHashMap);
            }
        }, this.nThreads);
    }

    /**
     * Ranks the matched hits (SSC indices) pairwise according to the following: <br>
     * 1. the number of set assignments (matched signals) (highest) <br>
     * 2. the Tanimoto coefficient regarding the query spectrum (highest) <br>
     * 3. the match factor regarding the query spectrum (lowest) <br>
     * 4. the total substructure size (highest)
     */
    private void rank() {
        this.rankedSSCIndices.clear();
        // use indices of SSC of the input SSC library which are valid hits (not null)
        this.rankedSSCIndices.addAll(this.hits.keySet());

        this.rankedSSCIndices.sort((indexSSC1, indexSSC2) -> {
            // ranking by number of overlapping signals
            final int setAssignmentsCountComp = -1 * Integer.compare(
                    ((Assignment) hits.get(indexSSC1)[SSCRanker.ASSIGNMENT_IDX]).getSetAssignmentsCount(0),
                    ((Assignment) hits.get(indexSSC2)[SSCRanker.ASSIGNMENT_IDX]).getSetAssignmentsCount(0));
            if (setAssignmentsCountComp != 0) {
                return setAssignmentsCountComp;
            }
            // ranking by tanimoto coefficient
            final int tanimotoCoefficientComp = -1 * Float.compare(
                    ((Float) hits.get(indexSSC1)[SSCRanker.TANIMOTO_COEFFICIENT_IDX]),
                    ((Float) hits.get(indexSSC2)[SSCRanker.TANIMOTO_COEFFICIENT_IDX]));
            if (tanimotoCoefficientComp != 0) {
                return tanimotoCoefficientComp;
            }
            // ranking by match factor if available
            if((hits.get(indexSSC1)[SSCRanker.MATCHFACTOR_IDX] != null) && (hits.get(indexSSC2)[SSCRanker.MATCHFACTOR_IDX] != null)){
                final int matchFactorComp = Double.compare(
                        ((Double) hits.get(indexSSC1)[SSCRanker.MATCHFACTOR_IDX]),
                        ((Double) hits.get(indexSSC2)[SSCRanker.MATCHFACTOR_IDX]));
                if (matchFactorComp != 0) {
                    return matchFactorComp;
                }
            }
            // ranking by total substructure size
            final int substructureSizeComp = -1 * Integer.compare(
                    sscLibrary.getSSC(indexSSC1).getAtomCount(),
                    sscLibrary.getSSC(indexSSC2).getAtomCount());
            if(substructureSizeComp != 0){
                return substructureSizeComp;
            }

            return 0;
        });
    }

    private void buildRankedSSCList() throws Exception {
        this.rankedSSCList.clear();
        for (final long rankedSSCIndex : this.rankedSSCIndices) {
            this.rankedSSCList.add(this.getSSCLibrary().getSSC(rankedSSCIndex).getClone(true));
        }
    }

}
