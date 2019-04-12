/*
 * The MIT License
 *
 * Copyright 2019 Michael Wenk [https://github.com/michaelwenk].
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
package start;

import assembly.Assembly;
import casekit.NMR.DB;
import casekit.NMR.Utils;
import casekit.NMR.model.Spectrum;
import com.google.gson.Gson;
import com.mongodb.MongoClient;
import com.mongodb.client.FindIterable;
import com.mongodb.client.MongoCollection;
import com.mongodb.client.model.Filters;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.logging.Level;
import java.util.logging.Logger;
import match.Match;
import model.SSC;
import model.SSCLibrary;
import org.bson.Document;
import org.bson.conversions.Bson;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import search.MultiplicitySectionsBuilder;
import search.SSCRanker;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class ProcessQueriesMongoDB {
    
    private final SSCLibrary sscLibrary;
    private final String pathToQueriesFile, pathToOutputsFolder, mongoUser, mongoPassword, mongoAuthDB, mongoDBName, mongoDBCollection;
    private final int nThreads, nStarts, minMatchingSphere;
    private final SmilesParser smilesParser;
    private final BufferedReader br;
    private MongoClient mongo;
    private MongoCollection<Document> collection;
    private final TimeMeasurement tm;
    private final double shiftTol, matchFactorThrs;

    
    public ProcessQueriesMongoDB(final String pathToQueriesFile, final String pathToOutputsFolder, 
            final String mongoUser, final String mongoPassword, final String mongoAuthDB, final String mongoDBName, final String mongoDBCollection,
            final int nThreads, final int nStarts, final double shiftTol, final double matchFactorThrs, final int minMatchingSphere) throws FileNotFoundException, CDKException {
        
        this.shiftTol = shiftTol;
        this.matchFactorThrs = matchFactorThrs;
        this.minMatchingSphere = minMatchingSphere;
        this.pathToQueriesFile = pathToQueriesFile;
        this.pathToOutputsFolder = pathToOutputsFolder;
        this.mongoUser = mongoUser;
        this.mongoPassword = mongoPassword;
        this.mongoAuthDB = mongoAuthDB;
        this.mongoDBName = mongoDBName;
        this.mongoDBCollection = mongoDBCollection;
        this.nThreads = nThreads;
        this.sscLibrary = new SSCLibrary(nThreads);
        this.nStarts = nStarts;
        this.smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        this.br = new BufferedReader(new FileReader(this.pathToQueriesFile));
        this.tm = new TimeMeasurement();        
    }

    public void process() throws FileNotFoundException, InterruptedException, IOException, InvalidSmilesException, CDKException, CloneNotSupportedException {                                     
        
        long nStartSSCs;
//        HashMap<String, SSC> solutions;
        SSCRanker sscRanker = new SSCRanker(this.sscLibrary, this.nThreads);
        final MultiplicitySectionsBuilder multiplicitySectionsBuilder = new MultiplicitySectionsBuilder();
        HashMap<String, ArrayList<Integer>> multiplicitySectionsQuerySpectrum;
        String filterExprMultSectionQSize, filterExprMultSectionTSize, filterExprMultSectionDSize, filterExprMultSectionSSize;
        final Bson[] filters = new Bson[4];

        IAtomContainer solutionAtomContainer;
        int querySpectrumCounter = 0;        
        BufferedWriter bw;
        final HashMap<String, Float> tanimotoCoefficients = new HashMap<>();
        ArrayList<Double> sortedQuerySpectrumShifts, sortedSolutionSpectrumShifts;

        this.mongo = DB.login(this.mongoUser, this.mongoPassword, this.mongoAuthDB);
        this.collection = DB.getCollection(this.mongo, mongoDBName, mongoDBCollection);
        final Gson gson = new Gson();
        final List<String> keys = new ArrayList<>();
        keys.add("subspectrum");
        final long collectionSize = this.collection.countDocuments();
        final int rangesCount = this.nThreads;
        final long rangesSize = collectionSize / rangesCount;
        
        System.out.println("\n\n-> processing query file: \"" + this.pathToQueriesFile + "\" ...");

        Iterator<String> it = this.br.lines().iterator();
        String line;
        while (it.hasNext()) {
            line = it.next();
            if (line.startsWith("//")) {
                continue;
            }
            if (line.startsWith("#")) {
                System.out.println("\n\nskip query: " + querySpectrumCounter + " -> " + line + "\n");
                querySpectrumCounter++;
                continue;
            }
            System.out.println("\n\nnow processing query: " + querySpectrumCounter + " -> " + line + "\n");

            final Spectrum querySpectrum = DB.NMRShiftDBSpectrumToSpectrum(line, "13C");
            System.out.println("\nquery spectrum:\t" + querySpectrum.getShifts(0)
                    + "\nequivalents:\t" + querySpectrum.getEquivalences()
                    + "\nmultiplicities:\t" + querySpectrum.getMultiplicities()
                    + "\nequivalent signals classes: " + querySpectrum.getEquivalentSignalClasses());
            System.out.println("multiplicities and shift sections:");
            multiplicitySectionsQuerySpectrum = multiplicitySectionsBuilder.getMultiplicitySections(querySpectrum);
            for (final String mult : multiplicitySectionsQuerySpectrum.keySet()) {
                System.out.println("-> " + mult + ":\t" + multiplicitySectionsQuerySpectrum.get(mult));
            }

            filterExprMultSectionQSize = "{$lte: [{$size: \"$multSections.Q\"}, " + multiplicitySectionsQuerySpectrum.get("Q").size() + "]}";
            filterExprMultSectionTSize = "{$lte: [{$size: \"$multSections.T\"}, " + multiplicitySectionsQuerySpectrum.get("T").size() + "]}";
            filterExprMultSectionDSize = "{$lte: [{$size: \"$multSections.D\"}, " + multiplicitySectionsQuerySpectrum.get("D").size() + "]}";
            filterExprMultSectionSSize = "{$lte: [{$size: \"$multSections.S\"}, " + multiplicitySectionsQuerySpectrum.get("S").size() + "]}";
            filters[0] = Filters.expr(Document.parse(filterExprMultSectionQSize));
            filters[1] = Filters.expr(Document.parse(filterExprMultSectionTSize));
            filters[2] = Filters.expr(Document.parse(filterExprMultSectionDSize));
            filters[3] = Filters.expr(Document.parse(filterExprMultSectionSSize));
            
            // empty SSC library
            this.sscLibrary.removeAll();
            // initialize an executor for parallelization
            final ExecutorService executor = Utils.initExecuter(this.nThreads);
            final ArrayList<Callable<ArrayList<Document>>> callables = new ArrayList<>();
            
            System.out.println("\nextending from pre-search result...");
            this.tm.start();
            // add all task to do
            for (int n = 0; n < rangesCount; n++) {
                // specify ranges to use for pre-searches
                final Bson filterRange;
                if(n == (rangesCount - 1)){
                    filterRange = Document.parse("{$and: [{\"index\": {$gte: " + (n * rangesSize) + "}}, {\"index\": {$lte: " + collectionSize + "}}]}");
//                    System.out.println("n: " + n + " -> begin: " + (n * rangeSizes) + ", end: " + this.collection.countDocuments() + " -> " + this.collection.countDocuments());
                } else {
                    filterRange = Document.parse("{$and: [{\"index\": {$gte: " + (n * rangesSize) + "}}, {\"index\": {$lte: " + (((n + 1) * rangesSize) - 1) + "}}]}");
//                    System.out.println("n: " + n + " -> begin: " + (n * rangeSizes) + ", end: " + (((n + 1) * rangeSizes) - 1) + " -> " + this.collection.countDocuments());
                }
                
                callables.add((Callable<ArrayList<Document>>) () -> {
//                    System.out.println("pre-search 1...");
//                    this.tm.start();
                    FindIterable<Document> resultPreSearch1 = this.collection.find(Filters.and(filters)).filter(filterRange);
                    ArrayList<Document> resultPreSearch2 = new ArrayList<>();                    
//                    System.out.println("pre-search 1 done!!!");
//                    this.tm.stop();
//                    System.out.println("--> time needed: " + this.tm.getResult() + " s");

                    
//                    System.out.println("pre-search 2...");
//                    this.tm.start();
                    // bottleneck: the iteration over result (FindIterable object) is very slow, even for small datasets 
                    // -> split whole collection over all given threads
                    Spectrum subspectrum;
                    Double matchFactor;
                    for (final Document document : resultPreSearch1) {   
                        subspectrum = gson.fromJson(document.getEmbedded(keys, Document.class).toJson(), Spectrum.class);
                        matchFactor = Match.calculateMatchFactor(subspectrum, querySpectrum, this.shiftTol);                        
                        if ((matchFactor != null) && (matchFactor <= this.matchFactorThrs) && Match.matchSpectra(subspectrum, querySpectrum, this.shiftTol).isFullyAssigned(0)) {
                            resultPreSearch2.add(document);
                        }
                    }
//                    System.out.println("pre-search 2 done!!!");
//                    this.tm.stop();
//                    System.out.println("--> time needed: " + this.tm.getResult() + " s");
                    return resultPreSearch2;
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
                    .forEach((resultPreSearches) -> {
                        try {
                            this.sscLibrary.extend(resultPreSearches);
                        } catch (CDKException | CloneNotSupportedException | InterruptedException ex) {
                            try {
                                throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": insertion of SSCs failed");
                            } catch (CDKException e) {
                                Logger.getLogger(SSCLibrary.class.getName()).log(Level.SEVERE, null, e);
                            }
                        }
                    });
            // shut down the executor service
            Utils.stopExecuter(executor, 5);
            
            System.out.println("extension done!!!");
            this.tm.stop();
            System.out.println("--> time needed: " + this.tm.getResult() + " s");
            System.out.println("-> #SSCs in SSC library: " + this.sscLibrary.getSSCCount());
                       
            sscRanker.rank(querySpectrum, this.shiftTol);
            System.out.println("\n\nno. of matches: " + sscRanker.getRankedSSCsCount());
            System.out.println("match factors: " + sscRanker.getMatchFactors());
            System.out.println("ranked SSC indices: " + sscRanker.getRankedSSCLibrary().getSSCIndices());
            final LinkedHashSet<Double> rankedMatchedFactors = new LinkedHashSet<>();
            for (final long sscIndex : sscRanker.getRankedSSCLibrary().getSSCIndices()) {
                rankedMatchedFactors.add(sscRanker.getMatchFactors().get(sscIndex));
            }
            System.out.println("ranked match factors: " + rankedMatchedFactors);            
            
            
            
            
            
            
            
            
            final SSCLibrary rankedSSCLibrary = sscRanker.getRankedSSCLibrary();                                 
            
//            final IAtomContainer ac = new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles("CC(=CCC1=C(C=C(C(=C1)C(=O)C(CO)C2=CC=C(C=C2)OC)O)OC)C");
//            Utils.setAromaticitiesInAtomContainer(ac);
//            
//            HashMap<String, ArrayList<Long>> uniqueSSCIndices = new HashMap<>();
//            SSC ssc;
//            String smilesTemp;
//            SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Absolute);
//            final SSCLibrary rankedSSCLibraryTemp = sscRanker.getRankedSSCLibrary();
//            final SSCLibrary rankedSSCLibrary = new SSCLibrary(this.nThreads);
//            final ArrayList<Double> rankedMatchFactors = new ArrayList<>();
//            for (final long sscIndex : rankedSSCLibraryTemp.getSSCIndices()) {
//                ssc = rankedSSCLibraryTemp.getSSC(sscIndex).getClone();
//                smilesTemp = smilesGenerator.create(ssc.getSubstructure());
//                Pattern pattern = DfPattern.findSubstructure(ssc.getSubstructure());                
//                if(!pattern.matches(ac)){
//                    continue;
//                }
//                
//                if(!uniqueSSCIndices.containsKey(smilesTemp)){
//                    uniqueSSCIndices.put(smilesTemp, new ArrayList<>());
//                    rankedMatchFactors.add(sscRanker.getMatchFactors().get(sscIndex));
//                    ssc.setIndex(rankedSSCLibrary.getSSCCount());
//                    rankedSSCLibrary.insert(ssc);
////                    System.out.println("matches for ssc " + ssc.getIndex() + " : " + Arrays.toString(pattern.match(ac)));
//                }
//                uniqueSSCIndices.get(smilesTemp).add(sscIndex);
//            }
//            System.out.println("-> unique substructures count in hit list: " + uniqueSSCIndices.size());
//            System.out.println("-> ranked SSC indices: " + rankedSSCLibrary.getSSCIndices());
//            System.out.println("-> ranked match factors: " + rankedMatchFactors);   
//     
////            final LinkedHashSet<Long> wantedSSCIndices = new LinkedHashSet<>();
////            wantedSSCIndices.add(3L);
////            wantedSSCIndices.add(9L);
////            wantedSSCIndices.add(39L);
////            wantedSSCIndices.add(63L);
////            wantedSSCIndices.add(104L);
////            wantedSSCIndices.add(114L);
////            wantedSSCIndices.add(126L);
////            wantedSSCIndices.add(131L);    
////            rankedSSCLibraryTemp.removeAll();
////            for (final Long wantedSSCIndex : wantedSSCIndices) {
////                rankedSSCLibraryTemp.insert(rankedSSCLibrary.getSSC(wantedSSCIndex));
////            }            
////            rankedSSCLibrary.removeAll();
////            rankedSSCLibrary.extend(rankedSSCLibraryTemp);
//            
//            System.out.println("--> ranked SSC indices: " + rankedSSCLibrary.getSSCIndices());

            
            
            





            
            if ((this.nStarts > 0) && (this.nStarts < rankedSSCLibrary.getSSCCount()/*sscRanker.getRankedSSCsCount()*/)) {
                nStartSSCs = this.nStarts;
            } else {
                nStartSSCs = rankedSSCLibrary.getSSCCount();//sscRanker.getRankedSSCsCount();
            }
            System.out.println("\nnumber of start SSCs for query " + querySpectrumCounter + ":\t" + nStartSSCs);

            final HashMap<String, SSC> solutions = Assembly.assemble(nStartSSCs, this.nThreads, rankedSSCLibrary, this.minMatchingSphere, querySpectrum, this.matchFactorThrs, this.shiftTol);

            System.out.println("\nsolutions for query " + querySpectrumCounter + " (" + line + "):\t" + solutions.size());

            final String[] solutionsSMILESToSort = new String[solutions.size()];
            int solutionsCounter = 0;
            for (final String smiles : solutions.keySet()) {
                solutionsSMILESToSort[solutionsCounter] = smiles;
                solutionsCounter++;
            }
            tanimotoCoefficients.clear();
            for (final String smiles : solutionsSMILESToSort) {
                try {
                    tanimotoCoefficients.put(smiles, Match.calculateTanimotoCoefficient(querySpectrum, solutions.get(smiles).getSubspectrum(), 0));
                } catch (CDKException ex) {
                    tanimotoCoefficients.put(smiles, (float) 0.0);
                }
            }
            
            Arrays.parallelSort(solutionsSMILESToSort, new Comparator<String>() {
                @Override
                public int compare(final String smiles1, final String smiles2) {                    
                    return -1 * Float.compare(tanimotoCoefficients.get(smiles1), tanimotoCoefficients.get(smiles2));
                }
            });

            bw = new BufferedWriter(new FileWriter(this.pathToOutputsFolder + "/results_" + querySpectrumCounter + ".smiles"));
            String smiles;
            sortedQuerySpectrumShifts = querySpectrum.getShifts(0);
            Collections.sort(sortedQuerySpectrumShifts);
            for (int i = 0; i < solutionsSMILESToSort.length; i++) {
                smiles = solutionsSMILESToSort[i];
                solutionAtomContainer = this.smilesParser.parseSmiles(smiles);
                sortedSolutionSpectrumShifts = solutions.get(smiles).getSubspectrum().getShifts(0);
                Collections.sort(sortedSolutionSpectrumShifts);

                System.out.println("query spectrum   :\t" + sortedQuerySpectrumShifts);
                System.out.println("solution spectrum:\t" + sortedSolutionSpectrumShifts);
                System.out.println("SMILES:          :\t" + smiles);
                System.out.println("              --> \t" + "atoms: " + solutionAtomContainer.getAtomCount() + ", bonds: " + solutionAtomContainer.getBondCount());
                System.out.println("              --> \t" + "tanimoto: " + tanimotoCoefficients.get(smiles));
                
                
                bw.append(smiles + " " + i + " " + tanimotoCoefficients.get(smiles));
                bw.newLine();
                bw.flush();
            }
            
            
            
            bw.close();
            querySpectrumCounter++;            
        }
        this.br.close();
        DB.logout(this.mongo);
    }
}
