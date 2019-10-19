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

package start;

import analysis.MultiplicitySectionsBuilder;
import assembly.Assembly;
import casekit.NMR.dbservice.MongoDB;
import casekit.NMR.match.Matcher;
import casekit.NMR.model.Signal;
import casekit.NMR.model.Spectrum;
import com.google.gson.Gson;
import com.mongodb.MongoClient;
import com.mongodb.client.FindIterable;
import com.mongodb.client.MongoCollection;
import com.mongodb.client.model.Filters;
import model.SSC;
import model.SSCLibrary;
import org.bson.Document;
import org.bson.conversions.Bson;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import search.SSCRanker;

import java.io.*;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

public class ProcessQueries {

    private SSCLibrary sscLibrary;
    private final String pathToQueriesFile, pathToOutputsFolder;
    private final int nThreads, nStarts, minMatchingSphere;
    private final SmilesParser smilesParser;
    private final BufferedReader br;
    private final TimeMeasurement tm;
    private final double shiftTol, matchFactorThrs;

    private String mongoUser, mongoPassword, mongoAuthDB, mongoDBName, mongoDBCollection;
    private MongoClient mongo;
    private MongoCollection<Document> collection;
    private boolean useMongoDB;


    public ProcessQueries(final SSCLibrary sscLibrary, final String pathToQueriesFile, final String pathToOutputsFolder, final int nThreads, final int nStarts, final double shiftTol, final double matchFactorThrs, final int minMatchingSphere) throws FileNotFoundException {
        this.pathToQueriesFile = pathToQueriesFile;
        this.pathToOutputsFolder = pathToOutputsFolder;
        this.nThreads = nThreads;
        this.sscLibrary = sscLibrary;
        this.sscLibrary.setNThreads(this.nThreads);
        this.nStarts = nStarts;
        this.shiftTol = shiftTol;
        this.matchFactorThrs = matchFactorThrs;
        this.minMatchingSphere = minMatchingSphere;

        this.smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        this.br = new BufferedReader(new FileReader(this.pathToQueriesFile));
        this.tm = new TimeMeasurement();

        this.useMongoDB = false;
    }


    public void initMongoDBProcessing(final String mongoUser, final String mongoPassword, final String mongoAuthDB, final String mongoDBName, final String mongoDBCollection) throws CDKException {
        this.mongoUser = mongoUser;
        this.mongoPassword = mongoPassword;
        this.mongoAuthDB = mongoAuthDB;
        this.mongoDBName = mongoDBName;
        this.mongoDBCollection = mongoDBCollection;

        this.mongo = MongoDB.login(this.mongoUser, this.mongoPassword, this.mongoAuthDB);
        this.collection = MongoDB.getCollection(this.mongo, mongoDBName, mongoDBCollection);

        this.sscLibrary = new SSCLibrary(this.nThreads);
        this.useMongoDB = true;
    }

    public void process() throws Exception {

        System.out.println("\n\n-> processing query file: \"" + this.pathToQueriesFile + "\" ...");

        final SSCRanker sscRanker = new SSCRanker(this.sscLibrary);

        int querySpectrumCounter = 0;
        Spectrum querySpectrum = new Spectrum(new String[]{Start.SIGNAL_NUCLEUS});
        final Iterator<String> it = this.br.lines().iterator();
        String line;
        String[] signalProperties;
        // process every query spectrum in queries file
        while (it.hasNext()) {
            line = it.next();
            if (line.trim().startsWith("//")) {
                // create new query spectrum with description
                querySpectrum = new Spectrum(new String[]{Start.SIGNAL_NUCLEUS});
                querySpectrum.setSpecDescription(line.split("//")[1].trim());
                while (it.hasNext()){
                    line = it.next();
                    if(line.trim().isEmpty()){
                        break;
                    }
                    signalProperties = line.trim().split(",");
                    if(!signalProperties[0].trim().equals(Start.SIGNAL_NUCLEUS)){
                        System.out.println("Query spectrum " + querySpectrumCounter + " contains a signal with different nucleus!!!");
                        continue;
                    }
                    querySpectrum.addSignal(new Signal(querySpectrum.getNuclei(), new Double[]{Double.parseDouble(signalProperties[1].trim())}, signalProperties[2].trim(), Double.parseDouble(signalProperties[3].trim())));
                }
                querySpectrum.detectEquivalences();
            }
            System.out.println("\n\nnow processing query: " + querySpectrumCounter + " -> " + querySpectrum.getSpecDescription() + "\n");
            System.out.println("\nquery spectrum: " + querySpectrum.getShifts(0)
                    + "\nmultiplicities: " + querySpectrum.getMultiplicities()
                    + "\nintensities   : " + querySpectrum.getIntensities()
                    + "\nequivalents   : " + querySpectrum.getEquivalences()
                    + "\nequivalent signals classes: " + querySpectrum.getEquivalentSignalClasses());

            // if MongoDB is used then a fast presearch can be done to reduce the amount of SSC to check
            if(this.useMongoDB){
                this.presearch(querySpectrum);
            }
            // main processing of the query spectrum against the SSC library
            this.processCore(sscRanker, querySpectrum, querySpectrumCounter);

            querySpectrumCounter++;
        }
        this.br.close();
    }

    public void processCore(final SSCRanker sscRanker, final Spectrum querySpectrum, final long querySpectrumCounter) throws Exception {

        sscRanker.findHits(querySpectrum, this.shiftTol);
        System.out.println("\n\nno. of matches:    " + sscRanker.getHitsCount());
        System.out.println("ranked SSC indices:    " + sscRanker.getRankedSSCIndices());
        System.out.println("ranked match factors:  " + sscRanker.getRankedMatchFactors());
        System.out.println("ranked match tanimoto: " + sscRanker.getRankedTanimotoCoefficients() + "\n");

        final ArrayList<SSC> rankedSSCList = sscRanker.getHits();

        long nStartSSCs;
        if ((this.nStarts > 0) && (this.nStarts < rankedSSCList.size())) {
            nStartSSCs = this.nStarts;
        } else {
            nStartSSCs = rankedSSCList.size();
        }
        System.out.println("\nnumber of start SSC for query " + querySpectrumCounter + ":\t" + nStartSSCs);

        final ConcurrentHashMap<String, SSC> solutions = Assembly.assemble(nStartSSCs, sscRanker.getNThreads(), rankedSSCList, this.minMatchingSphere, querySpectrum, this.matchFactorThrs, this.shiftTol, this.pathToOutputsFolder, querySpectrumCounter);

        System.out.println("\nsolutions for query " + querySpectrumCounter + " (" + querySpectrum.getSpecDescription() + "):\t" + solutions.size());

        final String[] solutionsSMILESToSort = new String[solutions.size()];
        int solutionsCounter = 0;
        for (final String smiles : solutions.keySet()) {
            solutionsSMILESToSort[solutionsCounter] = smiles;
            solutionsCounter++;
        }
        final HashMap<String, Float> tanimotoCoefficients = new HashMap<>();
        for (final String smiles : solutionsSMILESToSort) {
            try {
                tanimotoCoefficients.put(smiles, Matcher.calculateTanimotoCoefficient(querySpectrum, solutions.get(smiles).getSubspectrum(), 0, 0));
            } catch (CDKException ex) {
                tanimotoCoefficients.put(smiles, (float) 0.0);
            }
        }

        Arrays.parallelSort(solutionsSMILESToSort, (smiles1, smiles2) -> -1 * Float.compare(tanimotoCoefficients.get(smiles1), tanimotoCoefficients.get(smiles2)));

        final BufferedWriter bw = new BufferedWriter(new FileWriter(this.pathToOutputsFolder + "/results_" + querySpectrumCounter + ".smiles"));
        final ArrayList<Double> sortedQuerySpectrumShifts = querySpectrum.getShifts(0);
        Collections.sort(sortedQuerySpectrumShifts);

        String smiles;
        ArrayList<Double> sortedSolutionSpectrumShifts;
        for (int i = 0; i < solutionsSMILESToSort.length; i++) {
            smiles = solutionsSMILESToSort[i];
            sortedSolutionSpectrumShifts = new ArrayList<>(solutions.get(smiles).getSubspectrum().getShifts(0));
            Collections.sort(sortedSolutionSpectrumShifts);
            System.out.println("query spectrum   :\t" + sortedQuerySpectrumShifts);
            System.out.println("solution spectrum:\t" + sortedSolutionSpectrumShifts);
            System.out.println("SMILES:          :\t" + smiles);
            System.out.println("              --> \t" + "atoms: " + solutions.get(smiles).getSubstructure().getAtomCount() + ", bonds: " + solutions.get(smiles).getSubstructure().getBondCount());
            System.out.println("              --> \t" + "tanimoto: " + tanimotoCoefficients.get(smiles));


            bw.append(smiles + " " + i + " " + tanimotoCoefficients.get(smiles));
            bw.newLine();
            bw.flush();
        }

        bw.close();
    }

    private void presearch(final Spectrum querySpectrum) throws CDKException {

        final MultiplicitySectionsBuilder multiplicitySectionsBuilder = new MultiplicitySectionsBuilder();
        HashMap<String, ArrayList<Integer>> multiplicitySectionsQuerySpectrum;
        String filterExprMultSectionQSize, filterExprMultSectionTSize, filterExprMultSectionDSize, filterExprMultSectionSSize;
        final Bson[] filters = new Bson[4];

        final Gson gson = new Gson();
        final List<String> keys = new ArrayList<>();
        keys.add("subspectrum");
        final long collectionSize = this.collection.countDocuments();
        final int rangesCount = this.nThreads;
        final long rangesSize = collectionSize / rangesCount;

        System.out.println("multiplicities and shift sections:");
        multiplicitySectionsQuerySpectrum = multiplicitySectionsBuilder.buildMultiplicitySections(querySpectrum);
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

        final ArrayList<Document> resultPresearch = new ArrayList<>();
        FindIterable<Document> resultPresearch1;
        ArrayList<Document> resultPresearch2;
        Spectrum subspectrum;
        Double matchFactor;
        Bson filterRange;
        System.out.println("\nextending from pre-search result...");
        this.tm.start();
        for (int n = 0; n < rangesCount; n++) {
            // specify ranges to use for pre-searches
            if(n == (rangesCount - 1)){
                filterRange = Document.parse("{$and: [{\"index\": {$gte: " + (n * rangesSize) + "}}, {\"index\": {$lte: " + collectionSize + "}}]}");
//                    System.out.println("n: " + n + " -> begin: " + (n * rangeSizes) + ", end: " + this.collection.countDocuments() + " -> " + this.collection.countDocuments());
            } else {
                filterRange = Document.parse("{$and: [{\"index\": {$gte: " + (n * rangesSize) + "}}, {\"index\": {$lte: " + (((n + 1) * rangesSize) - 1) + "}}]}");
//                    System.out.println("n: " + n + " -> begin: " + (n * rangeSizes) + ", end: " + (((n + 1) * rangeSizes) - 1) + " -> " + this.collection.countDocuments());
            }
            resultPresearch1 = this.collection.find(Filters.and(filters)).filter(filterRange);
            resultPresearch2 = new ArrayList<>();
            // bottleneck: the iteration over result (FindIterable object) is very slow, even for small datasets
            // -> split whole collection over all given threads
            for (final Document document : resultPresearch1) {
                subspectrum = gson.fromJson(document.getEmbedded(keys, Document.class).toJson(), Spectrum.class);
                matchFactor = Matcher.calculateAverageDeviation(subspectrum, querySpectrum, 0, 0, this.shiftTol);
                if ((matchFactor != null) && (matchFactor <= this.matchFactorThrs) && Matcher.matchSpectra(subspectrum, querySpectrum, 0, 0, this.shiftTol).isFullyAssigned(0)) {
                    resultPresearch2.add(document);
                }
            }
            resultPresearch.addAll(resultPresearch2);
        }
        this.sscLibrary.extend(resultPresearch);
        System.out.println("extension done!!!");
        this.tm.stop();
        System.out.println("--> time needed: " + this.tm.getResult() + " s");
        System.out.println("-> #SSC in SSC library: " + this.sscLibrary.getSSCCount());
    }
}
