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
import com.mongodb.MongoClientOptions;
import com.mongodb.MongoCredential;
import com.mongodb.ServerAddress;
import com.mongodb.client.FindIterable;
import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoDatabase;
import com.mongodb.client.model.Filters;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
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
    private final int nThreads, nStarts;
    private final SmilesParser smilesParser;
    private final BufferedReader br;
    private MongoClient mongo;
    private MongoCollection<Document> collection;
    private final TimeMeasurement tm;

    
    public ProcessQueriesMongoDB(final String pathToQueriesFile, final String pathToOutputsFolder, final String mongoUser, final String mongoPassword,
            final String mongoAuthDB, final String mongoDBName, final String mongoDBCollection,
            final int nThreads, final int nStarts) throws FileNotFoundException, CDKException {
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

    private void login() throws CDKException {
        try {
            // Creating a Mongo client   
            this.tm.start();
            this.mongo = new MongoClient(
                    new ServerAddress("127.0.0.1", 27017),
                    MongoCredential.createCredential(
                            this.mongoUser,
                            this.mongoAuthDB,
                            this.mongoPassword.toCharArray()),
                    MongoClientOptions.builder().build());
            System.out.println("Login to MongoDB was successfull");
            // Accessing the database 
            final MongoDatabase database = this.mongo.getDatabase(this.mongoDBName);
            System.out.println("Access to database \"" + this.mongoDBName + "\" was successfull");
            // Retrieving a collection
            this.collection = database.getCollection(this.mongoDBCollection);
            System.out.println("Retrieval of collection \"" + this.mongoDBCollection + "\" was successfull -> size: " + this.collection.countDocuments());
            this.tm.stop();
            System.out.println("--> time needed: " + this.tm.getResult() + " s");
        } catch (Exception e) {
            throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": could not connect to database \"" + this.mongoDBName + "\" or collection \"" + this.mongoDBCollection + "\"");
        }
    }
    
    private void logout(){
        this.mongo.close();
    }
    
    public void process() throws FileNotFoundException, InterruptedException, IOException, InvalidSmilesException, CDKException, CloneNotSupportedException {
        int nStartSSCs;
        HashMap<String, SSC> solutions;
        SSCRanker sscRanker = new SSCRanker(this.sscLibrary, this.nThreads);
        final MultiplicitySectionsBuilder multiplicitySectionsBuilder = new MultiplicitySectionsBuilder();
        HashMap<String, ArrayList<Integer>> multiplicitySectionsQuerySpectrum;
        String filterExprMultSectionQSize, filterExprMultSectionTSize, filterExprMultSectionDSize, filterExprMultSectionSSize;
        final Bson[] filters = new Bson[4];

        IAtomContainer solutionAtomContainer;
        ArrayList<Double> shiftsQuerySpectrum, shiftsSolutionSpectrum;
        int querySpectrumCounter = 0;        
        BufferedWriter bw;

        this.login();
        
        final Gson gson = new Gson();
        final List<String> keys = new ArrayList<>();
        keys.add("ssc");
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

            final Spectrum querySpectrum = DB.NMRShiftDBSpectrumToSpectrum(line, Start.SPECTRUM_PROPERTY_ATOMTYPE);
            System.out.println("\nquery spectrum:\t" + querySpectrum.getShifts(0)
                    + "\nequivalents:\t" + querySpectrum.getEquivalences()
                    + "\nmultiplicities:\t" + querySpectrum.getMultiplicities()
                    + "\nequivalent signals classes: " + querySpectrum.getEquivalentSignalClasses());
            System.out.println("multiplicities and shift sections:");
            multiplicitySectionsQuerySpectrum = multiplicitySectionsBuilder.getMultiplicitySections(querySpectrum);
            for (final String mult : multiplicitySectionsQuerySpectrum.keySet()) {
                System.out.println("-> " + mult + ":\t" + multiplicitySectionsQuerySpectrum.get(mult));
            }

            filterExprMultSectionQSize = "{$lte: [{$size: \"$ssc.multSections.Q\"}, " + multiplicitySectionsQuerySpectrum.get("Q").size() + "]}";
            filterExprMultSectionTSize = "{$lte: [{$size: \"$ssc.multSections.T\"}, " + multiplicitySectionsQuerySpectrum.get("T").size() + "]}";
            filterExprMultSectionDSize = "{$lte: [{$size: \"$ssc.multSections.D\"}, " + multiplicitySectionsQuerySpectrum.get("D").size() + "]}";
            filterExprMultSectionSSize = "{$lte: [{$size: \"$ssc.multSections.S\"}, " + multiplicitySectionsQuerySpectrum.get("S").size() + "]}";
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
                    filterRange = Document.parse("{$and: [{\"_id\": {$gte: " + (n * rangesSize) + "}}, {\"_id\": {$lte: " + collectionSize + "}}]}");
//                    System.out.println("n: " + n + " -> begin: " + (n * rangeSizes) + ", end: " + this.collection.countDocuments() + " -> " + this.collection.countDocuments());
                } else {
                    filterRange = Document.parse("{$and: [{\"_id\": {$gte: " + (n * rangesSize) + "}}, {\"_id\": {$lte: " + (((n + 1) * rangesSize) - 1) + "}}]}");
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
                    for (Document doc : resultPreSearch1) {                        
                        if (Match.matchSpectra(gson.fromJson(doc.getEmbedded(keys, Document.class).toJson(), Spectrum.class), querySpectrum, Start.PICK_PRECISION).isFullyAssigned(0)) {
                            resultPreSearch2.add(doc);
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
                       
            sscRanker.rank(querySpectrum, Start.PICK_PRECISION);
            System.out.println("\n\nno. of matches: " + sscRanker.getMatchFactors().size());
            System.out.println("match factors: " + sscRanker.getMatchFactors());
            System.out.println("ranked SSC indices: " + sscRanker.getRankedSSCIndices());
            System.out.println("ranked match factors: " + sscRanker.getRankedMatchFactors());
            System.out.println("ranked SSC library indices: " + sscRanker.getRankedSSCLibrary().getSSCIndices() + "\n");

            if ((this.nStarts > 0) && (this.nStarts < sscRanker.getMatchFactors().size())) {
                nStartSSCs = this.nStarts;
            } else {
                nStartSSCs = sscRanker.getMatchFactors().size();
            }
            System.out.println("\nnumber of start SSCs for query " + querySpectrumCounter + ":\t" + nStartSSCs);

            solutions = Assembly.assemble(nStartSSCs, this.nThreads, sscRanker.getRankedSSCLibrary(), Start.MIN_MATCHING_SPHERE_COUNT, querySpectrum, Start.SHIFT_TOL, Start.MATCH_FACTOR_THRS, Start.PICK_PRECISION);

            System.out.println("\nsolutions for query " + querySpectrumCounter + " (" + line + "):\t" + solutions.size());

            shiftsQuerySpectrum = querySpectrum.getShifts(0);
            Collections.sort(shiftsQuerySpectrum);

            bw = new BufferedWriter(new FileWriter(this.pathToOutputsFolder + "/results_" + querySpectrumCounter + ".smiles"));
            for (final String smiles : solutions.keySet()) {
                solutionAtomContainer = this.smilesParser.parseSmiles(smiles);

                shiftsSolutionSpectrum = solutions.get(smiles).getSubspectrum().getShifts(0);
                Collections.sort(shiftsSolutionSpectrum);
                System.out.println("query spectrum   :\t" + shiftsQuerySpectrum);
                System.out.println("solution spectrum:\t" + shiftsSolutionSpectrum);
                System.out.println("SMILES:          :\t" + smiles);
                System.out.println("              --> \t" + "atoms: " + solutionAtomContainer.getAtomCount() + ", bonds: " + solutionAtomContainer.getBondCount());

                bw.append(smiles);
                bw.newLine();
                bw.flush();
            }
            bw.close();
            querySpectrumCounter++;
        }
        this.br.close();
        this.logout();
    }
}
