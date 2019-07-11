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
package analysis;

import casekit.NMR.dbservice.MongoDB;
import casekit.NMR.Utils;
import casekit.NMR.dbservice.NMRShiftDB;
import com.mongodb.MongoClient;
import com.mongodb.client.MongoCollection;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.bson.Document;
import org.openscience.cdk.exception.CDKException;

/**
 * Class to build a solvent deviations lookup table and to store it into MongoDB collection.
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class SolventDeviationsLookupTableBuilder {
    
    public static boolean buildSolventDeviationsLookupTable(final String pathToNMRShiftDB, final String nucleus, final int nThreads,
            final String mongoUser, final String mongoPassword, final String mongoAuthDB, final String mongoDBName, final String mongoDBCollection) throws InterruptedException {

        try {
            final MongoClient mongo = MongoDB.login(mongoUser, mongoPassword, mongoAuthDB);
            final MongoCollection<Document> collection = MongoDB.getCollection(mongo, mongoDBName, mongoDBCollection);
            collection.drop();
            System.out.println("Build and export lookup table...");
            final HashMap<String, ArrayList<Double>> deviationsMap = NMRShiftDB.getSolventDeviations(pathToNMRShiftDB, nucleus);
            
            
            // initialize an executor for parallelization
            final ExecutorService executor = Utils.initExecuter(nThreads);
            final ArrayList<Callable<Document>> callables = new ArrayList<>();
            // add all task to do
            for (final String combiKey : deviationsMap.keySet()) {
                callables.add((Callable<Document>) () -> {
                    final Document document = new Document();
                    final ArrayList<Double> deviations = deviationsMap.get(combiKey);
                    document.append("solvent1", combiKey.split("_")[0]);
                    document.append("solvent2", combiKey.split("_")[1]);
                    document.append("count", deviations.size());
                    document.append("min", (deviations.size() > 0) ? Collections.min(deviations) : null);
                    document.append("max", (deviations.size() > 0) ? Collections.max(deviations) : null);
                    document.append("rms", Utils.getRMS(deviations));
                    document.append("mean", Utils.getMean(deviations));
                    document.append("median", Utils.getMedian(deviations));
                    document.append("var", Utils.getVariance(deviations));
                    document.append("sd", Utils.getStandardDeviation(deviations));
                    document.append("deviations", deviations);

                    return document;
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
                    .forEach((document) -> {
                        try {
                            collection.insertOne(document);
                        } catch (Exception e) {
                            System.err.println("export for key \"" + document.get("_id") + "\" failed: " + e.getMessage());
                            System.out.println("-> document: \n" + document.toJson());
                        }
                    });
            // shut down the executor service
            Utils.stopExecuter(executor, 5);
            
            System.out.println("done...");

            MongoDB.logout(mongo);
        } catch (CDKException | FileNotFoundException ex) {
            Logger.getLogger(HOSECodeLookupTableBuilder.class.getName()).log(Level.SEVERE, null, ex);

            return false;
        } 

        return true;
    }

    private static Options setupOptions(String[] args) {
        Options options = new Options();
        Option nmrshiftdbOption = Option.builder("nmrdb")
                .required(true)
                .hasArg()
                .longOpt("nmrshiftdb")
                .desc("Path to NMRShiftDB (SDF).")
                .build();
        options.addOption(nmrshiftdbOption);
        Option spectrumPropertyOption = Option.builder("nu")
                .required(true)
                .hasArg()
                .longOpt("nucleus")
                .desc("Nucleus for spectra to use from NMRShiftDB, e.g. for nucleus \"13C\" use \"Spectrum 13C 0\" and \"Spectrum 13C 1\" etc.")
                .build();
        options.addOption(spectrumPropertyOption);
        Option mongoUserOption = Option.builder("u")
                .required(true)
                .hasArg()
                .longOpt("user")
                .desc("User name to use for login into MongoDB.")
                .build();
        options.addOption(mongoUserOption);
        Option mongoPasswordOption = Option.builder("p")
                .required(true)
                .hasArg()
                .longOpt("password")
                .desc("User password to use for login into MongoDB.")
                .build();
        options.addOption(mongoPasswordOption);
        Option mongoAuthDBOption = Option.builder("a")
                .required(true)
                .hasArg()
                .longOpt("authDB")
                .desc("Authentication database name to use for login into MongoDB.")
                .build();
        options.addOption(mongoAuthDBOption);
        Option mongoDatabaseNameOption = Option.builder("db")
                .required(true)
                .hasArg()
                .longOpt("database")
                .desc("Database name to use for operations in MongoDB.")
                .build();
        options.addOption(mongoDatabaseNameOption);
        Option mongoCollectionNameOption = Option.builder("c")
                .required(true)
                .hasArg()
                .longOpt("collection")
                .desc("Collection name to fetch from selected database in MongoDB.")
                .build();
        options.addOption(mongoCollectionNameOption);
        Option nthreadsOption = Option.builder("nt")
                .required(false)
                .hasArg()
                .longOpt("nthreads")
                .desc("Number of threads to use for parallelization. The default is set to 1.")
                .build();
        options.addOption(nthreadsOption);

        return options;
    }

    public static void main(String[] args) throws ParseException, CDKException {
        final Options options = setupOptions(args);
        CommandLineParser parser = new DefaultParser();
        try {
            final CommandLine cmd = parser.parse(options, args);
            final String pathToNMRShiftDB = cmd.getOptionValue("nmrdb");
            final String nucleus = cmd.getOptionValue("nucleus");
            final String mongoUser = cmd.getOptionValue("u");
            final String mongoPassword = cmd.getOptionValue("p");
            final String mongoAuthDB = cmd.getOptionValue("a");
            final String mongoDBName = cmd.getOptionValue("db");
            final String mongoDBCollection = cmd.getOptionValue("c");
            final int nThreads = Integer.parseInt(cmd.getOptionValue("nt", "1"));

            SolventDeviationsLookupTableBuilder.buildSolventDeviationsLookupTable(pathToNMRShiftDB, nucleus, nThreads, mongoUser, mongoPassword, mongoAuthDB, mongoDBName, mongoDBCollection);

        } catch (org.apache.commons.cli.ParseException e) {
            // TODO Auto-generated catch block
            HelpFormatter formatter = new HelpFormatter();
            formatter.setOptionComparator(null);
            String header = "Space for a header.\n\n";
            String footer = "\nSpace for a footer.";
            formatter.printHelp("java -jar something", header, options, footer, true);
            throw new org.apache.commons.cli.ParseException("Problem parsing command line");
        } catch (InterruptedException ex) {
            Logger.getLogger(SolventDeviationsLookupTableBuilder.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
}
