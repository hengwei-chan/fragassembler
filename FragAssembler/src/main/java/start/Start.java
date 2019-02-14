package start;

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

import assembly.Assembly;
import model.SSC;
import casekit.NMR.DB;
import casekit.NMR.model.Spectrum;
import com.mongodb.BasicDBObject;
import com.mongodb.MongoClient;
import com.mongodb.MongoClientOptions;
import com.mongodb.MongoCredential;
import com.mongodb.ServerAddress;
import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoDatabase;
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
import model.SSCLibrary;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.bson.Document;
import org.json.simple.parser.ParseException;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import rank.SSCRanker;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Start {  
    
    private String pathToNMRShiftDB, mongoUser, mongoPassword, mongoAuthDB, mongoDBName, mongoDBCollection, pathToQueryFile;
    private int nThreads, nStarts, maxSphere;
    private boolean importFromNMRShiftDB;
    private SSCLibrary sscLibrary;   
    private SSCRanker sscRanker;
    
    private void parseArgs(String[] args) throws ParseException, org.apache.commons.cli.ParseException {
        Options options = setupOptions(args);
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine cmd = parser.parse(options, args);
            if(cmd.hasOption("import")){
                this.importFromNMRShiftDB = true;
                if(!cmd.hasOption("nmrshiftdb") || !cmd.hasOption("maxsphere")){
                    this.pathToNMRShiftDB = "";
                    this.maxSphere = -1;
                } else {
                    this.pathToNMRShiftDB = cmd.getOptionValue("nmrshiftdb");
                    this.maxSphere = Integer.parseInt(cmd.getOptionValue("maxsphere"));                
                }
                
            } else {
                this.importFromNMRShiftDB = false;
            }     
            this.nThreads = Integer.parseInt(cmd.getOptionValue("nthreads"));
            if(cmd.hasOption("nstarts")){
                this.nStarts = Integer.parseInt(cmd.getOptionValue("nstarts"));
            } else {
                this.nStarts = -1;
            }
            
            this.mongoUser = cmd.getOptionValue("user");
            this.mongoPassword = cmd.getOptionValue("password");
            this.mongoAuthDB = cmd.getOptionValue("authDB");
            this.mongoDBName = cmd.getOptionValue("database");
            this.mongoDBCollection = cmd.getOptionValue("collection");
            this.pathToQueryFile = cmd.getOptionValue("query");
            
        } catch (org.apache.commons.cli.ParseException e) {
            // TODO Auto-generated catch block
            HelpFormatter formatter = new HelpFormatter();
            formatter.setOptionComparator(null);
            String header = "Space for a header.\n\n";
            String footer = "\nSpace for a footer.";
            formatter.printHelp("java -jar something", header, options, footer, true);
            throw new org.apache.commons.cli.ParseException("Problem parsing command line");
        }
    }

    private Options setupOptions(String[] args) {
        Options options = new Options();
        Option importFromNMRShiftDBOption = Option.builder("i")
                .required(false)
                .longOpt("import")
                .desc("Indicates that a NMRShiftDB file (SDF) will be used to build a SSC library from that and to overwrite all entries within a MongoDB collection. The parameters \"nmrshiftdb\" and \"maxsphere\" must be set too.")
                .build();
        options.addOption(importFromNMRShiftDBOption);
        Option nmrshiftdb = Option.builder("nmrshiftdb")
                .required(false)
                .hasArg()
                .desc("Path to NMRShiftDB (SDF).")
                .build();
        options.addOption(nmrshiftdb);
        Option maxsphereOption = Option.builder("maxsphere")
                .required(false)
                .hasArg()
                .desc("Maximum sphere limit for SSC creation. Minimum value is 2. SSCs with spherical limit (>= 2) will be created until \"maxsphere\" is reached.")
                .build();
        options.addOption(maxsphereOption);
//        Option jsonLibrary = Option.builder("lib")
//                .required(true)
//                .hasArg()
//                .longOpt("library")
//                .desc("path to SSC Library as JSON file (required)")
//                .build();
//        options.addOption(jsonLibrary); 
        Option nthreadsOption = Option.builder("nt")
                .required(true)
                .hasArg()
                .longOpt("nthreads")
                .desc("Number of threads to use for parallelization.")
                .build();
        options.addOption(nthreadsOption); 
        Option nstartsOption = Option.builder("ns")
                .required(false)
                .hasArg()
                .longOpt("nstarts")
                .desc("Specified number of ranked SSCs to use for assembly process. The default is to use all matched SSC given a query spectrum.")
                .build();
        options.addOption(nstartsOption);
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
                .desc("Collection name to take from \"database\".")
                .build();
        options.addOption(mongoCollectionNameOption);
        Option pathToQueryFileOption = Option.builder("q")
                .required(true)
                .hasArg()
                .longOpt("query")
                .desc("Path to a file containing the query spectra in NMRShiftDB format. One spectrum for each line.")
                .build();
        options.addOption(pathToQueryFileOption);
        
        
        
        return options;
    }

    public static void main(String[] args) {
        // TODO Auto-generated method stub
        final Start start = new Start();
        try {
            start.parseArgs(args);
            start.start();
        } catch (IOException | CloneNotSupportedException | InterruptedException | org.apache.commons.cli.ParseException | ParseException | CDKException e) {
            System.err.println(e.getMessage());
        }

    }
    
    
    
    public void start() throws FileNotFoundException, CDKException, CloneNotSupportedException, InterruptedException, IOException {
        
        final MongoClient mongo;
        final MongoCollection<Document> collection;
        try {
            // Creating a Mongo client             
            mongo = new MongoClient(new ServerAddress(),
                    MongoCredential.createCredential(this.mongoUser, this.mongoAuthDB,
                            this.mongoPassword.toCharArray()),
                    MongoClientOptions.builder().build());
            System.out.println("Login to MongoDB was successfull");
            // Accessing the database 
            final MongoDatabase database = mongo.getDatabase(this.mongoDBName);
            System.out.println("Access to database \"" + this.mongoDBName + "\" was successfull");
            // Retrieving a collection
            collection = database.getCollection(this.mongoDBCollection);
            System.out.println("Retrieving to collection \"" + this.mongoDBCollection + "\" was successfull");
        } catch (Exception e) {
            throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": could not find database \"" + this.mongoDBName + "\" or collection \"" + this.mongoDBCollection + "\"");
        }                        
        
        this.sscLibrary = new SSCLibrary(this.nThreads);
        if(this.importFromNMRShiftDB){   
            if(this.pathToNMRShiftDB.isEmpty()){
                throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": invalid path to NMRShiftDB file: \"" + this.pathToNMRShiftDB + "\"");
            }
            if (this.maxSphere < 2) {
                throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": invalid number of maximum sphere: \"" + this.maxSphere + "\"");
            }
            
            // empty MongoDB collection
            collection.deleteMany(new BasicDBObject());
            System.out.println("Collection \"" + this.mongoDBCollection + "\" is now empty");
            
            int offset = 0;
            // create SSC library for a specific max sphere and insert into into the MongoDB collection
            for (int m = 2; m <= this.maxSphere; m++) {
                System.out.println("Building SSCs for " + m + "-spheres...");                
                this.sscLibrary.extend(this.pathToNMRShiftDB, "Spectrum 13C 0", "C", m, offset);                
                System.out.println("SSCs for " + m + "-spheres build!!!");
                System.out.println("-> #SSCs in SSC library: " + this.sscLibrary.getSSCCount());
                
                offset = Collections.max(sscLibrary.getSSCIndices()) + 1;
            }    
            this.sscLibrary.exportToMongoDB(collection);
            System.out.println("-> SSC library stored into MongoDB in database \"" + this.mongoDBName + "\" in collection \"" + this.mongoDBCollection + "\" -> collection size: " + collection.countDocuments() + "!!!");            
        } else {            
            this.sscLibrary.importFromMongoDB(collection);
            System.out.println("-> SSC library imported from MongoDB!!!");
            System.out.println("--> SSC library size:\t" + this.sscLibrary.getSSCCount());
        }
        mongo.close();                        
        
        this.processQueries();
    }
    
    private void processQueries() throws FileNotFoundException, InterruptedException, IOException, InvalidSmilesException, CDKException, CloneNotSupportedException{
        final SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        final BufferedReader br = new BufferedReader(new FileReader(this.pathToQueryFile));
        BufferedWriter bw;
        
        final double shiftTol = 3, thrsMatchFactor = 5, pickPrecision = 0.3;
        final int minMatchingSphereCount = 0;
        int nStartSSCs;
        Spectrum querySpectrum;
        HashMap<String, SSC> solutions;
        this.sscRanker = new SSCRanker(this.sscLibrary, this.nThreads);
        
        IAtomContainer solutionAtomContainer;
        ArrayList<Double> shiftsQuerySpectrum, shiftsSolutionSpectrum;
        int querySpectrumCounter = 0;
        
        System.out.println("\n\n-> processing query file \"" + this.pathToQueryFile + "\" ...");        
        
        Iterator<String> it = br.lines().iterator();
        while (it.hasNext()) {
            System.out.println("\n\now processing query " + querySpectrumCounter + "\n");
            
            querySpectrum = DB.NMRShiftDBSpectrumToSpectrum(it.next(), "C");
            this.sscRanker.rank(querySpectrum, pickPrecision);
            System.out.println("no. of matches: " + this.sscRanker.getMatchFactors().size());
            System.out.println("match factors: " + this.sscRanker.getMatchFactors());
            System.out.println("ranked SSC indices: " + this.sscRanker.getRankedSSCIndices());
            System.out.println("ranked match factors: " + this.sscRanker.getRankedMatchFactors());
            System.out.println("ranked SSC library indices: " + this.sscRanker.getRankedSSCLibrary().getSSCIndices() + "\n");
            
            if((this.nStarts > 0) && (this.nStarts < this.sscRanker.getMatchFactors().size())){
                nStartSSCs = this.nStarts;
            } else {
                nStartSSCs = this.sscRanker.getMatchFactors().size();
            }
            System.out.println("\nnumber of start SSCs for query " + querySpectrumCounter + ":\t" + nStartSSCs);
            
            solutions = Assembly.assemble(nStartSSCs, nThreads, sscRanker.getRankedSSCLibrary(), minMatchingSphereCount, querySpectrum, shiftTol, thrsMatchFactor, pickPrecision);
            
            System.out.println("\nsolutions for query " + querySpectrumCounter + ":\t" + solutions.size());
            
            shiftsQuerySpectrum = querySpectrum.getShifts(0);
            Collections.sort(shiftsQuerySpectrum);
            
            bw = new BufferedWriter(new FileWriter("results_" + querySpectrumCounter + ".smiles"));
            for (final String smiles : solutions.keySet()) {
                solutionAtomContainer = smilesParser.parseSmiles(smiles);
                
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
        br.close();
    }
        
}
