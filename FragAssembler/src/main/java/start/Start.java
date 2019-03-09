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

import java.io.FileNotFoundException;
import java.io.IOException;
import model.SSCLibrary;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.json.simple.parser.ParseException;
import org.openscience.cdk.exception.CDKException;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Start {  
    
    private String pathToNMRShiftDB, mongoUser, mongoPassword, mongoAuthDB, mongoDBName, mongoDBCollection, pathToQueriesFile, pathToOutputsFolder, pathToJSON, format;
    private int nThreads, nStarts, maxSphere;
    private boolean importFromNMRShiftDB, extendFromNMRShiftDB, useMongoDB, useJSON;
    private SSCLibrary sscLibrary;   
    private ProcessQueriesJSON processQueriesJSON;
    private ProcessQueriesMongoDB processQueriesMongoDB;
    public static final String SPECTRUM_PROPERTY = "Spectrum 13C 0";
    public static final String SPECTRUM_PROPERTY_ATOMTYPE = "C";
    public static final double PICK_PRECISION = 0.3, MATCH_FACTOR_THRS = 3.0, SHIFT_TOL = 3.0;
    public static final int MIN_MATCHING_SPHERE_COUNT = 1;
    
    public void start() throws FileNotFoundException, CDKException, CloneNotSupportedException, InterruptedException, IOException {
        if (this.useMongoDB) {
            final PrepareMongoDB prepareMongoDB = new PrepareMongoDB(this.nThreads, this.importFromNMRShiftDB, this.extendFromNMRShiftDB, this.pathToNMRShiftDB, this.maxSphere);
            prepareMongoDB.prepare(this.mongoUser, this.mongoPassword, this.mongoAuthDB, this.mongoDBName, this.mongoDBCollection);
            this.processQueriesMongoDB = new ProcessQueriesMongoDB(this.pathToQueriesFile, this.pathToOutputsFolder, this.mongoUser, this.mongoPassword, this.mongoAuthDB, this.mongoDBName, this.mongoDBCollection, this.nThreads, this.nStarts);
            this.processQueriesMongoDB.process();
        } else if (this.useJSON) {
            final PrepareJSONFile prepareJSONFile = new PrepareJSONFile(this.nThreads, this.importFromNMRShiftDB, this.extendFromNMRShiftDB, this.pathToNMRShiftDB, this.maxSphere);
            this.sscLibrary = prepareJSONFile.prepare(this.pathToJSON);
            this.processQueriesJSON = new ProcessQueriesJSON(this.sscLibrary, this.pathToQueriesFile, this.pathToOutputsFolder, this.nThreads, this.nStarts);
            this.processQueriesJSON.process();
        } else {
            throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": invalid format: \"" + this.format + "\"");
        }
        System.gc();
    }
    
    
    private void parseArgs(String[] args) throws ParseException, org.apache.commons.cli.ParseException, CDKException {
        Options options = setupOptions(args);
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine cmd = parser.parse(options, args);
            this.format = cmd.getOptionValue("format");            
            switch (this.format) {
                case "m":
                    this.useMongoDB = true;
                    this.useJSON = false;
                    
                    this.mongoUser = cmd.getOptionValue("user");
                    this.mongoPassword = cmd.getOptionValue("password");
                    this.mongoAuthDB = cmd.getOptionValue("authDB");
                    this.mongoDBName = cmd.getOptionValue("database");
                    this.mongoDBCollection = cmd.getOptionValue("collection");
                    if ((this.mongoUser == null)
                            || (this.mongoPassword == null)
                            || (this.mongoAuthDB == null)
                            || (this.mongoDBName == null)
                            || (this.mongoDBCollection == null)) {
                        throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": at least one the parameters \"u\", \"p\", \"a\", \"db\" and \"c\" is missing");
                    }
                    break;
                case "j":
                    this.useMongoDB = false;
                    this.useJSON = true;
                    
                    this.pathToJSON = cmd.getOptionValue("json");
                    if (this.pathToJSON == null) {
                        throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": parameter \"j\" is missing");
                    }
                    break;
                default:
                    this.useMongoDB = false;
                    this.useJSON = false;
                    throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": invalid format: \"" + this.format + "\"");
            }            
            
            if(cmd.hasOption("import") || cmd.hasOption("extend")){
                if (cmd.hasOption("import")) {
                    this.importFromNMRShiftDB = true;
                    this.extendFromNMRShiftDB = false;                    
                } else {//if (cmd.hasOption("extend")) {
                    this.importFromNMRShiftDB = false;
                    this.extendFromNMRShiftDB = true;                    
                }                  
                if (!cmd.hasOption("nmrshiftdb") || !cmd.hasOption("maxsphere")) {
                    throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": parameter \"nmrshiftdb\" and/or \"maxsphere\" is missing\"");
                } else {
                    this.pathToNMRShiftDB = cmd.getOptionValue("nmrshiftdb");   
                    this.maxSphere = Integer.parseInt(cmd.getOptionValue("maxsphere"));
                    if (this.maxSphere < 2) {
                        throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": invalid number of maximum sphere: \"" + this.maxSphere + "\" < 2");
                    }
                    System.out.println("-nmrshiftdb: " + this.pathToNMRShiftDB);
                    System.out.println("-maxsphere: " + this.maxSphere);
                }
            } else {
                this.importFromNMRShiftDB = false;
                this.extendFromNMRShiftDB = false;
            }     
            
            System.out.println("-importFromNMRShiftDB: " + this.importFromNMRShiftDB);
            System.out.println("-extendFromNMRShiftDB: " + this.extendFromNMRShiftDB);
            this.nThreads = Integer.parseInt(cmd.getOptionValue("nthreads", "1"));            
            this.nStarts = Integer.parseInt(cmd.getOptionValue("nstarts", "-1"));                        
            
            this.pathToQueriesFile = cmd.getOptionValue("query");
            this.pathToOutputsFolder = cmd.getOptionValue("output", ".");
            this.pathToJSON = cmd.getOptionValue("json");
            
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
        Option dbFormatOption = Option.builder("f")
                .required(true)
                .hasArg()
                .longOpt("format")
                .desc("Format to use: " 
                        + "\ncase 1: \"j\" for JSON. The parameter \"j\" has to be set." 
                        + "\ncase 2: \"m\" for MongoDB. The parameters \"u\", \"p\", \"a\", \"db\" and \"c\" have to be set.")
                .build();
        options.addOption(dbFormatOption);
        Option pathToQueryFileOption = Option.builder("q")
                .required(true)
                .hasArg()
                .longOpt("query")
                .desc("Path to a file containing the query spectra in NMRShiftDB format. One spectrum for each line.")
                .build();
        options.addOption(pathToQueryFileOption);
        Option maxsphereOption = Option.builder("maxsphere")
                .required(false)
                .hasArg()
                .desc("Maximum sphere limit for SSC creation. Minimum value is 2."
                        + "\nIf the \"import\" option is selected: SSCs with spherical limit (>= 2) will be created until \"maxsphere\" is reached."
                        + "\nIf the \"extend\" option is selected: only SSCs with exactly spherical limit \"maxsphere\" will be created." 
                        + "\nIf both parameters are given: the \"import\" will be done.")
                .build();
        options.addOption(maxsphereOption);
        Option pathToOutputsFolderOption = Option.builder("o")
                .required(false)
                .hasArg()
                .longOpt("output")
                .desc("Path to a output directory for results. The default is set to \".\".")
                .build();
        options.addOption(pathToOutputsFolderOption);
        Option nthreadsOption = Option.builder("nt")
                .required(false)
                .hasArg()
                .longOpt("nthreads")
                .desc("Number of threads to use for parallelization. The default is set to 1.")
                .build();
        options.addOption(nthreadsOption); 
        Option nstartsOption = Option.builder("ns")
                .required(false)
                .hasArg()
                .longOpt("nstarts")
                .desc("Specified number of ranked SSCs to use for assembly process. The default is set to use all matched SSC given a query spectrum.")
                .build();
        options.addOption(nstartsOption);
        Option importFromNMRShiftDBOption = Option.builder("import")
                .required(false)
                .desc("Indicates that a NMRShiftDB file (SDF) will be used to build a SSC library from that and to overwrite all entries within a MongoDB collection or JSON file. The parameters \"nmrshiftdb\" and \"maxsphere\" must be set too.")
                .build();
        options.addOption(importFromNMRShiftDBOption);
        Option extendFromNMRShiftDBOption = Option.builder("extend")
                .required(false)
                .desc("Indicates that a NMRShiftDB file (SDF) will be used to build a SSC library from that and to extend a MongoDB collection or JSON file. The parameters \"nmrshiftdb\" and \"maxsphere\" must be set too.")
                .build();
        options.addOption(extendFromNMRShiftDBOption);
        Option nmrshiftdb = Option.builder("nmrshiftdb")
                .required(false)
                .hasArg()
                .desc("Path to NMRShiftDB (SDF).")
                .build();
        options.addOption(nmrshiftdb);
        Option mongoUserOption = Option.builder("u")
                .required(false)
                .hasArg()
                .longOpt("user")
                .desc("User name to use for login into MongoDB.")
                .build();
        options.addOption(mongoUserOption); 
        Option mongoPasswordOption = Option.builder("p")
                .required(false)
                .hasArg()
                .longOpt("password")
                .desc("User password to use for login into MongoDB.")
                .build();
        options.addOption(mongoPasswordOption);
        Option mongoAuthDBOption = Option.builder("a")
                .required(false)
                .hasArg()
                .longOpt("authDB")
                .desc("Authentication database name to use for login into MongoDB.")
                .build();
        options.addOption(mongoAuthDBOption);
        Option mongoDatabaseNameOption = Option.builder("db")
                .required(false)
                .hasArg()
                .longOpt("database")
                .desc("Database name to use for operations in MongoDB.")
                .build();
        options.addOption(mongoDatabaseNameOption);
        Option mongoCollectionNameOption = Option.builder("c")
                .required(false)
                .hasArg()
                .longOpt("collection")
                .desc("Collection name to fetch from selected database in MongoDB.")
                .build();
        options.addOption(mongoCollectionNameOption);
        Option pathToJSONLibraryOption = Option.builder("j")
                .required(false)
                .hasArg()
                .longOpt("json")
                .desc("Path to SSC library in JSON format.")
                .build();
        options.addOption(pathToJSONLibraryOption);
        
        
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
    
}
