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
package analysis;

import casekit.NMR.dbservice.MongoDB;
import com.mongodb.MongoClient;
import com.mongodb.client.MongoCollection;
import model.SSCLibrary;
import org.apache.commons.cli.*;
import org.bson.Document;
import org.openscience.cdk.exception.CDKException;
import start.Start;

import java.io.FileNotFoundException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Class to build a HOSE code lookup table and to store it into MongoDB collection.
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class HOSECodeLookupTableBuilder {
       
    public static boolean buildHOSECodeShiftsLookupTable(final String pathToNMRShiftDB, final String NMRShiftDBSpectrumProperty, final int maxSphere, final int nThreads, 
            final String mongoUser, final String mongoPassword, final String mongoAuthDB, final String mongoDBName, final String mongoDBCollection, final boolean removeDuplicates){
        
        try {            
            final SSCLibrary sscLibrary = new SSCLibrary(nThreads);
            long offset = 0;
            for (int m = 0; m <= maxSphere; m++) {
                System.out.println("Build fragments for maxsphere: " + m);
                sscLibrary.extend(pathToNMRShiftDB, NMRShiftDBSpectrumProperty, m, offset);
                if (removeDuplicates) {
                    sscLibrary.removeDuplicates(Start.DUPLICATES_SHIFT_TOL);
                }
                System.out.println("done...");
                offset = sscLibrary.getLastSSCIndex() + 1;
            }            
            
            final MongoClient mongo = MongoDB.login(mongoUser, mongoPassword, mongoAuthDB);
            final MongoCollection<Document> collection = MongoDB.getCollection(mongo, mongoDBName, mongoDBCollection);
            collection.drop();
            System.out.println("Build and export lookup table...");
            sscLibrary.exportHOSECodeLookupTable(collection);
            System.out.println("done...");

            MongoDB.logout(mongo);
        } catch (FileNotFoundException | CDKException | InterruptedException | CloneNotSupportedException ex) {
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
        Option spectrumPropertyOption = Option.builder("prop")
                .required(true)
                .hasArg()
                .longOpt("property")
                .desc("Property string of spectrum type to use from NMRShiftDB, e.g. \"Spectrum 13C 0\".")
                .build();
        options.addOption(spectrumPropertyOption);        
        Option maxsphereOption = Option.builder("max")
                .required(true)
                .hasArg()
                .longOpt("maxsphere")
                .desc("Maximum sphere limit for SSC creation. Minimum value is 0.")
                .build();
        options.addOption(maxsphereOption);        
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
        Option removeDuplicatesOption = Option.builder("nd")
                .required(false)
                .longOpt("noduplicates")
                .desc("If given, the SSC library to build/extend will contain no structural duplicates.")
                .build();
        options.addOption(removeDuplicatesOption);

        return options;
    }
    
    public static void main(String[] args) throws ParseException, CDKException {
        final Options options = setupOptions(args);
        CommandLineParser parser = new DefaultParser();
        try {
            final CommandLine cmd = parser.parse(options, args);
            final String pathToNMRShiftDB = cmd.getOptionValue("nmrdb");
            final String NMRShiftDBSpectrumProperty = cmd.getOptionValue("prop");
            final int maxSphere = Integer.parseInt(cmd.getOptionValue("max"));
            if(maxSphere < 0){
                throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": invalid number of maximum sphere: \"" + maxSphere + "\" < 1");
            }
            final String mongoUser = cmd.getOptionValue("u");
            final String mongoPassword = cmd.getOptionValue("p");
            final String mongoAuthDB = cmd.getOptionValue("a");
            final String mongoDBName = cmd.getOptionValue("db");
            final String mongoDBCollection = cmd.getOptionValue("c");
            final int nThreads = Integer.parseInt(cmd.getOptionValue("nt", "1"));    
            boolean removeDuplicates = false;
            if(cmd.hasOption("noduplicates")){
                removeDuplicates = true;
            }
            
            HOSECodeLookupTableBuilder.buildHOSECodeShiftsLookupTable(pathToNMRShiftDB, NMRShiftDBSpectrumProperty, maxSphere, nThreads, mongoUser, mongoPassword, mongoAuthDB, mongoDBName, mongoDBCollection, removeDuplicates);
            
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
}
