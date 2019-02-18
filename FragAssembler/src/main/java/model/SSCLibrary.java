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
package model;

import casekit.NMR.DB;
import casekit.NMR.Utils;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.stream.JsonWriter;
import com.mongodb.client.MongoCollection;
import com.mongodb.util.JSON;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import org.openscience.cdk.exception.CDKException;
import fragmentation.Fragmentation;
import java.util.Collections;
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
import org.bson.Document;
import org.json.simple.parser.ParseException;

/**
 * Class to search for matches between a SSC library and query spectra.
 * 
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public final class SSCLibrary {

    private final HashMap<Integer, SSC> map;
    private final HashMap<String, ArrayList<Double>> HOSECodeLookupTable;
    private int nThreads;
    
    /**
     * Instanciates a new object of this class.
     *
     */
    public SSCLibrary(){
        this(1);
    }
    
    /**
     * Instanciates a new object of this class.
     *
     * @param nThreads number of threads to use for parallelization
     */
    public SSCLibrary(final int nThreads){
        this.map = new HashMap<>();
        this.HOSECodeLookupTable = new HashMap<>();
        this.nThreads = nThreads;
    }
    
    /**
     * Creates a lookup table for generated HOSE codes regarding the given 
     * SSC library. Each HOSE code entry contains a list of belonging shift 
     * values. 
     *
     * @throws java.lang.InterruptedException
     */
    public void buildHOSELookupTable() throws InterruptedException{  
        this.HOSECodeLookupTable.clear();
        for (final SSC ssc : this.map.values()) {
            if(ssc != null){
                Utils.combineHashMaps(this.HOSECodeLookupTable, ssc.getHOSECodeLookupShifts());                          
            }            
        }
    }
    
    /**
     * Returns a HashMap of the shift lists in HOSE code lookup table. 
     * Before usage of this method, the HOSE code lookup table has to been 
     * built via {@link #buildHOSELookupTable() }.
     *
     * @return HashMap with HOSE codes as keys and lists of chemical shifts as
     * values
     *
     * @see #buildHOSELookupTable() 
     */
    public HashMap<String, ArrayList<Double>> getHOSECodeLookupTable() {
        return this.HOSECodeLookupTable;
    }
     
    public int getNThreads(){
        return this.nThreads;
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
     * Writes this SSC library into a MongoDB collection.
     *
     * @param collection
     */
    public void exportToMongoDB(final MongoCollection<Document> collection){
        final Gson gson = new Gson();
        Document document;
        for (SSC ssc : this.getSSCs()) {
            document = new Document();
            document.append("_id", ssc.getIndex());
            document.append("ssc", JSON.parse(gson.toJson(gson.toJsonTree(new SSCInJSON(ssc), SSCInJSON.class))));
            try {
                collection.insertOne(document);
//                    System.out.println("Document inserted successfully");
            } catch (Exception e) {
                System.err.println("export for key \"" + ssc.getIndex() + "\" failed: " + e.getMessage());
                System.out.println("-> document: \n" + document.toJson());
            }
        }
    }
    
    public void exportToJSONFile(final String pathToJSON) throws IOException, InterruptedException {                        
        final Gson gson = new GsonBuilder().setPrettyPrinting().create();
        final BufferedWriter bw = new BufferedWriter(new FileWriter(pathToJSON));
        final JsonObject jsonObject = new JsonObject();
        for (int sscIndex : this.map.keySet()) {
            jsonObject.add(String.valueOf(sscIndex), gson.toJsonTree(new SSCInJSON(this.map.get(sscIndex)), SSCInJSON.class));
        }
        
       
        
//        // initialize an executor for parallelization
//        final ExecutorService executor = Utils.initExecuter(this.nThreads);
//        final ArrayList<Callable<JsonElement>> callables = new ArrayList<>();
//        // add all task to do
//        for (final int sscIndex : this.map.keySet()) {
//            callables.add((Callable<JsonElement>) () -> {
//                return gson.toJsonTree(new SSCInJSON(this.map.get(sscIndex)), SSCInJSON.class);
//            });
//        }
//        // execute all task in parallel
//        executor.invokeAll(callables)
//                .stream()
//                .map(future -> {
//                    try {
//                        return future.get();
//                    } catch (InterruptedException | ExecutionException e) {
//                        throw new IllegalStateException(e);
//                    }
//                })
//                .forEach((jsonElement) -> {
//                    if (jsonElement != null) {
////                        System.out.println("index: " + jsonElement.getAsJsonObject().get("index").toString());
//                        jsonObject.add(jsonElement.getAsJsonObject().get("index").toString(), jsonElement);
//                    }
//                });
//        // shut down the executor service
//        Utils.stopExecuter(executor, 5);
        
        
        
        
        JsonWriter writer = new JsonWriter(bw);
        writer.setIndent("\t"); // <-- without this the file is written without line-breaks and indents        
//        writer.setLenient(true);
        gson.toJson(jsonObject, writer);
        bw.close();
    }
    
    /**
     * Imports a SSC library from a JSON file. All current SSCs will be 
     * deleted.
     *
     * @param pathToJSON JSON file containing ths SSC library
     * @param offset offset value for start index (index of first SSC)      
     * @throws java.lang.InterruptedException      
     * @throws java.io.FileNotFoundException      
     * @throws org.openscience.cdk.exception.CDKException      
     * @throws java.lang.CloneNotSupportedException      
     */
    public void importFromJSONFile(final String pathToJSON, final int offset) throws InterruptedException, FileNotFoundException, CDKException, CloneNotSupportedException  {
        this.removeAll();
        this.extend(pathToJSON, offset);
    }
    
    /**
     * Imports a SSC library by documents from a MongoDB collection 
     * containing the SSC information.  All current SSCs will be 
     * deleted.
     *
     * @param collection
     * @throws org.openscience.cdk.exception.CDKException
     * @throws java.lang.CloneNotSupportedException
     * @throws java.lang.InterruptedException
     * 
     */
    public void importFromMongoDB(final MongoCollection<Document> collection) throws CDKException, CloneNotSupportedException, InterruptedException  {
        this.removeAll();
        this.extend(collection);
    }
    
    /**
     * Creates a SSC library from a given NMRShiftDB SDF via fragmentation of
     * each structure and each of its atoms as individual starting point 
     * (central atom). All current SSCs will be deleted.
     *
     * @param pathToNMRShiftDB path to NMRShiftDB SDF
     * @param property property string of spectrum to use, 
     * e.g. "Spectrum 13C 0"
     * @param atomType atom type of spectrum to use, e.g. "C"
     * @param maxSphere maximum number of spheres in fragmention
     * @param offset offset number for next SSC indices to use as keys in SSC 
     * library
     * @throws FileNotFoundException
     * @throws CDKException
     * @throws CloneNotSupportedException
     * @throws java.lang.InterruptedException
     * 
     * @see #removeAll() 
     * @see #extend(java.lang.String, java.lang.String, java.lang.String, int, int)  
     */
    public void importFromNMRShiftDB(final String pathToNMRShiftDB, final String property, final String atomType, final int maxSphere, final int offset) throws FileNotFoundException, CDKException, CloneNotSupportedException, InterruptedException {
        this.removeAll();
        this.extend(pathToNMRShiftDB, property, atomType, maxSphere, offset);
    }
    
    /**
     * Extends this SSC library by a given second one. All SSCs in this 
     * library object whose indices also exist in the given SSC library will be 
     * replaced.
     *
     * @param sscLibrary second SSC library to add
     * @see #containsSSC(int) 
     */
    public void extend(final SSCLibrary sscLibrary){
        this.map.putAll(sscLibrary.getMap());
    }
    
    /**
     * Extends this SSC library by documents from a MongoDB collection 
     * containing the SSC information. 
     * All SSCs in this library object whose indices also exist in the 
     * given map will be replaced.
     *
     * @param collection 
     * @throws org.openscience.cdk.exception.CDKException  
     * @throws java.lang.CloneNotSupportedException  
     * @throws java.lang.InterruptedException  
     * 
     * @see #containsSSC(int)
     */
    public void extend(final MongoCollection<Document> collection) throws CDKException, CloneNotSupportedException, InterruptedException  {
        final Gson gson = new Gson();        
        SSC ssc;
        for (final Document sscDocDB : collection.find()) {
            ssc = gson.fromJson(((Document) sscDocDB.get("ssc")).toJson(), SSCInJSON.class).toSSC();
            if (!this.insert(ssc)) {
                throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": insertion SSC with index " + ssc.getIndex() + " failed");
            }
        }
        
//        // initialize an executor for parallelization
//        final ExecutorService executor = Utils.initExecuter(this.nThreads);
//        final ArrayList<Callable<SSC>> callables = new ArrayList<>();
//        // add all task to do
//        for (final Document sscDocDB : collection.find()) {
//            callables.add((Callable<SSC>) () -> {
//                return gson.fromJson(((Document) sscDocDB.get("ssc")).toJson(), SSCInJSON.class).toSSC();
//            });
//        }
//        // execute all task in parallel
//        executor.invokeAll(callables)
//                .stream()
//                .map(future -> {
//                    try {
//                        return future.get();
//                    } catch (InterruptedException | ExecutionException e) {
//                        throw new IllegalStateException(e);
//                    }
//                })
//                .forEach((ssc) -> {
//                    if (ssc != null) {
//                        ssc.setIndex(ssc.getIndex());
//                        if (!this.insert(ssc)) {
//                            try {
//                                throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": insertion SSC with index " + ssc.getIndex() + " failed");
//                            } catch (CDKException ex) {
//                                Logger.getLogger(SSCLibrary.class.getName()).log(Level.SEVERE, null, ex);
//                            }
//                        }
//                    }
//                });
//        // shut down the executor service
//        Utils.stopExecuter(executor, 5);
    }
    
    /**
     * Extends this SSC library by created SSCs from a NMRShiftDB file. 
     * All SSCs in this library object whose indices also exist in the 
     * given map will be replaced.
     *
     * @param pathToNMRShiftDB path to NMRShiftDB SDF
     * @param property property string of spectrum to use, 
     * e.g. "Spectrum 13C 0"
     * @param atomType atom type of spectrum to use, e.g. "C"
     * @param maxSphere maximum number of spheres in fragmention
     * @param offset offset number for next SSC indices to use as keys in SSC 
     * library
     * @throws java.io.FileNotFoundException
     * @throws org.openscience.cdk.exception.CDKException
     * @throws java.lang.InterruptedException
     * @throws java.lang.CloneNotSupportedException
     * 
     * @see Fragmentation#buildSSCLibrary(java.util.HashMap, int, int, int)  
     * @see #containsSSC(int) 
     */
    public void extend(final String pathToNMRShiftDB, final String property, final String atomType, final int maxSphere, final int offset) throws FileNotFoundException, CDKException, InterruptedException, CloneNotSupportedException {
        final HashMap<Integer, Object[]> SSCComponentsSet = DB.getSSCComponentsFromNMRShiftDB(pathToNMRShiftDB, property, atomType);
        this.extend(Fragmentation.buildSSCLibrary(SSCComponentsSet, maxSphere, this.nThreads, offset));
    }
    
    /**
     * Extends this SSC library by SSCs stored in a JSON file. All SSCs in this 
     * library object whose indices also exist in the given file will be 
     * replaced.
     *
     * @param pathToJSON path to JSON file containing SSCs
     * @param offset offset number for next SSC indices to use as keys in SSC 
     * library
     * @throws java.lang.InterruptedException
     * @throws java.io.FileNotFoundException
     * @throws org.openscience.cdk.exception.CDKException
     * @throws java.lang.CloneNotSupportedException
     * @see #containsSSC(int) 
     */
    public void extend(final String pathToJSON, final int offset) throws InterruptedException, FileNotFoundException, CDKException, CloneNotSupportedException{
        final Gson gson = new Gson();
//        SSC ssc;
        final BufferedReader br = new BufferedReader(new FileReader(pathToJSON));
        final JsonObject jsonObject = new JsonParser().parse(br).getAsJsonObject();
//        for (final String sscIndex: jsonObject.keySet()) {
//            ssc = gson.fromJson(jsonObject.get(sscIndex), SSCInJSON.class).toSSC();
//            ssc.setIndex(offset + ssc.getIndex());
//            if(!this.insert(ssc)){
//                throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": insertion SSC with index " + (offset + ssc.getIndex()) + "(incl. offset: " + offset + ") failed");
//            }
//        }
        
        
        
        // initialize an executor for parallelization
        final ExecutorService executor = Utils.initExecuter(this.nThreads);
        final ArrayList<Callable<SSC>> callables = new ArrayList<>();
        // add all task to do
        for (final String sscIndex: jsonObject.keySet()) {
            callables.add((Callable<SSC>) () -> {
                return gson.fromJson(jsonObject.get(sscIndex), SSCInJSON.class).toSSC();
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
                .forEach((ssc) -> {
                    if (ssc != null) {
                        ssc.setIndex(offset + ssc.getIndex());
                        if(!this.insert(ssc)){
                            try {
                                throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": insertion SSC with index " + (offset + ssc.getIndex()) + "(incl. offset: " + offset + ") failed");
                            } catch (CDKException ex) {
                                Logger.getLogger(SSCLibrary.class.getName()).log(Level.SEVERE, null, ex);
                            }
                        }                        
                    }
                });
        // shut down the executor service
        Utils.stopExecuter(executor, 5);
        
    }
    
    public int getSSCCount(){
        return this.map.size();
    }
    
    public SSC getSSC(final int sscIndex){
        return this.map.get(sscIndex);
    } 
    
    public boolean containsSSC(final int sscIndex){
        return this.map.containsKey(sscIndex);
    }    
    
    public HashMap<Integer, SSC> getMap(){
        return this.map;
    }
    
    public HashSet<Integer> getSSCIndices(){
        return new HashSet<>(this.map.keySet());
    }
    
    public Collection<SSC> getSSCs(){
        return this.map.values();
    } 
    
    /**
     * Inserts a new SSC to this SSC library.
     *
     * @param ssc SSC to add to this library
     * @return false if this SSC library already contains the index of {@code ssc} 
     * or {@code ssc.getIndex() < 0}; otherwise true
     */
    public boolean insert(final SSC ssc) {
        if(this.containsSSC(ssc.getIndex()) || ssc.getIndex() < 0){
            return false;
        }
        this.map.put(ssc.getIndex(), ssc);        
        
        return true;
    }
    
    public boolean remove(final int sscIndex){
        if(this.containsSSC(sscIndex)){
            return false;
        }
        this.map.remove(sscIndex);
        
        return true;
    }
    
    public SSCLibrary getClone() throws CDKException, CloneNotSupportedException {
        final SSCLibrary sscLibrary = new SSCLibrary();
        sscLibrary.extend(this);
        
        return sscLibrary;
    } 
    
    public void removeAll(){
        this.map.clear();
    }
    
    // ################################
    
    public static void main(String[] args) {
        final SSCLibrary sscLibrary = new SSCLibrary();
        try {
            final Object[] argValues = sscLibrary.parseArgs(args);
            final String pathToNMRShiftDB = (String) argValues[0];
            final String pathToJSON = (String) argValues[1];            
            final int maxsphere = Integer.parseInt((String) argValues[2]);
            final int nthreads = Integer.parseInt((String) argValues[3]);
            final int offset = Integer.parseInt((String) argValues[4]);
            final String property = (String) argValues[5];
            final String atomType = (String) argValues[6];
            System.out.println("-> importing from: " + pathToNMRShiftDB + " ...");
            sscLibrary.setNThreads(nthreads);
            sscLibrary.importFromNMRShiftDB(pathToNMRShiftDB, property, atomType, maxsphere, offset);
            System.out.println(" -> imported! -> #SSC: " + sscLibrary.getSSCCount() + ", next offset: " + (Collections.max(sscLibrary.getSSCIndices()) + 1));
            System.out.println("-> exporting to: " + pathToJSON + " ...");
            sscLibrary.exportToJSONFile(pathToJSON);
            System.out.println(" -> exported!!!");
            
            
        } catch (IOException | CloneNotSupportedException | InterruptedException | org.apache.commons.cli.ParseException | ParseException | CDKException e) {
            System.err.println(e.getMessage());
        }
    }
    
    
    private Object[] parseArgs(String[] args) throws ParseException, org.apache.commons.cli.ParseException {
        final Options options = setupOptions();
        final CommandLineParser parser = new DefaultParser();
        final Object[] argValues = new Object[options.getOptions().size()];
        try {
            final CommandLine cmd = parser.parse(options, args);
            argValues[0] = cmd.getOptionValue("nmrshiftdb");
            argValues[1] = cmd.getOptionValue("library");
            argValues[2] = cmd.getOptionValue("maxsphere");            
            argValues[3] = cmd.getOptionValue("nthreads", "1");            
            argValues[4] = cmd.getOptionValue("offset", "0");            
            argValues[5] = cmd.getOptionValue("property", "Spectrum 13C 0");            
            argValues[6] = cmd.getOptionValue("atomtype", "C");            

        } catch (org.apache.commons.cli.ParseException e) {
            // TODO Auto-generated catch block
            HelpFormatter formatter = new HelpFormatter();
            formatter.setOptionComparator(null);
            String header = "Space for a header.\n\n";
            String footer = "\nSpace for a footer.";
            formatter.printHelp("java -jar something", header, options, footer, true);
            throw new org.apache.commons.cli.ParseException("Problem parsing command line");
        }
        
        return argValues;
    }

    private Options setupOptions() {
        Options options = new Options();
        Option nmrshiftdb = Option.builder("db")
                .required(true)
                .hasArg()
                .longOpt("nmrshiftdb")
                .desc("path to NMRShiftDB sdf (required)")
                .build();
        options.addOption(nmrshiftdb);
        Option jsonLibrary = Option.builder("lib")
                .required(true)
                .hasArg()
                .longOpt("library")
                .desc("path to SSC Library as JSON file (required)")
                .build();
        options.addOption(jsonLibrary);
        Option maxsphere = Option.builder("m")
                .required(true)
                .hasArg()
                .longOpt("maxsphere")
                .desc("maximum sphere limit for SSC library fragments (required)")
                .build();
        options.addOption(maxsphere);
        Option nthreads = Option.builder("n")
                .required(false)
                .hasArg()
                .longOpt("nthreads")
                .desc("number of threads to use for parallelization (optional); this is set to \"1\" by default")
                .build();
        options.addOption(nthreads);
        Option offset = Option.builder("o")
                .required(false)
                .hasArg()
                .longOpt("offset")
                .desc("offset value for SSC indexing (optional); this is set to \"0\" by default")
                .build();
        options.addOption(offset);
        Option property = Option.builder("p")
                .required(false)
                .hasArg()
                .longOpt("property")
                .desc("property string for spectrum to use in NMRShiftDB (optional); this is set to \"Spectrum 13C 0\" by default")
                .build();
        options.addOption(property);
        Option atomType = Option.builder("t")
                .required(false)
                .hasArg()
                .longOpt("atomtype")
                .desc("atom type (element symbol) for spectrum to use in NMRShiftDB (optional); this is set to \"C\" by default")
                .build();
        options.addOption(atomType);

        return options;
    }


}
