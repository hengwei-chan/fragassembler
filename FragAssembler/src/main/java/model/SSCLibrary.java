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
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.stream.JsonWriter;
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
    public boolean setNthreads(final int nThreads){
        if(nThreads < 1){
            return false;
        }
        this.nThreads = nThreads;
        
        return true;
    }
    
    public void exportToJSON(final String pathToJSON) throws IOException, InterruptedException {                        
        final Gson gson = new GsonBuilder().setPrettyPrinting().create();
        final BufferedWriter bw = new BufferedWriter(new FileWriter(pathToJSON));
        JsonObject jsonObject = new JsonObject();
//        for (final int sscIndex : this.map.keySet()) {
//            jsonObject.add(String.valueOf(sscIndex), gson.toJsonTree(new SSCInJSON(this.map.get(sscIndex)), SSCInJSON.class));
//        }
        
       
        
        // initialize an executor for parallelization
        final ExecutorService executor = Utils.initExecuter(this.nThreads);
        final ArrayList<Callable<JsonElement>> callables = new ArrayList<>();
        // add all task to do
        for (final int sscIndex : this.map.keySet()) {
            callables.add((Callable<JsonElement>) () -> {
                return gson.toJsonTree(new SSCInJSON(this.map.get(sscIndex)), SSCInJSON.class);
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
                .forEach((jsonElement) -> {
                    if (jsonElement != null) {
//                        System.out.println("index: " + jsonElement.getAsJsonObject().get("index").toString());
                        jsonObject.add(jsonElement.getAsJsonObject().get("index").toString(), jsonElement);
                    }
                });
        // shut down the executor service
        Utils.stopExecuter(executor, 5);
        
        
        
        
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
    public void importFromJSON(final String pathToJSON, final int offset) throws InterruptedException, FileNotFoundException, CDKException, CloneNotSupportedException  {
        this.removeAll();
        this.extend(pathToJSON, offset);
    }
    
    /**
     * Creates a SSC library from a given NMRShiftDB SDF via fragmentation of
     * each structure and each of its atoms as individual starting point 
     * (central atom). All current SSCs will be deleted.
     *
     * @param pathToNMRShiftDB path to NMRShiftDB SDF
     * @param maxSphere maximum number of spheres in fragmention
     * @param offset offset number for next SSC indices to use as keys in SSC 
     * library
     * @throws FileNotFoundException
     * @throws CDKException
     * @throws CloneNotSupportedException
     * @throws java.lang.InterruptedException
     * 
     * @see Fragmentation#buildSSCLibrary(java.util.HashMap, int, int, int) 
     */
    public void importFromNMRShiftDB(final String pathToNMRShiftDB, final int maxSphere, final int offset) throws FileNotFoundException, CDKException, CloneNotSupportedException, InterruptedException {
        this.removeAll();
        this.extend(pathToNMRShiftDB, maxSphere, offset);
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
     * Extends this SSC library by a given map. All SSCs in this 
     * library object whose indices also exist in the given map will be 
     * replaced.
     *
     * @param map map of indices (keys) and SSCs (values) to add
     * @see #containsSSC(int) 
     */
    public void extend(final HashMap<Integer, SSC> map){
        this.map.putAll(map);
    }
    
    /**
     * Extends this SSC library by created SSCs from a NMRShiftDB file. 
     * All SSCs in this library object whose indices also exist in the 
     * given map will be replaced.
     *
     * @param pathToNMRShiftDB path to NMRShiftDB SDF
     * @param maxSphere maximum number of spheres in fragmention
     * @param offset offset number for next SSC indices to use as keys in SSC 
     * library
     * @throws java.io.FileNotFoundException
     * @throws org.openscience.cdk.exception.CDKException
     * @throws java.lang.InterruptedException
     * @throws java.lang.CloneNotSupportedException
     * @see #containsSSC(int)
     */
    public void extend(final String pathToNMRShiftDB, final int maxSphere, final int offset) throws FileNotFoundException, CDKException, InterruptedException, CloneNotSupportedException {
        final String atomType = "C";
        final String spectrumProp = "Spectrum 13C 0";
        final HashMap<Integer, Object[]> SSCComponentsSet = DB.getSSCComponentsFromNMRShiftDB(pathToNMRShiftDB, spectrumProp, atomType);
        this.extend(Fragmentation.buildSSCLibrary(SSCComponentsSet, maxSphere, nThreads, offset));
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
        SSCInJSON sscInJSON;
        final BufferedReader br = new BufferedReader(new FileReader(pathToJSON));
        final JsonObject jsonObject = new JsonParser().parse(br).getAsJsonObject();
        System.out.println("-> parsing done!!! -> #SSCs: " + jsonObject.size());
        System.out.println("-> now inserting ... ");
        for (final String sscIndex: jsonObject.keySet()) {
            sscInJSON = gson.fromJson(jsonObject.get(sscIndex), SSCInJSON.class);
            if(!this.insert(sscInJSON.toSSC(), offset + sscInJSON.getIndex())){
                throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": insertion SSC with index " + (offset + sscInJSON.getIndex()) + " failed");
            }
        }
        
        
        
//        // initialize an executor for parallelization
//        final ExecutorService executor = Utils.initExecuter(this.nThreads);
//        final ArrayList<Callable<SSCInJSON>> callables = new ArrayList<>();
//        // add all task to do
//        for (final String sscIndex: jsonObject.keySet()) {
//            callables.add((Callable<SSCInJSON>) () -> {
//                return gson.fromJson(jsonObject.get(sscIndex), SSCInJSON.class);
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
//                .forEach((sscInJSON) -> {
//                    if (sscInJSON != null) {
//                        try {
//                            this.insert(sscInJSON.toSSC(), offset + sscInJSON.getIndex());
//                        } catch (CDKException | CloneNotSupportedException ex) {
//                            Logger.getLogger(SSCLibrary.class.getName()).log(Level.SEVERE, null, ex);
//                        }
//                    }
//                });
//        // shut down the executor service
//        Utils.stopExecuter(executor, 5);
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
    
    /**
     * Inserts a new SSC to this SSC library with given SSC index. The index 
     * will then be stored in the new added SSC object.
     *
     * @param ssc SSC to add to this library
     * @param sscIndex SSC index to use/store
     * @return false if this SSC library already contains {@code sscIndex} 
     * or {@code sscIndex < 0}; otherwise true
     */
    public boolean insert(final SSC ssc, final int sscIndex){
        if(this.containsSSC(sscIndex) || (sscIndex < 0)){
            return false;
        }
        ssc.setIndex(sscIndex);
        this.map.put(sscIndex, ssc);        
        
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
            final int nthreads = Integer.parseInt((String) argValues[2]);
            final int maxsphere = Integer.parseInt((String) argValues[3]);
            final int offset = Integer.parseInt((String) argValues[4]);
            System.out.println("-> importing from: " + pathToNMRShiftDB + " ...");
            sscLibrary.setNthreads(nthreads);
            sscLibrary.importFromNMRShiftDB(pathToNMRShiftDB, maxsphere, offset);
            System.out.println(" -> imported! -> #SSC: " + sscLibrary.getSSCCount() + ", next offset: " + (Collections.max(sscLibrary.getSSCIndices()) + 1));
            System.out.println("-> exporting to: " + pathToJSON + " ...");
            sscLibrary.exportToJSON(pathToJSON);
            System.out.println(" -> exported!!!");
            
            
        } catch (IOException | CloneNotSupportedException | InterruptedException | org.apache.commons.cli.ParseException | ParseException | CDKException e) {
            System.err.println(e.getMessage());
        }
    }
    
    
    private Object[] parseArgs(String[] args) throws ParseException, org.apache.commons.cli.ParseException {
        final Options options = setupOptions(args);
        final CommandLineParser parser = new DefaultParser();
        final Object[] argValues = new Object[options.getOptions().size()];
        try {
            final CommandLine cmd = parser.parse(options, args);
            argValues[0] = cmd.getOptionValue("nmrshiftdb");
            argValues[1] = cmd.getOptionValue("library");
            argValues[2] = cmd.getOptionValue("nthreads");
            argValues[3] = cmd.getOptionValue("maxsphere");
            argValues[4] = cmd.getOptionValue("offset");

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

    private Options setupOptions(String[] args) {
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
        Option nthreads = Option.builder("nt")
                .required(true)
                .hasArg()
                .longOpt("nthreads")
                .desc("number of threads to use for parallelization (required)")
                .build();
        options.addOption(nthreads);
        Option maxsphere = Option.builder("m")
                .required(true)
                .hasArg()
                .longOpt("maxsphere")
                .desc("maximum sphere limit for SSC library fragments (required)")
                .build();
        options.addOption(maxsphere);
        Option offset = Option.builder("o")
                .required(true)
                .hasArg()
                .longOpt("offset")
                .desc("offset value for SSC indexing (required)")
                .build();
        options.addOption(offset);

        return options;
    }


}
