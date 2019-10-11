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
package model;

import casekit.NMR.Utils;
import casekit.NMR.dbservice.NMRShiftDB;
import com.google.gson.JsonParser;
import com.mongodb.client.MongoCollection;
import fragmentation.Fragmentation;
import org.bson.Document;
import org.openscience.cdk.exception.CDKException;
import start.TimeMeasurement;

import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Class to maintain a library of SSC objects.
 * 
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public final class SSCLibrary {

    private final ConcurrentHashMap<Long, SSC> map;
    private final ConcurrentHashMap<String, ArrayList<Long>> HOSECodeLookupTable;
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
        this.map = new ConcurrentHashMap<>();
        this.HOSECodeLookupTable = new ConcurrentHashMap<>();
        this.nThreads = nThreads;
    }


    
    /**
     * Returns a HashMap of the SSC indices for each HOSE code of that SSC 
     * library.
     *
     * @return HashMap with HOSE codes as keys and lists of SSC indices as
     * values
     *
     */
    public ConcurrentHashMap<String, ArrayList<Long>> getHOSECodeLookupTable() {
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
    
    public void exportToJSONFile(final String pathToJSON) throws IOException {
        final TimeMeasurement tm = new TimeMeasurement();
        System.out.println("Now storing SSC library into JSON file \"" + pathToJSON + "\"...");
        tm.start();

        final BufferedWriter bw = new BufferedWriter(new FileWriter(pathToJSON));
        long sscCounter = 0;
        String json;
        bw.write("{");
        bw.newLine();
        for (final SSC ssc : this.getSSCs()) {
            json = SSCConverter.SSCToJSON(ssc, sscCounter);
            bw.write(json.substring(1, json.length() - 1));
            if(sscCounter < this.getSSCCount() - 1){
                bw.write(",");
            }
            bw.newLine();
            bw.flush();

            sscCounter++;
        }

        bw.write("}");
        bw.close();

        tm.stop();
        System.out.println("--> time needed: " + tm.getResult() + " s");
        System.out.println("-> SSC library stored into JSON file");
    }
    
    /**
     * Imports a SSC library from a JSON file. All current SSCs will be 
     * deleted.
     *
     * @param pathToJSON JSON file containing ths SSC library
     *
     * @throws java.lang.InterruptedException
     * @throws java.io.FileNotFoundException
     * 
     * @see #removeAll() 
     * @see #extend(String)
     */
    public void buildFromJSON(final String pathToJSON) throws InterruptedException, FileNotFoundException  {
        this.removeAll();

        final TimeMeasurement tm = new TimeMeasurement();
        System.out.println("-> importing SSC library from JSON file...");
        tm.start();

        this.extend(pathToJSON);

        System.out.println("-> SSC library imported from JSON file!!!");
        tm.stop();
        System.out.println("--> time needed: " + tm.getResult() + " s");
        System.out.println("--> SSC library size:\t" + this.getSSCCount());
    }
    
    /**
     * Creates a SSC library from a given NMRShiftDB SDF via fragmentation of
     * each structure and each of its atoms as individual starting point 
     * (central atom). All current SSCs will be deleted.
     * Signals shifts outlier in each SSC will be removed, see {@link SSC#removeSignalShiftsOutlier()}.
     *
     * @param pathToNMRShiftDB path to NMRShiftDB SDF
     * @param property property string of spectrum to use, 
     * e.g. "Spectrum 13C 0"
     * @param maxSphere maximum number of spheres in fragmention
     *
     * @throws FileNotFoundException
     * @throws CDKException
     * @throws java.lang.InterruptedException
     * 
     * @see #removeAll() 
     * @see Fragmentation#buildSSCLibrary(HashMap, int, int)
     * @see #extend(SSCLibrary)
     * @see #removeSignalShiftsOutlier()
     */
    public void buildFromNMRShiftDB(final String pathToNMRShiftDB, final String property, final int maxSphere) throws FileNotFoundException, CDKException, InterruptedException {
        this.removeAll();

        final TimeMeasurement tm = new TimeMeasurement();
        for (int m = 2; m <= maxSphere; m++) {
            System.out.println("Building SSCs for " + m + "-spheres...");
            tm.start();
            this.extend(Fragmentation.buildSSCLibrary(NMRShiftDB.getSSCComponentsFromNMRShiftDB(pathToNMRShiftDB, property), m, this.nThreads));
            System.out.println("SSCs for " + m + "-spheres build!!!");
            tm.stop();
            System.out.println("--> time needed: " + tm.getResult() + " s");
            System.out.println("-> #SSCs in SSC library: " + this.getSSCCount());
        }

        System.out.println("Building SSCs done!!!");
        System.out.println("Removing signal shifts oulier");
        tm.start();
        this.removeSignalShiftsOutlier();
        tm.stop();
        System.out.println("--> time needed: " + tm.getResult() + " s");
    }
    
    /**
     * Extends this SSC library by all SSCs from given second one.
     *
     * @param sscLibrary SSC library with SSCs to add to this library
     *
     */
    public void extend(final SSCLibrary sscLibrary) {
        this.extend(sscLibrary.getSSCs());
    }

    public void extend(final Collection<SSC> sscCollection){
        for (final SSC ssc : sscCollection) {
            if (ssc != null) {
                if (!this.insert(ssc)) {
                    try {
                        throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": insertion SSC with index " + ssc.getIndex() + " failed");
                    } catch (CDKException ex) {
                        Logger.getLogger(SSCLibrary.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
            }
        }
    }
    
    /**
     * Extends this SSC library by SSCs stored in a JSON file.
     *
     * @param pathToJSON path to JSON file containing SSCs
     * library
     *
     * @throws java.io.FileNotFoundException
     */
    public void extend(final String pathToJSON) throws FileNotFoundException, InterruptedException {
        final ConcurrentLinkedQueue<SSC> convertedSSCs = new ConcurrentLinkedQueue<>();
        final BufferedReader br = new BufferedReader(new FileReader(pathToJSON));
        final JsonParser jsonParser = new JsonParser();
        // initialize an executor
        final ExecutorService executor = Utils.initExecuter(nThreads);
        final ArrayList<Callable<SSC>> callables = new ArrayList<>();
        // add all task to do
        br.lines().forEach( line -> {
            if((line.trim().length() > 1) || (!line.trim().startsWith("{") && !line.trim().endsWith("}"))){
                final StringBuilder sscInJSON = new StringBuilder();
                if(line.endsWith(",")){
                    sscInJSON.append(line, 0, line.length() - 1);
                } else {
                    sscInJSON.append(line);
                }
                try {
                    callables.add(() -> SSCConverter.JSONToSSC(jsonParser.parse(sscInJSON.substring(sscInJSON.toString().indexOf("{"))).getAsJsonObject().toString()));
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        });
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
                .forEach(convertedSSCs::add);
        // shut down the executor service
        Utils.stopExecuter(executor, 5);

        this.extend(convertedSSCs);
    }

    /**
     * Checks whether a SSC with same structural properties already exists.
     * This is based on HOSE code comparisons and hydrogen counts in last sphere. </br>
     * If the attached hydrogens in last sphere differ, then two SSCs are not considered to be the same.
     *
     * @param ssc
     * @return true if there already exists a SSC with same structural properties
     *
     */
    public boolean containsSSC(final SSC ssc) {
        return this.findSSC(ssc) != null;
    }

    /**
     * Returns a SSC with same structural properties which already exists.
     * This is based on HOSE code comparisons and hydrogen counts in last sphere. </br>
     * If the attached hydrogens in last sphere differ, then two SSCs are not considered to be the same.
     *
     * @param ssc
     * @return index in this SSC library if there is already a SSC entry, otherwise null
     *
     */
    public Long findSSC(final SSC ssc) {
        if(!this.HOSECodeLookupTable.containsKey(ssc.getAsHOSECode())){
            return null;
        }
        for (final long sscIndexInLibrary : this.getHOSECodeLookupTable().get(ssc.getAsHOSECode())){
            if(Arrays.equals(this.getSSC(sscIndexInLibrary).getAttachedHydrogensInOuterSphere(), ssc.getAttachedHydrogensInOuterSphere())){
                return sscIndexInLibrary;
            }
        }

        return null;
    }
    
    public boolean isEmpty(){
        return this.getSSCCount() == 0;
    }
    
    /**
     * Returns the last (numerical highest) SSC index of this SSC library. 
     * If this SSC library is empty then null will be returned.
     *
     * @return
     * 
     * @see #isEmpty() 
     */
    public Long getLastSSCIndex(){
        if(this.isEmpty()){
            return null;
        }
        
        return Collections.max(this.getSSCIndices());
    }
    
    public long getSSCCount(){
        return this.map.mappingCount();
    }
    
    public SSC getSSC(final long sscIndex){
        return this.map.get(sscIndex);
    } 
    
    public boolean containsSSC(final long sscIndex){
        return this.map.containsKey(sscIndex);
    }    
    
    public ConcurrentHashMap<Long, SSC> getMap(){
        return this.map;
    }
    
    /**
     * Returns the SSC indices of this SSC library.
     *
     * @return
     */
    public LinkedHashSet<Long> getSSCIndices(){
        return new LinkedHashSet<>(this.map.keySet());
    }
    
    /**
     * Returns the SSCs of this SSC library in same order as the SSCs 
     * were inserted.
     *
     * @return
     */
    public Collection<SSC> getSSCs(){
        return this.map.values();
    }
    
    /**
     * Inserts a new SSC to this SSC library.
     *
     * @param ssc SSC to add to this library
     * @return false if the given SSC could not be cloned successfully
     * 
     * @see #findSSC(SSC)
     */
    public boolean insert(final SSC ssc) {
        if(ssc == null){
            return false;
        }
        final Long sscIndexInLibrary = this.findSSC(ssc);
        if(sscIndexInLibrary != null){
            this.getSSC(sscIndexInLibrary).addShiftsFromSubspectrum(ssc.getSubspectrum());

            return true;
        }
        final String HOSECode = ssc.getAsHOSECode();
        if (!this.HOSECodeLookupTable.containsKey(HOSECode)) {
            this.HOSECodeLookupTable.put(HOSECode, new ArrayList<>());
        }
        final SSC sscClone;
        try {
            sscClone = ssc.getClone();
        } catch (Exception e) {
            return false;
        }
        sscClone.setIndex(this.getSSCCount());
        this.HOSECodeLookupTable.get(HOSECode).add(sscClone.getIndex());
        this.map.put(sscClone.getIndex(), sscClone);
        
        return true;
    }

    /**
     * Removes shifts for each signal of a SSC subspectrum which are considered as outlier ({@link SSC#removeSignalShiftsOutlier()}).
     */
    public void removeSignalShiftsOutlier(){
        for (final SSC ssc : this.getSSCs()){
            ssc.removeSignalShiftsOutlier();
        }
    }
    
    public boolean remove(final long sscIndex){
        if(!this.containsSSC(sscIndex)){
            return false;
        }

        this.HOSECodeLookupTable.get(this.getSSC(sscIndex).getAsHOSECode()).remove(sscIndex);
        return this.map.remove(sscIndex, this.getSSC(sscIndex));
    }

    public void removeAll(){
        this.HOSECodeLookupTable.clear();
        this.map.clear();
    }

    /**
     * Returns a deep clone if this SSC library and its SSCs object. The number of threads is
     * set too the number of threads of this class object.
     *
     * @return
     * @throws CDKException
     * @throws CloneNotSupportedException
     * @throws InterruptedException
     *
     */
    public SSCLibrary getClone() throws Exception {
        final SSCLibrary sscLibrary = new SSCLibrary(this.nThreads);
        for (final long sscIndex : this.map.keySet()) {
            sscLibrary.insert(this.getSSC(sscIndex).getClone());
        }
        
        return sscLibrary;
    }





    /**
     * Writes this SSC library into a MongoDB collection.
     *
     * @param collection
     * @throws java.lang.InterruptedException
     *
     */
    public void exportToMongoDB(final MongoCollection<Document> collection) throws InterruptedException{

        // initialize an executor for parallelization
        final ExecutorService executor = Utils.initExecuter(this.nThreads);
        final ArrayList<Callable<Document>> callables = new ArrayList<>();
        // add all task to do
        long sscCounter = 0;
        for (final SSC ssc : this.getSSCs()) {
            final long sscCounterCopy = sscCounter;
            callables.add(() -> {
                return SSCConverter.SSCToDocument(ssc, sscCounterCopy);
            });
            sscCounter++;
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
     * @deprecated
     */
    public void buildFromMongoDB(final MongoCollection<Document> collection) throws CDKException, CloneNotSupportedException, InterruptedException  {
        this.removeAll();
        this.extend(collection);
    }

    /**
     * Imports a SSC library by documents from a MongoDB
     * containing the SSC information.  All current SSCs will be
     * deleted.
     *
     * @param documents
     * @throws org.openscience.cdk.exception.CDKException
     * @throws java.lang.CloneNotSupportedException
     * @throws java.lang.InterruptedException
     *
     * @deprecated
     */
    public void buildFromMongoDB(final ArrayList<Document> documents) throws CDKException, CloneNotSupportedException, InterruptedException  {
        this.removeAll();
        this.extend(documents);
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
     * @deprecated
     */
    public void extend(final MongoCollection<Document> collection) throws CDKException, CloneNotSupportedException, InterruptedException  {
        // initialize an executor for parallelization
        final ExecutorService executor = Utils.initExecuter(this.nThreads);
        final ArrayList<Callable<SSC>> callables = new ArrayList<>();
        // add all task to do
        for (final Document sscDocument : collection.find()) {
            callables.add(() -> SSCConverter.DocumentToSSC(sscDocument));
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
                        ssc.setIndex(ssc.getIndex());
                        if (!this.insert(ssc)) {
                            try {
                                throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": insertion SSC with index " + ssc.getIndex() + " failed");
                            } catch (CDKException ex) {
                                Logger.getLogger(SSCLibrary.class.getName()).log(Level.SEVERE, null, ex);
                            }
                        }
                    }
                });
        // shut down the executor service
        Utils.stopExecuter(executor, 5);
    }

    /**
     * Extends this SSC library by documents from a MongoDB
     * containing the SSC information.
     * All SSCs in this library object whose indices also exist in the
     * given map will be replaced.
     *
     * @param documents
     * @throws org.openscience.cdk.exception.CDKException
     * @throws java.lang.CloneNotSupportedException
     * @throws java.lang.InterruptedException
     *
     */
    public void extend(final ArrayList<Document> documents) throws CDKException, CloneNotSupportedException, InterruptedException  {
        // initialize an executor for parallelization
        final ExecutorService executor = Utils.initExecuter(this.nThreads);
        final ArrayList<Callable<SSC>> callables = new ArrayList<>();
        // add all task to do
        for (final Document sscDocument : documents) {
            callables.add(() -> SSCConverter.DocumentToSSC(sscDocument));
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
                        if (!this.insert(ssc)) {
                            try {
                                throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": insertion SSC with index " + ssc.getIndex() + " failed");
                            } catch (CDKException ex) {
                                Logger.getLogger(SSCLibrary.class.getName()).log(Level.SEVERE, null, ex);
                            }
                        }
                    }
                });
        // shut down the executor service
        Utils.stopExecuter(executor, 5);
    }

    //    /**
//     * Builds both HOSE code lookup tables of this SSC library containing a list
//     * of shifts as well as a list of SSC indices for each HOSE code.
//     * Then it will be exported into a given MongoDB collection.
//     *
//     * @param collection MongoDB collection to store in
//     * @throws java.lang.InterruptedException
//     *
//     * @see #buildHOSELookupTables()
//     */
//    public void exportHOSECodeLookupTable(final MongoCollection<Document> collection) throws InterruptedException {
//
//        this.buildHOSELookupTables();
//        // initialize an executor for parallelization
//        final ExecutorService executor = Utils.initExecuter(this.nThreads);
//        final ArrayList<Callable<Document>> callables = new ArrayList<>();
//        // add all task to do
//        for (final String HOSECode : this.HOSECodeLookupTableShifts.keySet()) {
//            callables.add((Callable<Document>) () -> {
//                final Document document = new Document();
//                final ArrayList<Double> shifts = this.HOSECodeLookupTableShifts.get(HOSECode);
//                document.append("HOSECode", HOSECode);
//                document.append("spheres", hose.Utils.getSpheresCount(HOSECode));
//                final Document documentShifts = new Document();
//                documentShifts.append("count", shifts.size());
//                documentShifts.append("min", (shifts.size() > 0) ? Collections.min(shifts) : null);
//                documentShifts.append("max", (shifts.size() > 0) ? Collections.max(shifts) : null);
//                documentShifts.append("rms", Utils.getRMS(shifts));
//                documentShifts.append("mean", Utils.getMean(shifts));
//                documentShifts.append("median", Utils.getMedian(shifts));
//                documentShifts.append("var", Utils.getVariance(shifts));
//                documentShifts.append("sd", Utils.getStandardDeviation(shifts));
//                documentShifts.append("shifts", shifts);
//                final Document documentSSCIndices = new Document();
//                documentSSCIndices.append("count", this.HOSECodeLookupTable.size());
//                documentSSCIndices.append("indices", this.HOSECodeLookupTable.get(HOSECode));
//
//                document.append("shift", documentShifts);
//                document.append("ssc", documentSSCIndices);
//
//                return document;
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
//                .forEach((document) -> {
//                    try {
//                        collection.insertOne(document);
//                    } catch (Exception e) {
//                        System.err.println("export for key \"" + document.get("_id") + "\" failed: " + e.getMessage());
//                        System.out.println("-> document: \n" + document.toJson());
//                    }
//                });
//        // shut down the executor service
//        Utils.stopExecuter(executor, 5);
//    }
}
