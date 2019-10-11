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
import com.mongodb.client.FindIterable;
import com.mongodb.client.MongoCollection;
import fragmentation.Fragmentation;
import org.bson.Document;
import org.openscience.cdk.exception.CDKException;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
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
     * Creates/Updates lookup tables for generated HOSE codes regarding the 
     * given SSC library. 
     **
     * @see #getHOSECodeLookupTable()
     *
     * @deprecated
     */
    public void buildHOSELookupTables() throws CDKException {
        this.HOSECodeLookupTable.clear();
        String HOSECode;
        for (final SSC ssc : this.map.values()) {
            if(ssc != null){
                HOSECode = ssc.getAsHOSECode();
                if (!this.HOSECodeLookupTable.containsKey(HOSECode)) {
                    this.HOSECodeLookupTable.put(HOSECode, new ArrayList<>());
                }
                this.HOSECodeLookupTable.get(HOSECode).add(ssc.getIndex());
            }            
        }
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
    
    /**
     * Returns a HashMap of the SSC indices for each HOSE code of that SSC 
     * library. 
     * Before usage of this method, the HOSE code lookup table has to been 
     * built via {@link #buildHOSELookupTables()}.
     *
     * @return HashMap with HOSE codes as keys and lists of SSC indices as
     * values
     *
     * @see #buildHOSELookupTables() 
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
    
    /**
     * Writes this SSC library into a MongoDB collection.
     *
     * @param collection
     * @throws java.lang.InterruptedException
     *
     * @deprecated
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
    
    public void exportToJSONFile(final String pathToJSON) throws IOException, InterruptedException {                        
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
    }
    
    /**
     * Imports a SSC library from a JSON file. All current SSCs will be 
     * deleted.
     *
     * @param pathToJSON JSON file containing ths SSC library
     * @throws java.lang.InterruptedException
     * @throws java.io.FileNotFoundException
     */
    public void importFromJSONFile(final String pathToJSON) throws InterruptedException, FileNotFoundException  {
        this.removeAll();
        this.extend(pathToJSON);
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
     * Imports a SSC library by documents from a MongoDB query result 
     * containing the SSC information.  All current SSCs will be 
     * deleted.
     *
     * @param queryResult 
     * @throws org.openscience.cdk.exception.CDKException
     * @throws java.lang.CloneNotSupportedException
     * @throws java.lang.InterruptedException
     * 
     */
    public void importFromMongoDB(final FindIterable<Document> queryResult) throws CDKException, CloneNotSupportedException, InterruptedException  {
        this.removeAll();
        this.extend(queryResult);
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
     */
    public void importFromMongoDB(final ArrayList<Document> documents) throws CDKException, CloneNotSupportedException, InterruptedException  {
        this.removeAll();
        this.extend(documents);
    }
    
    /**
     * Creates a SSC library from a given NMRShiftDB SDF via fragmentation of
     * each structure and each of its atoms as individual starting point 
     * (central atom).All current SSCs will be deleted.
     *
     * @param pathToNMRShiftDB path to NMRShiftDB SDF
     * @param property property string of spectrum to use, 
     * e.g. "Spectrum 13C 0"
     * @param maxSphere maximum number of spheres in fragmention
     *
     * @throws FileNotFoundException
     * @throws CDKException
     * @throws CloneNotSupportedException
     * @throws java.lang.InterruptedException
     * 
     * @see #removeAll() 
     * @see #extend(String, String, int)
     */
    public void importFromNMRShiftDB(final String pathToNMRShiftDB, final String property, final int maxSphere) throws FileNotFoundException, CDKException, CloneNotSupportedException, InterruptedException {
        this.removeAll();
        this.extend(pathToNMRShiftDB, property, maxSphere);
    }
    
    /**
     * Extends this SSC library by a given second one.
     *
     * @param sscLibrary SSC library to add to this library
     *
     * @throws java.lang.InterruptedException
     */
    public void extend(final SSCLibrary sscLibrary) throws InterruptedException{

//        // initialize an executor for parallelization
//        final ExecutorService executor = Utils.initExecuter(this.nThreads);
//        final ArrayList<Callable<SSC>> callables = new ArrayList<>();
//        // add all task to do
//        for (final SSC ssc : sscLibrary.getSSCs()) {
//            callables.add((Callable<SSC>) () -> {
//                return ssc;
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


        for (final SSC ssc : sscLibrary.getSSCs()) {
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

//        sscLibrary.getSSCs().parallelStream().forEach(ssc -> {
//            if (!this.insert(ssc)) {
//                try {
//                    throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": insertion SSC with index " + ssc.getIndex() + " failed");
//                } catch (CDKException ex) {
//                    Logger.getLogger(SSCLibrary.class.getName()).log(Level.SEVERE, null, ex);
//                }
//            }
//        });
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
            callables.add((Callable<SSC>) () -> {                
                return SSCConverter.DocumentToSSC(sscDocument);
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
     * Extends this SSC library by documents from a MongoDB query result 
     * containing the SSC information. 
     * All SSCs in this library object whose indices also exist in the 
     * given map will be replaced.
     *
     * @param queryResult 
     * @throws org.openscience.cdk.exception.CDKException  
     * @throws java.lang.CloneNotSupportedException  
     * @throws java.lang.InterruptedException  
     *
     * @deprecated
     */
    public void extend(final FindIterable<Document> queryResult) throws CDKException, CloneNotSupportedException, InterruptedException  {
        // initialize an executor for parallelization
        final ExecutorService executor = Utils.initExecuter(this.nThreads);
        final ArrayList<Callable<SSC>> callables = new ArrayList<>();
        // add all task to do
        for (final Document sscDocument : queryResult) {
            callables.add((Callable<SSC>) () -> {
                return SSCConverter.DocumentToSSC(sscDocument);
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
     * @deprecated
     */
    public void extend(final ArrayList<Document> documents) throws CDKException, CloneNotSupportedException, InterruptedException  {
        // initialize an executor for parallelization
        final ExecutorService executor = Utils.initExecuter(this.nThreads);
        final ArrayList<Callable<SSC>> callables = new ArrayList<>();
        // add all task to do
        for (final Document sscDocument : documents) {
            callables.add((Callable<SSC>) () -> {
                return SSCConverter.DocumentToSSC(sscDocument);
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
     * Extends this SSC library by created SSCs from a NMRShiftDB file.
     *
     * @param pathToNMRShiftDB path to NMRShiftDB SDF
     * @param property property string of spectrum to use, 
     * e.g. "Spectrum 13C 0"
     * @param maxSphere maximum number of spheres in fragmention
     * library
     * @throws java.io.FileNotFoundException
     * @throws org.openscience.cdk.exception.CDKException
     * @throws java.lang.InterruptedException
     * @throws java.lang.CloneNotSupportedException
     * 
     * @see Fragmentation#buildSSCLibrary(HashMap, int, int, boolean)
     * @see #extend(String, String, int)
     */
    public void extend(final String pathToNMRShiftDB, final String property, final int maxSphere) throws FileNotFoundException, CDKException, InterruptedException, CloneNotSupportedException {
        this.extend(Fragmentation.buildSSCLibrary(NMRShiftDB.getSSCComponentsFromNMRShiftDB(pathToNMRShiftDB, property), maxSphere, this.nThreads, true));
    }
    
    /**
     * Extends this SSC library by SSCs stored in a JSON file.
     *
     * @param pathToJSON path to JSON file containing SSCs
     * library
     *
     * @throws java.io.FileNotFoundException
     */
    public void extend(final String pathToJSON) throws FileNotFoundException {
        final BufferedReader br = new BufferedReader(new FileReader(pathToJSON));

        final JsonParser jsonParser = new JsonParser();
//        br.lines().parallel().forEach( line -> {
        br.lines().forEach( line -> {
            if(!line.startsWith("{") && !line.endsWith("}")){
                String asJSON = "";
                if(line.endsWith(",")){
                    asJSON += line.substring(0,line.length() - 1);
                } else {
                    asJSON += line;
                }
                asJSON = asJSON.substring(asJSON.indexOf("{"));
                try {
                    SSC ssc = SSCConverter.JSONToSSC(jsonParser.parse(asJSON).getAsJsonObject().toString());
                    if(!this.insert(ssc)){
                        try {
                            throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": insertion SSC with index " + ssc.getIndex() + " failed");
                        } catch (CDKException ex) {
                            Logger.getLogger(SSCLibrary.class.getName()).log(Level.SEVERE, null, ex);
                        }
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        });
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
     * Returns a deep clone if this SSC library and its SSCs object. The set number of threads is 
     * set too.
     *
     * @return
     * @throws CDKException
     * @throws CloneNotSupportedException
     * @throws InterruptedException
     * 
     * @see #setNThreads(int) 
     */
    public SSCLibrary getClone() throws Exception {
        final SSCLibrary sscLibrary = new SSCLibrary(this.nThreads);
        for (final long sscIndex : this.map.keySet()) {
            sscLibrary.insert(this.getSSC(sscIndex).getClone());
        }
        
        return sscLibrary;
    }

}
