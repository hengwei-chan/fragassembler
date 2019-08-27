/*
 * The MIT License
 *
 * Copyright (c) 2019. Michael Wenk [https://github.com/michaelwenk]
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
import casekit.NMR.model.Signal;
import com.google.gson.Gson;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.mongodb.client.FindIterable;
import com.mongodb.client.MongoCollection;
import fragmentation.Fragmentation;
import hose.model.ConnectionTree;
import org.bson.Document;
import org.openscience.cdk.exception.CDKException;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;
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

    private final LinkedHashMap<Long, SSC> map;
    private final HashMap<String, ArrayList<Double>> HOSECodeLookupTableShifts;
    private final HashMap<String, ArrayList<Long>> HOSECodeLookupTableSSCIndices;    
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
        this.map = new LinkedHashMap<>();
        this.HOSECodeLookupTableShifts = new HashMap<>();
        this.HOSECodeLookupTableSSCIndices = new HashMap<>();
        this.nThreads = nThreads;
    }
    
    /**
     * Creates/Updates lookup tables for generated HOSE codes regarding the 
     * given SSC library. 
     *
     * @throws java.lang.InterruptedException
     * 
     * @see #getHOSECodeLookupTableSSCIndices()
     * @see #getHOSECodeLookupTableShifts() 
     */
    public void buildHOSELookupTables() throws InterruptedException {  
        this.HOSECodeLookupTableShifts.clear();
        this.HOSECodeLookupTableSSCIndices.clear();
        String HOSECode;
        Signal signal; 
        for (final SSC ssc : this.map.values()) {
            if(ssc != null){
                HOSECode = ssc.getHOSECode(ssc.getRootAtomIndex());
                if (!this.HOSECodeLookupTableSSCIndices.containsKey(HOSECode)) {
                    this.HOSECodeLookupTableSSCIndices.put(HOSECode, new ArrayList<>());
                    this.HOSECodeLookupTableShifts.put(HOSECode, new ArrayList<>());
                }
                this.HOSECodeLookupTableSSCIndices.get(HOSECode).add(ssc.getIndex());
                signal = ssc.getSubspectrum().getSignal(ssc.getAssignments().getIndex(0, ssc.getRootAtomIndex()));
                if(signal != null){
                    this.HOSECodeLookupTableShifts.get(HOSECode).add(signal.getShift(0));
                }                
            }            
        }
    }
    
    /**
     * Returns a HashMap of the shifts for each HOSE code of that SSC 
     * library. 
     * Before usage of this method, the HOSE code lookup table has to been 
     * built via {@link #buildHOSELookupTables()}.
     *
     * @return HashMap with HOSE codes as keys and lists of chemical shifts as
     * values
     *
     * @see #buildHOSELookupTables() 
     */
    public HashMap<String, ArrayList<Double>> getHOSECodeLookupTableShifts() {
        return this.HOSECodeLookupTableShifts;
    }
    
    /**
     * Builds both HOSE code lookup tables of this SSC library containing a list 
     * of shifts as well as a list of SSC indices for each HOSE code. 
     * Then it will be exported into a given MongoDB collection.      
     *
     * @param collection MongoDB collection to store in
     * @throws java.lang.InterruptedException 
     *
     * @see #buildHOSELookupTables() 
     */
    public void exportHOSECodeLookupTable(final MongoCollection<Document> collection) throws InterruptedException {
        
        this.buildHOSELookupTables();
        // initialize an executor for parallelization
        final ExecutorService executor = Utils.initExecuter(this.nThreads);
        final ArrayList<Callable<Document>> callables = new ArrayList<>();
        // add all task to do
        for (final String HOSECode : this.HOSECodeLookupTableShifts.keySet()) {
            callables.add((Callable<Document>) () -> {
                final Document document = new Document();  
                final ArrayList<Double> shifts = this.HOSECodeLookupTableShifts.get(HOSECode);
                document.append("HOSECode", HOSECode);
                document.append("spheres", hose.Utils.getSpheresCount(HOSECode));
                final Document documentShifts = new Document();
                documentShifts.append("count", shifts.size());
                documentShifts.append("min", (shifts.size() > 0) ? Collections.min(shifts) : null);
                documentShifts.append("max", (shifts.size() > 0) ? Collections.max(shifts) : null);
                documentShifts.append("rms", Utils.getRMS(shifts));
                documentShifts.append("mean", Utils.getMean(shifts));
                documentShifts.append("median", Utils.getMedian(shifts));
                documentShifts.append("var", Utils.getVariance(shifts));
                documentShifts.append("sd", Utils.getStandardDeviation(shifts));
                documentShifts.append("shifts", shifts);
                final Document documentSSCIndices = new Document();
                documentSSCIndices.append("count", this.HOSECodeLookupTableSSCIndices.size());
                documentSSCIndices.append("indices", this.HOSECodeLookupTableSSCIndices.get(HOSECode));
                
                document.append("shift", documentShifts);
                document.append("ssc", documentSSCIndices);

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
    }
    
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
    public HashMap<String, ArrayList<Long>> getHOSECodeLookupTableSSCIndices() {
        return this.HOSECodeLookupTableSSCIndices;
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
        // initialize an executor for parallelization
        final ExecutorService executor = Utils.initExecuter(this.nThreads);
        final ArrayList<Callable<Document>> callables = new ArrayList<>();
        // add all task to do
        long sscCounter = 0;
        for (final SSC ssc : this.getSSCs()) {
            final long sscCounterCopy = sscCounter;
            callables.add((Callable<Document>) () -> {                        
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
                        bw.write(document.toJson());
                        bw.newLine();
                        bw.flush();
                    } catch (IOException ex) {
                        Logger.getLogger(SSCLibrary.class.getName()).log(Level.SEVERE, null, ex);
                    }
                });
        // shut down the executor service
        Utils.stopExecuter(executor, 5);
        
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
    public void importFromJSONFile(final String pathToJSON, final long offset) throws InterruptedException, FileNotFoundException, CDKException, CloneNotSupportedException  {
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
     * @param offset offset number for next SSC indices to use as keys in SSC
     * 
     * @throws FileNotFoundException
     * @throws CDKException
     * @throws CloneNotSupportedException
     * @throws java.lang.InterruptedException
     * 
     * @see #removeAll() 
     * @see #extend(String, String, int, long)
     */
    public void importFromNMRShiftDB(final String pathToNMRShiftDB, final String property, final int maxSphere, final long offset) throws FileNotFoundException, CDKException, CloneNotSupportedException, InterruptedException {
        this.removeAll();
        this.extend(pathToNMRShiftDB, property, maxSphere, offset);
    }
    
    /**
     * Extends this SSC library by a given second one. If a SSC in the given 
     * SSC library exists which index also exists in this SSC library then an 
     * exception will be thrown.
     *
     * @param sscLibrary SSC library to add to this library
     * @throws java.lang.InterruptedException
     * @see #containsSSC(long)
     */
    public void extend(final SSCLibrary sscLibrary) throws InterruptedException{
//        this.map.putAll(sscLibrary.getMap());
        
        // initialize an executor for parallelization
        final ExecutorService executor = Utils.initExecuter(this.nThreads);
        final ArrayList<Callable<SSC>> callables = new ArrayList<>();
        // add all task to do
        for (final SSC ssc : sscLibrary.getSSCs()) {
            callables.add((Callable<SSC>) () -> {
                return ssc;
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
     * @see #containsSSC(long)
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
     * @see #containsSSC(long)
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
     * @see #containsSSC(long)
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
     * Extends this SSC library by created SSCs from a NMRShiftDB file.All SSCs in this library object whose indices also exist in the 
 given map will be replaced.
     *
     * @param pathToNMRShiftDB path to NMRShiftDB SDF
     * @param property property string of spectrum to use, 
     * e.g. "Spectrum 13C 0"
     * @param maxSphere maximum number of spheres in fragmention
     * @param offset offset number for next SSC indices to use as keys in SSC 
     * library
     * @throws java.io.FileNotFoundException
     * @throws org.openscience.cdk.exception.CDKException
     * @throws java.lang.InterruptedException
     * @throws java.lang.CloneNotSupportedException
     * 
     * @see Fragmentation#buildSSCLibrary(HashMap, int, int, long)
     * @see #extend(String, String, int, long)
     */
    public void extend(final String pathToNMRShiftDB, final String property, final int maxSphere, final long offset) throws FileNotFoundException, CDKException, InterruptedException, CloneNotSupportedException {
        this.extend(Fragmentation.buildSSCLibrary(NMRShiftDB.getSSCComponentsFromNMRShiftDB(pathToNMRShiftDB, property), maxSphere, this.nThreads, offset));
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
     * @see #containsSSC(long)
     */
    public void extend(final String pathToJSON, final long offset) throws InterruptedException, FileNotFoundException, CDKException, CloneNotSupportedException{
        final Gson gson = new Gson();
        final BufferedReader br = new BufferedReader(new FileReader(pathToJSON));
        final JsonObject jsonObject = new JsonParser().parse(br).getAsJsonObject();
        // initialize an executor for parallelization
        final ExecutorService executor = Utils.initExecuter(this.nThreads);
        final ArrayList<Callable<SSC>> callables = new ArrayList<>();
        // add all task to do
        for (final String sscIndex: jsonObject.keySet()) {
            callables.add((Callable<SSC>) () -> {
                return SSCConverter.JsonObjectToSSC(jsonObject.getAsJsonObject(sscIndex));
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

    /**
     * Removes all duplicated SSCs based on HOSE code comparisons.
     * If the multiplicities in last sphere differ or shift deviations are higher than a given tolerance value,
     * then two SSCs are not considered to be the same.
     *
     * @param shiftTol shift tolerance value [ppm] in which chemical shifts are considered as the same
     *
     * @throws InterruptedException
     */
    public void removeDuplicates(final double shiftTol) throws InterruptedException {
        // build/update HOSE code lookup tables, i.e. the one with SSC indices
        this.buildHOSELookupTables();
        
        final HashSet<Long> duplicatesSSCIndices = new HashSet<>();   
        SSC ssc1, ssc2;
        Signal signalSSC1, signalSSC2;
        long sscIndex1, sscIndex2;
        ConnectionTree connectionTreeSSC1;
        boolean isDuplicate;
        for (final String HOSECode : this.getHOSECodeLookupTableSSCIndices().keySet()) {
            for (int i = 0; i < this.getHOSECodeLookupTableSSCIndices().get(HOSECode).size(); i++) {
                sscIndex1 = this.getHOSECodeLookupTableSSCIndices().get(HOSECode).get(i);
                ssc1 = this.getSSC(sscIndex1);
                connectionTreeSSC1 = ssc1.getConnectionTree(ssc1.getRootAtomIndex());                
                for (int j = i + 1; j < this.getHOSECodeLookupTableSSCIndices().get(HOSECode).size(); j++) {
                    sscIndex2 = this.getHOSECodeLookupTableSSCIndices().get(HOSECode).get(j);
                    ssc2 = this.getSSC(sscIndex2);
                    isDuplicate = true;
                    // check for same hydrogen count in last sphere, because there it could be different
                    for (final int nodeKey : connectionTreeSSC1.getNodeKeysInSphere(connectionTreeSSC1.getMaxSphere())) {
                        if(connectionTreeSSC1.getNode(nodeKey).isRingClosureNode()){
                            continue;
                        }
                        signalSSC1 = ssc1.getSubspectrum().getSignal(nodeKey);
                        signalSSC2 = ssc2.getSubspectrum().getSignal(nodeKey);
                        if((ssc1.getSubstructure().getAtom(nodeKey).getImplicitHydrogenCount() != ssc2.getSubstructure().getAtom(nodeKey).getImplicitHydrogenCount())
                                || ((signalSSC1 != null) && (signalSSC2 != null) && (Math.abs(signalSSC1.getShift(0) - signalSSC2.getShift(0)) > shiftTol))){
                            isDuplicate = false;
                            break;
                        }
                    }
                    if(isDuplicate){
                        duplicatesSSCIndices.add(sscIndex2);
                    }
                }
            }
        }
        // remove all duplicates
        for (final long duplicatesSSCIndex : duplicatesSSCIndices) {            
            this.remove(duplicatesSSCIndex);
        }
        // update lookup tables
        this.buildHOSELookupTables();
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
        return this.map.size();
    }
    
    public SSC getSSC(final long sscIndex){
        return this.map.get(sscIndex);
    } 
    
    public boolean containsSSC(final long sscIndex){
        return this.map.containsKey(sscIndex);
    }    
    
    public LinkedHashMap<Long, SSC> getMap(){
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
     * @return false if this SSC library already contains {@code ssc.getIndex()} 
     * or {@code ssc.getIndex() < 0}; otherwise true
     */
    public boolean insert(final SSC ssc) {
        if(this.containsSSC(ssc.getIndex()) || ssc.getIndex() < 0){
            return false;
        }
        this.map.put(ssc.getIndex(), ssc);        
        
        return true;
    }
    
    public boolean remove(final long sscIndex){
        if(!this.containsSSC(sscIndex)){
            return false;
        }
        
        return this.map.remove(sscIndex, this.getSSC(sscIndex));
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
    
    public void removeAll(){
        this.map.clear();
    }    

}
