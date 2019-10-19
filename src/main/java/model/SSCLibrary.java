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

import casekit.NMR.dbservice.NMRShiftDB;
import casekit.NMR.model.Signal;
import com.mongodb.client.MongoCollection;
import conversion.Converter;
import fragmentation.Fragmentation;
import org.bson.Document;
import org.openscience.cdk.exception.CDKException;
import start.TimeMeasurement;

import java.io.FileNotFoundException;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Class to maintain a library of SSC objects.
 * 
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public final class SSCLibrary {

    /**
     * Holds all SSC entries as LinkedHashMap to keep them in input order. <br>
     * That is important for holding ranked SSCs in order in SSCRanker class.
     */
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
     * Instantiates a new object of this class.
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
    public HashMap<String, ArrayList<Long>> getHOSECodeLookupTable() {
        return new HashMap<>(this.HOSECodeLookupTable);
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
     * Stores this SSC library in JSON format into a file.
     *
     * @param pathToJSONFile
     *
     * @see Converter#SSCLibraryToJSONFile(SSCLibrary, String)
     */
    public boolean exportToJSONFile(final String pathToJSONFile) {
        final TimeMeasurement tm = new TimeMeasurement();
        System.out.println("Now storing SSC library into JSON file \"" + pathToJSONFile + "\"...");
        tm.start();
        if(!Converter.SSCLibraryToJSONFile(this, pathToJSONFile)){
            return false;
        }
        tm.stop();
        System.out.println("--> time needed: " + tm.getResult() + " s");
        System.out.println("-> SSC library stored into JSON file");

        return true;

//        return Converter.SSCLibraryToJSONFile(this, pathToJSONFile);
    }


    
    /**
     * Imports a SSC library from a JSON file. All current SSC will be
     * deleted.
     *
     * @param pathToJSONFile JSON file containing ths SSC library
     *
     * @see #removeAll() 
     * @see #extend(String)
     */
    public boolean importFromJSON(final String pathToJSONFile) {
        this.removeAll();

        final TimeMeasurement tm = new TimeMeasurement();
        System.out.println("-> importing SSC library from JSON file...");
        tm.start();
        if(!this.extend(pathToJSONFile)){
            tm.stop();
            return false;
        }
        System.out.println("-> SSC library imported from JSON file!!!");
        tm.stop();
        System.out.println("--> time needed: " + tm.getResult() + " s");
        System.out.println("--> SSC library size:\t" + this.getSSCCount());

        return true;
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
     * @param maxSphere maximum number of spheres in fragmentation
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
            System.out.println("Building SSC for " + m + "-spheres...");
            tm.start();

            this.extend(Fragmentation.buildSSCLibrary(NMRShiftDB.getSSCComponentsFromNMRShiftDB(pathToNMRShiftDB, property), m, this.nThreads));

            System.out.println("SSC for " + m + "-spheres build!!!");
            tm.stop();
            System.out.println("--> time needed: " + tm.getResult() + " s");
            System.out.println("-> #SSC in SSC library: " + this.getSSCCount());
        }
        System.out.println("Building SSC done!!!");
        System.out.println("Removing signal shifts outlier");
        tm.start();

        this.removeSignalShiftsOutlier();

        tm.stop();
        System.out.println("--> time needed: " + tm.getResult() + " s");
    }
    
    /**
     * Extends this SSC library by all SSC from given second one.
     *
     * @param sscLibrary SSC library with SSCs to add to this library
     *
     * @see #getSSCs()
     * @see #extend(Collection)
     */
    public boolean extend(final SSCLibrary sscLibrary) {
        return this.extend(sscLibrary.getSSCs());
    }

    /**
     * Extends this SSC library by all SSC from given second one.
     *
     * @param collection collection with SSCs to add to this library
     * @return true if all SSCs in collection could be added
     * 
     * @see #insert(SSC)
     */
    private boolean extendBySSCs(final Collection<SSC> collection){
        // no further parallelization needed because insert() method will block anyway
        for (final SSC ssc : collection) {
            if (ssc != null) {
                if (!this.insert(ssc)) {
                    System.err.println(Thread.currentThread().getStackTrace()[1].getMethodName() + ": insertion SSC with index " + ssc.getIndex() + " failed");
                    try {
                        throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": insertion SSC with index " + ssc.getIndex() + " failed");
                    } catch (CDKException ex) {
                        Logger.getLogger(SSCLibrary.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
            }
        }

        return true;
    }
    
    /**
     * Extends this SSC library by SSC stored in a JSON file.
     *
     * @param pathToJSONFile path to JSON file containing SSC
     * library
     *
     * @return
     */
    public boolean extend(final String pathToJSONFile) {
        try {
            return this.extend(Converter.JSONFileToSSCLibrary(pathToJSONFile, this.nThreads));
        } catch (InterruptedException | FileNotFoundException e) {
            e.printStackTrace();
            return false;
        }
    }

    /**
     * Checks whether a SSC with same structural properties already exists.
     *
     * @param ssc
     * @return true if there already exists a SSC with same structural properties
     *
     * @see #findSSC(SSC)
     */
    public boolean containsSSC(final SSC ssc) {
        return this.findSSC(ssc) != null;
    }

    /**
     * Returns a SSC with same structural properties which already exists.
     * This is based on HOSE code comparisons and hydrogen counts in last sphere as well as multiplicities. </br>
     * If the attached hydrogens in last sphere or the multiplicities in both subspectra differ,
     * then two SSC are not considered to be the same.
     *
     * @param ssc
     * @return index in this SSC library if there is already a SSC entry, otherwise null
     *
     */
    public Long findSSC(final SSC ssc) {
        if(!this.HOSECodeLookupTable.containsKey(ssc.getAsHOSECode())){
            return null;
        }
        SSC sscInLibrary;
        Signal signalSSC, signalSSCInLibrary;
        boolean valid;
        for (final long sscIndexInLibrary : this.getHOSECodeLookupTable().get(ssc.getAsHOSECode())){
            sscInLibrary = this.getSSC(sscIndexInLibrary);
            // as pre-filter
            if(Arrays.equals(sscInLibrary.getAttachedHydrogensCountsInOuterSphere(), ssc.getAttachedHydrogensCountsInOuterSphere())
                    && Arrays.equals(sscInLibrary.getMultiplicitiesInOuterSphere(), ssc.getMultiplicitiesInOuterSphere())){

                valid = true;
                for (int i = 0; i < ssc.getAtomCount(); i++) {
                    signalSSC = ssc.getSubspectrum().getSignal(ssc.getAssignments().getIndex(0, i));
                    signalSSCInLibrary = sscInLibrary.getSubspectrum().getSignal(sscInLibrary.getAssignments().getIndex(0, i));
                    // if both signals exist and the multiplicities are not equal
                    if((signalSSC != null) && (signalSSCInLibrary != null)
                            && !signalSSC.getMultiplicity().equals(signalSSCInLibrary.getMultiplicity())){
                        valid = false;
                        break;
                    } else {
                        // if one signal exists and the other not then something is wrong
                        if(((signalSSC == null) && (signalSSCInLibrary != null)) || ((signalSSC != null) && (signalSSCInLibrary == null))){
                            valid = false;
                            break;
                        }
                        // else: atoms which have no signal
                    }
                }
                if(valid){
                    return sscIndexInLibrary;
                }
            }
        }

        return null;
    }
    
    public boolean isEmpty(){
        return this.getSSCCount() == 0;
    }
    
    public long getSSCCount(){
        return this.map.mappingCount();
    }

    public SSC getSSC(final long sscIndex){
        return this.getSSC(sscIndex);
    }
    
    public boolean containsSSC(final long sscIndex){
        return this.map.containsKey(sscIndex);
    }    

    /**
     * Returns the SSC indices of this SSC library (in order).
     *
     * @return
     */
    public LinkedHashSet<Long> getSSCIndices(){
        return new LinkedHashSet<>(this.map.keySet());
    }
    
    /**
     * Returns the SSC of this SSC library in same order as the SSC
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
//     * @param clone indicates whether SSCs to insert should be cloned beforehand ({@link model.SSC#getClone(boolean)})
//     * @param reduceShifts used when {code clone == true} ; indicates whether the shift lists of each signal to be
//     *                     reduced to the representative signal shift ({@link model.SSC#getClone(boolean)})
//     *
//     * @return false if the given SSC could not be cloned successfully
     * @return false if given SSC is null
     * 
     * @see #findSSC(SSC)
     * @see SSC#getClone(boolean)
     */
    public boolean insert(final SSC ssc){//}, final boolean clone, final boolean reduceShifts) {
        if(ssc == null){
            return false;
        }
        final Long sscIndexInLibrary = this.findSSC(ssc);
        if(sscIndexInLibrary != null){
            if(!this.getSSC(sscIndexInLibrary).addShiftsFromSubspectrum(ssc.getSubspectrum())){
                System.err.println(" -> could not add subspectrum to already existing SSC with HOSE code: " + ssc.getAsHOSECode());
                return false;
            }

            return true;
        }
        final String HOSECode = ssc.getAsHOSECode();
        this.HOSECodeLookupTable.put(HOSECode, new ArrayList<>());

//        final SSC sscToInsert;
//        if(clone){
//            try {
//                sscToInsert = ssc.getClone(reduceShifts);
//            } catch (Exception e) {
//                return false;
//            }
//        } else {
//            sscToInsert = ssc;
//        }
//        sscToInsert.setIndex(this.getSSCCount());
//        this.HOSECodeLookupTable.get(HOSECode).add(sscToInsert.getIndex());
//        this.map.put(sscToInsert.getIndex(), sscToInsert);

        ssc.setIndex(this.getSSCCount());
        this.HOSECodeLookupTable.get(HOSECode).add(ssc.getIndex());
        this.map.put(ssc.getIndex(), ssc);

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
     * Returns a deep clone if this SSC library and its SSC object. The number of threads is
     * set too the number of threads of this class object.
     *
     * @return
     */
    public SSCLibrary getClone() throws Exception {
        final SSCLibrary sscLibrary = new SSCLibrary(this.nThreads);
        for (final long sscIndex : this.map.keySet()) {
            sscLibrary.insert(this.getSSC(sscIndex).getClone(false));
        }
        
        return sscLibrary;
    }

    /**
     * Inserts this SSC library into a MongoDB collection.
     *
     * @param collection
     *
     */
    public boolean exportToMongoDB(final MongoCollection<Document> collection) {
        final List<Document> convertedDocuments = Collections.synchronizedList(new ArrayList<>());
        // add all task to do
        final ArrayList<Callable<Document>> callables = new ArrayList<>();
        for (final SSC ssc : this.getSSCs()) {
            callables.add(() -> Converter.SSCToDocument(ssc));
        }
        try {
            Converter.convertSSCsToDocuments(callables, convertedDocuments, this.nThreads);
            collection.insertMany(convertedDocuments);

            return true;
        } catch (InterruptedException e) {
            e.printStackTrace();
            return false;
        }
    }

    /**
     * Imports a SSC library by documents from a MongoDB collection
     * containing the SSC information.  All current SSC will be
     * deleted.
     *
     * @param collection
     *
     */
    public boolean importFromMongoDB(final MongoCollection<Document> collection) {
        this.removeAll();

        final TimeMeasurement tm = new TimeMeasurement();
        System.out.println("-> importing SSC library from JSON file...");
        tm.start();
        if(!this.extend(collection)){
            tm.stop();
            return false;
        }
        System.out.println("-> SSC library imported from JSON file!!!");
        tm.stop();
        System.out.println("--> time needed: " + tm.getResult() + " s");
        System.out.println("--> SSC library size:\t" + this.getSSCCount());

        return true;
    }

    /**
     * Extends this SSC library by documents from a MongoDB collection
     * containing the SSC information.
     *
     * @param collection
     *
     */
    public boolean extend(final MongoCollection<Document> collection) {
        return this.extend(Converter.convertMongoDBCollectionToDocuments(collection));
    }

    /**
     * Extends this SSC library by documents which contain the SSC information.
     *
     * @param documents
     *
     */
    private boolean extendByDocuments(final Collection<Document> documents) {
        final ConcurrentLinkedQueue<SSC> convertedSSCs = new ConcurrentLinkedQueue<>();
        final ArrayList<Callable<SSC>> callables = new ArrayList<>();
        // add all task to do
        for (final Document document : documents) {
            callables.add(() -> Converter.DocumentToSSC(document));
        }
        try {
            Converter.convertDocumentsToSSCs(callables, convertedSSCs, this.nThreads);
            return this.extend(convertedSSCs);
        } catch (InterruptedException e) {
            e.printStackTrace();
            return false;
        }
    }

    /**
     * Extends this SSC library by a collection containing either SSCs ({@link model.SSC}) or Documents ({@link org.bson.Document}) with SSC information.
     *
     * @param collection
     *
     * @see model.SSC
     * @see org.bson.Document
     */
    public boolean extend(final Collection<?> collection) {
        if (collection.isEmpty()){
            return true;
        }
        final Object firstElement = collection.iterator().next();
        if(firstElement instanceof SSC){
            return this.extendBySSCs((Collection<SSC>) collection);
        } else if(firstElement instanceof Document){
            return this.extendByDocuments((Collection<Document>) collection);
        }
        System.err.println("object type in collection is not supported!!!");
        return false;
    }

//    /**
//     * Builds both HOSE code lookup tables of this SSC library containing a list
//     * of shifts for each HOSE code.
//     * Then it will be exported as JSON into a file.
//     *
//     * @param pathToJSONFile path to JSON file to store in
//     *
//     */
//    public void exportHOSECodeLookupTable(final String pathToJSONFile) {

//        this.buildHOSELookupTables();
//        // initialize an executor for parallelization
//        final ExecutorService executor = Utils.initExecuter(this.nThreads);
//        final ArrayList<Callable<withDocument>> callables = new ArrayList<>();
//        // add all task to do
//        for (final String HOSECode : this.HOSECodeLookupTableShifts.keySet()) {
//            callables.add((Callable<withDocument>) () -> {
//                final withDocument document = new withDocument();
//                final ArrayList<Double> shifts = this.HOSECodeLookupTableShifts.get(HOSECode);
//                document.append("HOSECode", HOSECode);
//                document.append("spheres", hose.Utils.getSpheresCount(HOSECode));
//                final withDocument documentShifts = new withDocument();
//                documentShifts.append("count", shifts.size());
//                documentShifts.append("min", (shifts.size() > 0) ? Collections.min(shifts) : null);
//                documentShifts.append("max", (shifts.size() > 0) ? Collections.max(shifts) : null);
//                documentShifts.append("rms", Utils.getRMS(shifts));
//                documentShifts.append("mean", Utils.getMean(shifts));
//                documentShifts.append("median", Utils.getMedian(shifts));
//                documentShifts.append("var", Utils.getVariance(shifts));
//                documentShifts.append("sd", Utils.getStandardDeviation(shifts));
//                documentShifts.append("shifts", shifts);
//                final withDocument documentSSCIndices = new withDocument();
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
    }
}
