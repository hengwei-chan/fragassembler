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
package conversion;

import casekit.NMR.model.Assignment;
import casekit.NMR.model.Spectrum;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.mongodb.client.MongoCollection;
import model.ExtendedConnectionMatrix;
import model.SSC;
import model.SSCLibrary;
import org.bson.Document;
import parallel.ParallelTasks;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.function.Consumer;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Converter {

    private final static Gson GSON = new GsonBuilder().setLenient().create();
    
    public static SSC DocumentToSSC(final Document sscDocument) throws Exception {
        final SSC ssc = new SSC(
                GSON.fromJson(((Document) sscDocument.get("subspectrum")).toJson(), Spectrum.class),
                GSON.fromJson(((Document) sscDocument.get("assignment")).toJson(), Assignment.class),
                GSON.fromJson(((Document) sscDocument.get("substructure")).toJson(), ExtendedConnectionMatrix.class).toAtomContainer(),
                sscDocument.getInteger("rootAtomIndex"),
                sscDocument.getInteger("maxSphere"),
//                GSON.fromJson(((withDocument) sscDocument.get("HOSECodes")).toJson(), HashMap.class),
//                GSON.fromJson(((withDocument) sscDocument.get("connectionTrees")).toJson(), HashMap.class),
//                GSON.fromJson(((withDocument) sscDocument.get("attachedHydrogensInOuterSphere")).toJson(), ArrayList.class),
                GSON.fromJson(((Document) sscDocument.get("shifts")).toJson(), HashMap.class),
                GSON.fromJson(((Document) sscDocument.get("shiftsRanges")).toJson(), HashMap.class),
                sscDocument.get("unsaturatedAtomIndices", ArrayList.class),
                GSON.fromJson(((Document) sscDocument.get("multiplicitySections")).toJson(), HashMap.class)
                );
        ssc.setIndex(sscDocument.getLong("index"));

        return ssc;
    }
    
    public static Document SSCToDocument(final SSC ssc){
        final Document document = new Document();
        document.append("HOSECode", ssc.getHOSECode(ssc.getRootAtomIndex()));
        document.append("substructure", Document.parse(GSON.toJson(GSON.toJsonTree(new ExtendedConnectionMatrix(ssc.getSubstructure()), ExtendedConnectionMatrix.class))));
        document.append("subspectrum", Document.parse(GSON.toJson(GSON.toJsonTree(ssc.getSubspectrum(), Spectrum.class))));
        document.append("assignment", Document.parse(GSON.toJson(GSON.toJsonTree(ssc.getAssignments(), Assignment.class))));
        document.append("maxSphere", ssc.getMaxSphere());
        document.append("rootAtomIndex", ssc.getRootAtomIndex());
        document.append("index", ssc.getIndex());
//        document.append("HOSECodes", withDocument.parse(GSON.toJson(GSON.toJsonTree(ssc.getHOSECodes(), HashMap.class))));
//        document.append("connectionTrees", withDocument.parse(GSON.toJson(GSON.toJsonTree(ssc.getConnectionTrees(), HashMap.class))));
        document.append("shifts", Document.parse(GSON.toJson(GSON.toJsonTree(ssc.getShifts(), HashMap.class))));
        document.append("shiftsRanges", Document.parse(GSON.toJson(GSON.toJsonTree(ssc.getShiftsRanges(), HashMap.class))));
        document.append("unsaturatedAtomIndices", ssc.getUnsaturatedAtomIndices());
        document.append("multiplicitySections", ssc.getMultiplicitySections());

        return document;
    }

    public static SSC JSONToSSC(final String json) throws Exception {
        final JsonObject jsonObject = new JsonParser().parse(json).getAsJsonObject();
        final SSC ssc = new SSC(
                GSON.fromJson(jsonObject.get("subspectrum"), Spectrum.class),
                GSON.fromJson(jsonObject.get("assignment"), Assignment.class),
                GSON.fromJson(jsonObject.get("substructure"), ExtendedConnectionMatrix.class).toAtomContainer(),
                jsonObject.get("rootAtomIndex").getAsInt(),
                jsonObject.get("maxSphere").getAsInt(),
//                GSON.fromJson(jsonObject.get("HOSECodes"), HashMap.class),
//                GSON.fromJson(jsonObject.get("connectionTrees"), HashMap.class),
//                GSON.fromJson(jsonObject.get("attachedHydrogensInOuterSphere"), ArrayList.class),
                GSON.fromJson(jsonObject.get("shifts"), HashMap.class),
                GSON.fromJson(jsonObject.get("shiftsRanges"), HashMap.class),
                GSON.fromJson(jsonObject.get("unsaturatedAtomIndices"), ArrayList.class),
                GSON.fromJson(jsonObject.get("multiplicitySections"), HashMap.class)

        );
        ssc.setIndex(jsonObject.get("index").getAsJsonObject().get("$numberLong").getAsLong());

        return ssc;
    }

    public static String SSCToJSON(final SSC ssc){
        return new Document().append(String.valueOf(ssc.getIndex()), Converter.SSCToDocument(ssc)).toJson();
    }

    /**
     * Converts SSCs to documents in parallel and adds them to the given collection.
     *
     * @param callables
     * @param collection
     * @param nThreads
     * @throws InterruptedException
     *
     * @see ParallelTasks#processTasks(Collection, Consumer, int)
     * @see Collection#add(Object)
     */
    public static void convertSSCsToDocuments(final Collection<Callable<Document>> callables, final Collection<Document> collection, final int nThreads) throws InterruptedException {
        ParallelTasks.processTasks(callables, collection::add, nThreads);
    }

    /**
     * Converts documents to SSCs in parallel and adds them to the given collection.
     *
     * @param callables
     * @param collection
     * @param nThreads
     * @throws InterruptedException
     *
     * @see ParallelTasks#processTasks(Collection, Consumer, int)
     * @see Collection#add(Object)
     */
    public static void convertDocumentsToSSCs(final Collection<Callable<SSC>> callables, final Collection<SSC> collection, final int nThreads) throws InterruptedException {
        ParallelTasks.processTasks(callables, collection::add, nThreads);
    }

    /**
     * Adds documents from a MongoDB collection into a collection and returns it.
     *
     * @param collection
     *
     */
    public static Collection<Document> convertMongoDBCollectionToDocuments(final MongoCollection<Document> collection) {
        final Collection<Document> documents = new ArrayList<>();
        for (final Document document : collection.find()) {
            documents.add(document);
        }
        return documents;
    }

    /**
     * Stores this SSC library into a string in JSON format.
     *
     * @return this SSC library as JSON string
     */
    public static String SSCLibraryToJSON(final SSCLibrary sscLibrary){
        final StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append("{");
        stringBuilder.append("\n");
        String json;
        long sscCounter = 0;
        for (final SSC ssc : sscLibrary.getSSCs()) {
            json = Converter.SSCToJSON(ssc);
            stringBuilder.append(json.substring(1, json.length() - 1));
            if(sscCounter < sscLibrary.getSSCCount() - 1){
                stringBuilder.append(",");
            }
            stringBuilder.append("\n");

            sscCounter++;
        }
        stringBuilder.append("}");

        return stringBuilder.toString();
    }

    public static boolean SSCLibraryToJSONFile(final SSCLibrary sscLibrary, final String pathToJSONFile) {
        try {
            final BufferedWriter bw = new BufferedWriter(new FileWriter(pathToJSONFile));
            bw.write(Converter.SSCLibraryToJSON(sscLibrary));
            bw.close();

            return true;
        } catch (IOException e) {
            e.printStackTrace();
            return false;
        }
    }

    public static SSCLibrary JSONToSSCLibrary(final String json, final int nThreads) throws InterruptedException, FileNotFoundException {
        return Converter.JSONToSSCLibrary(new BufferedReader(new StringReader(json)), nThreads);
    }

    public static SSCLibrary JSONFileToSSCLibrary(final String pathToJSONFile, final int nThreads) throws InterruptedException, FileNotFoundException {
        return Converter.JSONToSSCLibrary(new BufferedReader(new FileReader(pathToJSONFile)), nThreads);
    }

    private static SSCLibrary JSONToSSCLibrary(final BufferedReader br, final int nThreads) throws InterruptedException, FileNotFoundException {
        final ConcurrentLinkedQueue<SSC> convertedSSCs = new ConcurrentLinkedQueue<>();
        final JsonParser jsonParser = new JsonParser();
        final List<Callable<SSC>> callables = Collections.synchronizedList(new ArrayList<>());
        // add all task to do
        br.lines().parallel().forEach( line -> {
            if((line.trim().length() > 1) || (!line.trim().startsWith("{") && !line.trim().endsWith("}"))){
                final StringBuilder sscInJSON = new StringBuilder();
                if(line.endsWith(",")){
                    sscInJSON.append(line, 0, line.length() - 1);
                } else {
                    sscInJSON.append(line);
                }
                callables.add(() -> Converter.JSONToSSC(jsonParser.parse(sscInJSON.substring(sscInJSON.toString().indexOf("{"))).getAsJsonObject().toString()));
            }
        });
        Converter.convertDocumentsToSSCs(callables, convertedSSCs, nThreads);

        final SSCLibrary sscLibrary = new SSCLibrary(nThreads);
        sscLibrary.extend(convertedSSCs);

        return sscLibrary;
    }
}
