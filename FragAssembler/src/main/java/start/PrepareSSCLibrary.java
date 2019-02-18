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
package start;

import com.mongodb.BasicDBObject;
import com.mongodb.MongoClient;
import com.mongodb.MongoClientOptions;
import com.mongodb.MongoCredential;
import com.mongodb.ServerAddress;
import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoDatabase;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Collections;
import model.SSCLibrary;
import org.bson.Document;
import org.openscience.cdk.exception.CDKException;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class PrepareSSCLibrary {
    
    private final boolean importFromNMRShiftDB, extendFromNMRShiftDB;
    private final int nThreads, maxSphere;
    private final String pathToNMRShiftDB;
    private SSCLibrary sscLibrary;   
    
    public PrepareSSCLibrary(final int nThreads){
        this(nThreads, false, false, "", -1);
    }    
    
    public PrepareSSCLibrary(final int nThreads, final boolean importFromNMRShiftDB, final boolean extendFromNMRShiftDB, final String pathToNMRShiftDB, final int maxSphere){
        this.importFromNMRShiftDB = importFromNMRShiftDB;
        this.extendFromNMRShiftDB = extendFromNMRShiftDB;
        this.pathToNMRShiftDB = pathToNMRShiftDB;
        this.maxSphere = maxSphere;
        this.nThreads = nThreads;
    }
    
    public SSCLibrary processJSONFile(final String pathToJSON) throws CDKException, FileNotFoundException, InterruptedException, CloneNotSupportedException, IOException {

        this.sscLibrary = new SSCLibrary(this.nThreads);
        if (this.importFromNMRShiftDB) {
            int offset = 0;
            // create SSC library for a specific max sphere and insert into the JSON file
            for (int m = 2; m <= this.maxSphere; m++) {
                System.out.println("Building SSCs for " + m + "-spheres...");
                this.sscLibrary.extend(this.pathToNMRShiftDB, Start.SPECTRUM_PROPERTY, Start.SPECTRUM_PROPERTY_ATOMTYPE, m, offset);
                System.out.println("SSCs for " + m + "-spheres build!!!");
                System.out.println("-> #SSCs in SSC library: " + this.sscLibrary.getSSCCount());

                offset = Collections.max(this.sscLibrary.getSSCIndices()) + 1;
            }
            System.out.println("now storing SSC library into JSON file \"" + pathToJSON + "\"...");
            this.sscLibrary.exportToJSONFile(pathToJSON);
            System.out.println("-> SSC library stored into JSON file");
        } else {
            System.out.println("-> importing SSC library from JSON file...");
            this.sscLibrary.importFromJSONFile(pathToJSON, 0);
            System.out.println("-> SSC library imported from JSON file!!!");
            System.out.println("--> SSC library size:\t" + this.sscLibrary.getSSCCount());
            if (this.extendFromNMRShiftDB) {
                int offset = Collections.max(this.sscLibrary.getSSCIndices()) + 1;
                System.out.println("offset: " + offset);
                System.out.println("Building SSCs for " + this.maxSphere + "-spheres...");
                this.sscLibrary.extend(pathToNMRShiftDB, Start.SPECTRUM_PROPERTY, Start.SPECTRUM_PROPERTY_ATOMTYPE, maxSphere, offset);
                System.out.println("SSCs for " + this.maxSphere + "-spheres build and added!!!");
                System.out.println("-> #SSCs in SSC library: " + this.sscLibrary.getSSCCount());
                this.sscLibrary.exportToJSONFile(pathToJSON);
                System.out.println("-> SSC library stored into JSON file \"" + pathToJSON + "\"");
            }
        }
        
        return this.sscLibrary;
    }

    public SSCLibrary processMongoDB(final String mongoUser, final String mongoPassword, 
            final String mongoAuthDB, final String mongoDBName, final String mongoDBCollection) throws CDKException, FileNotFoundException, InterruptedException, CloneNotSupportedException {
        final MongoClient mongo;
        final MongoCollection<Document> collection;
        try {
            // Creating a Mongo client   
//            final List<ServerAddress> seeds = new ArrayList<>();
//            final List<MongoCredential> credentials = new ArrayList<>();
//            for (int i = 0; i < this.nThreads; i++) {
////                seeds.add(new ServerAddress("127.0.0.1", 27017));
//                credentials.add(MongoCredential.createCredential(
//                        this.mongoUser, 
//                        this.mongoAuthDB,
//                        this.mongoPassword.toCharArray()));
//            }
            MongoClientOptions.Builder builder = MongoClientOptions.builder();
//            builder.connectionsPerHost(this.nThreads);
//            builder.threadsAllowedToBlockForConnectionMultiplier(this.nThreads);
            mongo = new MongoClient(
                    new ServerAddress("127.0.0.1", 27017),
                    MongoCredential.createCredential(
                            mongoUser,
                            mongoAuthDB,
                            mongoPassword.toCharArray()),
                    builder.build());
            System.out.println("Login to MongoDB was successfull");
            // Accessing the database 
            final MongoDatabase database = mongo.getDatabase(mongoDBName);
            System.out.println("Access to database \"" + mongoDBName + "\" was successfull");
            // Retrieving a collection
            collection = database.getCollection(mongoDBCollection);
            System.out.println("Retrieval of collection \"" + mongoDBCollection + "\" was successfull");
        } catch (Exception e) {
            throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": could not connect to database \"" + mongoDBName + "\" or collection \"" + mongoDBCollection + "\"");
        }

        this.sscLibrary = new SSCLibrary(this.nThreads);
        if (this.importFromNMRShiftDB) {
            // empty MongoDB collection
            collection.deleteMany(new BasicDBObject());
            System.out.println("Collection \"" + mongoDBCollection + "\" is now empty");

            int offset = 0;
            // create SSC library for a specific max sphere and insert into into the MongoDB collection
            for (int m = 2; m <= this.maxSphere; m++) {
                System.out.println("Building SSCs for " + m + "-spheres...");
                this.sscLibrary.extend(this.pathToNMRShiftDB, Start.SPECTRUM_PROPERTY, Start.SPECTRUM_PROPERTY_ATOMTYPE, m, offset);
                System.out.println("SSCs for " + m + "-spheres build!!!");
                System.out.println("-> #SSCs in SSC library: " + this.sscLibrary.getSSCCount());

                offset = Collections.max(this.sscLibrary.getSSCIndices()) + 1;
            }
            System.out.println("storing SSC library into MongoDB in database \"" + mongoDBName + "\" in collection \"" + mongoDBCollection + "\"...");
            this.sscLibrary.exportToMongoDB(collection);
            System.out.println("-> SSC library stored into MongoDB in database \"" + mongoDBName + "\" in collection \"" + mongoDBCollection + "\" -> collection size: " + collection.countDocuments() + "!!!");
        } else {
            System.out.println("-> importing SSC library from MongoDB...");
            this.sscLibrary.importFromMongoDB(collection);
            System.out.println("-> SSC library imported from MongoDB!!!");
            System.out.println("--> SSC library size:\t" + this.sscLibrary.getSSCCount());
            if (this.extendFromNMRShiftDB) {
                // db.data.find({}, {_id: 1}).sort({_id: -1}).limit(1)
                int offset = Collections.max(this.sscLibrary.getSSCIndices()) + 1;//((int) collection.find().sort(new BasicDBObject("_id", -1)).limit(1).first().get("_id")) + 1;
                SSCLibrary sscLibraryTemp = new SSCLibrary(this.nThreads);
                System.out.println("offset: " + offset);
                System.out.println("Building SSCs for " + this.maxSphere + "-spheres...");
                sscLibraryTemp.extend(this.pathToNMRShiftDB, Start.SPECTRUM_PROPERTY, Start.SPECTRUM_PROPERTY_ATOMTYPE, this.maxSphere, offset);
                System.out.println("SSCs for " + this.maxSphere + "-spheres build!!!");
                System.out.println("-> #SSCs in SSC library: " + sscLibraryTemp.getSSCCount());
                sscLibraryTemp.exportToMongoDB(collection);
                System.out.println("-> SSC library stored into MongoDB in database \"" + mongoDBName + "\" in collection \"" + mongoDBCollection + "\" -> collection size: " + collection.countDocuments() + "!!!");
                this.sscLibrary.extend(sscLibraryTemp);
                System.out.println("SSCs for " + this.maxSphere + "-spheres added!!!");
                System.out.println("-> #SSCs in SSC library: " + this.sscLibrary.getSSCCount());
            }
        }
        mongo.close();
        
        return this.sscLibrary;
    }
}
