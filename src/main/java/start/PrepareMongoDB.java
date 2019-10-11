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
package start;

import casekit.NMR.dbservice.MongoDB;
import com.mongodb.BasicDBObject;
import com.mongodb.MongoClient;
import com.mongodb.client.MongoCollection;
import model.SSCLibrary;
import org.bson.Document;
import org.openscience.cdk.exception.CDKException;

import java.io.FileNotFoundException;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class PrepareMongoDB {
    
    private final boolean importFromNMRShiftDB, extendFromNMRShiftDB;
    private final int nThreads, maxSphere;
    private final String pathToNMRShiftDB;
    private SSCLibrary sscLibrary;   
    private MongoClient mongo;
    private MongoCollection<Document> collection;
    private final TimeMeasurement tm;
    
    public PrepareMongoDB(final int nThreads) throws CDKException {
        this(nThreads, false, false, "", -1);
    }    
    
    public PrepareMongoDB(final int nThreads, final boolean importFromNMRShiftDB, final boolean extendFromNMRShiftDB, final String pathToNMRShiftDB, final int maxSphere) throws CDKException {
        this.importFromNMRShiftDB = importFromNMRShiftDB;
        this.extendFromNMRShiftDB = extendFromNMRShiftDB;
        this.pathToNMRShiftDB = pathToNMRShiftDB;        
        this.nThreads = nThreads;
        this.maxSphere = maxSphere;
        if ((this.importFromNMRShiftDB || this.extendFromNMRShiftDB) && (this.maxSphere < 2)) {
            throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": invalid number of maximum sphere: \"" + this.maxSphere + "\" < 2");
        }
        this.tm = new TimeMeasurement();        
    }  
    
    public void prepare(final String mongoUser, final String mongoPassword, 
            final String mongoAuthDB, final String mongoDBName, final String mongoDBCollection) throws CDKException, FileNotFoundException, InterruptedException, CloneNotSupportedException {
        
        this.mongo = MongoDB.login(mongoUser, mongoPassword, mongoAuthDB);
        this.collection = MongoDB.getCollection(this.mongo, mongoDBName, mongoDBCollection);
        this.sscLibrary = new SSCLibrary(this.nThreads);
        if (this.importFromNMRShiftDB) {
            // empty MongoDB collection
            this.tm.start();
            this.collection.deleteMany(new BasicDBObject());
            System.out.println("Collection \"" + mongoDBCollection + "\" is now empty");
            this.tm.stop();
            System.out.println("--> time needed: " + this.tm.getResult() + " s");
            
            // create SSC library for a specific max sphere and insert into into the MongoDB collection
            for (int m = 2; m <= this.maxSphere; m++) {
                System.out.println("Building SSCs for " + m + "-spheres...");
                this.tm.start();
                this.sscLibrary.extend(this.pathToNMRShiftDB, Start.SPECTRUM_PROPERTY, m);
                System.out.println("SSCs for " + m + "-spheres build!!!");
                this.tm.stop();
                System.out.println("--> time needed: " + this.tm.getResult() + " s");
                System.out.println("-> #SSCs in SSC library: " + this.sscLibrary.getSSCCount());
            }
            System.out.println("storing SSC library into MongoDB in database \"" + mongoDBName + "\" in collection \"" + mongoDBCollection + "\"...");
            this.tm.start();
            this.sscLibrary.exportToMongoDB(this.collection);
            System.out.println("-> SSC library stored into MongoDB in database \"" + mongoDBName + "\" in collection \"" + mongoDBCollection + "\" -> collection size: " + this.collection.countDocuments() + "!!!");
            this.tm.stop();
            System.out.println("--> time needed: " + tm.getResult() + " s");
        } else if (this.extendFromNMRShiftDB) {
            System.out.println("Building SSCs for " + this.maxSphere + "-spheres...");
            this.tm.start();
            this.sscLibrary.extend(this.pathToNMRShiftDB, Start.SPECTRUM_PROPERTY, this.maxSphere);
            System.out.println("SSCs for " + this.maxSphere + "-spheres build!!!");
            this.tm.stop();
            System.out.println("--> time needed: " + this.tm.getResult() + " s");
            System.out.println("-> #SSCs in SSC library: " + this.sscLibrary.getSSCCount());
            System.out.println("Extension of MongoDB collection \"" + mongoDBCollection + "\"...");
            this.tm.start();
            this.sscLibrary.exportToMongoDB(this.collection);
            System.out.println("-> MongoDB extended in database \"" + mongoDBName + "\" in collection \"" + mongoDBCollection + "\" by " + this.sscLibrary.getSSCCount() + " SSCs -> new collection size: " + this.collection.countDocuments() + " !!!");
            this.tm.stop();
            System.out.println("--> time needed: " + this.tm.getResult() + " s");
        }           
        MongoDB.logout(this.mongo);
        
    }
}
