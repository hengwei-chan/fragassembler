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
import java.io.IOException;

public class Prepare {
    private final boolean buildFromNMRShiftDB;
    private final int nThreads, maxSphere;
    private final String pathToNMRShiftDB;
    private final TimeMeasurement tm;

    public Prepare(final int nThreads, final boolean importFromNMRShiftDB, final String pathToNMRShiftDB, final int maxSphere) throws CDKException {
        this.buildFromNMRShiftDB = importFromNMRShiftDB;
        this.pathToNMRShiftDB = pathToNMRShiftDB;
        this.nThreads = nThreads;
        this.maxSphere = maxSphere;
        if (this.buildFromNMRShiftDB && (this.maxSphere < 2)) {
            throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": invalid number of maximum sphere: \"" + this.maxSphere + "\" < 2");
        }
        this.tm = new TimeMeasurement();
    }


    public SSCLibrary prepareJSON(final String pathToJSON) throws CDKException, InterruptedException, IOException {
        final SSCLibrary sscLibrary = new SSCLibrary(this.nThreads);
        if (this.buildFromNMRShiftDB) {
            sscLibrary.buildFromNMRShiftDB(this.pathToNMRShiftDB, Start.SPECTRUM_PROPERTY, this.maxSphere);
            sscLibrary.exportToJSONFile(pathToJSON);
        } else {
            sscLibrary.buildFromJSON(pathToJSON);
        }

        return sscLibrary;
    }

    public void prepareMongoDB(final String mongoUser, final String mongoPassword, final String mongoAuthDB, final String mongoDBName, final String mongoDBCollection) throws CDKException, FileNotFoundException, InterruptedException {
        final SSCLibrary sscLibrary = new SSCLibrary(this.nThreads);
        final MongoClient mongo = MongoDB.login(mongoUser, mongoPassword, mongoAuthDB);
        final MongoCollection<Document> collection = MongoDB.getCollection(mongo, mongoDBName, mongoDBCollection);
        if (this.buildFromNMRShiftDB) {
            // empty MongoDB collection
            this.tm.start();
            collection.deleteMany(new BasicDBObject());
            System.out.println("Collection \"" + mongoDBCollection + "\" is now empty");
            this.tm.stop();
            System.out.println("--> time needed: " + this.tm.getResult() + " s");

            // create SSC library for a specific max sphere and insert it into the MongoDB collection
            sscLibrary.buildFromNMRShiftDB(this.pathToNMRShiftDB, Start.SPECTRUM_PROPERTY, this.maxSphere);

            System.out.println("storing SSC library into MongoDB in database \"" + mongoDBName + "\" in collection \"" + mongoDBCollection + "\"...");
            this.tm.start();
            sscLibrary.exportToMongoDB(collection);
            System.out.println("-> SSC library stored into MongoDB in database \"" + mongoDBName + "\" in collection \"" + mongoDBCollection + "\" -> collection size: " + collection.countDocuments() + "!!!");
            this.tm.stop();
            System.out.println("--> time needed: " + tm.getResult() + " s");
        }

        MongoDB.logout(mongo);
    }
}
