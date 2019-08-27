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
package start;

import model.SSCLibrary;
import org.openscience.cdk.exception.CDKException;

import java.io.FileNotFoundException;
import java.io.IOException;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class PrepareJSONFile {
    
    private final boolean importFromNMRShiftDB, extendFromNMRShiftDB;
    private final int nThreads, maxSphere;
    private final String pathToNMRShiftDB;
    private SSCLibrary sscLibrary;
    private final TimeMeasurement tm;

    public PrepareJSONFile(final int nThreads) throws CDKException {
        this(nThreads, false, false, "", -1);
    }

    public PrepareJSONFile(final int nThreads, final boolean importFromNMRShiftDB, final boolean extendFromNMRShiftDB, final String pathToNMRShiftDB, final int maxSphere) throws CDKException {
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

    public SSCLibrary prepare(final String pathToJSON, final boolean removeDuplicates) throws CDKException, FileNotFoundException, InterruptedException, CloneNotSupportedException, IOException {
        this.sscLibrary = new SSCLibrary(this.nThreads);
        if (this.importFromNMRShiftDB) {
            long offset = 0;
            // create SSC library for a specific max sphere and insert into the JSON file
            for (int m = 2; m <= this.maxSphere; m++) {
                System.out.println("Building SSCs for " + m + "-spheres...");
                this.tm.start();
                this.sscLibrary.extend(this.pathToNMRShiftDB, Start.SPECTRUM_PROPERTY, m, offset);
                if (removeDuplicates) {
                    this.sscLibrary.removeDuplicates(Start.DUPLICATES_SHIFT_TOL);
                }
                System.out.println("SSCs for " + m + "-spheres build!!!");
                this.tm.stop();
                System.out.println("--> time needed: " + this.tm.getResult() + " s");
                System.out.println("-> #SSCs in SSC library: " + this.sscLibrary.getSSCCount());

                offset = this.sscLibrary.getLastSSCIndex() + 1;
            }
            System.out.println("now storing SSC library into JSON file \"" + pathToJSON + "\"...");
            this.tm.start();
            this.sscLibrary.exportToJSONFile(pathToJSON);
            this.tm.stop();
            System.out.println("--> time needed: " + this.tm.getResult() + " s");
            System.out.println("-> SSC library stored into JSON file");
        } else {
            System.out.println("-> importing SSC library from JSON file...");
            this.tm.start();
            this.sscLibrary.importFromJSONFile(pathToJSON, 0);
            System.out.println("-> SSC library imported from JSON file!!!");
            this.tm.stop();
            System.out.println("--> time needed: " + this.tm.getResult() + " s");
            System.out.println("--> SSC library size:\t" + this.sscLibrary.getSSCCount());
            if (this.extendFromNMRShiftDB) {
                long offset = this.sscLibrary.getLastSSCIndex() + 1;
                System.out.println("offset: " + offset);
                System.out.println("Building SSCs for " + this.maxSphere + "-spheres...");
                this.tm.start();
                this.sscLibrary.extend(this.pathToNMRShiftDB, Start.SPECTRUM_PROPERTY, maxSphere, offset);
                if (removeDuplicates) {
                    this.sscLibrary.removeDuplicates(Start.DUPLICATES_SHIFT_TOL);
                }
                System.out.println("SSCs for " + this.maxSphere + "-spheres build and added!!!");
                this.tm.stop();
                System.out.println("--> time needed: " + this.tm.getResult() + " s");
                System.out.println("-> #SSCs in SSC library: " + this.sscLibrary.getSSCCount());
                this.tm.start();
                this.sscLibrary.exportToJSONFile(pathToJSON);
                System.out.println("-> SSC library stored into JSON file \"" + pathToJSON + "\"");
                this.tm.stop();
                System.out.println("--> time needed: " + this.tm.getResult() + " s");
            }
        }

        return this.sscLibrary;
    }
}
