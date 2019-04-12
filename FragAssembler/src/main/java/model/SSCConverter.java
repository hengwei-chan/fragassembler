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

import casekit.NMR.model.Assignment;
import casekit.NMR.model.Spectrum;
import com.google.gson.Gson;
import com.google.gson.JsonObject;
import com.mongodb.util.JSON;
import org.bson.Document;
import org.openscience.cdk.exception.CDKException;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class SSCConverter {
    
    public static SSC JsonObjectToSSC(final JsonObject sscJsonObject) throws CDKException, CloneNotSupportedException {
        final Gson gson = new Gson();    
        
        System.out.println("\nsubspectrum: ");
        System.out.println(gson.fromJson(sscJsonObject.get("subspectrum"), Spectrum.class).getShifts(0));
        
        final SSC ssc = new SSC(                                
                gson.fromJson(sscJsonObject.get("subspectrum"), Spectrum.class),
                gson.fromJson(sscJsonObject.get("assignment"), Assignment.class),
                gson.fromJson(sscJsonObject.get("substructure"), ExtendedConnectionMatrix.class).toAtomContainer(),
                sscJsonObject.get("rootAtomIndex").getAsInt(),
                sscJsonObject.get("maxSphere").getAsInt()
        );
        ssc.setIndex(sscJsonObject.get("index").getAsLong());

        return ssc;
    }
    
    public static SSC DocumentToSSC(final Document sscDocument) throws CDKException, CloneNotSupportedException{
        final Gson gson = new Gson();
        final SSC ssc = new SSC(
                gson.fromJson(((Document) sscDocument.get("subspectrum")).toJson(), Spectrum.class),
                gson.fromJson(((Document) sscDocument.get("assignment")).toJson(), Assignment.class),
                gson.fromJson(((Document) sscDocument.get("substructure")).toJson(), ExtendedConnectionMatrix.class).toAtomContainer(),
                sscDocument.getInteger("rootAtomIndex"),
                sscDocument.getInteger("maxSphere")
                );
        ssc.setIndex(sscDocument.getLong("index"));
        
        return ssc;
    }
    
    public static Document SSCToDocument(final SSC ssc, final Long sscindex){
        final Gson gson = new Gson();
        final Document document = new Document();        
        document.append("substructure", JSON.parse(gson.toJson(gson.toJsonTree(new ExtendedConnectionMatrix(ssc.getSubstructure()), ExtendedConnectionMatrix.class))));
        document.append("subspectrum", JSON.parse(gson.toJson(gson.toJsonTree(ssc.getSubspectrum(), Spectrum.class))));
        document.append("assignment", JSON.parse(gson.toJson(gson.toJsonTree(ssc.getAssignments(), Assignment.class))));
        document.append("maxSphere", ssc.getMaxSphere());
        document.append("rootAtomIndex", ssc.getRootAtomIndex());
        document.append("index", sscindex);
        document.append("multSections", ssc.getMultiplicitySections());
        
        return document;
    }
}
