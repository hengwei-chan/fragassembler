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
import org.bson.Document;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class SSCConverter {

    private final static Gson GSON = new Gson();

    public static SSC JsonObjectToSSC(final JsonObject sscJsonObject) throws Exception {
        final SSC ssc = new SSC(
                GSON.fromJson(sscJsonObject.get("subspectrum"), Spectrum.class),
                GSON.fromJson(sscJsonObject.get("assignment"), Assignment.class),
                GSON.fromJson(sscJsonObject.get("substructure"), ExtendedConnectionMatrix.class).toAtomContainer(),
                sscJsonObject.get("rootAtomIndex").getAsInt(),
                sscJsonObject.get("maxSphere").getAsInt()
        );
        ssc.setIndex(sscJsonObject.get("index").getAsLong());

        return ssc;
    }
    
    public static SSC DocumentToSSC(final Document sscDocument) throws Exception {
        final SSC ssc = new SSC(
                GSON.fromJson(((Document) sscDocument.get("subspectrum")).toJson(), Spectrum.class),
                GSON.fromJson(((Document) sscDocument.get("assignment")).toJson(), Assignment.class),
                GSON.fromJson(((Document) sscDocument.get("substructure")).toJson(), ExtendedConnectionMatrix.class).toAtomContainer(),
                sscDocument.getInteger("rootAtomIndex"),
                sscDocument.getInteger("maxSphere")
                );
        ssc.setIndex(sscDocument.getLong("index"));

        return ssc;
    }
    
    public static Document SSCToDocument(final SSC ssc, final Long sscIndex){
        final Document document = new Document();
        document.append("substructure", Document.parse(GSON.toJson(GSON.toJsonTree(new ExtendedConnectionMatrix(ssc.getSubstructure()), ExtendedConnectionMatrix.class))));
        document.append("subspectrum", Document.parse(GSON.toJson(GSON.toJsonTree(ssc.getSubspectrum(), Spectrum.class))));
        document.append("assignment", Document.parse(GSON.toJson(GSON.toJsonTree(ssc.getAssignments(), Assignment.class))));
        document.append("maxSphere", ssc.getMaxSphere());
        document.append("rootAtomIndex", ssc.getRootAtomIndex());
        document.append("index", sscIndex);
        document.append("multSections", ssc.getMultiplicitySections());
        
        return document;
    }
}
