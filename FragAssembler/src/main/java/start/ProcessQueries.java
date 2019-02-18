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

import assembly.Assembly;
import casekit.NMR.DB;
import casekit.NMR.model.Spectrum;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import model.SSC;
import model.SSCLibrary;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import rank.SSCRanker;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class ProcessQueries {
    
    private final SSCLibrary sscLibrary;
    private final String pathToQueriesFile, pathToOutputsFolder;
    private final int nThreads, nStarts;
    private final SmilesParser smilesParser;
    private final BufferedReader br;
    private final int minMatchingSphereCount;
    
    public ProcessQueries(final SSCLibrary sscLibrary, final String pathToQueriesFile, final String pathToOutputsFolder) throws FileNotFoundException {
        this(sscLibrary, pathToQueriesFile, pathToOutputsFolder, 1, -1);
    }

    public ProcessQueries(final SSCLibrary sscLibrary, final String pathToQueriesFile, final String pathToOutputsFolder, final int nThreads, final int nStarts) throws FileNotFoundException {
        this.minMatchingSphereCount = 0;
        this.pathToQueriesFile = pathToQueriesFile;
        this.pathToOutputsFolder = pathToOutputsFolder;
        this.nThreads = nThreads;        
        this.sscLibrary = sscLibrary;
        this.sscLibrary.setNThreads(this.nThreads);
        this.nStarts = nStarts;
        this.smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        this.br = new BufferedReader(new FileReader(this.pathToQueriesFile));
    }    
    
    
    public void process() throws FileNotFoundException, InterruptedException, IOException, InvalidSmilesException, CDKException, CloneNotSupportedException {
        int nStartSSCs;
        Spectrum querySpectrum;
        HashMap<String, SSC> solutions;
        SSCRanker sscRanker = new SSCRanker(this.sscLibrary, this.nThreads);

        IAtomContainer solutionAtomContainer;
        ArrayList<Double> shiftsQuerySpectrum, shiftsSolutionSpectrum;
        int querySpectrumCounter = 0;
        BufferedWriter bw;

        System.out.println("\n\n-> processing query file: \"" + this.pathToQueriesFile + "\" ...");

        Iterator<String> it = this.br.lines().iterator();
        String line;
        while (it.hasNext()) {
            line = it.next();
            if(line.startsWith("#")){
                System.out.println("\n\nskip query: " + querySpectrumCounter + " -> " + line + "\n");
                querySpectrumCounter++;
                continue;
            }
            System.out.println("\n\nnow processing query: " + querySpectrumCounter + " -> " + line + "\n");

            querySpectrum = DB.NMRShiftDBSpectrumToSpectrum(line, "C");
            sscRanker.rank(querySpectrum, Start.PICK_PRECISION);
            System.out.println("no. of matches: " + sscRanker.getMatchFactors().size());
            System.out.println("match factors: " + sscRanker.getMatchFactors());
            System.out.println("ranked SSC indices: " + sscRanker.getRankedSSCIndices());
            System.out.println("ranked match factors: " + sscRanker.getRankedMatchFactors());
            System.out.println("ranked SSC library indices: " + sscRanker.getRankedSSCLibrary().getSSCIndices() + "\n");

            if ((this.nStarts > 0) && (this.nStarts < sscRanker.getMatchFactors().size())) {
                nStartSSCs = this.nStarts;
            } else {
                nStartSSCs = sscRanker.getMatchFactors().size();
            }
            System.out.println("\nnumber of start SSCs for query " + querySpectrumCounter + ":\t" + nStartSSCs);

            solutions = Assembly.assemble(nStartSSCs, this.nThreads, sscRanker.getRankedSSCLibrary(), minMatchingSphereCount, querySpectrum, Start.SHIFT_TOL, Start.MATCH_FACTOR_THRS, Start.PICK_PRECISION);

            System.out.println("\nsolutions for query " + querySpectrumCounter + " (" + line + "):\t" + solutions.size());

            shiftsQuerySpectrum = querySpectrum.getShifts(0);
            Collections.sort(shiftsQuerySpectrum);

            bw = new BufferedWriter(new FileWriter(this.pathToOutputsFolder + "/results_" + querySpectrumCounter + ".smiles"));
            for (final String smiles : solutions.keySet()) {
                solutionAtomContainer = this.smilesParser.parseSmiles(smiles);

                shiftsSolutionSpectrum = solutions.get(smiles).getSubspectrum().getShifts(0);
                Collections.sort(shiftsSolutionSpectrum);
                System.out.println("query spectrum   :\t" + shiftsQuerySpectrum);
                System.out.println("solution spectrum:\t" + shiftsSolutionSpectrum);
                System.out.println("SMILES:          :\t" + smiles);
                System.out.println("              --> \t" + "atoms: " + solutionAtomContainer.getAtomCount() + ", bonds: " + solutionAtomContainer.getBondCount());

                bw.append(smiles);
                bw.newLine();
                bw.flush();
            }
            bw.close();
            querySpectrumCounter++;
        }
        this.br.close();
    }
    
}
