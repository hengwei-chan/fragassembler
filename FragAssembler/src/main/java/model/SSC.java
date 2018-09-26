/*
 * The MIT License
 *
 * Copyright 2018 Michael Wenk [https://github.com/michaelwenk].
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
import java.util.ArrayList;
import java.util.HashMap;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Class for representing a subspectrum-substructure-correlation.
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class SSC {
    
    private final Spectrum subspectrum;
    private final Assignment assignment;
    private final IAtomContainer substructure;
    private final HashMap<String, ArrayList<Double>> hoseShifts;
    private final HashMap<String, ArrayList<Integer>> hoseIndices;

    public SSC(final Spectrum subspectrum, final Assignment assignment, final IAtomContainer substructure, 
            final HashMap<String, ArrayList<Double>> hoseShifts, final HashMap<String, ArrayList<Integer>> hoseIndices) {
        this.subspectrum = subspectrum;
        this.assignment = assignment;
        this.substructure = substructure;
        this.hoseShifts = hoseShifts;
        this.hoseIndices = hoseIndices;
    }
    
    public Spectrum getSubspectrum(){
        return this.subspectrum;
    }
    
    public Assignment getAssignments(){
        return this.assignment;
    }
    
    public IAtomContainer getSubstructure(){
        return this.substructure;
    }
    
    public HashMap<String, ArrayList<Double>> getHOSEShifts(){
        return this.hoseShifts;
    }
    
    public HashMap<String, ArrayList<Integer>> getHOSEIndices(){
        return this.hoseIndices;
    }
}
