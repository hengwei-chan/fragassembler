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

import casekit.NMR.Utils;
import casekit.NMR.model.Assignment;
import casekit.NMR.model.Spectrum;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.matrix.ConnectionMatrix;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.Atom;
import org.openscience.cdk.silent.Bond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class SSCInJSON {
    
    final double[][] connectionMatrix;
    final String[] atomTypes;
    final Integer[] hydrogenCounts, valencies;
    final Double[] charges;
    final boolean[] isInRingAtoms, isAromaticAtoms, isInRingBonds, isAromaticBonds;
    final Integer[][] bondIDs;
    final Spectrum subspectrum;
    final Assignment assignment;
    final int index, rootAtomIndex, maxSphere;
    
    public SSCInJSON(final SSC ssc){
        this.connectionMatrix = ConnectionMatrix.getMatrix(ssc.getSubstructure());
        this.atomTypes = new String[ssc.getAtomCount()];
        this.hydrogenCounts = new Integer[ssc.getAtomCount()];
        this.isInRingAtoms = new boolean[ssc.getAtomCount()];
        this.isAromaticAtoms = new boolean[ssc.getAtomCount()];
        this.valencies = new Integer[ssc.getAtomCount()];
        this.charges = new Double[ssc.getAtomCount()];
        this.bondIDs = new Integer[ssc.getBondCount()][2];        
        this.isInRingBonds = new boolean[ssc.getSubstructure().getBondCount()];
        this.isAromaticBonds = new boolean[ssc.getSubstructure().getBondCount()];
                
        this.subspectrum = ssc.getSubspectrum();
        this.assignment = ssc.getAssignments();
        this.index = ssc.getIndex();
        this.rootAtomIndex = ssc.getRootAtomIndex();
        this.maxSphere = ssc.getMaxSphere();
        
        this.initAtomsProperties(ssc.getSubstructure());
        this.initBondsProperties(ssc.getSubstructure());
    }
    
    private void initAtomsProperties(final IAtomContainer structure){
        for (final IAtom atom : structure.atoms()) {
            this.atomTypes[atom.getIndex()] = atom.getSymbol();
            this.hydrogenCounts[atom.getIndex()] = atom.getImplicitHydrogenCount();
            this.isInRingAtoms[atom.getIndex()] = atom.isInRing();
            this.isAromaticAtoms[atom.getIndex()] = atom.isAromatic();
            this.valencies[atom.getIndex()] = atom.getValency();
            this.charges[atom.getIndex()] = atom.getCharge();
        }
    }
    
    private void initBondsProperties(final IAtomContainer structure){
        for (final IBond bond : structure.bonds()) {
            this.bondIDs[bond.getIndex()][0] = bond.getAtom(0).getIndex();
            this.bondIDs[bond.getIndex()][1] = bond.getAtom(1).getIndex();
            this.isInRingBonds[bond.getIndex()] = bond.isInRing();
            this.isAromaticBonds[bond.getIndex()] = bond.isAromatic();
        }
    }
    
    public int getIndex(){
        return this.index;
    }
    
    public SSC toSSC() throws CDKException, CloneNotSupportedException {
        final SSC ssc = new SSC(this.subspectrum, this.assignment, this.buildAtomContainer(), this.rootAtomIndex, this.maxSphere);
        ssc.setIndex(this.index);
        return ssc;
    }
    
    private IAtomContainer buildAtomContainer(){
        final IAtomContainer substructure = SilentChemObjectBuilder.getInstance().newAtomContainer();
        IAtom atom;
        for (int i = 0; i < this.connectionMatrix.length; i++) {
            atom = new Atom(this.atomTypes[i]);
            atom.setImplicitHydrogenCount(this.hydrogenCounts[i]);
            atom.setIsInRing(this.isInRingAtoms[i]);
            atom.setIsAromatic(this.isAromaticAtoms[i]);
            atom.setValency(this.valencies[i]);
            atom.setCharge(this.charges[i]);
            
            substructure.addAtom(atom);
        }
        int atomIndex1, atomIndex2;
        IBond bond;
        for (int i = 0; i < this.bondIDs.length; i++) {
            atomIndex1 = this.bondIDs[i][0];
            atomIndex2 = this.bondIDs[i][1];
            bond = new Bond(substructure.getAtom(atomIndex1), substructure.getAtom(atomIndex2), Utils.getBondOrder((int) connectionMatrix[atomIndex1][atomIndex2]));
            bond.setIsInRing(this.isInRingBonds[i]);
            bond.setIsAromatic(this.isAromaticBonds[i]);

            substructure.addBond(bond);
        }               
        
        return substructure;
    }
    
}
