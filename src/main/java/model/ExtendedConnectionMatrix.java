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
package model;

import casekit.NMR.Utils;
import org.openscience.cdk.graph.matrix.ConnectionMatrix;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType.Hybridization;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.Atom;
import org.openscience.cdk.silent.Bond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class ExtendedConnectionMatrix {

    private final double[][] connectionMatrix;
    private final String[] atomTypes;
    private final Integer[][] atomPropertiesNumeric;// hydrogenCounts, valencies, formalCharges;
    private final Hybridization[] hybridizations;
    private final Boolean[][] atomPropertiesBoolean;// isInRingAtoms, isAromaticAtoms;
    private final Boolean[][][] bondProperties;

    public ExtendedConnectionMatrix(final IAtomContainer ac){
        this.connectionMatrix = ConnectionMatrix.getMatrix(ac);
        this.atomTypes = new String[this.connectionMatrix.length];
        this.hybridizations = new Hybridization[this.connectionMatrix.length];
        this.atomPropertiesNumeric = new Integer[this.connectionMatrix.length][];
        this.atomPropertiesBoolean = new Boolean[this.connectionMatrix.length][];
        this.bondProperties = new Boolean[this.connectionMatrix.length][][];

        this.init(ac);
    }
    
    private void init(final IAtomContainer ac){

        IAtom atom1, atom2;
        IBond bond;
        for (int i = 0; i < this.connectionMatrix.length; i++) {
            atom1 = ac.getAtom(i);
            this.atomTypes[i] = atom1.getSymbol();
            this.atomPropertiesNumeric[i] = new Integer[3];
            this.atomPropertiesNumeric[i][0] = atom1.getImplicitHydrogenCount();
            this.atomPropertiesNumeric[i][1] = atom1.getValency();
            this.atomPropertiesNumeric[i][2] = atom1.getFormalCharge();
            this.atomPropertiesBoolean[i] = new Boolean[2];
            this.atomPropertiesBoolean[i][0] = atom1.isInRing();
            this.atomPropertiesBoolean[i][1] = atom1.isAromatic();
            this.hybridizations[i] = atom1.getHybridization();

            this.bondProperties[i] = new Boolean[this.connectionMatrix.length][2];
            for (int k = 0; k < this.connectionMatrix.length; k++) {
                atom2 = ac.getAtom(k);
                bond = ac.getBond(atom1, atom2);
                if(bond != null){
                    this.bondProperties[i][k][0] = bond.isInRing();
                    this.bondProperties[i][k][1] = bond.isAromatic();
                }
            }
        }

    }
    
    public IAtomContainer toAtomContainer(){
        final IAtomContainer ac = SilentChemObjectBuilder.getInstance().newAtomContainer();
        IAtom atom;
        for (int i = 0; i < this.connectionMatrix.length; i++) {
            atom = new Atom(this.atomTypes[i]);
            atom.setImplicitHydrogenCount(this.atomPropertiesNumeric[i][0]);
            atom.setValency(this.atomPropertiesNumeric[i][1]);
            atom.setFormalCharge(this.atomPropertiesNumeric[i][2]);
            atom.setHybridization(this.hybridizations[i]);
            atom.setIsInRing(this.atomPropertiesBoolean[i][0]);
            atom.setIsAromatic(this.atomPropertiesBoolean[i][1]);

            ac.addAtom(atom);
        }
        IBond bond;
        for (int i = 0; i < this.bondProperties.length; i++) {
            for (int k = i + 1; k < this.bondProperties.length; k++) {
                if(this.connectionMatrix[i][k] > 0.0){
                    bond = new Bond(ac.getAtom(i), ac.getAtom(k), Utils.getBondOrder((int) this.connectionMatrix[i][k]));
                    bond.setIsInRing(this.bondProperties[i][k][0]);
                    bond.setIsAromatic(this.bondProperties[i][k][1]);
                    ac.addBond(bond);
                }
            }
        }
        
        return ac;
    }
    
}
