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

import hose.HOSECodeBuilder;
import hose.model.ConnectionTree;
import hose.model.ConnectionTreeNode;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class ExtendedConnectionMatrixHOSECode extends ExtendedConnectionMatrix {

    private final int rooAtomIndex;
    private String[] HOSECodes;
    private int[][] atomIndicesInHOSECodeOrder; // atom index, atom indices in HOSE code order
    private int[][] spheresInHOSECodeOrder; // atom index, atom indices in HOSE code order
    private HashMap<Integer, ArrayList<Integer[]>> ringClosures;
    private Integer maxSphere;
    private Integer[] attachedHydrogensCountsInOuterSphere;

    public ExtendedConnectionMatrixHOSECode(final IAtomContainer ac, final int rootAtomIndex, final Integer maxSphere) {
        super(ac);
        this.rooAtomIndex = rootAtomIndex;
        this.maxSphere = maxSphere;
        this.update();
    }

    private void update() {
        // initialize all
        this.HOSECodes = new String[this.getAtomCount()];
        this.atomIndicesInHOSECodeOrder = new int[this.getAtomCount()][this.getAtomCount()];
        this.spheresInHOSECodeOrder = new int[this.getAtomCount()][this.getAtomCount()];
        this.ringClosures = new HashMap<>();

        for (int i = 0; i < this.getAtomCount(); i++) {
            for (int j = 0; j < this.getAtomCount(); j++) {
                this.atomIndicesInHOSECodeOrder[i][j] = -1;
                this.spheresInHOSECodeOrder[i][j] = -1;
                this.ringClosures.put(i, new ArrayList<>());
            }
        }
        // update HOSE code information
        // @TODO replace IAtomContainer (ac) usage in BFS by usage of ExtendedConnectionMatrix
        final IAtomContainer ac = this.toAtomContainer();
        ConnectionTree connectionTree;
        ArrayList<ConnectionTreeNode> nodes;
        ConnectionTreeNode node;
        // go through all atoms
        for (int i = 0; i < this.getAtomCount(); i++) {
            connectionTree = HOSECodeBuilder.buildConnectionTree(ac, i, this.maxSphere);
            // update max. sphere to last sphere in conn. tree which contains non-ring closure nodes
            if(i == this.getRootAtomIndex()){
                for (int s = 0; s <= connectionTree.getMaxSphere(); s++) {
                    if(connectionTree.getNodesInSphere(s, false).isEmpty()){
                        break;
                    }
                    this.maxSphere = s;
                }
            }
            // store HOSE codes
            try {
                this.HOSECodes[i] = HOSECodeBuilder.buildHOSECode(connectionTree, false);
            } catch (CDKException e) {
                this.HOSECodes[i] = null;
                e.printStackTrace();
            }
            // store atom indices in HOSE code order as well as ring closure and sphere information
            nodes = connectionTree.getNodes(true);
            int indexCounter = 0;
            for (int j = 0; j < nodes.size(); j++) {
                node = nodes.get(j);
                if(node.isRingClosureNode()){
                    this.ringClosures.get(i).add(new Integer[]{node.getParent().getKey(), node.getRingClosureParent().getKey()});
                    continue;
                }
                this.atomIndicesInHOSECodeOrder[i][indexCounter] = node.getKey();
                this.spheresInHOSECodeOrder[i][indexCounter] = node.getSphere();
                indexCounter++;
            }
        }

        this.updateAttachedHydrogensCountsInOuterSphere();
    }

    public int getRootAtomIndex() {
        return this.rooAtomIndex;
    }

    public int getMaxSphere(){
        return this.maxSphere;
    }


    public String getHOSECode(final int atomIndex) {
        if(!this.hasAtom(atomIndex)){
            return null;
        }

        return this.HOSECodes[atomIndex];
    }

    public int[] getAtomIndicesInHOSECodeOrder(final int rooAtomIndex) {
        if(!this.hasAtom(rooAtomIndex)){
            return null;
        }

        return this.atomIndicesInHOSECodeOrder[rooAtomIndex];
    }

    public int[] getSpheresInHOSECodeOrder(final int rooAtomIndex) {
        if(!this.hasAtom(rooAtomIndex)){
            return null;
        }

        return this.spheresInHOSECodeOrder[rooAtomIndex];
    }

    public ArrayList<Integer[]> getRingClosures(final int rooAtomIndex) {
        if(!this.hasAtom(rooAtomIndex)){
            return null;
        }

        return this.ringClosures.get(rooAtomIndex);
    }

    public Integer[] getAttachedHydrogensCountsInOuterSphere(){
        return this.attachedHydrogensCountsInOuterSphere;
    }

    public Integer getSphere(final int rootAtomIndex, final int atomIndex){
        if(!this.hasAtom(rootAtomIndex) || !this.hasAtom(atomIndex)){
            return null;
        }

        return this.getSpheresInHOSECodeOrder(rootAtomIndex)[atomIndex];
    }

    public Integer getParent(final int rooAtomIndex, final int atomIndex) {
        if(!this.hasAtom(rooAtomIndex) || !this.hasAtom(atomIndex)){
            return null;
        }

        final int sphere = this.getSphere(rooAtomIndex, atomIndex);
        for (int s = 0; s < this.getSpheresInHOSECodeOrder(rooAtomIndex).length; s++) {
            // out of maximum spherical limit
            if(this.getSpheresInHOSECodeOrder(rooAtomIndex)[s] == -1){
                return -1;
            }
            if(this.getSpheresInHOSECodeOrder(rooAtomIndex)[s] == (sphere - 1)){
                if(this.hasBond(this.getAtomIndicesInHOSECodeOrder(rooAtomIndex)[s], atomIndex)){
                    return this.getAtomIndicesInHOSECodeOrder(rooAtomIndex)[s];
                }
            }
        }

        return -1;
    }

    public Boolean hasParent(final int rooAtomIndex, final int parentAtomIndex, final int childAtomIndex) {
        if(!this.hasAtom(rooAtomIndex) || !this.hasAtom(parentAtomIndex) || !this.hasAtom(childAtomIndex)){
            return null;
        }

        return this.getParent(rooAtomIndex, childAtomIndex) == parentAtomIndex;
    }

    public Integer[] getChildren(final int rooAtomIndex, final int atomIndex){
        if(!this.hasAtom(rooAtomIndex) || !this.hasAtom(atomIndex)){
            return null;
        }

        final int sphere = this.getSphere(rooAtomIndex, atomIndex);
        final ArrayList<Integer> childrenAsList = new ArrayList<>();
        for (int s = 0; s < this.getSpheresInHOSECodeOrder(rooAtomIndex).length; s++) {
            // out of maximum spherical limit
            if(this.getSpheresInHOSECodeOrder(rooAtomIndex)[s] == -1){
                break;
            }
            if(this.getSpheresInHOSECodeOrder(rooAtomIndex)[s] == (sphere + 1)){
                if(this.hasBond(this.getAtomIndicesInHOSECodeOrder(rooAtomIndex)[s], atomIndex)){
                    childrenAsList.add(this.getAtomIndicesInHOSECodeOrder(rooAtomIndex)[s]);
                }
            }
        }
//        final int[] children = new int[childrenAsList.size()];
//        for (int k = 0; k < children.length; k++) {
//            children[k] = childrenAsList.get(k);
//        }

        return childrenAsList.toArray(new Integer[]{});
    }

    public Boolean hasChild(final int rooAtomIndex, final int parentAtomIndex, final int childAtomIndex) {
        if(!this.hasAtom(rooAtomIndex) || !this.hasAtom(parentAtomIndex) || !this.hasAtom(childAtomIndex)){
            return null;
        }
        final Integer[] children = this.getChildren(rooAtomIndex, parentAtomIndex);
        for (int k = 0; k < children.length; k++) {
            if(children[k] == childAtomIndex){
                return true;
            }
        }

        return false;
    }

    private void updateAttachedHydrogensCountsInOuterSphere(){
        final ArrayList<Integer> atomIndicesInOuterSphere = new ArrayList<>();
        for (int i = 0; i < this.getAtomCount(); i++) {
            if(this.getSphere(this.getRootAtomIndex(), i) == this.getMaxSphere()){
                atomIndicesInOuterSphere.add(i);
            }
        }
        this.attachedHydrogensCountsInOuterSphere = new Integer[atomIndicesInOuterSphere.size()];
        for (int k = 0; k < atomIndicesInOuterSphere.size(); k++) {
            this.attachedHydrogensCountsInOuterSphere[k] = this.getHydrogenCount(atomIndicesInOuterSphere.get(k));
        }
    }

    public Boolean formRingClosure(final int rootAtomIndex, final int ringClosureAtomIndex1, final int ringClosureAtomIndex2){
        if(!this.hasAtom(rooAtomIndex) || !this.hasAtom(ringClosureAtomIndex1) || !this.hasAtom(ringClosureAtomIndex2)){
            return null;
        }
        for (final Integer[] ringClosurePartner : (this.getRingClosures(rootAtomIndex))){
            if((ringClosurePartner[0] == ringClosureAtomIndex1) && (ringClosurePartner[1] == ringClosureAtomIndex2)
                    || (ringClosurePartner[1] == ringClosureAtomIndex1) && (ringClosurePartner[0] == ringClosureAtomIndex2)
            ){
                return true;
            }
        }

        return false;
    }

    @Override
    public void addAtom(final String atomType, final Integer implicitHydrogenCount, final Integer valency, final Integer formalCharge, final Boolean isInRing, final Boolean isAromatic, final IAtomType.Hybridization hybridization){
        super.addAtom(atomType, implicitHydrogenCount, valency, formalCharge, isInRing, isAromatic, hybridization);
        this.update();
    }

    @Override
    public boolean addBond(final int atomIndex1, final int atomIndex2, final double order, final Boolean isInRing, final Boolean isAromatic){
        if(!super.addBond(atomIndex1, atomIndex2, order,isInRing, isAromatic)){
            return false;
        }
        this.update();

        return true;
    }

    @Override
    public String toString(){
        final StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append("atom count: ").append(this.getAtomCount()).append(", bond count: ").append(this.getBondCount()).append("\n");
        stringBuilder.append("max. sphere: " + this.getMaxSphere()).append("\n");
        for (int i = 0; i < this.getAtomCount(); i++) {
            stringBuilder.append("for atom: ").append(i).append("\n");
            stringBuilder.append(" ->").append(this.HOSECodes[i]).append("\n");
            stringBuilder.append(" ->").append(Arrays.toString(this.atomIndicesInHOSECodeOrder[i])).append("\n");
            stringBuilder.append(" ->").append(Arrays.toString(this.spheresInHOSECodeOrder[i])).append("\n");
            for (final Integer[] ringClosure : this.ringClosures.get(i)){
                stringBuilder.append(" -->").append(Arrays.toString(ringClosure)).append("\n");
            }
        }

        return stringBuilder.toString();
    }

    /**
     * Returns the SSC substructure as HOSE code with the beginning at the SSC's root atom.
     *
     * @return HOSE code of SSC's root atom
     *
     */
    public String getAsHOSECode() {
        return this.getHOSECode(this.getRootAtomIndex());
    }

    public ExtendedConnectionMatrixHOSECode getClone(){
        return new ExtendedConnectionMatrixHOSECode(this.toAtomContainer(), this.getRootAtomIndex(), this.getMaxSphere());
    }
}

