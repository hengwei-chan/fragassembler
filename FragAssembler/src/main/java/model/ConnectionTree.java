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

import hose.HOSECodeBuilder;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class ConnectionTree implements Cloneable {
    private final ConnectionTreeNode root;
    private final LinkedHashSet<Integer> keySet;
    private int maxSphere;

    public ConnectionTree(final IAtom rootAtom, final int key, final int startSphere) {
        this.root = new ConnectionTreeNode(rootAtom, key, startSphere, false); 
        this.keySet = new LinkedHashSet<>();
        this.keySet.add(this.root.getKey());
        this.maxSphere = 0;
    }
    
    public ConnectionTreeNode getRootNode(){
        return this.root;
    }
    
    public boolean addNode(final IAtom newNodeData, final int newNodeKey, final int parentNodeKey, final IBond bond, final int sphere, final boolean isRingClosure){ 
        if(this.containsKey(newNodeKey) && (!isRingClosure)){
            return false;
        }      
        final ConnectionTreeNode parentNode = this.getNode(parentNodeKey);
        final ConnectionTreeNode newNode = new ConnectionTreeNode(newNodeData, newNodeKey, sphere, isRingClosure);
        newNode.addParentNode(parentNode, bond);
        parentNode.addChildNode(newNode, bond);
        this.keySet.add(newNodeKey);
        if (sphere > this.maxSphere) {
            this.maxSphere = sphere;
        }
        
        return true;
    }
    
    public int getMaxSphere(){
        return this.maxSphere;
    }
    
    public int getNodesCount(){
        return this.keySet.size();
    }
    
    public int getNodesCountInSphere(final int sphere){
        return this.getNodesInSphere(sphere).size();
    }
    
    public LinkedHashSet<Integer> getKeys(final boolean withoutRingClosureNodes){
        final LinkedHashSet<Integer> keys = new LinkedHashSet<>();
        for (int s = 0; s <= this.getMaxSphere(); s++) {
            for (final ConnectionTreeNode nodeInSphere : this.getNodesInSphere(s)) {
                if (withoutRingClosureNodes) {
                    if (!nodeInSphere.isRingClosureNode()) {
                        keys.add(nodeInSphere.getKey());
                    }
                } else {
                    keys.add(nodeInSphere.getKey());
                }
            }
        }
        
        return keys;
    }
   
    public boolean containsKey(final int key){
        return this.keySet.contains(key);
    }
    
    public ConnectionTreeNode getNode(final int key){
        if(!this.containsKey(key)){
            return null;
        }
        
        return this.findNode(key, this.root);
    }
    
    private ConnectionTreeNode findNode(final int key, final ConnectionTreeNode currentNode){
        if(currentNode.getKey()== key){
            return currentNode;
        }
        ConnectionTreeNode result = null;
        for (final ConnectionTreeNode childNode : currentNode.getChildNodes()) {
            result = this.findNode(key, childNode);
            if((result != null) && (result.getKey() == key)){
                break;
            }
        }
            
        return result;
    }
    
    public int getNodeIndexInSphere(final ConnectionTreeNode node, final int sphere){
        
        return this.getNodesInSphere(sphere).indexOf(node);
    }
    
    public ArrayList<Integer> getNodeKeysInSphere(final int sphere){
        final ArrayList<Integer> keys = new ArrayList<>();
        for (final ConnectionTreeNode treeNode : this.getNodesInSphere(sphere)) {
            keys.add(treeNode.getKey());
        }
        
        return keys;
    }
        
    public ArrayList<ConnectionTreeNode> getNodesInSphere(final int sphere){        
        return this.findNodesInSphere(sphere, this.root, new ArrayList<>());
    }
    
    private ArrayList<ConnectionTreeNode> findNodesInSphere(final int sphere, final ConnectionTreeNode currentNode, final ArrayList<ConnectionTreeNode> indicesInSphere){
        if(currentNode.getSphere() == sphere){
            indicesInSphere.add(currentNode);
            return indicesInSphere;
        }
        for (final ConnectionTreeNode childNode : currentNode.getChildNodes()) {
            this.findNodesInSphere(sphere, childNode, indicesInSphere);
        }
            
        return indicesInSphere;
    }    
    
    public IBond getBond(final int parentKey, final int childKey){        
        if(!this.containsKey(parentKey) || !this.containsKey(childKey)){
            return null;
        }
        final ConnectionTreeNode parentNode = this.getNode(parentKey);
        for (int i = 0; i < parentNode.getChildNodes().size(); i++) {
            if(parentNode.getChildNodes().get(i).getKey() == childKey){
                return parentNode.getBondsToChildren().get(i);
            }
        }
        
        return null;
    }
    
    public boolean hasParent(final int key, final int parentKey) { 
        if (!this.containsKey(key) || !this.containsKey(parentKey)) {
            return false;
        }
        
        return this.getNode(key).hasParent(parentKey);
    }

    public boolean hasChild(final int key, final int childKey) {
        if (!this.containsKey(key) || !this.containsKey(childKey)) {
            return false;
        }
        
        return this.getNode(key).hasChild(childKey);
    } 
    
    @Override
    public String toString(){
        String treeString = "";
        for (int s = 0; s <= this.maxSphere; s++) {
            treeString += "[" + s  + "]";
            for (final ConnectionTreeNode nodeInSphere : this.getNodesInSphere(s)) {
                treeString += " ";
                treeString += nodeInSphere.getKey() + ": ";
                if (nodeInSphere.getParentNodes().size() >= 1) {
                    treeString += HOSECodeBuilder.getSymbolForBond(getBond(nodeInSphere.getParentNodes().get(0).getKey(), nodeInSphere.getKey()));
                }
                if(nodeInSphere.getAtom() != null){
                    treeString += nodeInSphere.getAtom().getSymbol();
                } else if(nodeInSphere.isRingClosureNode()){
                    treeString += "&";
                }
                treeString += " {";
                
                if (nodeInSphere.getParentNodes().size() > 1) {
                    for (final ConnectionTreeNode parentNode : nodeInSphere.getParentNodes()) {                        
                        treeString += parentNode.getKey() + " ";
                    }
                } else if(nodeInSphere.getParentNodes().size() == 1){
                    treeString += nodeInSphere.getParentNodes().get(0).getKey();
                }
                treeString += "} ";
            }
        }
        
        return treeString;
    }
    
    @Override
    public ConnectionTree clone() throws CloneNotSupportedException{
        return (ConnectionTree) super.clone();
    }
}
