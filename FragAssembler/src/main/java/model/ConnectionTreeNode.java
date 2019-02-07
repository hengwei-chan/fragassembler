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

import java.util.ArrayList;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class ConnectionTreeNode {

    private final IAtom atom;    
    private final ArrayList<ConnectionTreeNode> parents;
    private final ArrayList<IBond> bondsToParents;
    private final ArrayList<ConnectionTreeNode> children;
    private final ArrayList<IBond> bondsToChildren;
    private final int key;
    private final int sphere;
    private final boolean isRingClosure;

    public ConnectionTreeNode(final IAtom atom, final int key, final int sphere, final boolean isRingClosure) {        
        this.atom = atom;                
        this.key = key;
        this.sphere = sphere;
        this.isRingClosure = isRingClosure;
        this.parents = new ArrayList<>();
        this.children = new ArrayList<>();
        this.bondsToParents = new ArrayList<>();
        this.bondsToChildren = new ArrayList<>();
    }

    public IAtom getAtom(){
        return this.atom;
    }
    
    public ArrayList<ConnectionTreeNode> getParentNodes(){
        return this.parents;
    }
    
    public ArrayList<IBond> getBondsToParents(){
        return this.bondsToParents;
    }
    
    public ArrayList<ConnectionTreeNode> getChildNodes(){
        return this.children;
    }
    
    public ArrayList<IBond> getBondsToChildren(){
        return this.bondsToChildren;
    }
    
    public int getKey(){
        return this.key;
    }
    
    public int getSphere(){
        return this.sphere;
    }
    
    public boolean isRingClosureNode(){
        return this.isRingClosure;
    }
    
//    public boolean isPlaceholderNode(){
//        return (this.atom == null) && !this.isRingClosure;
//    }
    
    public void addParentNode(final ConnectionTreeNode parentNode, final IBond bondToParent){
        this.getParentNodes().add(parentNode);
        this.getBondsToParents().add(bondToParent);
    }
    
    public void addChildNode(final ConnectionTreeNode childNode, final IBond bondToChild){
        this.getChildNodes().add(childNode);
        this.getBondsToChildren().add(bondToChild);
    }
    
    public boolean hasParent(final int parentKey){
        for (final ConnectionTreeNode parentNode : this.parents) {
            if(parentNode.getKey() == parentKey){
                return true;
            }
        }
        
        return false;
    }
    
    public boolean hasParents(){
        return !this.parents.isEmpty();
    }
    
    public boolean hasChild(final int childKey){
        for (final ConnectionTreeNode childNode : this.children) {
            if(childNode.getKey() == childKey){
                return true;
            }
        }
        
        return false;
    }
    
    public boolean hasChildren(){
        return !this.children.isEmpty();
    }
}
