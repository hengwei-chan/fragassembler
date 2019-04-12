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
package hose;

import casekit.NMR.Utils;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Queue;
import model.ConnectionTree;
import model.ConnectionTreeNode;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.Atom;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

/**
 *
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class HOSECodeBuilder {
    
    /**
     * Returns the number of non-empty spheres. 
     * For example: C-3;() -> 1, C-3;=N(C/) -> 3
     *
     * @param HOSECode
     * @return
     */
    public static int getSpheresCount(final String HOSECode){
        int spheresCount = 0;
        for (final String sphere : HOSECodeBuilder.splitHOSECodeIntoSpheres(HOSECode)) {
            if (!sphere.trim().isEmpty()) {
                spheresCount++;
            }
        }
        return spheresCount;
    }
    
    /**
     * Splits a HOSE code into a list of spheres as strings.
     *
     * @param HOSECode HOSE code
     * @return ArrayList of all sphere strings
     */
    public static ArrayList<String> splitHOSECodeIntoSpheres(final String HOSECode) {
        final ArrayList<String> HOSECodeSpheres = new ArrayList<>();
        final String[] splitSpheres_0_1 = HOSECode.split(";");
        final String[] splitSpheres_1_2 = splitSpheres_0_1[1].split("\\(");
        final String[] splitSpheres_2_n = splitSpheres_1_2[1].substring(0, splitSpheres_1_2[1].length() - 1).split("\\/");        
        HOSECodeSpheres.add(splitSpheres_0_1[0]);
        HOSECodeSpheres.add(splitSpheres_1_2[0]);
        for (int s = 0; s < splitSpheres_2_n.length; s++) {
            HOSECodeSpheres.add(splitSpheres_2_n[s]);
        }
       
        return HOSECodeSpheres;
    }

    /**
     * Splits a HOSE code sphere into its positions. Each position includes all
     * its elements.
     * Example: {@code /CC,*N&,C/} results in: {@code {0: [C,C], 1: [*N,&], 2: [C]}}
     *
     * @param HOSECodeSphere  HOSE code sphere
     * @return HashMap of ArrayLists containing elements for each position of that HOSE code sphere
     */
    public static HashMap<Integer, ArrayList<String>> splitHOSECodeSphereIntoPositions(final String HOSECodeSphere) {
        final HashMap<Integer, ArrayList<String>> positions = new HashMap<>();
        // zeroth sphere
        if (HOSECodeSphere.contains("-")) {
            positions.put(0, new ArrayList<>());
            positions.get(0).add(HOSECodeSphere.split("-")[0]);
            positions.get(0).add(HOSECodeSphere.split("-")[1]);
            return positions;
        }
        // higher spheres
        char c;
        String elem = "";
        int positionCounter = 0;
        for (int i = 0; i < HOSECodeSphere.length(); i++) {
            c = HOSECodeSphere.charAt(i);
            if ((c == '=') || (c == '%') || (c == '*')) {
                if (!elem.isEmpty()) {
                    if (!positions.containsKey(positionCounter)) {
                        positions.put(positionCounter, new ArrayList<>());
                    }
                    positions.get(positionCounter).add(elem);
                    elem = "";
                }
                elem += c;
            } else if (Character.isUpperCase(c) || (c == '&')) {
                if (!elem.isEmpty() && (Character.isLetter(elem.charAt(elem.length() - 1)) || (elem.charAt(elem.length() - 1) == '&'))) {
                    if (!positions.containsKey(positionCounter)) {
                        positions.put(positionCounter, new ArrayList<>());
                    }
                    positions.get(positionCounter).add(elem);
                    elem = "";
                    elem += c;
                } else {
                    elem += c;
                }
            } else if (Character.isLowerCase(c)) {
                elem += c;
            } else if (c == ',') {
                if (!positions.containsKey(positionCounter)) {
                    positions.put(positionCounter, new ArrayList<>());
                }
                if (elem.isEmpty()) {
                    positions.get(positionCounter).add(null);
                } else {
                    positions.get(positionCounter).add(elem);
                    elem = "";
                }
                positionCounter++;
            } 
        }
        // add last element
        if (!elem.isEmpty()) {
            if (!positions.containsKey(positionCounter)) {
                positions.put(positionCounter, new ArrayList<>());
            }
            positions.get(positionCounter).add(elem);
        } else if (HOSECodeSphere.endsWith(",")) {
            if (!positions.containsKey(positionCounter)) {
                positions.put(positionCounter, new ArrayList<>());
            }
            positions.get(positionCounter).add(null);
        }
        
        return positions;
    }

    /**
     * Returns the weight/cost for an HOSE code symbol regarding its priority.
     *
     * @param symbol HOSE code symbol
     * @return weight/cost for the symbol
     */
    private static int getSymbolPriorityWeight(final String symbol) {
        switch (symbol) {            
            case "%":
                return 20;
            case "=":
                return 19;
            case "*":
                return 18;
            case "C":
                return 17;
            case "O":
                return 16;
            case "N":
                return 15;
            case "S":
                return 14;
            case "P":
                return 13;
            case "Si":
            case "Q":
                return 12;
            case "B":
                return 11;
            case "F":
                return 10;
            case "Cl":
            case "X":
                return 9;
            case "Br":
            case "Y":
                return 8;
            case ";":
                return 7;
            case "I":
                return 6;
            case "#":
                return 5;
            case "&":
                return 4;
            case ",":
                return 3;
        }

        return 0;
    }

    /**
     * Converts an element symbol into notation as shown in origin article by 
     * Bremser.
     * That includes: Si -> Q, Cl -> X, Br -> Y
     *
     * @param element
     * @return HOSE code symbol as in origin article by Bremser
     */
    public static String toHOSECodeSymbol(final String element) {
        if (element.equals("Si")) {
            return "Q";
        }
        if (element.equals("Cl")) {
            return "X";
        }
        if (element.equals("Br")) {
            return "Y";
        }

        return element;
    }
    
    /**
     * Converts an HOSE code symbol as shown in origin article by 
     * Bremser into default element notation.
     * That includes: Q -> Si, X -> Cl, Y -> Br
     *
     * @param HOSECodeSymbol
     * @return default element notation
     */
    public static String toElementSymbol(final String HOSECodeSymbol) {
        if (HOSECodeSymbol.equals("Q")) {
            return "Si";
        }
        if (HOSECodeSymbol.equals("X")) {
            return "Cl";
        }
        if (HOSECodeSymbol.equals("Y")) {
            return "Br";
        }

        return HOSECodeSymbol;
    }
    
    /**
     * Returns the notation of bond information used in HOSE code.
     * The bond has to contain its bond order and aromaticity information.
     *
     * @param bond bond containing bond order and aromaticity information
     * @return HOSE code symbol for a bond
     */
    public static String getSymbolForBond(final IBond bond){
        if(bond != null){
            if (bond.isAromatic()) {
                return "*";
            }
            switch (bond.getOrder()) {
                case SINGLE:
                    return "";
                case DOUBLE:
                    return "=";
                case TRIPLE:
                    return "%";
            }  
        }        
        
        return null;
    }    
    
    /**
     * Returns the bond order from an HOSE code bond symbol.
     * One has to consider that in this direction the aromatic HOSE code symbol 
     * (*) is ambiguous. It means either a single or a double aromatic bond. 
     * For that case, a single bond will always be returned.  
     *
     * @param symbol HOSE code bond symbol
     * @return bond order for a bond symbol
     */
    public static IBond.Order getBondOrderForSymbol(final String symbol) {
        if (symbol == null) {
            return null;
        }
        if(symbol.isEmpty() || symbol.equals("*")){ 
            return IBond.Order.SINGLE;
        }
        
        return Utils.getBondOrderFromString(symbol);
    }
    
    /**
     * Returns the summed subtree weight starting at a specific node in a connection 
     * tree. The weight of starting node is included here.
     *
     * @param connectionTree
     * @param node
     * @return
     *
     */
    public static int calculateSubtreeWeight(final ConnectionTree connectionTree, final ConnectionTreeNode node) {
        return HOSECodeBuilder.getSubtreeWeight(connectionTree, node, null);
    }

    /**
     *
     * @param connectionTree
     * @param node
     * @param parentNode
     * @return
     *
     */
    private static int getSubtreeWeight(final ConnectionTree connectionTree, final ConnectionTreeNode node, final ConnectionTreeNode parentNode) {
        int weight = HOSECodeBuilder.getNodeWeight(connectionTree, node, parentNode);
        for (final ConnectionTreeNode childNode : node.getChildNodes()) {
            weight += HOSECodeBuilder.getSubtreeWeight(connectionTree, childNode, node);
        }

        return weight;
    }
    
    /**
     * Returns the weight for a node and its connection to a parent node 
     * (optional).
     *
     * @param connectionTree connection tree 
     * @param node node to get the weight from
     * @param parentNode parent node of node or null
     * @return the priority weight for node; plus the weight of 
     * the bond to its parent node if the parent node is not null
     * 
     * @see #getSymbolPriorityWeight(java.lang.String) 
     */
    public static Integer getNodeWeight(final ConnectionTree connectionTree, final ConnectionTreeNode node, final ConnectionTreeNode parentNode){
        int weight = 0;
        if (parentNode != null) {
            String bondSymbol = HOSECodeBuilder.getSymbolForBond(connectionTree.getBond(parentNode.getKey(), node.getKey()));
            if(bondSymbol == null){
                return null;
            }
            // add weight for bond type priority 
            if (!bondSymbol.isEmpty()) {
                weight += HOSECodeBuilder.getSymbolPriorityWeight(bondSymbol);
            }
            
        }
        // add weight for further symbol priority    
        if (node.isRingClosureNode()) {
            weight += HOSECodeBuilder.getSymbolPriorityWeight("&");
        } else {
            weight += HOSECodeBuilder.getSymbolPriorityWeight(node.getAtom().getSymbol());
        }
        
        return weight;
    }
    
    /**
     * Returns an ArrayList of ranked child node indices for a tree node.
     *
     * @param connectionTree
     * @param node
     * @return
     * 
     * @see #getNodeWeight(model.ConnectionTree, model.ConnectionTreeNode, model.ConnectionTreeNode) 
     */
    private static ArrayList<Integer> getRankedChildNodesIndices(final ConnectionTree connectionTree, final ConnectionTreeNode node){
        final ArrayList<ConnectionTreeNode> childNodes = node.getChildNodes();
        final ArrayList<Integer> rankedChildNodesIndices = new ArrayList<>();
        for (int i = 0; i < childNodes.size(); i++) {
            rankedChildNodesIndices.add(i);
        }
        rankedChildNodesIndices.sort(new Comparator<Integer>() {
            @Override
            public int compare(final Integer childNodeIndex1, final Integer childNodeIndex2) {      
                
                int nodeWeightsDifference = -1 * Integer.compare(HOSECodeBuilder.getNodeWeight(connectionTree, childNodes.get(childNodeIndex1), node), 
                        HOSECodeBuilder.getNodeWeight(connectionTree, childNodes.get(childNodeIndex2), node));
                if(nodeWeightsDifference != 0){
                    return nodeWeightsDifference;
                }
//                return -1 * Integer.compare(HOSECodeBuilder.getNodeWeight(connectionTree, childNodes.get(childNodeIndex1), node), 
//                        HOSECodeBuilder.getNodeWeight(connectionTree, childNodes.get(childNodeIndex2), node));
                return -1 * Integer.compare(HOSECodeBuilder.calculateSubtreeWeight(connectionTree, childNodes.get(childNodeIndex1)),
                        HOSECodeBuilder.calculateSubtreeWeight(connectionTree, childNodes.get(childNodeIndex2)));
            }
        });
        
        return rankedChildNodesIndices;
    }
    
    /**
     * Sorts the child nodes of a node by HOSE code priority and weight.
     *
     * @param connectionTree connection tree containing the node
     * @param node node with child nodes to rank
     * 
     * @see #getNodeWeight(model.ConnectionTree, model.ConnectionTreeNode, model.ConnectionTreeNode) 
     */
    public static void rankChildNodes(final ConnectionTree connectionTree, final ConnectionTreeNode node){
        final ArrayList<Integer> rankedChildNodesIndices = getRankedChildNodesIndices(connectionTree, node);
        final ArrayList<ConnectionTreeNode> rankedChildNodes = new ArrayList<>();
        final ArrayList<IBond> rankedChildNodeBonds = new ArrayList<>();
        for (int i = 0; i < rankedChildNodesIndices.size(); i++) {
            rankedChildNodes.add(node.getChildNodes().get(rankedChildNodesIndices.get(i)));
            rankedChildNodeBonds.add(node.getBondsToChildren().get(rankedChildNodesIndices.get(i)));
        }
        node.getChildNodes().clear();
        node.getBondsToChildren().clear();
        node.getChildNodes().addAll(rankedChildNodes);
        node.getBondsToChildren().addAll(rankedChildNodeBonds);        
    }
    
    /**
     * Sorts the child nodes of each node in the connection tree by HOSE code 
     * priority and weight.
     *
     * @param connectionTree connection tree where to rank the child nodes of 
     * each node.
     * 
     * @see #rankChildNodes(model.ConnectionTree, model.ConnectionTreeNode) 
     */
    public static void rankChildNodes(final ConnectionTree connectionTree) {
        ArrayList<ConnectionTreeNode> nodesInSphere;
        for (int sphere = 0; sphere < connectionTree.getMaxSphere(); sphere++) {
            nodesInSphere = connectionTree.getNodesInSphere(sphere);
            // for all nodes in sphere
            for (int i = 0; i < nodesInSphere.size(); i++) {
                // rank all child nodes of that node 
                if (nodesInSphere.get(i).hasChildren()) {
                    HOSECodeBuilder.rankChildNodes(connectionTree, nodesInSphere.get(i));
                } 
            }
        }
    }  
    
    /**
     * Returns the parent node list index for node2 in the parents list of 
     * node1.
     *
     * @param node1 first node which should contain node2 as parent in its
     * parent node list
     * @param node2 second node which should contain node1 as parent in its 
     * parent node list
     * @return index in parent nodes list of node1; -1 if node2 is not parent 
     * of node1 and/or vice versa
     */
    public static int getParentNodesListIndex(final ConnectionTreeNode node1, final ConnectionTreeNode node2){
        ConnectionTreeNode parentNode1, parentNode2;
        // in conn. tree, a ring closures is if there are nodes with more than one parent
        // and since only the parent information is stored and not the child information, it is possible to distinguish that
        if (HOSECodeBuilder.isAtRingClosure(node1)) {
            for (int i = 0; i < node1.getParentNodes().size(); i++) {
                parentNode1 = node1.getParentNodes().get(i);
                // if one of the parents of node1 is node2
                if(parentNode1.getKey() == node2.getKey()){
                    for (int j = 0; j < parentNode1.getParentNodes().size(); j++) {
                        parentNode2 = parentNode1.getParentNodes().get(j);
                        // cross validation: if the both know each other as parent; 
                        if ((parentNode2.getParentNodes().size() > 1) && (parentNode2.getKey() == node1.getKey())) {
                            return i;
                        }
                    }
                }                
            }
        }
        
        return -1;
    }
    
    public static boolean nodesFormRingClosure(final ConnectionTreeNode node1, final ConnectionTreeNode node2){                
        
        return HOSECodeBuilder.getParentNodesListIndex(node1, node2) != -1;
    }
    
    public static boolean isAtRingClosure(final ConnectionTreeNode node){        
        
        return node.getParentNodes().size() > 1;
    }
    
//    private static String buildRingClosureString(final ConnectionTreeNode node){
//        ConnectionTreeNode parentNode2;
//        IBond bond;
//        // in conn. tree, a ring closures is if there are nodes with more than one parent
//        // and since only the parent information is stored and not the child information, it is possible to distinguish that
//        if (node.getParentNodes().size() > 1) {
//            for (final ConnectionTreeNode parentNode : node.getParentNodes()) {
//                for (int j = 0; j < parentNode.getParentNodes().size(); j++) {
//                    parentNode2 = parentNode.getParentNodes().get(j);
//                    // cross validation: if the both know each other as parent; now take the bond information between them
//                    if ((parentNode2.getParentNodes().size() > 1) && (parentNode2.getKey() == node.getKey())) {
//                        bond = parentNode.getBondsToParents().get(j);
//                        return HOSECodeBuilder.getSymbolForBond(bond) + "&";// + node.getKey() + "-" + parentNode.getKey();
//                    }
//                }
//            }
//        }
//        
//        return "";
//    }

    private static ArrayList<String> buildPositionsInSphere(final ConnectionTreeNode nodeInPrevSphere, final boolean useBremserElementNotation) throws CDKException {
        final ArrayList<ConnectionTreeNode> nodesInSphere = nodeInPrevSphere.getChildNodes();
        final ArrayList<String> positionsInSphere = new ArrayList<>();
        ConnectionTreeNode nodeInSphere;
        IBond bond;
        String position;
        for (int j = 0; j < nodesInSphere.size(); j++) {
            nodeInSphere = nodesInSphere.get(j);
            bond = nodeInPrevSphere.getBondsToChildren().get(j);
            position = "";
            if (HOSECodeBuilder.getSymbolForBond(bond) == null) {
                throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": no bond information");
            }            
            position += HOSECodeBuilder.getSymbolForBond(bond);
            if(nodeInSphere.isRingClosureNode()){
                position += "&";// + nodeInPrevSphere.getKey() + "-" + nodeInPrevSphere.getParentNodes().get(1).getKey();
            } else {
                if (useBremserElementNotation) {
                    position += HOSECodeBuilder.toHOSECodeSymbol(nodeInSphere.getAtom().getSymbol());
                } else {
                    position += nodeInSphere.getAtom().getSymbol();
                }
            }
            positionsInSphere.add(position);
        }
        
        return positionsInSphere;
    }
    
    private static String buildSphereString(final ConnectionTree connectionTree, final int sphere, final String delimiter, final boolean useBremserElementNotation) throws CDKException{
        String sphereString = "";
        final ArrayList<ConnectionTreeNode> nodesInPrevSphere = connectionTree.getNodesInSphere(sphere - 1);
        ConnectionTreeNode nodeInPrevSphere;
        // for all nodes in previous sphere
        for (int i = 0; i < nodesInPrevSphere.size(); i++) {
            nodeInPrevSphere = nodesInPrevSphere.get(i);
            // skip ring closure nodes
            if(nodeInPrevSphere.isRingClosureNode()){
                if((i == nodesInPrevSphere.size() - 1) && sphereString.endsWith(",")){
                    sphereString = sphereString.substring(0, sphereString.length() - 1);
                }
                continue;
            }
            // for all child nodes in the actual requested sphere
            if(nodeInPrevSphere.hasChildren()){                
                for (final String position : HOSECodeBuilder.buildPositionsInSphere(nodeInPrevSphere, useBremserElementNotation)) {
                    sphereString += position;
                }
            } 
            // check for ring closure and add it
//            sphereString += HOSECodeBuilder.buildRingClosureString(nodeInPrevSphere);
                        
            // add delimiter
            if(i < nodesInPrevSphere.size() - 1){                
                sphereString += delimiter;
            } 
        }
        
        return sphereString;
    }
    
    private static String buildHOSECodeString(final ConnectionTree connectionTree, final boolean useBremserElementNotation) throws CDKException{
        final IAtom rootAtom = connectionTree.getRootNode().getAtom();
        final int maxSphere = connectionTree.getMaxSphere();
        // zeroth sphere
        String HOSECode = rootAtom.getSymbol() + "-" + (connectionTree.getRootNode().getChildNodes().size());// + rootAtom.getImplicitHydrogenCount());
        HOSECode += ";";
        String delimiter;
        // go through each sphere of the connection tree
        for (int s = 1; s <= maxSphere; s++) { 
            if (s == 1) {
                delimiter = "";
            } else {
                delimiter = ",";
            }
            // create sphere string and add it to HOSE code string
            HOSECode += HOSECodeBuilder.buildSphereString(connectionTree, s, delimiter, useBremserElementNotation);
            if (s == 1) {
                HOSECode += "(";
            } 
            if (s > 1 && s < maxSphere) {
                HOSECode += "/";
            }
        }
        if (maxSphere == 0) {
            HOSECode += "(";
        }
        HOSECode += ")";

        return HOSECode;
    }
    
    public static String buildHOSECode(final ConnectionTree connectionTree, final boolean useBremserElementNotation) throws CDKException {
        return HOSECodeBuilder.buildHOSECodeString(connectionTree, useBremserElementNotation);
    }
    
    public static String buildHOSECode(final IAtomContainer ac, final int rootAtomIndex, final Integer maxSphere, final boolean useBremserElementNotation) throws CDKException {
        return HOSECodeBuilder.buildHOSECodeString(HOSECodeBuilder.buildConnectionTree(ac, rootAtomIndex, maxSphere), useBremserElementNotation);
    }
    
    /**
     * Builds a connection tree of an atom container with specific start atom
     * and maximum number of spheres. 
     * If the atoms in the atom container are not fully connected, then the 
     * connection tree will be built until the last atom of all connected atoms 
     * to the start atom is reached.
     *
     * @param ac atom container
     * @param rootAtomIndex starting atom
     * @param maxSphere if this is set to null, then the connection tree of whole 
     * structure will be created
     * @return
     * 
     * @see model.ConnectionTree
     * @throws org.openscience.cdk.exception.CDKException
     */
    public static ConnectionTree buildConnectionTree(final IAtomContainer ac, final int rootAtomIndex, final Integer maxSphere) throws CDKException {            
        return HOSECodeBuilder.buildConnectionTree(ac, rootAtomIndex, maxSphere, new HashSet<>());
    }  
    
    /**
     * Builds a connection tree of an atom container with specific start atom
     * and maximum number of spheres. 
     * If the atoms in the atom container are not fully connected, then the 
     * connection tree will be built until the last atom of all connected atoms 
     * to the start atom is reached.
     *
     * @param ac atom container
     * @param rootAtomIndex starting atom
     * @param maxSphere if this is set to null, then the connection tree of whole 
     * structure will be created
     * @param visited certain atom indices can be given here to ignore atoms 
     * in BFS; they are then seen as already visited and not included in 
     * the connection tree
     * @return
     * 
     * @see model.ConnectionTree
     * @throws org.openscience.cdk.exception.CDKException
     */
    public static ConnectionTree buildConnectionTree(final IAtomContainer ac, final int rootAtomIndex, final Integer maxSphere, final HashSet<Integer> visited) throws CDKException {
        // create queue for BFS and add root atom index
        final Queue<Integer> queue = new LinkedList<>();
        queue.add(rootAtomIndex);
        final ConnectionTree connectionTree = new ConnectionTree(ac.getAtom(rootAtomIndex), rootAtomIndex, 0);
        HOSECodeBuilder.BFS(ac, connectionTree, queue, visited, maxSphere);

        HOSECodeBuilder.rankChildNodes(connectionTree);

        return connectionTree;
    }
    
    private static void BFS(final IAtomContainer ac, ConnectionTree connectionTree, final Queue<Integer> queue, final HashSet<Integer> visited, final Integer maxSphere){
        // all nodes visited?
        if(queue.isEmpty()){
            return;
        }
        final int atomIndex = queue.remove();
        final IAtom atom = ac.getAtom(atomIndex);  
        final ConnectionTreeNode node = connectionTree.getNode(atomIndex);
        final int sphere = node.getSphere();
        // check whether the current sphere is to high, if maxSphere parameter is set
        if((maxSphere != null) && (sphere > maxSphere)){
            return;
        }
        // mark atom as visited
        visited.add(atomIndex);     
                 
        IBond bond;
        if((maxSphere != null) && (sphere == maxSphere)){
            // set connections (parent nodes) in last sphere nodes which have to be connected -> ring closures
            // only parent nodes will be set to detect those ring closures again
            for (final ConnectionTreeNode nodeInLastSphere : connectionTree.getNodesInSphere(maxSphere)) {
                if ((ac.getBond(atom, nodeInLastSphere.getAtom()) != null)
                        && (!node.hasParent(nodeInLastSphere.getKey()))
                        && (!nodeInLastSphere.hasParent(node.getKey()))
                        ) {
                    bond = ac.getBond(node.getAtom(), nodeInLastSphere.getAtom());
                    node.addParentNode(nodeInLastSphere, bond);
                    nodeInLastSphere.addParentNode(node, bond);
                }
            }
        } else { 
            // add nodes and bonds in lower spheres
            // go to all child nodes
            int connectedAtomIndex;
            for (final IAtom connectedAtom : ac.getConnectedAtomsList(atom)) {
                connectedAtomIndex = ac.indexOf(connectedAtom);
                bond = ac.getBond(atom, connectedAtom);
                // add children to queue if not already visited
                if (!visited.contains(connectedAtomIndex)) {
                    // and not already waiting in queue
                    if (!queue.contains(connectedAtomIndex)) {
                        queue.add(connectedAtomIndex);
                        connectionTree.addNode(connectedAtom, connectedAtomIndex, node.getKey(), bond, sphere + 1, false);
                    } 
                    else {
                        // node already exists in tree; add a further parent to connected atom (for ring closures)
                        connectionTree.getNode(connectedAtomIndex).addParentNode(node, bond);
                        if(maxSphere != null){
                            if((sphere + 1 + 1) <= maxSphere){
                                connectionTree.addNode(null, -1 * connectedAtomIndex, connectedAtomIndex, bond, sphere + 1 + 1, true);
                            }                            
                        } else {
                            connectionTree.addNode(null, -1 * connectedAtomIndex, connectedAtomIndex, bond, sphere + 1 + 1, true);
                        }
                        node.addParentNode(connectionTree.getNode(connectedAtomIndex), bond);
                        connectionTree.addNode(null, -1 * node.getKey(), node.getKey(), bond, sphere + 1, true);
                    }
                }                 
            }
        }                   
        // further extension of connectivity tree
        BFS(ac, connectionTree, queue, visited, maxSphere);
    }
    
    public static IAtomContainer buildAtomContainer(final ConnectionTree connectionTree) throws CDKException{
        ConnectionTreeNode nodeInSphere, nodeInPrevSphere;
        // create new atom container and add root atom
        final IAtomContainer ac = SilentChemObjectBuilder.getInstance().newAtomContainer();
        ac.addAtom(connectionTree.getRootNode().getAtom());     
        // for each sphere: add the atom which is stored as node to atom container and set bonds between parent nodes                
        for (int s = 1; s <= connectionTree.getMaxSphere(); s++) {            
            // some of the parent atoms will be added to container later
            // so first add all atoms to structure
            for (int i = 0; i < connectionTree.getNodesInSphere(s).size(); i++) {
                nodeInSphere = connectionTree.getNodesInSphere(s).get(i);     
                if(!nodeInSphere.isRingClosureNode()){
                    ac.addAtom(nodeInSphere.getAtom());
                }                
            }
        }
        for (int s = 1; s <= connectionTree.getMaxSphere(); s++) {
            // and as second add the bonds to structure
            for (int i = 0; i < connectionTree.getNodesInSphere(s).size(); i++) {
                nodeInSphere = connectionTree.getNodesInSphere(s).get(i);
                if(nodeInSphere.isRingClosureNode()){
                   continue;
                }
                for (int j = 0; j < nodeInSphere.getParentNodes().size(); j++) {
                    nodeInPrevSphere = nodeInSphere.getParentNodes().get(j);
                    if(ac.getBond(nodeInSphere.getAtom(), nodeInPrevSphere.getAtom()) == null){
                        ac.addBond(ac.indexOf(nodeInSphere.getAtom()), ac.indexOf(nodeInPrevSphere.getAtom()), nodeInSphere.getBondsToParents().get(j).getOrder());
                    }                    
                }
            }
        }
        //AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(ac);
                
        return ac;
    }
    
    public static ConnectionTree buildConnectionTree(final String HOSECode, final boolean useBremserElementNotation) throws CDKException {
        final ArrayList<String> sphereStrings = HOSECodeBuilder.splitHOSECodeIntoSpheres(HOSECode);
        IAtom atom; IBond bond;
        int sphere, maxSphere;
        String bondTypeString, atomTypeString;
        HashMap<Integer, ArrayList<String>> positionsInSphere;
        // set maxSphere
        maxSphere = sphereStrings.size() - 1;
        // zeroth sphere
        sphere = 0;
        positionsInSphere = HOSECodeBuilder.splitHOSECodeSphereIntoPositions(sphereStrings.get(sphere));
        // charge: positionsInSphere.get(0).get(1)
        atom = new Atom(positionsInSphere.get(0).get(0));
        final ConnectionTree connectionTree = new ConnectionTree(atom, 0, 0);        
        // higher spheres
        for(sphere = 1; sphere <= maxSphere; sphere++){
            // get positions (sections separated by comma) of current sphere
            positionsInSphere = HOSECodeBuilder.splitHOSECodeSphereIntoPositions(sphereStrings.get(sphere)); 
            // for all positions
            for (final int positionIndex : positionsInSphere.keySet()) {
                // for each child elements (symbols) in position
                for (final String childElement : positionsInSphere.get(positionIndex)) {
                    // ignore children containing null value from previous nodes; previous node has no further (unvisited in BFS) connected atoms
                    if(childElement == null){
                        continue;
                    }
                    bondTypeString = "";
                    atomTypeString = "";
                    // add new node and set bond to parent node or set a ring closure
                    if (childElement.contains("&")) { // ring closure
                        if (childElement.length() == 2) {
                            bondTypeString = String.valueOf(childElement.charAt(0));
                        }
                        // the parent node/atom in previous sphere and its key we already have of a ring closure; 
                        // the bond information we already have too (see below)
                        ConnectionTreeNode parentNodeInPrevSphere = connectionTree.getNodesInSphere(sphere - 1).get(positionIndex);
                        // check whether the node in previous sphere is already involved in a ring closure which should be not valid
                        if(HOSECodeBuilder.isAtRingClosure(parentNodeInPrevSphere)){                            
                            continue;
                        }
                        
                        // what we still not can detect for sure is the correct second node/atom of a ring closure 
                        // -> still open
                        ConnectionTreeNode parentNodeInSphere = null; // null is just a dummy value and should be replaced by the correct ConnectionTreeNode object                        
                        
                        // check that the detected node in sphere is not null; could be removed after the implementation of detection of that node
                        if(parentNodeInSphere == null){                        
                            continue;
                        }
                        // after that both node detections, check if that ring closure was already set beforehand by the the reversed node order case
                        if(HOSECodeBuilder.nodesFormRingClosure(parentNodeInPrevSphere, parentNodeInSphere)){
                            continue;
                        }
                        // otherwise build a new bond and fill it with 
                        bond = SilentChemObjectBuilder.getInstance().newBond();
                        bond.setAtom(parentNodeInPrevSphere.getAtom(), 0);
                        bond.setAtom(parentNodeInSphere.getAtom(), 1);
                        bond.setOrder(HOSECodeBuilder.getBondOrderForSymbol(bondTypeString));
                        if (bondTypeString.equals("*")) {
                            bond.setIsAromatic(true);
                        } else {
                            bond.setIsAromatic(false);
                        }
                        // set parent nodes as parents to each other, that one can detect them as ring closure afterwards
                        parentNodeInPrevSphere.addParentNode(parentNodeInSphere, bond);
                        connectionTree.addNode(null, -1 * parentNodeInPrevSphere.getKey(), parentNodeInPrevSphere.getKey(), bond, sphere + 1, true);
                        parentNodeInSphere.addParentNode(parentNodeInPrevSphere, bond);
                        connectionTree.addNode(null, -1 * parentNodeInSphere.getKey(), parentNodeInSphere.getKey(), bond, sphere + 1, true);
                        
                    } else if(Utils.countElements(childElement) == 1){ // each position contains either ring closures (&) or one element (e.g. C, Br), plus the bond information
                        if(childElement.length() == 3){ // in case of bond type and an element with two letters, e.g. *Cl or =Br
                            bondTypeString = String.valueOf(childElement.charAt(0));
                            atomTypeString = String.valueOf(childElement.charAt(1));
                            atomTypeString += String.valueOf(childElement.charAt(2));
                        } else if(childElement.length() == 2){ // in case of bond type and an element with one letter or an element with two letters, e.g. Cl or =N
                            if(Character.isLetter(childElement.charAt(0))){
                                atomTypeString = String.valueOf(childElement.charAt(0));
                            } else {
                                bondTypeString = String.valueOf(childElement.charAt(0));
                            }
                            atomTypeString += String.valueOf(childElement.charAt(1));
                        } else if(childElement.length() == 1){ // in case of an element with only one letter
                            atomTypeString = String.valueOf(childElement.charAt(0));
                        } 
                        // there has to be some information (at least an element)
                        if(atomTypeString.isEmpty()){ 
                            throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": no atom information in child element");
                        }
                        // otherwise set a new bond
                        bond = SilentChemObjectBuilder.getInstance().newBond();
                        if(useBremserElementNotation){
                            atomTypeString = HOSECodeBuilder.toElementSymbol(atomTypeString);
                        }
                        atom = new Atom(atomTypeString);
                        bond.setAtom(atom, 0);
                        bond.setAtom(connectionTree.getNodesInSphere(sphere - 1).get(positionIndex).getAtom(), 1);
                        bond.setOrder(HOSECodeBuilder.getBondOrderForSymbol(bondTypeString));
                        if(bondTypeString.equals("*")){
                            bond.setIsAromatic(true);
                        } else {
                            bond.setIsAromatic(false);
                        }
                        // create a new node with the new build bond information to its parent node in the connection tree
                        connectionTree.addNode(atom, connectionTree.getNodesCount(),
                                connectionTree.getNodeKeysInSphere(sphere - 1).get(positionIndex), bond, sphere, false);
                    } else {
                        throw new CDKException(Thread.currentThread().getStackTrace()[1].getMethodName() + ": no valid components in child element");
                    }                  
                }
            }
        }
        
        HOSECodeBuilder.rankChildNodes(connectionTree);
        
        return connectionTree;        
    }
    
    public static IAtomContainer buildAtomContainer(final String HOSECode, final boolean useBremserElementNotation) throws CDKException{
        return HOSECodeBuilder.buildAtomContainer(HOSECodeBuilder.buildConnectionTree(HOSECode, useBremserElementNotation));
    }
}