/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package utils.taxTree;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

/**
 *
 * @author bagci
 */
public class TaxTree {

	private ArrayList<TaxNode> leaves;
	private TaxNode root;
	private HashMap<Integer, TaxNode> id2node;

	public TaxTree() {
		this.id2node = new HashMap<>();
	}

	public TaxNode getRoot() {
		return root;
	}

	public void addNode(int id, TaxNode n) {
		id2node.put(id, n);
	}

	public TaxNode getNode(int id) {
		return id2node.get(id);
	}

	public boolean hasNode(int id) {
		return id2node.containsKey(id);
	}

	public int getNumNodes() {
		return id2node.size();
	}

	public Collection<TaxNode> getNodes() {
		return id2node.values();
	}

	public ArrayList<TaxNode> getLeaves() {
		if (leaves == null) {
			leaves = new ArrayList<TaxNode>();
			cmpLeafSetRec(root);
		}
		return leaves;
	}

	private void cmpLeafSetRec(TaxNode v) {
		if (v.getChildren().isEmpty())
			leaves.add(v);
		else {
			for (TaxNode c : v.getChildren())
				cmpLeafSetRec(c);
		}
	}

	public void setRoot(TaxNode node) {
		this.root = node;
		root.setRank("root");
		root.setName("root");
	}

}
