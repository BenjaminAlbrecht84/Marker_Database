/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mairaDatabase.utils.taxTree;

import java.io.*;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author bagci
 */
public class TaxDump {

    public TaxTree parse(String src) {

        long time = System.currentTimeMillis();
        System.out.println(">Parsing taxonomic tree");

        String namesPath = src + File.separator + "taxdump" + File.separator + "names.dmp";
        String nodesPath = src + File.separator + "taxdump" + File.separator + "nodes.dmp";
        HashMap<Integer, String> id2Name = parseNamesFile(namesPath);
        TaxTree tree = new TaxTree();
        try {
            BufferedReader br = new BufferedReader(new FileReader(new File(nodesPath)));
            String line = br.readLine();
            while ((line = br.readLine()) != null) {
                String[] lineSplit = line.split("\t");
                int id = Integer.parseInt(lineSplit[0]);
                int parentId = Integer.parseInt(lineSplit[2]);
                String rank = lineSplit[4];
                TaxNode n = null;
                if (tree.hasNode(id)) {
                    n = tree.getNode(id);
                    n.setRank(rank);
                    if (tree.hasNode(parentId)) {
                        n.setParent(tree.getNode(parentId));
                        tree.getNode(parentId).addChild(n);
                    } else {
                        TaxNode parent = new TaxNode(parentId);
                        parent.addChild(n);
                        parent.setName(id2Name.get(parentId));
                        n.setParent(parent);
                        tree.addNode(parentId, parent);
                    }
                } else {
                    n = new TaxNode(id);
                    n.setRank(rank);
                    if (tree.hasNode(parentId)) {
                        n.setParent(tree.getNode(parentId));
                        tree.getNode(parentId).addChild(n);
                    } else {
                        TaxNode parent = new TaxNode(parentId);
                        parent.addChild(n);
                        parent.setName(id2Name.get(parentId));
                        tree.addNode(parentId, parent);
                        n.setParent(parent);
                    }
                    n.setName(id2Name.get(id));
                    tree.addNode(id, n);
                }
            }
            tree.setRoot(tree.getNode(1));
            br.close();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(TaxDump.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(TaxDump.class.getName()).log(Level.SEVERE, null, ex);
        }
        TaxNode unknown = new TaxNode(-1);
        unknown.setParent(tree.getNode(1));
        tree.getNode(1).addChild(unknown);
        tree.addNode(-1, unknown);

        System.out.println("Runtime: " + ((System.currentTimeMillis() - time) / 1000) + "s");

        return tree;
    }

    private HashMap<Integer, String> parseNamesFile(String namesPath) {
        HashMap<Integer, String> id2Name = new HashMap<>();
        try {

            BufferedReader br = new BufferedReader(new FileReader(new File(namesPath)));
            String line;
            while ((line = br.readLine()) != null) {
                String[] lineSplit = line.split("\t");
                // System.out.println(Arrays.toString(lineSplit));
                if (lineSplit[6].equals("scientific name")) {
                    id2Name.put(Integer.parseInt(lineSplit[0]), lineSplit[2]);
                }
            }
            br.close();
            return id2Name;

        } catch (FileNotFoundException ex) {
            Logger.getLogger(TaxDump.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(TaxDump.class.getName()).log(Level.SEVERE, null, ex);
        }
        return null;
    }

}
