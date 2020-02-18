package refseq.step2_selecting;

import refseq.RefseqManager;
import refseq.utils.aliHelper.SQLAlignmentDatabase;
import utils.FastaReader;
import utils.SQLMappingDatabase;
import utils.SparseString;
import utils.Statistics;
import utils.taxTree.TaxTree;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;

public class Filtering {

    private int ID_THRESHOLD, COV_THRESHOLD;
    private SQLMappingDatabase mappingDatabase;
    private SQLAlignmentDatabase alignmentDatabase;
    private String table;

    public void run(File faaFile, String genus, BufferedWriter factorWriter, BufferedWriter markerWriter, TaxTree taxTree, SQLMappingDatabase mappingDatabase, SQLAlignmentDatabase alignmentDatabase, int NUM_OF_PROTEINS, int MIN_ID) {

        this.table = genus + "_clusterTable";
        this.COV_THRESHOLD = MIN_ID;
        this.ID_THRESHOLD = MIN_ID;
        this.mappingDatabase = mappingDatabase;
        this.alignmentDatabase = alignmentDatabase;
        long time = System.currentTimeMillis();

        ArrayList<Object[]> markerProteins = FastaReader.read(faaFile);
        ArrayList<MarkerNode> markerNodes = new ArrayList<>(markerProteins.size());
        HashMap<String, MarkerNode> acc2node = new HashMap<>(markerProteins.size());
        for (Object[] protein : markerProteins) {
            SparseString acc = (SparseString) protein[0];
            int count = alignmentDatabase.getAlignmentCount(acc.toString(), table);
            MarkerNode v = new MarkerNode(acc, count);
            markerNodes.add(v);
            acc2node.put(acc.toString(), v);
        }
        Collections.sort(markerNodes);

        int selectedNodes = markerNodes.size();
        HashMap<String, Integer> gcf2Counts = new HashMap<>();
        try {
            for (MarkerNode v : markerNodes) {
                HashMap<String, Integer> localCounts = new HashMap<>();
                for (String gcf : getCoveredGenomes(v, acc2node)) {
                    gcf2Counts.putIfAbsent(gcf, 0);
                    localCounts.putIfAbsent(gcf, 0);
                    localCounts.put(gcf, localCounts.get(gcf) + 1);
                    if (gcf2Counts.get(gcf) < NUM_OF_PROTEINS)
                        v.setSelected(true);
                }
                if (v.isSelected) {
                    for (String gcf : localCounts.keySet()) {
                        gcf2Counts.putIfAbsent(gcf, 0);
                        gcf2Counts.put(gcf, gcf2Counts.get(gcf) + localCounts.get(gcf));
                    }
                } else
                    selectedNodes--;
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        // writing out marker proteins
        for (Object[] protein : markerProteins) {
            String acc = protein[0].toString();
            if (acc2node.get(acc).isSelected) {
                try {
                    markerWriter.write(">" + acc + "\n" + protein[1].toString() + "\n");
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

        // writing out marker protein weights
        for (MarkerNode v : acc2node.values()) {
            if (v.isSelected) {
                ArrayList<Double> gcfFactors = new ArrayList<>();
                for (String gcf : getCoveredGenomes(v, acc2node))
                    gcfFactors.add(1. / (double) gcf2Counts.get(gcf));
                double mean = Statistics.getMean(gcfFactors);
                try {
                    factorWriter.write(v.getAcc() + "\t" + mean);
                    factorWriter.newLine();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }

        long runtime = (System.currentTimeMillis() - time) / 1000;
        System.err.println(genus + ": " + selectedNodes + "/" + markerProteins.size() + " marker proteins selected (" + runtime + "s)");

    }

    private HashSet<String> getCoveredGenomes(MarkerNode v, HashMap<String, MarkerNode> acc2node) {
        ArrayList<SQLAlignmentDatabase.AlignmentInfo> alis = alignmentDatabase.getAlignments(v.getAcc(), table);
        HashSet<String> coveredGenomes = new HashSet<>();
        for (SQLAlignmentDatabase.AlignmentInfo ali : alis) {
            double qLen = ali.getQlen(), slen = ali.getSlen();
            MarkerNode w = acc2node.get(ali.getRef());
            if (w == null || qLen < RefseqManager.MIN_LENGTH || slen < RefseqManager.MIN_LENGTH)
                continue;
            if (ali.getIdentity() > ID_THRESHOLD && ali.getQueryCoverage() > COV_THRESHOLD) {
                for (String gcf : mappingDatabase.getGCFByAcc(ali.getRef()))
                    coveredGenomes.add(gcf);
            }
        }
        return coveredGenomes;
    }

    public class MarkerNode implements Comparable<MarkerNode> {

        private int outDegree;
        private SparseString acc;
        private boolean isSelected = false;

        public MarkerNode(SparseString acc, int outDegree) {
            this.acc = acc;
            this.outDegree = outDegree;
        }

        @Override
        public int compareTo(MarkerNode v) {
            return Integer.compare(v.outDegree, outDegree);
        }

        public String getAcc() {
            return acc.toString();
        }

        public void setSelected(boolean selected) {
            isSelected = selected;
        }

    }

}
