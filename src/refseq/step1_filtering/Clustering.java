package refseq.step1_filtering;

import refseq.RefseqManager;
import refseq.utils.aliHelper.SQLAlignmentDatabase;
import utils.FastaReader;
import utils.SparseString;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class Clustering {

    private int ID_THRESHOLD, COV_THRESHOLD;
    private final static int MIN_PROTEINS_CLUSTER = 1000;

    public void run(String genus, SQLAlignmentDatabase alignmentDatabase, File faaFile, File outFile, int MIN_COV, int MIN_ID) {

        String table = genus + "_clusterTable";
        COV_THRESHOLD = MIN_COV;
        ID_THRESHOLD = MIN_ID;
        ArrayList<Object[]> genusProteins = FastaReader.read(faaFile);
        long time = System.currentTimeMillis();

        ArrayList<ClusterNode> clusterNodes = new ArrayList<>(genusProteins.size());
        for (Object[] protein : genusProteins) {
            SparseString acc = (SparseString) protein[0];
            int count = alignmentDatabase.getAlignmentCount(acc.toString(), table);
            clusterNodes.add(new ClusterNode(acc, count));
        }
        HashMap<String, ClusterNode> acc2node = new HashMap<>(clusterNodes.size());
        clusterNodes.stream().forEach(v -> acc2node.put(v.getAcc(), v));
        Collections.sort(clusterNodes);

        int selectedNodes = clusterNodes.size();
        for (ClusterNode v : clusterNodes) {
            if (!v.isDominated) {
                ArrayList<SQLAlignmentDatabase.AlignmentInfo> alis = alignmentDatabase.getAlignments(v.getAcc(), table);
                for (SQLAlignmentDatabase.AlignmentInfo ali : alis) {
                    double qLen = ali.getQlen(), slen = ali.getSlen();
                    ClusterNode w = acc2node.get(ali.getRef());
                    if (w == null || w.isDominated || qLen < RefseqManager.MIN_LENGTH || slen < RefseqManager.MIN_LENGTH)
                        continue;
                    if (!ali.getQuery().equals(ali.getRef()) && ali.getIdentity() > ID_THRESHOLD && ali.getQueryCoverage() > COV_THRESHOLD) {
                        w.setDominated(true);
                        selectedNodes--;
                    }
                }
            }
            if (selectedNodes < MIN_PROTEINS_CLUSTER)
                break;
        }

        long written = 0;
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
            try {
                for (Object[] protein : genusProteins) {
                    String acc = protein[0].toString();
                    String seq = protein[1].toString();
                    if (!acc2node.get(acc).isDominated) {
                        writer.write(">" + acc + "\n" + seq + "\n");
                        written++;
                    }
                }
            } finally {
                writer.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        long runtime = (System.currentTimeMillis() - time) / 1000;
        System.err.println(genus + ": " + written + " proteins reported (" + runtime + "s)");

    }

    public class ClusterNode implements Comparable<ClusterNode> {

        private int outDegree;
        private SparseString acc;
        private boolean isDominated = false;

        public ClusterNode(SparseString acc, int outDegree) {
            this.acc = acc;
            this.outDegree = outDegree;
        }

        @Override
        public int compareTo(ClusterNode v) {
            return Integer.compare(v.outDegree, outDegree);
        }

        public String getAcc() {
            return acc.toString();
        }

        public void setDominated(boolean dominated) {
            isDominated = dominated;
        }

    }

}
