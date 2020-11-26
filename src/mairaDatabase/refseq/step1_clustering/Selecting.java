package mairaDatabase.refseq.step1_clustering;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;

import mairaDatabase.refseq.RefseqManager;
import mairaDatabase.refseq.utils.Formatter;
import mairaDatabase.refseq.utils.aliHelper.SQLAlignmentDatabase;
import mairaDatabase.utils.FastaReader;
import mairaDatabase.utils.SQLMappingDatabase;
import mairaDatabase.utils.FastaReader.FastaEntry;
import mairaDatabase.utils.taxTree.TaxNode;
import mairaDatabase.utils.taxTree.TaxTree;

public class Selecting {

    private int ID_THRESHOLD, COV_THRESHOLD;
    private TaxTree taxTree;
    private final static int MIN_PROTEINS_SELECT = 100;

    public void run(String genus, TaxTree taxTree, SQLMappingDatabase mappingDatabase, SQLAlignmentDatabase alignmentDatabase, File faaFile, File outFile, int MIN_COV, int MIN_ID) {

        this.COV_THRESHOLD = MIN_COV;
        this.ID_THRESHOLD = MIN_ID;
        this.taxTree = taxTree;
        String table = genus + "_markerTable";
        ArrayList<FastaEntry> clusteringProteins = FastaReader.read(faaFile);
        long time = System.currentTimeMillis();

        int selectedNodes = clusteringProteins.size();
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
            try {
                for (FastaEntry protein : clusteringProteins) {
                    String acc = protein.getName();
                    ArrayList<SQLAlignmentDatabase.AlignmentInfo> alis = alignmentDatabase.getAlignments(acc, table);
                    boolean isUnique = true;
                    for (SQLAlignmentDatabase.AlignmentInfo ali : alis) {
                        double qLen = ali.getQlen(), sLen = ali.getSlen();
                        if (qLen < RefseqManager.MIN_LENGTH || sLen < RefseqManager.MIN_LENGTH)
                            continue;
                        String refGenus = getRank(mappingDatabase.getTaxIdByAcc(ali.getRef()), "genus");
                        if (refGenus != null && !refGenus.equals(genus) && ali.getIdentity() > ID_THRESHOLD && ali.getQueryCoverage() > COV_THRESHOLD) {
                            isUnique = false;
                            selectedNodes--;
                            break;
                        }
                    }
                    if (isUnique)
                        writer.write(">" + acc + "\n" + protein.getSequence() + "\n");
                    if (selectedNodes < MIN_PROTEINS_SELECT)
                        break;
                }
            } finally {
                writer.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        long runtime = (System.currentTimeMillis() - time) / 1000;
        System.err.println(genus + ": " + selectedNodes + " proteins reported (" + runtime + "s)");

    }

    private String getRank(List<Integer> taxids, String rank) {
        String name = null;
        for (int taxid : taxids) {
            TaxNode w = taxTree.getNode(taxid);
            while (w != null) {
                if (w.getRank().equals(rank)) {
                    if (name == null)
                        name = w.getName();
                    else if (!name.equals(w.getName()))
                        return null;
                    break;
                }
                w = w.getParent();
            }
        }
        if (name != null)
            return Formatter.removeNonAlphanumerics(name);
        return null;
    }

}
