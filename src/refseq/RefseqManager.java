package refseq;

import refseq.step0_downloading.NCBIDownloader;
import refseq.step0_downloading.ProteinDownloadManager;
import refseq.step1_filtering.SuppressorManager;
import refseq.step2_selecting.FilterManager;
import utils.SQLMappingDatabase;
import utils.taxTree.TaxDump;
import utils.taxTree.TaxTree;

import java.io.File;
import java.time.LocalDate;

public class RefseqManager {

    public final static int CLUSTER_ID = 90, MARKER_ID = 80;
    public final static int NUM_OF_PROTEINS = 1000;
    public final static int MIN_LENGTH = 100;

    public void run(File src, File tmp, String aliDir, int cores, int memory) {

        long time = System.currentTimeMillis();

        File srcFolder = new File(src.getAbsolutePath() + File.separator + "maira_refseq_" + LocalDate.now());
        srcFolder.mkdir();
        String srcPath = srcFolder.getAbsolutePath();
        File database = new File(srcPath + File.separator + "refseq_mapping.db");
        File tmpDir = tmp == null ? src : tmp;
        SQLMappingDatabase mappingDatabase = new SQLMappingDatabase(database.getAbsolutePath(), tmpDir);

        new NCBIDownloader().run(srcPath);
        TaxTree taxTree = new TaxDump().parse(srcPath);
        ProteinDownloadManager proteinDownloadManager = new ProteinDownloadManager();
        proteinDownloadManager.run(srcPath, mappingDatabase, taxTree, cores);

        String rank = "genus";
        new SuppressorManager().runClustering(rank, srcPath, aliDir, taxTree, mappingDatabase, cores, memory, CLUSTER_ID, tmp);
        new SuppressorManager().runMarker(rank, srcPath, aliDir, taxTree, mappingDatabase, cores, memory, MARKER_ID, tmp);
        new FilterManager().run(rank, srcPath, aliDir, tmpDir, taxTree, mappingDatabase, NUM_OF_PROTEINS, CLUSTER_ID, cores);

        long runtime = (System.currentTimeMillis() - time) / 1000;
        System.out.println("Total runtime: " + runtime + "ms");

    }

}
