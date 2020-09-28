package mairaDatabase.refseq;

import java.io.File;
import java.time.LocalDate;

import mairaDatabase.refseq.step0_downloading.NCBIDownloader;
import mairaDatabase.refseq.step0_downloading.ProteinDownloadManager;
import mairaDatabase.refseq.step1_clustering.ClusterManager;
import mairaDatabase.refseq.step1_clustering.MarkerManager;
import mairaDatabase.refseq.step2_filtering.FilterManager;
import mairaDatabase.refseq.utils.NewickTaxTreeWriter;
import mairaDatabase.utils.SQLMairaDatabase;
import mairaDatabase.utils.SQLMappingDatabase;
import mairaDatabase.utils.taxTree.TaxDump;
import mairaDatabase.utils.taxTree.TaxTree;

public class RefseqManager {

    public final static int CLUSTER_ID = 90, MARKER_ID = 80;
    public final static int MAX_PROTEINS_PER_GCF = 1000;
    public final static int MIN_LENGTH = 100;

    public void run(File src, File tmp, String aliDir, int cores, int memory, String[] genera) {

        long time = System.currentTimeMillis();

        File srcFolder = new File(src.getAbsolutePath() + File.separator + "maira_refseq_" + LocalDate.now());
        srcFolder.mkdir();
        String srcPath = srcFolder.getAbsolutePath();
        File database = new File(srcPath + File.separator + "refseq_mapping.db");
        File tmpDir = tmp == null ? src : tmp;
        SQLMappingDatabase mappingDatabase = new SQLMappingDatabase(database.getAbsolutePath(), tmpDir);

//        new NCBIDownloader().run(srcPath);
        TaxTree taxTree = new TaxDump().parse(srcPath);
        ProteinDownloadManager proteinDownloadManager = new ProteinDownloadManager();
        proteinDownloadManager.run(srcPath, mappingDatabase, taxTree, cores, genera);

        String rank = "genus";
        new ClusterManager().runClustering(rank, srcPath, aliDir, taxTree, mappingDatabase, cores, memory, CLUSTER_ID, tmp);
        new MarkerManager().runMarker(rank, srcPath, aliDir, taxTree, mappingDatabase, cores, memory, MARKER_ID, tmp);
        FilterManager filterManager = new FilterManager();
        filterManager.run(rank, srcPath, aliDir, tmpDir, taxTree, mappingDatabase, MAX_PROTEINS_PER_GCF, CLUSTER_ID, cores);
        
        File mairaDb = new File(srcPath + File.separator + "maira.db");
        SQLMairaDatabase mairaDatabase = new SQLMairaDatabase(mairaDb.getAbsolutePath(), tmpDir);
        mairaDatabase.createAcc2taxidTable(proteinDownloadManager.getMappingFile());
        mairaDatabase.createFactorsTable(filterManager.getWeightFile(), MAX_PROTEINS_PER_GCF);
        mairaDatabase.createProtCountsTable(proteinDownloadManager.getProteinCountsFile());
        
        NewickTaxTreeWriter.run(src);
        
        long runtime = (System.currentTimeMillis() - time) / 1000;
        System.out.println("Total runtime: " + runtime + "ms");

    }

}
