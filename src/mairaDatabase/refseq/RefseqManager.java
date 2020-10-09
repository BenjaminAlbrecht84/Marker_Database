package mairaDatabase.refseq;

import java.io.File;
import java.time.LocalDate;

import mairaDatabase.refseq.step0_downloading.NCBIDownloader;
import mairaDatabase.refseq.step0_downloading.ProteinDownloadManager;
import mairaDatabase.refseq.step1_clustering.ClusterManager;
import mairaDatabase.refseq.step1_clustering.MarkerManager;
import mairaDatabase.refseq.step2_filtering.FilterManager;
import mairaDatabase.refseq.utils.Cleaner;
import mairaDatabase.refseq.utils.NewickTaxTreeWriter;
import mairaDatabase.utils.SQLMairaDatabase;
import mairaDatabase.utils.SQLMappingDatabase;
import mairaDatabase.utils.taxTree.TaxDump;
import mairaDatabase.utils.taxTree.TaxTree;

public class RefseqManager {

	public final static int CLUSTER_MARKER_ID = 90, CLUSTER_GENUS_ID = 99, MARKER_ID = 80;
	public final static int MAX_PROTEINS_PER_GCF = 100;
	public final static int MIN_LENGTH = 100;

	public void run(File src, File tmp, String aliDir, int cores, int memory, String[] genera) {

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
		proteinDownloadManager.run(srcPath, mappingDatabase, taxTree, cores, genera);

		String rank = "genus";
		ClusterManager clusterManager = new ClusterManager();
		clusterManager.runClustering(rank, srcPath, aliDir, proteinDownloadManager.getProteinFolder(), taxTree,
				mappingDatabase, cores, memory, CLUSTER_MARKER_ID, CLUSTER_GENUS_ID, tmp);
		MarkerManager markerManager = new MarkerManager();
		markerManager.runMarker(rank, srcPath, aliDir, clusterManager.getMarkerClusterOutputFolder(), taxTree,
				mappingDatabase, cores, memory, MARKER_ID, tmp);
		FilterManager filterManager = new FilterManager();
		filterManager.run(rank, srcPath, aliDir, tmpDir, markerManager.getMarkerOutputFolder(), taxTree,
				mappingDatabase, MAX_PROTEINS_PER_GCF, CLUSTER_MARKER_ID, cores);

		File mairaDb = new File(srcPath + File.separator + "maira.db");
		mairaDb.delete();
		SQLMairaDatabase mairaDatabase = new SQLMairaDatabase(mairaDb.getAbsolutePath(), tmpDir);
		mairaDatabase.createAcc2taxidTable(proteinDownloadManager.getMappingFile());
		mairaDatabase.createFactorsTable(filterManager.getWeightFile(), MAX_PROTEINS_PER_GCF);
		mairaDatabase.createProtCountsTable(proteinDownloadManager.getProteinCountsFile());
		mairaDatabase.createAcc2DominatorsTable(clusterManager.getGenusDominationFile());

		NewickTaxTreeWriter.run(srcFolder, mairaDatabase);

//		Cleaner.apply(srcFolder, database);
		long runtime = (System.currentTimeMillis() - time) / 1000;
		System.out.println("Total runtime: " + runtime + "ms");

	}

}
