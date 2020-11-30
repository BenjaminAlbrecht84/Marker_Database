package mairaDatabase.refseq.step1_clustering;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import mairaDatabase.refseq.step1_clustering.Clustering.ClusteringMode;
import mairaDatabase.refseq.utils.DiamondRunner;
import mairaDatabase.refseq.utils.Formatter;
import mairaDatabase.refseq.utils.aliHelper.SQLAlignmentDatabase;
import mairaDatabase.utils.FastaReader;
import mairaDatabase.utils.FastaReader.FastaEntry;
import mairaDatabase.utils.FileUtils;
import mairaDatabase.utils.ResourceLoader;
import mairaDatabase.utils.SQLMappingDatabase;
import mairaDatabase.utils.taxTree.TaxTree;

public class ClusterManager {

	private ResourceLoader rL = new ResourceLoader();
	private File genusDominationFile;
	private File genusFolder;
	private File markerClusterOutputFolder;

	public void runClustering(String rank, String srcPath, String aliFolder, File proteinFolder, TaxTree taxTree,
			SQLMappingDatabase mappingDatabase, int cores, int memory, int markerIdentity, int genusIdentity,
			File tmpFile, String diamondBin) {

		markerClusterOutputFolder = new File(srcPath + File.separator + rank + "_marker_proteins_clustered");
		markerClusterOutputFolder.mkdir();
		genusFolder = new File(srcPath + File.separator + rank + "_dbs");
		genusFolder.mkdir();
		File[] initFaaFiles = proteinFolder
				.listFiles((dir, name) -> name.endsWith(".faa") && !name.endsWith("_new.faa"));
		ArrayList<File> faaFiles = new ArrayList<>(Arrays.asList(initFaaFiles));
		Collections.sort(faaFiles, Comparator.comparingLong(File::length).reversed());
		long totalFileLength = 0L;
		for (File f : faaFiles)
			totalFileLength += f.length();
		initFaaFiles = faaFiles.toArray(new File[faaFiles.size()]);
		int n = initFaaFiles.length;

		System.out.println(
				">Assessing new proteins for " + initFaaFiles.length + " protein " + ((n == 1) ? "set" : "sets"));
		List<Runnable> collectNewProteinsThreads = new ArrayList<>();
		for (File faaFile : initFaaFiles)
			collectNewProteinsThreads.add(new CollectNewProteinsThread(faaFile, tmpFile, aliFolder));
		rL.runThreads(cores, collectNewProteinsThreads, totalFileLength);

		System.out.println(">Running DIAMOND on " + initFaaFiles.length + " protein " + ((n == 1) ? "set" : "sets"));
		List<Runnable> alignProteinsThreads = new ArrayList<>();
		for (File faaFile : initFaaFiles)
			alignProteinsThreads.add(new AlignProteinsThread(faaFile, tmpFile, aliFolder, cores, memory, markerIdentity,
					mappingDatabase, diamondBin));
		rL.runThreads(1, alignProteinsThreads, totalFileLength);

		System.out.println(
				">Clustering for marker db " + initFaaFiles.length + " protein " + ((n == 1) ? "set" : "sets"));
		List<Runnable> clusterProteinsForMarkerDbThread = new ArrayList<>();
		for (File faaFile : initFaaFiles)
			clusterProteinsForMarkerDbThread.add(new ClusterProteinsThread(faaFile, tmpFile, markerClusterOutputFolder,
					null, aliFolder, markerIdentity, mappingDatabase, ClusteringMode.MARKER_DB));
		rL.runThreads(cores, clusterProteinsForMarkerDbThread, totalFileLength);

		try {
			if (rank.equals("genus")) {

				genusDominationFile = new File(srcPath + File.separator + "acc2dominator.tab");
				genusDominationFile.delete();

				System.out.println(
						">Clustering for genus db " + initFaaFiles.length + " protein " + ((n == 1) ? "set" : "sets"));
				List<Runnable> clusterProteinsForGenusDbThread = new ArrayList<>();
				BufferedWriter dominationWriter = new BufferedWriter(new FileWriter(genusDominationFile));
				try {
					for (File faaFile : initFaaFiles) {
						File genusOutputFolder = new File(genusFolder + File.separator + getGenus(faaFile));
						genusOutputFolder.mkdir();
						clusterProteinsForGenusDbThread
								.add(new ClusterProteinsThread(faaFile, tmpFile, genusOutputFolder, dominationWriter,
										aliFolder, genusIdentity, mappingDatabase, ClusteringMode.GENUS_DB));
					}
					rL.runThreads(cores, clusterProteinsForGenusDbThread, totalFileLength);
				} finally {
					dominationWriter.close();
				}

			}
		} catch (Exception e) {
			e.printStackTrace();
		}

		FileUtils.deleteDirectory(proteinFolder.getAbsolutePath());

	}

	public File getGenusFolder() {
		return genusFolder;
	}

	public String getGenus(File faaFile) {
		return faaFile.getName().split("\\.")[0];
	}

	public File getGenusDominationFile() {
		return genusDominationFile;
	}

	public File getMarkerClusterOutputFolder() {
		return markerClusterOutputFolder;
	}

	private class CollectNewProteinsThread implements Runnable {

		private File faaFile, tmpFile;
		private String aliFolder;

		public CollectNewProteinsThread(File faaFile, File tmpFile, String aliFolder) {
			this.faaFile = faaFile;
			this.tmpFile = tmpFile;
			this.aliFolder = aliFolder;
		}

		@Override
		public void run() {
			try {
				String genus = Formatter.removeNonAlphanumerics(faaFile.getName().replaceAll("\\.faa", ""));
				SQLAlignmentDatabase alignmentDatabase = createAlignmentDatabase(aliFolder, genus, tmpFile);
				File newFile = new File(faaFile.getAbsolutePath().replaceAll("\\.faa", "_new.faa"));
				newFile.delete();
				newFile.createNewFile();
				BufferedWriter writer = new BufferedWriter(new FileWriter(newFile));
				List<FastaEntry> tokens = FastaReader.read(faaFile);
				for (FastaEntry token : tokens) {
					String acc = token.getName();
					if (!alignmentDatabase.containsAcc(acc, genus + "_clusterTable"))
						writer.write(">" + acc + "\n" + token.getSequence() + "\n");
				}
				writer.close();
				alignmentDatabase.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
			rL.reportProgress(faaFile.length());
			rL.countDown();
		}

	}

	private class AlignProteinsThread implements Runnable {

		private File faaFile, tmpFile;
		private String aliFolder;
		private int cores, memory, identity;
		private SQLMappingDatabase mappingDatabase;
		private String diamondBin;

		public AlignProteinsThread(File faaFile, File tmpFile, String aliFolder, int cores, int memory, int identity,
				SQLMappingDatabase mappingDatabase, String diamondBin) {
			this.faaFile = faaFile;
			this.tmpFile = tmpFile;
			this.aliFolder = aliFolder;
			this.cores = cores;
			this.memory = memory;
			this.identity = identity;
			this.mappingDatabase = createMappingDatabase(mappingDatabase);
			this.diamondBin = diamondBin;
		}

		@Override
		public void run() {

			String genus = Formatter.removeNonAlphanumerics(faaFile.getName().replaceAll("\\.faa", ""));
			SQLAlignmentDatabase sqlAliDatabase = createAlignmentDatabase(aliFolder, genus, tmpFile);

			try {

				File newFile = new File(faaFile.getAbsolutePath().replaceAll("\\.faa", "_new.faa"));
				newFile.createNewFile();
				if (newFile.exists() && newFile.length() > 0) {
					File dbFile1 = DiamondRunner.makedb(faaFile, cores, diamondBin);
					File tabFile1 = DiamondRunner.blastp(dbFile1, newFile, tmpFile, identity, memory, cores,
							diamondBin);
					sqlAliDatabase.addAlignmentTable(genus + "_clusterTable", null, tabFile1, false);
					dbFile1.delete();
					tabFile1.delete();

					File dbFile2 = DiamondRunner.makedb(newFile, cores, diamondBin);
					File tabFile2 = DiamondRunner.blastp(dbFile2, faaFile, tmpFile, identity, memory, cores,
							diamondBin);
					sqlAliDatabase.addAlignmentTable(genus + "_clusterTable", null, tabFile2, false);
					dbFile2.delete();
					tabFile2.delete();
				}
				newFile.delete();

			} catch (Exception e) {
				e.printStackTrace();
			}

			mappingDatabase.close();
			sqlAliDatabase.close();
			rL.reportProgress(faaFile.length());
			rL.countDown();

		}
	}

	private class ClusterProteinsThread implements Runnable {

		private BufferedWriter dominationWriter;
		private String aliFolder, outFolder, genus;
		private File faaFile, tmpFile;
		private int identity;
		private SQLMappingDatabase mappingDatabase;
		private ClusteringMode mode;

		public ClusterProteinsThread(File faaFile, File tmpFile, File outFolder, BufferedWriter dominationWriter,
				String aliFolder, int identity, SQLMappingDatabase mappingDatabase, ClusteringMode mode) {
			this.outFolder = outFolder.getAbsolutePath();
			this.faaFile = faaFile;
			this.tmpFile = tmpFile;
			this.dominationWriter = dominationWriter;
			this.identity = identity;
			this.aliFolder = aliFolder;
			this.mappingDatabase = createMappingDatabase(mappingDatabase);
			this.genus = Formatter.removeNonAlphanumerics(faaFile.getName().replaceAll("\\.faa", ""));
			this.mode = mode;
		}

		@Override
		public void run() {
			SQLAlignmentDatabase aliDatabase = createAlignmentDatabase(aliFolder, genus, tmpFile);
			String fileName = mode == ClusteringMode.MARKER_DB
					? faaFile.getName().replaceAll("\\.faa", "_clustered.faa")
					: faaFile.getName();
			File proteinOutFile = new File(outFolder + File.separator + fileName);
			proteinOutFile.delete();
			new Clustering().run(genus, aliDatabase, faaFile, proteinOutFile, dominationWriter, identity, mode);
			mappingDatabase.close();
			aliDatabase.close();
			rL.reportProgress(faaFile.length());
			rL.countDown();
		}

	}

	private synchronized SQLMappingDatabase createMappingDatabase(SQLMappingDatabase mappingDatabase) {
		return new SQLMappingDatabase(mappingDatabase);
	}

	private synchronized SQLAlignmentDatabase createAlignmentDatabase(String aliFolder, String genus, File tmpFile) {
		return new SQLAlignmentDatabase(aliFolder, genus, tmpFile);
	}

}
