package mairaDatabase.refseq.step1_clustering;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.sql.SQLException;
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

	private List<File> faaFiles;
	private int faaFilePointer = 0;

	public void runClustering(String rank, String srcPath, String aliFolder, File proteinFolder, TaxTree taxTree,
			SQLMappingDatabase mappingDatabase, int cores, int memory, int markerIdentity, int genusIdentity,
			File tmpFile, String diamondBin) {

		try {

			markerClusterOutputFolder = new File(srcPath + File.separator + rank + "_marker_proteins_clustered");
			markerClusterOutputFolder.mkdir();
			genusFolder = new File(srcPath + File.separator + rank + "_dbs");
			genusFolder.mkdir();
			faaFiles = new ArrayList<>(Arrays.asList(
					proteinFolder.listFiles((dir, name) -> name.endsWith(".faa") && !name.endsWith("_new.faa"))));
			Collections.sort(faaFiles, Comparator.comparingLong(File::length).reversed());
			long totalFileLength = 0L;
			for (File f : faaFiles)
				totalFileLength += f.length();
			int n = faaFiles.size();

			System.out.println(
					">Assessing new proteins for " + faaFiles.size() + " protein " + ((n == 1) ? "set" : "sets"));
			List<Runnable> collectNewProteinsThreads = new ArrayList<>();
			faaFilePointer = 0;
			for (int i = 0; i < cores; i++)
				collectNewProteinsThreads.add(new CollectNewProteinsThread(tmpFile, aliFolder));
			rL.runThreads(cores, collectNewProteinsThreads, totalFileLength);

			System.out.println(">Running DIAMOND on " + faaFiles.size() + " protein " + ((n == 1) ? "set" : "sets"));
			List<Runnable> alignProteinsThreads = new ArrayList<>();
			faaFilePointer = 0;
			for (int i = 0; i < cores; i++)
				alignProteinsThreads.add(new AlignProteinsThread(tmpFile, aliFolder, cores, memory, markerIdentity,
						mappingDatabase, diamondBin));
			rL.runThreads(1, alignProteinsThreads, totalFileLength);

			System.out.println(
					">Clustering for marker db " + faaFiles.size() + " protein " + ((n == 1) ? "set" : "sets"));
			List<Runnable> clusterProteinsForMarkerDbThread = new ArrayList<>();
			faaFilePointer = 0;
			for (int i = 0; i < cores; i++)
				clusterProteinsForMarkerDbThread.add(new ClusterProteinsThread(tmpFile, markerClusterOutputFolder, null,
						aliFolder, markerIdentity, mappingDatabase, ClusteringMode.MARKER_DB));
			rL.runThreads(cores, clusterProteinsForMarkerDbThread, totalFileLength);

			if (rank.equals("genus")) {
				System.out.println(
						">Clustering for genus db " + faaFiles.size() + " protein " + ((n == 1) ? "set" : "sets"));
				genusDominationFile = new File(srcPath + File.separator + "acc2dominator.tab");
				genusDominationFile.delete();
				List<Runnable> clusterProteinsForGenusDbThread = new ArrayList<>();
				try (BufferedWriter dominationWriter = new BufferedWriter(new FileWriter(genusDominationFile))) {
					faaFilePointer = 0;
					for (int i = 0; i < cores; i++)
						clusterProteinsForGenusDbThread.add(new ClusterProteinsThread(tmpFile, null, dominationWriter,
								aliFolder, genusIdentity, mappingDatabase, ClusteringMode.GENUS_DB));
					rL.runThreads(cores, clusterProteinsForGenusDbThread, totalFileLength);
				}
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}

		FileUtils.deleteDirectory(proteinFolder.getAbsolutePath());

	}

	private synchronized File nextFaaFile() {
		if (faaFilePointer < faaFiles.size())
			return faaFiles.get(faaFilePointer++);
		return null;
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

		private File tmpFile;
		private String aliFolder;

		public CollectNewProteinsThread(File tmpFile, String aliFolder) {
			this.tmpFile = tmpFile;
			this.aliFolder = aliFolder;
		}

		@Override
		public void run() {
			File faaFile;
			while ((faaFile = nextFaaFile()) != null) {
				try {
					String genus = Formatter.removeNonAlphanumerics(faaFile.getName().replaceAll("\\.faa", ""));
					SQLAlignmentDatabase alignmentDatabase = createAlignmentDatabase(aliFolder, genus, tmpFile);
					File newFile = new File(faaFile.getAbsolutePath().replaceAll("\\.faa", "_new.faa"));
					newFile.delete();
					newFile.createNewFile();
					try (BufferedWriter writer = new BufferedWriter(new FileWriter(newFile))) {
						List<FastaEntry> tokens = FastaReader.read(faaFile);
						for (FastaEntry token : tokens) {
							String acc = token.getName();
							if (!alignmentDatabase.containsAcc(acc, genus + "_clusterTable"))
								writer.write(">" + acc + "\n" + token.getSequence() + "\n");
						}
					}
					alignmentDatabase.close();
				} catch (Exception e) {
					e.printStackTrace();
				}
				rL.reportProgress(faaFile.length());
			}

			rL.countDown();
		}

	}

	private class AlignProteinsThread implements Runnable {

		private File tmpFile;
		private String aliFolder;
		private int cores, memory, identity;
		private SQLMappingDatabase mappingDatabase;
		private String diamondBin;

		public AlignProteinsThread(File tmpFile, String aliFolder, int cores, int memory, int identity,
				SQLMappingDatabase mappingDatabase, String diamondBin) {
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

			File faaFile;
			while ((faaFile = nextFaaFile()) != null) {

				try {
					String genus = Formatter.removeNonAlphanumerics(faaFile.getName().replaceAll("\\.faa", ""));
					SQLAlignmentDatabase sqlAliDatabase = createAlignmentDatabase(aliFolder, genus, tmpFile);
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
					sqlAliDatabase.close();
				} catch (Exception e) {
					e.printStackTrace();
				}
				rL.reportProgress(faaFile.length());
			}

			mappingDatabase.close();
			rL.countDown();

		}
	}

	private class ClusterProteinsThread implements Runnable {

		private BufferedWriter dominationWriter;
		private String aliFolder, outFolder;
		private File tmpFile;
		private int identity;
		private SQLMappingDatabase mappingDatabase;
		private ClusteringMode mode;

		public ClusterProteinsThread(File tmpFile, File outFolder, BufferedWriter dominationWriter, String aliFolder,
				int identity, SQLMappingDatabase mappingDatabase, ClusteringMode mode) {
			this.outFolder = outFolder != null ? outFolder.getAbsolutePath() : null;
			this.tmpFile = tmpFile;
			this.dominationWriter = dominationWriter;
			this.identity = identity;
			this.aliFolder = aliFolder;
			this.mappingDatabase = createMappingDatabase(mappingDatabase);
			this.mode = mode;
		}

		@Override
		public void run() {
			File faaFile;
			while ((faaFile = nextFaaFile()) != null) {
				try {
					String genus = Formatter.removeNonAlphanumerics(faaFile.getName().replaceAll("\\.faa", ""));
					if (mode == ClusteringMode.GENUS_DB) {
						File genusOutFolder = new File(genusFolder + File.separator + getGenus(faaFile));
						genusOutFolder.mkdir();
						outFolder = genusOutFolder.getAbsolutePath();
					}
					SQLAlignmentDatabase aliDatabase = createAlignmentDatabase(aliFolder, genus, tmpFile);
					String fileName = mode == ClusteringMode.MARKER_DB
							? faaFile.getName().replaceAll("\\.faa", "_clustered.faa")
							: faaFile.getName();
					File proteinOutFile = new File(outFolder + File.separator + fileName);
					proteinOutFile.delete();
					new Clustering().run(genus, aliDatabase, faaFile, proteinOutFile, dominationWriter, identity, mode);
					mappingDatabase.close();
					aliDatabase.close();
				} catch (Exception e) {
					e.printStackTrace();
				}
				rL.reportProgress(faaFile.length());
			}
			rL.countDown();
		}

	}

	private synchronized SQLMappingDatabase createMappingDatabase(SQLMappingDatabase mappingDatabase) {
		return new SQLMappingDatabase(mappingDatabase);
	}

	private synchronized SQLAlignmentDatabase createAlignmentDatabase(String aliFolder, String genus, File tmpFile)
			throws ClassNotFoundException, SQLException {
		return new SQLAlignmentDatabase(aliFolder, genus, tmpFile);
	}

}
