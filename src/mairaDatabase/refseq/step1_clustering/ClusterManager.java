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
	private File genusDominationFile, speciesDisjoinFile;
	private File genusFolder;
	private File markerClusterOutputFolder;

	private List<File> faaFiles;
	private int faaFilePointer = 0;

	public void runClustering(String rank, String srcPath, String aliFolder, File proteinFolder, TaxTree taxTree,
			SQLMappingDatabase mappingDatabase, int cores, double blockSize, int minMarkerIdentity,
			int minGenusIdentity, int minDisjoinIdentity, int minDisjoinCoverage, File tmpFile, String diamondBin) {

		try {

			markerClusterOutputFolder = new File(srcPath + File.separator + rank + "_marker_proteins_clustered");
			markerClusterOutputFolder.mkdir();
			genusFolder = new File(srcPath + File.separator + rank + "_dbs");
			genusFolder.mkdir();
			faaFiles = new ArrayList<>(Arrays.asList(proteinFolder.listFiles(
					(dir, name) -> name.endsWith(".faa") && !name.endsWith("_new.faa") && !name.endsWith("_old.faa"))));
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
				alignProteinsThreads.add(new AlignProteinsThread(tmpFile, aliFolder, cores, blockSize,
						minMarkerIdentity, mappingDatabase, diamondBin));
			rL.runThreads(1, alignProteinsThreads, totalFileLength);
			
			if (rank.equals("genus")) {

				System.out.println(">Computing cross species genome disjoints for " + faaFiles.size() + " "
						+ ((n == 1) ? "genus" : "genera"));
				speciesDisjoinFile = new File(srcPath + File.separator + "species2disjoin.tab");
				speciesDisjoinFile.delete();
				try (BufferedWriter speciesDisjoinWriter = new BufferedWriter(new FileWriter(speciesDisjoinFile))) {
					faaFilePointer = 0;
					new SpeciesDisjoinComparatorThread(speciesDisjoinWriter, aliFolder, tmpFile, minDisjoinIdentity,
							minDisjoinCoverage, mappingDatabase, cores).run();
				}

			}

			System.out.println(
					">Clustering for marker db " + faaFiles.size() + " protein " + ((n == 1) ? "set" : "sets"));
			List<Runnable> clusterProteinsForMarkerDbThreads = new ArrayList<>();
			faaFilePointer = 0;
			for (int i = 0; i < cores; i++)
				clusterProteinsForMarkerDbThreads.add(new ClusterProteinsThread(tmpFile, markerClusterOutputFolder,
						null, aliFolder, minMarkerIdentity, mappingDatabase, taxTree, ClusteringMode.MARKER_DB));
			rL.runThreads(cores, clusterProteinsForMarkerDbThreads, totalFileLength);

			if (rank.equals("genus")) {

				System.out.println(
						">Clustering for genus db " + faaFiles.size() + " protein " + ((n == 1) ? "set" : "sets"));
				genusDominationFile = new File(srcPath + File.separator + "acc2dominator.tab");
				genusDominationFile.delete();
				List<Runnable> clusterProteinsForGenusDbThreads = new ArrayList<>();
				try (BufferedWriter dominationWriter = new BufferedWriter(new FileWriter(genusDominationFile))) {
					faaFilePointer = 0;
					for (int i = 0; i < cores; i++)
						clusterProteinsForGenusDbThreads.add(new ClusterProteinsThread(tmpFile, null, dominationWriter,
								aliFolder, minGenusIdentity, mappingDatabase, taxTree, ClusteringMode.GENUS_DB));
					rL.runThreads(cores, clusterProteinsForGenusDbThreads, totalFileLength);
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

	public File getSpeciesDisjoinFile() {
		return speciesDisjoinFile;
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
					String genus = Formatter
							.removeNonAlphanumerics(faaFile.getName().split("\\-")[1].replaceAll("\\.faa", ""));
					SQLAlignmentDatabase alignmentDatabase = createAlignmentDatabase(aliFolder, genus, tmpFile);
					try {
						File oldFile = new File(faaFile.getAbsolutePath().replaceAll("\\.faa", "_old.faa"));
						oldFile.delete();
						oldFile.createNewFile();
						File newFile = new File(faaFile.getAbsolutePath().replaceAll("\\.faa", "_new.faa"));
						newFile.delete();
						newFile.createNewFile();
						try (BufferedWriter newWriter = new BufferedWriter(new FileWriter(newFile));
								BufferedWriter oldWriter = new BufferedWriter(new FileWriter(oldFile));) {
							List<FastaEntry> tokens = FastaReader.read(faaFile);
							for (FastaEntry token : tokens) {
								String acc = token.getName();
								if (!alignmentDatabase.containsAcc(acc))
									newWriter.write(">" + acc + "\n" + token.getSequence() + "\n");
								else
									oldWriter.write(">" + acc + "\n" + token.getSequence() + "\n");
							}
						}
					} finally {
						alignmentDatabase.close();
					}
				} catch (Exception e) {
					e.printStackTrace();
				}
				rL.reportProgress(faaFile.length());
				try {
					Thread.sleep(1000);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}

			rL.countDown();
		}

	}

	private class AlignProteinsThread implements Runnable {

		private File tmpFile;
		private String aliFolder;
		private int cores, identity;
		private double blockSize;
		private SQLMappingDatabase mappingDatabase;
		private String diamondBin;

		public AlignProteinsThread(File tmpFile, String aliFolder, int cores, double blockSize, int identity,
				SQLMappingDatabase mappingDatabase, String diamondBin) {
			this.tmpFile = tmpFile;
			this.aliFolder = aliFolder;
			this.cores = cores;
			this.blockSize = blockSize;
			this.identity = identity;
			this.mappingDatabase = createMappingDatabase(mappingDatabase);
			this.diamondBin = diamondBin;
		}

		@Override
		public void run() {

			File faaFile;
			while ((faaFile = nextFaaFile()) != null) {

				try {
					String genus = Formatter
							.removeNonAlphanumerics(faaFile.getName().split("\\-")[1].replaceAll("\\.faa", ""));
					SQLAlignmentDatabase alignmentDatabase = createAlignmentDatabase(aliFolder, genus, tmpFile);
					File newFile = new File(faaFile.getAbsolutePath().replaceAll("\\.faa", "_new.faa"));
					newFile.createNewFile();
					File oldFile = new File(faaFile.getAbsolutePath().replaceAll("\\.faa", "_old.faa"));
					oldFile.createNewFile();
					try {
						if (newFile.exists() && newFile.length() > 0) {
							File db = DiamondRunner.makedb(faaFile, cores, diamondBin);
							File tabFile1 = DiamondRunner.blastp(db, newFile, tmpFile, identity, blockSize, cores,
									diamondBin);
							alignmentDatabase.addAlignmentTable(genus + "_clusterTable", null, tabFile1, false);
							db.delete();
							tabFile1.delete();

							File dbNew = DiamondRunner.makedb(newFile, cores, diamondBin);
							File tabFile2 = DiamondRunner.blastp(dbNew, oldFile, tmpFile, identity, blockSize, cores,
									diamondBin);
							alignmentDatabase.addAlignmentTable(genus + "_clusterTable", null, tabFile2, false);
							dbNew.delete();
							tabFile2.delete();
						}
					} finally {
						newFile.delete();
						oldFile.delete();
						alignmentDatabase.close();
					}
				} catch (Exception e) {
					e.printStackTrace();
				}
				rL.reportProgress(faaFile.length());
				try {
					Thread.sleep(1000);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
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
		private TaxTree taxTree;

		public ClusterProteinsThread(File tmpFile, File outFolder, BufferedWriter dominationWriter, String aliFolder,
				int identity, SQLMappingDatabase mappingDatabase, TaxTree taxTree, ClusteringMode mode) {
			this.outFolder = outFolder != null ? outFolder.getAbsolutePath() : null;
			this.tmpFile = tmpFile;
			this.dominationWriter = dominationWriter;
			this.identity = identity;
			this.aliFolder = aliFolder;
			this.mappingDatabase = createMappingDatabase(mappingDatabase);
			this.taxTree = taxTree;
			this.mode = mode;
		}

		@Override
		public void run() {
			File faaFile;
			while ((faaFile = nextFaaFile()) != null) {
				try {
					String genus = Formatter
							.removeNonAlphanumerics(faaFile.getName().split("\\-")[1].replaceAll("\\.faa", ""));
					int genusId = Integer.parseInt(faaFile.getName().split("\\-")[0]);
					if (mode == ClusteringMode.GENUS_DB) {
						File genusOutFolder = new File(genusFolder + File.separator + getGenus(faaFile));
						genusOutFolder.mkdir();
						outFolder = genusOutFolder.getAbsolutePath();
					}
					SQLAlignmentDatabase alignmentDatabase = createAlignmentDatabase(aliFolder, genus, tmpFile);
					try {
						String proteinFileName = mode == ClusteringMode.MARKER_DB
								? faaFile.getName().replaceAll("\\.faa", "_clustered.faa")
								: faaFile.getName();
						File proteinOutFile = new File(outFolder + File.separator + proteinFileName);
						proteinOutFile.delete();
						new Clustering().run(genusId, genus, alignmentDatabase, mappingDatabase, taxTree, faaFile,
								proteinOutFile, dominationWriter, identity, mode);
					} finally {
						alignmentDatabase.close();
					}
				} catch (Exception e) {
					e.printStackTrace();
				}
				rL.reportProgress(faaFile.length());
				try {
					Thread.sleep(1000);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			mappingDatabase.close();
			rL.countDown();
		}

	}

	private class SpeciesDisjoinComparatorThread {

		private BufferedWriter speciesDisjoinWriter;
		private String aliFolder;
		private File tmpFile;
		private int minIdentity, minCoverage;
		private SQLMappingDatabase mappingDatabase;
		private int cores;

		public SpeciesDisjoinComparatorThread(BufferedWriter speciesDisjoinWriter, String aliFolder, File tmpFile,
				int identity, int coverage, SQLMappingDatabase mappingDatabase, int cores) {
			this.speciesDisjoinWriter = speciesDisjoinWriter;
			this.aliFolder = aliFolder;
			this.tmpFile = tmpFile;
			this.minIdentity = identity;
			this.minCoverage = coverage;
			this.mappingDatabase = createMappingDatabase(mappingDatabase);
			this.cores = cores;
		}

		public void run() {
			File faaFile;
			while ((faaFile = nextFaaFile()) != null) {
				try {
					String genus = Formatter
							.removeNonAlphanumerics(faaFile.getName().split("\\-")[1].replaceAll("\\.faa", ""));
					int genusId = Integer.parseInt(faaFile.getName().split("\\-")[0]);
					SQLAlignmentDatabase alignmentDatabase = createAlignmentDatabase(aliFolder, genus, tmpFile);
					try {
						new SpeciesDisjoinComparator2().run(faaFile, genus, genusId, alignmentDatabase, mappingDatabase,
								speciesDisjoinWriter, minIdentity, minCoverage, cores);
					} finally {
						alignmentDatabase.close();
					}
				} catch (Exception e) {
					e.printStackTrace();
				}
				rL.reportProgress(faaFile.length());
				try {
					Thread.sleep(1000);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			mappingDatabase.close();
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
