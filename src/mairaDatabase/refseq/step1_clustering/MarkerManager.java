package mairaDatabase.refseq.step1_clustering;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import mairaDatabase.refseq.utils.DiamondRunner;
import mairaDatabase.refseq.utils.Formatter;
import mairaDatabase.refseq.utils.aliHelper.SQLAlignmentDatabase;
import mairaDatabase.utils.FastaReader;
import mairaDatabase.utils.FastaReader.FastaEntry;
import mairaDatabase.utils.FileUtils;
import mairaDatabase.utils.ResourceLoader;
import mairaDatabase.utils.SQLMappingDatabase;
import mairaDatabase.utils.taxTree.TaxTree;

public class MarkerManager {

	private ResourceLoader rL = new ResourceLoader();
	private File markerOutputFolder;
	private List<File> faaFiles;
	private int faaFilePointer = 0;
	private Set<String> oldAccessions = new HashSet<>();
	private Set<String> newAccessions = new HashSet<>();

	public void runMarker(String rank, String srcPath, String aliFolder, File markerClusterFolder, TaxTree taxTree,
			SQLMappingDatabase mappingDatabase, int cores, double blockSize, int identity, File tmpFile, String diamondBin) {

		try {

			SQLAlignmentDatabase genusAliDatabase = createAlignmentDatabase(aliFolder, "Genus", tmpFile);
			markerOutputFolder = new File(srcPath + File.separator + rank + "_marker_proteins");
			markerOutputFolder.mkdir();
			faaFiles = new ArrayList<>(Arrays.asList(
					markerClusterFolder.listFiles((dir, name) -> name.endsWith(".faa") && !name.endsWith("_new.faa"))));
			Collections.sort(faaFiles, Comparator.comparingLong(File::length).reversed());
			long totalFileLength = 0L;
			for (File f : faaFiles)
				totalFileLength += f.length();

			System.out.println(">Assessing newly computed representatives for " + faaFiles.size() + " protein sets");
			rL.setTime();
			List<Runnable> collectNewProteinsThreads = new ArrayList<>();
			faaFilePointer = 0;
			for (int i = 0; i < cores; i++)
				collectNewProteinsThreads
						.add(new CollectNewProteinsThread(tmpFile, markerOutputFolder.getAbsolutePath(), aliFolder));
			rL.runThreads(cores, collectNewProteinsThreads, totalFileLength);

			System.out.println(">Collecting old proteins");
			File oldProteins = new File(markerOutputFolder + File.separator + "oldProteins.faa");
			oldProteins.createNewFile();
			oldProteins.deleteOnExit();
			for (File f : markerClusterFolder.listFiles((dir, name) -> name.endsWith("_old.faa")))
				appendToFile(f, oldProteins, oldAccessions);
			oldAccessions = null;

			System.out.println(">Collecting new proteins");
			File newProteins = new File(markerOutputFolder + File.separator + "newProteins.faa");
			newProteins.createNewFile();
			newProteins.deleteOnExit();
			for (File f : markerOutputFolder.listFiles((dir, name) -> name.endsWith("_new.faa")))
				appendToFile(f, newProteins, newAccessions);
			newAccessions = null;

			System.out.println(">Aligning new vs old proteins using DIAMOND");
			rL.setTime();
			File oldDb = DiamondRunner.makedb(oldProteins, cores, diamondBin);
			oldDb.deleteOnExit();
			File tabFile1 = DiamondRunner.blastp(oldDb, newProteins, tmpFile, identity, blockSize, cores, diamondBin);
			genusAliDatabase.addAlignmentTable("Genus" + "_markerTable", null, tabFile1, false);
			tabFile1.delete();
			oldDb.delete();
			rL.getUptime();

			System.out.println(">Aligning old vs new proteins using DIAMOND");
			rL.setTime();
			File newDb = DiamondRunner.makedb(newProteins, cores, diamondBin);
			newDb.deleteOnExit();
			File tabFile2 = DiamondRunner.blastp(newDb, oldProteins, tmpFile, identity, blockSize, cores, diamondBin);
			genusAliDatabase.addAlignmentTable("Genus" + "_markerTable", null, tabFile2, false);
			tabFile2.delete();
			rL.getUptime();

			System.out.println(">Aligning new vs new proteins using DIAMOND");
			rL.setTime();
			File tabFile3 = DiamondRunner.blastp(newDb, newProteins, tmpFile, identity, blockSize, cores, diamondBin);
			genusAliDatabase.addAlignmentTable("Genus" + "_markerTable", null, tabFile3, false);
			tabFile3.delete();
			newDb.delete();
			rL.getUptime();

			oldProteins.delete();
			newProteins.delete();

			System.out.println(">Computing marker proteins for " + faaFiles.size() + " protein sets");
			List<Runnable> selectMarkersThreads = new ArrayList<>();
			faaFilePointer = 0;
			for (int i = 0; i < cores; i++)
				selectMarkersThreads.add(new SelectMarkersThread(tmpFile, identity, aliFolder, markerOutputFolder,
						taxTree, mappingDatabase));
			rL.runThreads(cores, selectMarkersThreads, totalFileLength);

			// removing new files
			for (Runnable t : collectNewProteinsThreads) {
				List<File> toDelete = ((CollectNewProteinsThread) t).newFiles;
				toDelete.addAll(((CollectNewProteinsThread) t).oldFiles);
				toDelete.stream().forEach(f -> f.delete());
			}

		} catch (Exception e) {
			e.printStackTrace();
		}

		FileUtils.deleteDirectory(markerClusterFolder.getAbsolutePath());

	}

	private synchronized File nextFaaFile() {
		if (faaFilePointer < faaFiles.size())
			return faaFiles.get(faaFilePointer++);
		return null;
	}

	public File getMarkerOutputFolder() {
		return markerOutputFolder;
	}

	private class SelectMarkersThread implements Runnable {

		private String outFolder, aliFolder;
		private File tmpFile;
		private int identity;
		private TaxTree taxTree;
		private SQLMappingDatabase mappingDatabase;

		public SelectMarkersThread(File tmpFile, int identity, String aliFolder, File outFolder, TaxTree taxTree,
				SQLMappingDatabase mappingDatabase) {
			this.tmpFile = tmpFile;
			this.identity = identity;
			this.aliFolder = aliFolder;
			this.outFolder = outFolder.getAbsolutePath();
			this.taxTree = taxTree;
			this.mappingDatabase = createMappingDatabase(mappingDatabase);
		}

		@Override
		public void run() {
			File faaFile;
			while ((faaFile = nextFaaFile()) != null) {
				try {
					String genus = Formatter
							.removeNonAlphanumerics(faaFile.getName().replaceAll("_clustered\\.faa", ""));
					SQLAlignmentDatabase alignmentDatabase = createAlignmentDatabase(aliFolder, "Genus", tmpFile);
					try {
						File outFile = new File(outFolder + File.separator + genus + "_marker.faa");
						new Selecting().run(genus, taxTree, mappingDatabase, alignmentDatabase, faaFile, outFile,
								identity, identity);
					} finally {
						alignmentDatabase.close();
					}
				} catch (Exception e) {
					e.printStackTrace();
				}
				rL.reportProgress(faaFile.length());
			}
			mappingDatabase.close();
			rL.countDown();
		}

	}

	private class CollectNewProteinsThread implements Runnable {

		private File tmpFile;
		private String outFolder, aliFolder;
		private List<File> newFiles = new ArrayList<>();
		private List<File> oldFiles = new ArrayList<>();

		public CollectNewProteinsThread(File tmpFile, String outFolder, String aliFolder) {
			this.tmpFile = tmpFile;
			this.outFolder = outFolder;
			this.aliFolder = aliFolder;
		}

		@Override
		public void run() {

			File faaFile;
			while ((faaFile = nextFaaFile()) != null) {
				try {
					String genus = Formatter
							.removeNonAlphanumerics(faaFile.getName().replaceAll("_clustered\\.faa", ""));
					SQLAlignmentDatabase sqlAliDatabase = createAlignmentDatabase(aliFolder, "Genus", tmpFile);
					try {
						File newFile = new File(outFolder + File.separator + genus + "_new.faa");
						newFile.createNewFile();
						newFiles.add(newFile);
						File oldFile = new File(outFolder + File.separator + genus + "_old.faa");
						oldFile.createNewFile();
						oldFiles.add(oldFile);
						try (BufferedWriter newWriter = new BufferedWriter(new FileWriter(newFile));
								BufferedWriter oldWriter = new BufferedWriter(new FileWriter(oldFile))) {
							ArrayList<FastaEntry> tokens = FastaReader.read(faaFile);
							for (FastaEntry o : tokens) {
								String acc = o.getName();
								String entry = ">" + acc + "\n" + o.getSequence() + "\n";
								if (!sqlAliDatabase.containsAcc(acc, "Genus_markerTable"))
									newWriter.write(entry);
								else
									oldWriter.write(entry);
							}
						}
					} finally {
						sqlAliDatabase.close();
					}
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

	private void appendToFile(File source, File target, Set<String> addedAccession) throws IOException {
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(target, true))) {
			for (FastaEntry token : FastaReader.read(source)) {
				if (!addedAccession.contains(token.getName())) {
					addedAccession.add(token.getName());
					writer.write(">" + token.getName() + "\n" + token.getSequence() + "\n");
				}
			}
		}
	}

}
