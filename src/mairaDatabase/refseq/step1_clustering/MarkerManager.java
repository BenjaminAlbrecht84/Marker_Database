package mairaDatabase.refseq.step1_clustering;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import mairaDatabase.refseq.utils.DiamondRunner;
import mairaDatabase.refseq.utils.Formatter;
import mairaDatabase.refseq.utils.aliHelper.SQLAlignmentDatabase;
import mairaDatabase.utils.FastaReader;
import mairaDatabase.utils.FileUtils;
import mairaDatabase.utils.ResourceLoader;
import mairaDatabase.utils.SQLMappingDatabase;
import mairaDatabase.utils.FastaReader.FastaEntry;
import mairaDatabase.utils.taxTree.TaxNode;
import mairaDatabase.utils.taxTree.TaxTree;

public class MarkerManager {

	private ResourceLoader rL = new ResourceLoader();
	private File markerOutputFolder;

	public void runMarker(String rank, String srcPath, String aliFolder, File markerClusterFolder, TaxTree taxTree,
			SQLMappingDatabase mappingDatabase, int cores, int memory, int identity, File tmpFile, String diamondBin) {

		try {

			markerOutputFolder = new File(srcPath + File.separator + rank + "_marker_proteins");
			markerOutputFolder.mkdir();
			List<File> faaFiles = new ArrayList<>(Arrays.asList(
					markerClusterFolder.listFiles((dir, name) -> name.endsWith(".faa") && !name.endsWith("_new.faa"))));
			Collections.sort(faaFiles, Comparator.comparingLong(File::length).reversed());
			long totalFileLength = 0L;
			for (File f : faaFiles)
				totalFileLength += f.length();

			System.out.println(">Assessing newly computed representatives for " + faaFiles.size() + " protein sets");
			rL.setTime();
			List<Runnable> collectNewProteinsThreads = new ArrayList<>();
			for (File faaFile : faaFiles)
				collectNewProteinsThreads.add(new CollectNewProteinsThread(faaFile, tmpFile,
						markerOutputFolder.getAbsolutePath(), aliFolder));
			rL.runThreads(cores, collectNewProteinsThreads, totalFileLength);

			System.out.println(">Aligning new vs all proteins using DIAMOND");
			rL.setTime();
			File clusteredProteins = new File(markerOutputFolder + File.separator + "cluster.faa");
			clusteredProteins.createNewFile();
			clusteredProteins.deleteOnExit();
			for (File f : markerClusterFolder.listFiles((dir, name) -> name.endsWith("_clustered.faa")))
				appendToFile(f, clusteredProteins);
			File clusterDb = DiamondRunner.makedb(clusteredProteins, cores, diamondBin);
			clusteredProteins.delete();
			clusterDb.deleteOnExit();
			File newProteins = new File(markerOutputFolder + File.separator + "newProteins.faa");
			newProteins.createNewFile();
			newProteins.deleteOnExit();
			for (File f : markerOutputFolder.listFiles((dir, name) -> name.endsWith("_new.faa")))
				appendToFile(newProteins, f);
			File tabFile1 = DiamondRunner.blastp(clusterDb, newProteins, tmpFile, identity, memory, cores, diamondBin);
			newProteins.delete();
			clusterDb.delete();

			System.out.println(">Aligning all vs new proteins using DIAMOND");
			File newClusteredProteins = new File(markerOutputFolder + File.separator + "newCluster.faa");
			newClusteredProteins.createNewFile();
			newClusteredProteins.deleteOnExit();
			for (File f : markerOutputFolder.listFiles((dir, name) -> name.endsWith("_new.faa")))
				appendToFile(f, newClusteredProteins);
			File newClusterDb = newClusteredProteins.exists() ? DiamondRunner.makedb(newClusteredProteins, cores, diamondBin)
					: null;
			newClusterDb.deleteOnExit();
			newClusteredProteins.delete();
			File allProteins = new File(markerOutputFolder + File.separator + "allProteins.faa");
			allProteins.createNewFile();
			allProteins.deleteOnExit();
			for (File f : faaFiles)
				appendToFile(allProteins, f);
			File tabFile2 = DiamondRunner.blastp(newClusterDb, allProteins, tmpFile, identity, memory, cores, diamondBin);
			allProteins.delete();
			newClusterDb.delete();
			System.out.println("Runtime: " + rL.getUptime());

			System.out.println(">Processing DIAMOND result for " + faaFiles.size() + " genera");
			List<Runnable> processingTabFileThreads = new ArrayList<>();
			File[] tabFiles = { tabFile1, tabFile2 };
			for (File faaFile : faaFiles)
				processingTabFileThreads.add(new ProcessingTabFileThread(faaFile, tabFiles, aliFolder, tmpFile, taxTree,
						mappingDatabase, markerOutputFolder));
			rL.runThreads(1, processingTabFileThreads, totalFileLength);

			System.out.println(">Computing marker proteins for " + faaFiles.size() + " protein sets");
			List<Runnable> selectMarkersThreads = new ArrayList<>();
			for (File faaFile : faaFiles)
				selectMarkersThreads.add(new SelectMarkersThread(faaFile, tmpFile, identity, aliFolder,
						markerOutputFolder, taxTree, mappingDatabase));
			rL.runThreads(cores, selectMarkersThreads, totalFileLength);

			// removing new files
			for (Runnable t : collectNewProteinsThreads)
				((CollectNewProteinsThread) t).newFile.delete();

		} catch (IOException e) {
			e.printStackTrace();
		}

		FileUtils.deleteDirectory(markerClusterFolder.getAbsolutePath());

	}

	public File getMarkerOutputFolder() {
		return markerOutputFolder;
	}

	private class SelectMarkersThread implements Runnable {

		private String outFolder, aliFolder;
		private File faaFile, tmpFile;
		private int identity;
		private TaxTree taxTree;
		private SQLMappingDatabase mappingDatabase;

		public SelectMarkersThread(File faaFile, File tmpFile, int identity, String aliFolder, File outFolder,
				TaxTree taxTree, SQLMappingDatabase mappingDatabase) {
			this.faaFile = faaFile;
			this.tmpFile = tmpFile;
			this.identity = identity;
			this.aliFolder = aliFolder;
			this.outFolder = outFolder.getAbsolutePath();
			this.taxTree = taxTree;
			this.mappingDatabase = createMappingDatabase(mappingDatabase);
		}

		@Override
		public void run() {
			String genus = Formatter.removeNonAlphanumerics(faaFile.getName().replaceAll("_clustered\\.faa", ""));
			SQLAlignmentDatabase aliDatabase = createAlignmentDatabase(aliFolder, genus, tmpFile);
			File outFile = new File(outFolder + File.separator + genus + "_marker.faa");
			new Selecting().run(genus, taxTree, mappingDatabase, aliDatabase, faaFile, outFile, identity, identity);
			aliDatabase.close();
			mappingDatabase.close();
			rL.reportProgress(faaFile.length());
			rL.countDown();
		}

	}

	private class ProcessingTabFileThread implements Runnable {

		private File[] tabFiles;
		private File faaFile, tmpFile, outFolder;
		private String aliFolder;
		private TaxTree taxTree;
		private SQLMappingDatabase mappingDatabase;

		public ProcessingTabFileThread(File faaFile, File[] tabFiles, String aliFolder, File tmpFile, TaxTree taxTree,
				SQLMappingDatabase mappingDatabase, File outFolder) {
			this.faaFile = faaFile;
			this.tabFiles = tabFiles;
			this.aliFolder = aliFolder;
			this.tmpFile = tmpFile;
			this.taxTree = taxTree;
			this.mappingDatabase = createMappingDatabase(mappingDatabase);
			this.outFolder = outFolder;
		}

		@Override
		public void run() {

			String genus = Formatter.removeNonAlphanumerics(faaFile.getName().replaceAll("_clustered\\.faa", ""));
			File tabFile = new File(outFolder + File.separator + genus + "_marker.tab");
			tabFile.deleteOnExit();
			try {
				BufferedWriter writer = new BufferedWriter(new FileWriter(tabFile));
				for (File tab : tabFiles) {
					if (tab == null || !tab.exists() || tab.length() == 0)
						continue;
					BufferedReader buf = new BufferedReader(new FileReader(tab));
					String line;
					while ((line = buf.readLine()) != null) {
						final String[] tokens = line.split("\t");
						final String qacc = tokens[0];
						final String g = getRank(mappingDatabase.getTaxIdByAcc(qacc), "genus");
						if (g != null && g.equals(genus))
							writer.write(line + "\n");
					}
					buf.close();
				}
				writer.close();
			} catch (Exception e) {
				e.printStackTrace();
			}

			SQLAlignmentDatabase aliDatabase = createAlignmentDatabase(aliFolder, genus, tmpFile);
			aliDatabase.addAlignmentTable(genus + "_markerTable", null, tabFile, false);
			tabFile.delete();
			aliDatabase.close();
			mappingDatabase.close();
			rL.reportProgress(faaFile.length());
			rL.countDown();

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

	private class CollectNewProteinsThread implements Runnable {

		private File faaFile, tmpFile;
		private String outFolder, aliFolder;
		private File newFile;

		public CollectNewProteinsThread(File faaFile, File tmpFile, String outFolder, String aliFolder) {
			this.faaFile = faaFile;
			this.tmpFile = tmpFile;
			this.outFolder = outFolder;
			this.aliFolder = aliFolder;
		}

		@Override
		public void run() {

			String genus = Formatter.removeNonAlphanumerics(faaFile.getName().replaceAll("_clustered\\.faa", ""));
			SQLAlignmentDatabase sqlAliDatabase = createAlignmentDatabase(aliFolder, genus, tmpFile);

			try {

				newFile = new File(outFolder + File.separator + genus + "_new.faa");
				newFile.createNewFile();
				BufferedWriter writer = new BufferedWriter(new FileWriter(newFile));
				ArrayList<FastaEntry> tokens = FastaReader.read(faaFile);
				for (FastaEntry o : tokens) {
					String acc = o.getName();
					if (!sqlAliDatabase.containsAcc(acc, genus + "_markerTable"))
						writer.write(">" + acc + "\n" + o.getSequence() + "\n");
				}
				writer.close();
			} catch (Exception e) {
				e.printStackTrace();
			}

			sqlAliDatabase.close();
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

	private void appendToFile(File source, File target) {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(target, true));
			for (FastaEntry token : FastaReader.read(source))
				writer.write(">" + token.getName() + "\n" + token.getSequence() + "\n");
			writer.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
