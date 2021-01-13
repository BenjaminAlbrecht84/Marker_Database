package mairaDatabase.refseq.step0_downloading;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import jloda.util.Pair;
import mairaDatabase.refseq.utils.AssemblyParser;
import mairaDatabase.refseq.utils.Downloader;
import mairaDatabase.utils.FastaReader;
import mairaDatabase.utils.FastaReader.FastaEntry;
import mairaDatabase.utils.ResourceLoader;
import mairaDatabase.utils.SQLMappingDatabase;
import mairaDatabase.utils.taxTree.TaxNode;
import mairaDatabase.utils.taxTree.TaxTree;

public class ProteinDownloadManager {

	private ResourceLoader rL = new ResourceLoader();
	private int downloadThreads = 1;
	private List<Runnable> threads;
	
	private File proteinFolder;
	private List<Pair<String, String>> ftpLinks;
	private int ftpPointer = 0;
	private File[] proteinDirs;
	private int dirPointer = 0;

	private AccessionMapper accessionMapper;

	public void run(String srcPath, SQLMappingDatabase mappingDatabase, TaxTree taxTree, int cores, String[] genera) {

		this.downloadThreads = cores * 2;

		try {

			File refseqProteins = new File(srcPath + File.separator + "refseq_proteins");
			List<String> genusList = genera == null ? null : Stream.of(genera).collect(Collectors.toList());

			// adding FTP address to leaves
			System.out.println(">Parsing assembly summary table ");
			rL.setTime();
			Map<Integer, ArrayList<String>> taxidToFTP = new AssemblyParser().getTaxidToFTP(srcPath);
			ftpLinks = new ArrayList<>();
			for (TaxNode v : taxTree.getNodes()) {
				String genus = getRank(v, "genus");
				if (taxidToFTP.containsKey(v.getTaxid()) && isBacterialOrArchaeaNode(v) && genus != null
						&& (genusList == null || genusList.contains(genus))) {
					if (taxidToFTP.containsKey(v.getTaxid())) {
						v.setInfo(taxidToFTP.get(v.getTaxid()));
						taxidToFTP.get(v.getTaxid()).stream()
								.forEach(ftp -> ftpLinks.add(new Pair<String, String>(ftp, genus)));
						v.reportFTPLeaf();
					}
				}
			}
			System.out.println("Runtime: " + rL.getUptime());

			// downloading protein data
			System.out.println(">Downloading " + ftpLinks.size() + " protein files");
			refseqProteins.mkdir();
			Map<String, File> proteinToFile = new ConcurrentHashMap<String, File>();
			storeExistingFilesRec(refseqProteins, proteinToFile);
			threads = new ArrayList<>();
			for (int i = 0; i < downloadThreads; i++)
				threads.add(new ProteinDownloadingThread(refseqProteins, proteinToFile));
			rL.runThreads(downloadThreads, threads, ftpLinks.size());

			// updating mapping database
			accessionMapper = new AccessionMapper();
			accessionMapper.run(srcPath, refseqProteins, cores, mappingDatabase, taxTree);

			// collecting protein data
			System.out
					.println(">Creating " + refseqProteins.listFiles(File::isDirectory).length + " genus protein sets");
			proteinFolder = new File(srcPath + File.separator + "genus_database_proteins");
			proteinFolder.mkdir();
			threads = new ArrayList<>();
			proteinDirs = refseqProteins.listFiles(File::isDirectory);
			dirPointer = 0;
			for (int i = 0; i < cores; i++)
				threads.add(new ProteinCollector(proteinFolder));
			int numOfFiles = 0;
			for (File dir : refseqProteins.listFiles(File::isDirectory))
				numOfFiles += dir.listFiles().length;
			rL.runThreads(cores, threads, numOfFiles);

		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	private synchronized File nextProteinDir() {
		if (dirPointer < proteinDirs.length)
			return proteinDirs[dirPointer++];
		return null;
	}

	public class ProteinCollector implements Runnable {

		private File proteinFolder;

		public ProteinCollector(File proteinFolder) {
			this.proteinFolder = proteinFolder;
		}

		@Override
		public void run() {

			File dir;
			while ((dir = nextProteinDir()) != null) {

				Map<String, Map<Integer, Set<String>>> seenAccessionMap = new HashMap<>();
				File out = new File(proteinFolder + "/" + dir.getName().replace("\\s+", "_") + ".faa");
				if (out.exists()) {
					ArrayList<FastaEntry> tokens = FastaReader.read(out);
					for (FastaEntry e : tokens) {
						if (!e.getSequence().isEmpty())
							addSeenAccession(seenAccessionMap, e.getName());
					}
				}

				try {
					BufferedWriter writer = new BufferedWriter(new FileWriter(out, true));
					File[] files = dir.listFiles();
					try {
						for (File f : files) {
							if (f.isFile()) {
								ArrayList<FastaEntry> tokens = FastaReader.read(f);
								for (FastaEntry e : tokens) {
									String acc = e.getName();
									String seq = e.getSequence();
									if (!isSeenAccession(seenAccessionMap, acc)) {
										addSeenAccession(seenAccessionMap, acc);
										writer.write(">" + acc + "\n" + seq + "\n");
									}
								}
							}
							rL.reportProgress(1);
						}
					} finally {
						writer.close();
					}
				} catch (Exception e) {
					e.printStackTrace();
				}

			}
			rL.countDown();
		}

		private void addSeenAccession(Map<String, Map<Integer, Set<String>>> seenAccessionMap, String acc) {
			String prefix = acc.split("_")[0];
			int id = Integer.parseInt(acc.split("_")[1].replaceAll("\\.", ""));
			int bucket = id % 10;
			seenAccessionMap.putIfAbsent(prefix, new HashMap<>());
			seenAccessionMap.get(prefix).putIfAbsent(bucket, new HashSet<>());
			seenAccessionMap.get(prefix).get(bucket).add(acc);
		}

		private boolean isSeenAccession(Map<String, Map<Integer, Set<String>>> seenAccessionMap, String acc) {
			String prefix = acc.split("_")[0];
			int id = acc.contains("_") ? Integer.parseInt(acc.split("_")[1].replaceAll("\\.", "")) : 0;
			int bucket = id % 10;
			if (!seenAccessionMap.containsKey(prefix) || !seenAccessionMap.get(prefix).containsKey(bucket)
					|| !seenAccessionMap.get(prefix).get(bucket).contains(acc))
				return false;
			return true;
		}

	}

	private void storeExistingFilesRec(File file, Map<String, File> proteinToFile) {
		if (file.isDirectory()) {
			for (File f : file.listFiles())
				storeExistingFilesRec(f, proteinToFile);
		} else
			proteinToFile.put(file.getName(), file);
	}

	private boolean isBacterialOrArchaeaNode(TaxNode v) {
		TaxNode p = v;
		while (p != null) {
			if (p.getName().equals("Bacteria") || p.getName().equals("Archaea"))
				return true;
			p = p.getParent();
		}
		return false;
	}

	private synchronized List<Pair<String, String>> nextFTPLinks() {
		if (ftpPointer < ftpLinks.size()) {
			int from = ftpPointer;
			int to = Math.min(ftpPointer + 100, ftpLinks.size());
			ftpPointer = to;
			return ftpLinks.subList(from, to);
		}
		return null;
	}

	public class ProteinDownloadingThread implements Runnable {

		private Map<String, File> proteinToFile;
		private File proteinFolder;

		public ProteinDownloadingThread(File proteinFolder, Map<String, File> proteinToFile) {
			this.proteinFolder = proteinFolder;
			this.proteinToFile = proteinToFile;
		}

		@Override
		public void run() {

			try {

				List<Pair<String, String>> ftpLinks;
				while ((ftpLinks = nextFTPLinks()) != null) {

					for (Pair<String, String> ftpLink : ftpLinks) {
						
						String ftp = ftpLink.get1();
						String genus = ftpLink.get2();

						boolean isReadable = false;
						while (!isReadable) {

							String[] split = ftp.split("/");
							String outFile = split[split.length - 1] + "_protein.ftp";
							String remoteFaaFile = split[split.length - 1] + "_protein.faa.gz";
							File localFolder = new File(proteinFolder + File.separator + genus.replaceAll("\\s+", "_"));
							if (!proteinToFile.containsKey(remoteFaaFile)) {
								localFolder.mkdir();
								Downloader.getFtpFile(ftp, remoteFaaFile, localFolder, outFile);
							}

							File localFaaFile = new File(localFolder + File.separator + remoteFaaFile);
							if (FastaReader.read(localFaaFile) == null) {
								proteinToFile.remove(remoteFaaFile);
								localFaaFile.delete();
							} else
								isReadable = true;

						}

						rL.reportProgress(1);

					}

				}

			} catch (Exception e) {
				rL.countDown();
				e.printStackTrace();
			}
			
			rL.countDown();

		}

	}
	
	public File getProteinFolder() {
		return proteinFolder;
	}

	private String getRank(TaxNode v, String rank) {
		TaxNode w = v;
		while (w != null && w.getRank() != null) {
			if (w.getRank().equals(rank))
				return w.getName();
			w = w.getParent();
		}
		return null;
	}

	public File getMappingFile() {
		return accessionMapper.getMappingFile();
	}

	public File getProteinCountsFile() {
		return accessionMapper.getProteinCountsFile();
	}

}
