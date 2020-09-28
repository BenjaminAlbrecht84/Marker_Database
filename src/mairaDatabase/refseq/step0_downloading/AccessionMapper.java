package mairaDatabase.refseq.step0_downloading;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.concurrent.ConcurrentHashMap;

import mairaDatabase.refseq.utils.AssemblyParser;
import mairaDatabase.utils.FastaReader;
import mairaDatabase.utils.ResourceLoader;
import mairaDatabase.utils.SQLMappingDatabase;
import mairaDatabase.utils.Statistics;
import mairaDatabase.utils.FastaReader.FastaEntry;
import mairaDatabase.utils.taxTree.TaxNode;
import mairaDatabase.utils.taxTree.TaxTree;

public class AccessionMapper {

	private ResourceLoader rL = new ResourceLoader();

	private List<File> faaFiles;
	private int filePointer = 0;
	private File mappingFile;
	private File proteinCountsFile;

	private Map<String, Integer> gcfToTaxID;
	private Map<Integer, List<Double>> proteinCounts;

	public void run(String src, File downloadFolder, int cores, SQLMappingDatabase mappingDatabase, TaxTree taxTree)
			throws IOException {

		mappingFile = new File(src + File.separator + "acc2gcfs2taxid.tab");
		proteinCounts = new ConcurrentHashMap<Integer, List<Double>>();

		gcfToTaxID = new AssemblyParser().getGcf2Taxid(src);
		if (mappingFile.exists())
			mappingFile.delete();

		try (BufferedWriter writer = new BufferedWriter(new FileWriter(mappingFile, true));) {
			System.out.println(">Creating mapping file");
			List<Runnable> mappers = new ArrayList<>();
			for (int i = 0; i < cores; i++)
				mappers.add(new Mapper(writer, taxTree));
			faaFiles = new ArrayList<>();
			for (File dir : Objects.requireNonNull(downloadFolder.listFiles(pathname -> pathname.isDirectory()))) {
				for (File f : Objects.requireNonNull(dir.listFiles((dir1, name) -> name.endsWith(".faa.gz"))))
					faaFiles.add(f);
			}
			rL.runThreads(cores, mappers, faaFiles.size());
		}

		proteinCountsFile = new File(src + File.separator + "taxid2proteins.tab");
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(proteinCountsFile, true));) {
			for (Integer taxid : proteinCounts.keySet()) {
				List<Double> counts = proteinCounts.get(taxid);
				int median = (int) Math.round(Statistics.getMedian(counts));
				writer.write(taxid + "\t" + median + "\n");
			}
		}
		proteinCounts.clear();

		mappingDatabase.createAcc2gcf2taxidTable(mappingFile);

	}

	private synchronized File nextFile() {
		if (filePointer < faaFiles.size())
			return faaFiles.get(filePointer++);
		return null;
	}

	public File getMappingFile() {
		return mappingFile;
	}

	public File getProteinCountsFile() {
		return proteinCountsFile;
	}

	private class Mapper implements Runnable {

		private BufferedWriter writer;
		private TaxTree taxTree;

		public Mapper(BufferedWriter writer, TaxTree taxTree) {
			this.writer = writer;
			this.taxTree = taxTree;
		}

		@Override
		public void run() {
			File faaFile;
			while ((faaFile = nextFile()) != null) {
				String searchID = "GCF_" + faaFile.getName().split("_")[1];
				Integer taxid = gcfToTaxID.get(searchID);
				if (taxid != null) {
					double proteinCounter = 0;
					for (FastaEntry token : FastaReader.read(faaFile)) {
						String acc = token.getName();
						proteinCounter++;
						try {
							writer.write(acc + "\t" + searchID + "\t" + taxid + "\n");
						} catch (IOException e) {
							e.printStackTrace();
						}
					}
					proteinCounts.putIfAbsent(taxid, new ArrayList<Double>());
					proteinCounts.get(taxid).add(proteinCounter);
				}
				rL.reportProgress(1);
			}
			rL.countDown();
		}

		private Integer getRankTaxid(TaxNode v, String rank) {
			TaxNode w = v;
			while (w != null && w.getRank() != null) {
				if (w.getRank().equals(rank))
					return w.getTaxid();
				w = w.getParent();
			}
			return null;
		}

	}

}
