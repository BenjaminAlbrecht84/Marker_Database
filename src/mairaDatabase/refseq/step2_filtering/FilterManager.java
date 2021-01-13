package mairaDatabase.refseq.step2_filtering;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import mairaDatabase.refseq.utils.Formatter;
import mairaDatabase.refseq.utils.aliHelper.SQLAlignmentDatabase;
import mairaDatabase.utils.FileUtils;
import mairaDatabase.utils.ResourceLoader;
import mairaDatabase.utils.SQLMappingDatabase;
import mairaDatabase.utils.taxTree.TaxTree;

public class FilterManager {

	private ResourceLoader rL = new ResourceLoader();

	private BufferedWriter factorWriter, markerWriter;
	private TaxTree taxTree;
	private int NUM_OF_PROTEINS, MIN_ID;
	private String aliFolder;
	private File tmpDir, weightFile, markerDatabase;

	private List<File> faaFiles;
	private int faaFilePointer = 0;

	public void run(String rank, String srcPath, String aliFolder, File tmpDir, File markerDir, TaxTree taxTree,
			SQLMappingDatabase mappingDatabase, int NUM_OF_PROTEINS, int MIN_ID, int cores) {

		this.taxTree = taxTree;
		this.NUM_OF_PROTEINS = NUM_OF_PROTEINS;
		this.MIN_ID = MIN_ID;
		this.aliFolder = aliFolder;
		this.tmpDir = tmpDir;
		rL.setTime();

		try {
			
			// creating marker writer
			markerDatabase = new File(srcPath + File.separator + "marker_db");
			markerDatabase.mkdir();
			File markerFile = new File(
					markerDatabase.getAbsolutePath() + File.separator + "genus_marker_db_" + NUM_OF_PROTEINS + ".faa");
			markerFile.delete();
			markerWriter = new BufferedWriter(new FileWriter(markerFile, true));
			
			// creating factor writer
			weightFile = new File(srcPath + File.separator + "genus_marker_db_" + NUM_OF_PROTEINS + "_weights.tab");
			weightFile.delete();
			factorWriter = new BufferedWriter(new FileWriter(weightFile, true));

			System.out.println(">Filtering marker proteins");
			try {
				faaFiles = new ArrayList<>(Arrays.asList(markerDir.listFiles()));
				List<Runnable> threads = new ArrayList<>();
				faaFilePointer = 0;
				for (int i = 0; i < cores; i++)
					threads.add(new FilterThread(mappingDatabase));
				rL.runThreads(cores, threads, threads.size());
			} finally {
				markerWriter.close();
				factorWriter.close();
			}

		} catch (Exception e) {
			e.printStackTrace();
		}

		FileUtils.deleteDirectory(markerDir.getAbsolutePath());

	}

	public File getMarkerDatabase() {
		return markerDatabase;
	}

	private synchronized File nextFaaFile() {
		if (faaFilePointer < faaFiles.size())
			return faaFiles.get(faaFilePointer++);
		return null;
	}

	private class FilterThread implements Runnable {

		private SQLMappingDatabase mappingDatabase;

		public FilterThread(SQLMappingDatabase mappingDatabase) {
			this.mappingDatabase = createMappingDatabase(mappingDatabase);
		}

		@Override
		public void run() {
			File faaFile;
			while ((faaFile = nextFaaFile()) != null) {
				try {
					String genus = Formatter.removeNonAlphanumerics(faaFile.getName().replaceAll("_marker\\.faa", ""));
					SQLAlignmentDatabase alignmentDatabase = createAlignmentDatabase(genus);
					try {
						new Filtering().run(faaFile, genus, factorWriter, markerWriter, taxTree, mappingDatabase,
								alignmentDatabase, NUM_OF_PROTEINS, MIN_ID);
					} finally {
						alignmentDatabase.close();
					}
				} catch (ClassNotFoundException | SQLException e) {
					e.printStackTrace();
				}
				rL.reportProgress(1);
			}
			mappingDatabase.close();
			rL.countDown();
		}
	}

	private synchronized SQLMappingDatabase createMappingDatabase(SQLMappingDatabase mappingDatabase) {
		return new SQLMappingDatabase(mappingDatabase);
	}

	private synchronized SQLAlignmentDatabase createAlignmentDatabase(String genus)
			throws ClassNotFoundException, SQLException {
		return new SQLAlignmentDatabase(aliFolder, genus, tmpDir);
	}

	public File getWeightFile() {
		return weightFile;
	}

}
