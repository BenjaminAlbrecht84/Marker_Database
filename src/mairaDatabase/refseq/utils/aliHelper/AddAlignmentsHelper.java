package mairaDatabase.refseq.utils.aliHelper;

import java.io.File;
import java.sql.SQLException;
import java.util.ArrayList;

import mairaDatabase.refseq.utils.DiamondRunner;
import mairaDatabase.refseq.utils.Formatter;
import mairaDatabase.utils.ResourceLoader;

public class AddAlignmentsHelper {

	private static File[] daaFiles;
	private static int daaPointer = 0;
	private static File tmpDir;

	private static ResourceLoader rL = new ResourceLoader();

	public static void main(String[] args) {

		File clusterAlignmentFolder = new File(args[0]);
		File markerAlignmentFolder = new File(args[1]);
		File databaseFolder = new File(args[2]);
		tmpDir = new File(args[3]);
		int cores = Integer.parseInt(args[4]);
		String diamondBin = args[5];

		daaFiles = markerAlignmentFolder.listFiles((dir, name) -> name.endsWith(".daa") && !name.endsWith("2015.daa"));
		ArrayList<Runnable> threads = new ArrayList<>();
		daaPointer = 0;
		for (int i = 0; i < cores; i++)
			threads.add(new AlignmentsAdder(databaseFolder.getAbsolutePath(), clusterAlignmentFolder.getAbsolutePath(),
					"markerTable", 1, diamondBin));
		rL.runThreads(cores, threads, daaFiles.length);

	}

	private static synchronized File nextDaaFile() {
		if (daaPointer < daaFiles.length)
			return daaFiles[daaPointer++];
		return null;
	}

	private static synchronized SQLAlignmentDatabase createDatabase(String databaseFile, String genus, File tmpDir)
			throws ClassNotFoundException, SQLException {
		return new SQLAlignmentDatabase(databaseFile, genus, tmpDir);
	}

	private static class AlignmentsAdder implements Runnable {

		private String dbFolder, srcFolder, dbType;
		private int cores;
		private String diamondBin;

		public AlignmentsAdder(String dbFolder, String srcFolder, String dbType, int cores, String diamondBin) {
			this.cores = cores;
			this.dbType = dbType;
			this.dbFolder = dbFolder;
			this.srcFolder = srcFolder;
		}

		@Override
		public void run() {
			File daa;
			while ((daa = nextDaaFile()) != null) {
				try {
					String genus = Formatter.removeNonAlphanumerics(
							daa.getName().replaceAll("_clustered", "").replaceAll("\\.daa", ""));
					SQLAlignmentDatabase aliDb = createDatabase(dbFolder, genus, tmpDir);
					File tab = DiamondRunner.view(daa, cores, diamondBin);
					File src = srcFolder != null
							? new File(srcFolder + File.separator + daa.getName().replaceAll("\\.daa", ".faa"))
							: null;
					aliDb.addAlignmentTable(genus + "_" + dbType, src, tab, true);
					tab.delete();
				} catch (Exception e) {
					e.printStackTrace();
				}
				rL.reportProgress(1);
			}
			rL.countDown();
		}
	}

}
