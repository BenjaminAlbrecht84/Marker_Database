package mairaDatabase.refseq.utils.aliHelper;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import mairaDatabase.refseq.utils.aliHelper.SQLAlignmentDatabase.AlignmentInfo;
import mairaDatabase.utils.ResourceLoader;

public class FilterAlignmentsHelper {

	private static File aliDatabaseFolder;
	private static List<File> dbFiles;
	private static int dbFilePointer = 0;

	private static ResourceLoader rL = new ResourceLoader();

	public static void main(String[] args) {

		aliDatabaseFolder = new File(args[0]);
		int cores = Integer.parseInt(args[1]);

		dbFiles = Arrays
				.asList(aliDatabaseFolder.listFiles((dir, name) -> name.endsWith(".db") && !name.endsWith("_tmp.db")));
		Collections.sort(dbFiles, (f1, f2) -> Long.compare(f1.length(), f2.length()));

		List<Runnable> threads = new ArrayList<>();
		dbFilePointer = 0;
		for (int i = 0; i < cores; i++)
			threads.add(new AlignmentFilter());
		rL.runThreads(cores, threads, dbFiles.size());

	}

	private static synchronized File nextDbFile() {
		if (dbFilePointer < dbFiles.size())
			return dbFiles.get(dbFilePointer++);
		return null;
	}

	private static synchronized SQLAlignmentDatabase createDatabase(String databaseFile, String genus, File tmpDir)
			throws ClassNotFoundException, SQLException {
		return new SQLAlignmentDatabase(databaseFile, genus, tmpDir);
	}

	private static class AlignmentFilter implements Runnable {

		@Override
		public void run() {
			File dbFile1, dbFile2;
			while ((dbFile1 = nextDbFile()) != null) {
				try {
					String genus = dbFile1.getName().split("\\.")[0];
					dbFile2 = new File(aliDatabaseFolder.getAbsolutePath() + File.separator + genus + "_tmp.db");
					dbFile2.delete();
					SQLAlignmentDatabase db1 = createDatabase(aliDatabaseFolder.getAbsolutePath(), genus,
							aliDatabaseFolder);
					SQLAlignmentDatabase db2 = createDatabase(aliDatabaseFolder.getAbsolutePath(), genus + "_tmp",
							aliDatabaseFolder);
					List<AlignmentInfo> distinctAlignments = db1.getAllDistinctAlignments(genus + "_clusterTable");
					File aliTab = new File(
							aliDatabaseFolder.getAbsolutePath() + File.separatorChar + genus + "_tmp.tab");
					try {
						try (BufferedWriter writer = new BufferedWriter(new FileWriter(aliTab))) {
							for (AlignmentInfo ali : distinctAlignments) {
								String line = ali.getQuery() + "\t" + ali.getRef() + "\t" + ali.getQueryStart() + "\t"
										+ (ali.getQueryStart() + ali.getQueryLen()) + "\t" + ali.getQueryLen() + "\t"
										+ ali.getSubjectStart() + "\t" + (ali.getSubjectStart() + ali.getSubjectLen())
										+ "\t" + ali.getSubjectLen() + "\t" + ali.getIdentity() + "\t" + ali.getBtop() + "\n";
								if (line.split("\t").length < 5)
									System.out.println("ERROR: " + line);
								writer.write(line);
							}
						}
						db2.addAlignmentTable(genus + "_clusterTable", null, aliTab, true);
					} finally {
						aliTab.delete();
						Files.move(dbFile2.toPath(), dbFile2.toPath().resolveSibling(dbFile1.getName()),
								StandardCopyOption.REPLACE_EXISTING);
					}

					rL.reportProgress(1);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
			rL.countDown();
		}

	}

}
