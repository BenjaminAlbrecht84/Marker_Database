package mairaDatabase.refseq.utils.aliHelper;

import java.io.File;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import mairaDatabase.utils.ResourceLoader;

public class ReduceAlignmentsHelper {

	private static File aliDatabaseFolder;
	private static List<File> dbFiles;
	private static int dbFilePointer = 0;

	private static ResourceLoader rL = new ResourceLoader();

	public static void main(String[] args) {

		aliDatabaseFolder = new File(args[0]);
		int cores = Integer.parseInt(args[1]);
		String genera = args.length >= 3 ? args[2] : null;

		dbFiles = Arrays
				.asList(aliDatabaseFolder.listFiles((dir, name) -> name.endsWith(".db") && !name.endsWith("_tmp.db")));
		if (genera != null)
			dbFiles = dbFiles.stream().filter(f -> genera.contains(f.getName().split("\\.")[0]))
					.collect(Collectors.toList());
		Collections.sort(dbFiles, (f1, f2) -> Long.compare(f1.length(), f2.length()));

		List<Runnable> threads = new ArrayList<>();
		dbFilePointer = 0;
		for (int i = 0; i < cores; i++)
			threads.add(new AlignmentReducer());
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

	private static class AlignmentReducer implements Runnable {

		@Override
		public void run() {
			File dbFile;
			while ((dbFile = nextDbFile()) != null) {
				try {

					long time = System.currentTimeMillis();
					String genus = dbFile.getName().split("\\.")[0];
					SQLAlignmentDatabase db = createDatabase(aliDatabaseFolder.getAbsolutePath(), genus,
							aliDatabaseFolder);

					db.getStmt().setFetchSize(1000000);
					Connection c = db.getConnection();

					createColumn(c, "qacc_id", "INTEGER", genus + "_clusterTable");
					createColumn(c, "racc_id", "INTEGER", genus + "_clusterTable");
					createColumn(c, "qacc_id", "INTEGER", genus + "_markerTable");
					createColumn(c, "racc_id", "INTEGER", genus + "_markerTable");
					createColumn(c, "qgcf_id", "INTEGER", genus + "_overlapTable");
					createColumn(c, "rgcf_id", "INTEGER", genus + "_overlapTable");
					
					for (String aliTable : Arrays.asList(genus + "_clusterTable", genus + "_markerTable")) {
						
						if (existsColumn(c, aliTable, "qacc") && existsColumn(c, aliTable, "racc")) {
							
							c.createStatement().execute(
									"CREATE INDEX IF NOT EXISTS " + aliTable + "_qaccIndex ON " + aliTable + " (qacc)");
							c.createStatement().execute(
									"CREATE INDEX IF NOT EXISTS " + aliTable + "_raccIndex ON " + aliTable + " (racc)");
							c.createStatement().execute("CREATE INDEX IF NOT EXISTS " + aliTable + "_pidentIndex ON "
									+ aliTable + " (pident)");
							
							c.setAutoCommit(false);
							PreparedStatement qAccStmt = c
									.prepareStatement("UPDATE " + aliTable + " SET qacc_id = ? WHERE qacc = ?");
							PreparedStatement rAccStmt = c
									.prepareStatement("UPDATE " + aliTable + " SET racc_id = ? WHERE racc = ?");
							try (ResultSet rs1 = c.createStatement().executeQuery("SELECT qacc, racc FROM " + aliTable
									+ " WHERE qacc_id IS NULL OR racc_id IS NULL")) {
								while (rs1.next()) {

									String qacc = rs1.getString(1);
									String racc = rs1.getString(2);
									int qaccId = db.getAccessionId(qacc);
									int raccId = db.getAccessionId(racc);

									qAccStmt.setInt(1, qaccId);
									qAccStmt.setString(2, qacc);
									qAccStmt.execute();
									rAccStmt.setInt(1, raccId);
									rAccStmt.setString(2, racc);
									rAccStmt.execute();

								}
							}
							c.commit();	
							c.setAutoCommit(true);
							
							c.createStatement().execute("DROP INDEX IF EXISTS " + aliTable + "_qaccIndex");
							c.createStatement().execute("DROP INDEX IF EXISTS " + aliTable + "_raccIndex");
							c.createStatement().execute("DROP INDEX IF EXISTS " + aliTable + "_pidentIndex");
							deleteColumns(c, new ArrayList<String>(Arrays.asList("qacc", "racc")),
									genus + "_clusterTable");
							c.createStatement().execute(
									"CREATE INDEX IF NOT EXISTS " + aliTable + "_qaccIndex ON " + aliTable + " (qacc_id)");
							c.createStatement().execute(
									"CREATE INDEX IF NOT EXISTS " + aliTable + "_raccIndex ON " + aliTable + " (racc_id)");
							c.createStatement().execute("CREATE INDEX IF NOT EXISTS " + aliTable + "_pidentIndex ON "
									+ aliTable + " (pident)");
		

						}
					}

					String overlapTable = genus + "_overlapTable";
					if (existsColumn(c, genus + "_overlapTable", "qgcf")
							&& existsColumn(c, genus + "_overlapTable", "rgcf")) {
						
						c.createStatement().execute("CREATE INDEX IF NOT EXISTS " + overlapTable + "_gcfIndex ON "
								+ overlapTable + " (qgcf, rgcf)");
						
						c.setAutoCommit(false);
						PreparedStatement qGenomeStmt = c
								.prepareStatement("UPDATE " + genus + "_overlapTable SET qgcf_id = ? WHERE qgcf = ?");
						PreparedStatement rGenomeStmt = c
								.prepareStatement("UPDATE " + genus + "_overlapTable SET rgcf_id = ? WHERE rgcf = ?");
						try (ResultSet rs3 = c.createStatement().executeQuery("SELECT qgcf, rgcf FROM " + genus
								+ "_overlapTable WHERE qgcf_id IS NULL OR rgcf_id IS NULL")) {
							while (rs3.next()) {

								String qgcf = rs3.getString(1);
								String rgcf = rs3.getString(2);
								int qgcfId = db.getGenomeId(qgcf);
								int rgcfId = db.getGenomeId(rgcf);

								qGenomeStmt.setInt(1, qgcfId);
								qGenomeStmt.setString(2, qgcf);
								qGenomeStmt.execute();
								rGenomeStmt.setInt(1, rgcfId);
								rGenomeStmt.setString(2, rgcf);
								rGenomeStmt.execute();

							}
						}
						c.commit();	
						c.setAutoCommit(true);
						
						c.createStatement().execute("DROP INDEX IF EXISTS " + overlapTable + "_gcfIndex");
						deleteColumns(c, new ArrayList<String>(Arrays.asList("qgcf", "rgcf")), overlapTable);
						c.createStatement().execute("CREATE INDEX IF NOT EXISTS " + overlapTable + "_gcfIndex ON "
								+ overlapTable + " (qgcf_id, rgcf_id)");

					}

					c.createStatement().execute("VACUUM");
					c.close();

					long runtime = (System.currentTimeMillis() - time) / 1000;
					System.out.println(genus + " " + runtime + "s");
					rL.reportProgress(1);

				} catch (Exception e) {
					System.err.println("ERROR: " + dbFile.getAbsolutePath());
					e.printStackTrace();
				}

			}

			rL.countDown();

		}

		private void createColumn(Connection c, String columnName, String columnType, String tableName)
				throws SQLException {
			if (!existsColumn(c, tableName, columnName))
				c.createStatement()
						.execute("ALTER TABLE " + tableName + " ADD COLUMN " + columnName + " " + columnType);
		}

		private void deleteColumns(Connection c, ArrayList<String> toDelete, String tableName) throws SQLException {
			try (ResultSet rs = c.createStatement().executeQuery("PRAGMA table_info(" + tableName + ")")) {
				List<String> toCopyCol = new ArrayList<>();
				List<String> toCopyDef = new ArrayList<>();
				int nDelColumns = 0;
				while (rs.next()) {
					String col = rs.getString(2);
					String type = rs.getString(3);
					if (!toDelete.contains(col)) {
						toCopyCol.add(col);
						toCopyDef.add(col + " " + type);
					} else
						nDelColumns++;
				}
				if (nDelColumns > 0) {
					String columnString = toCopyCol.stream().collect(Collectors.joining(", "));
					String defString = toCopyDef.stream().collect(Collectors.joining(", "));
					c.createStatement().execute("CREATE TABLE table_tmp(" + defString + ")");
					c.createStatement().execute("INSERT INTO table_tmp SELECT " + columnString + " FROM " + tableName);
					c.createStatement().execute("DROP TABLE " + tableName);
					c.createStatement().execute("CREATE TABLE " + tableName + "(" + defString + ")");
					c.createStatement()
							.execute("INSERT INTO " + tableName + " SELECT " + columnString + " FROM table_tmp");
					c.createStatement().execute("DROP TABLE table_tmp");
				}
			}
		}

		public boolean existsColumn(Connection c, String tableName, String columnName) throws SQLException {
			try (ResultSet rs = c.createStatement().executeQuery("PRAGMA table_info(" + tableName + ")")) {
				while (rs.next()) {
					if (rs.getString(2).equals(columnName))
						return true;
				}
			}
			return false;
		}

	}

}
