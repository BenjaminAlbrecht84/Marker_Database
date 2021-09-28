package mairaDatabase.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;

public class SQLMairaDatabase {

	private ResourceLoader rL = new ResourceLoader();

	private String databaseFile;
	private File tmpDir;
	private Connection c;
	private Statement stmt;

	public SQLMairaDatabase(String databaseFile, File tmpDir) {
		init(databaseFile, tmpDir);
	}

	private void init(String databaseFile, File tmpDir) {
		try {
			this.databaseFile = databaseFile;
			this.tmpDir = tmpDir;
			Class.forName("org.sqlite.JDBC");
			c = DriverManager.getConnection("jdbc:sqlite:" + this.databaseFile);
			stmt = c.createStatement();
			stmt.execute("PRAGMA temp_store_directory = '" + tmpDir.getAbsolutePath() + "'");
		} catch (ClassNotFoundException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, null, ex);
		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, null, ex);
		}
	}

	public void createProtCountsTable(File protCounts) {
		try {

			System.out.println(">Creating SQL Table prot_counts");
			rL.setTime();
			stmt.execute("CREATE TABLE IF NOT EXISTS prot_counts (taxid INTEGER, median INTEGER)");
			stmt.execute("DELETE FROM prot_counts");
			c.setAutoCommit(false);
			int count = 0;
			rL.setMaxProgress(protCounts.length());
			try (PreparedStatement insertStmd = c.prepareStatement("INSERT INTO prot_counts VALUES (?, ?);");
					BufferedReader buf = new BufferedReader(new FileReader(protCounts));) {
				String line;
				while ((line = buf.readLine()) != null) {
					final String[] tokens = line.split("\t");
					final int taxid = Integer.parseInt(tokens[0]);
					final int median = Integer.parseInt(tokens[1]);
					insertStmd.setInt(1, taxid);
					insertStmd.setInt(2, median);
					insertStmd.execute();
					count++;
					rL.reportProgress(line.length() + 1);
				}
			}
			c.commit();
			c.setAutoCommit(true);
			System.out.println(String.format("Table prot_counts: added %,d items", count));
			rL.reportFinish();
			rL.reportRuntime();

			System.out.println(">Indexing SQL Table prot_counts");
			rL.setTime();
			stmt.execute("DROP INDEX IF EXISTS prot_countsIndex");
			stmt.execute("CREATE INDEX prot_countsIndex ON prot_counts (taxid)");
			rL.reportFinish();
			rL.reportRuntime();

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void createAcc2taxidTable(File acc2gcf2taxid) {
		try {

			System.out.println(">Creating SQL Table acc2taxids");
			rL.setTime();
			stmt.execute("CREATE TABLE IF NOT EXISTS acc2taxids (acc TEXT, taxid integer)");
			stmt.execute("DELETE FROM acc2taxids");
			c.setAutoCommit(false);
			int count = 0;
			rL.setMaxProgress(acc2gcf2taxid.length());
			try (PreparedStatement insertStmd = c.prepareStatement("INSERT INTO acc2taxids VALUES (?, ?);");
					BufferedReader buf = new BufferedReader(new FileReader(acc2gcf2taxid));) {
				String line;
				while ((line = buf.readLine()) != null) {
					final String[] tokens = line.split("\t");
					final String accession = tokens[0];
					final int taxid = Integer.parseInt(tokens[2]);
					insertStmd.setString(1, accession);
					insertStmd.setInt(2, taxid);
					insertStmd.execute();
					count++;
					rL.reportProgress(line.length() + 1);
				}
			}
			c.commit();
			c.setAutoCommit(true);
			System.out.println(String.format("Table acc2taxids: added %,d items", count));
			rL.reportFinish();
			rL.reportRuntime();

			System.out.println(">Indexing SQL Table acc2taxids");
			rL.setTime();
			stmt.execute("DROP INDEX IF EXISTS acc2taxidIndex");
			stmt.execute("CREATE INDEX acc2taxidIndex ON acc2taxids (acc)");
			rL.reportFinish();
			rL.reportRuntime();

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void createSpecies2overlapTable(File species2overlap) {
		try {

			System.out.println(">Creating SQL Table species2overlap");
			rL.setTime();
			stmt.execute("CREATE TABLE IF NOT EXISTS species2overlap (source INTEGER, target INTEGER, max_overlap INTEGER)");
			stmt.execute("DELETE FROM species2overlap");
			c.setAutoCommit(false);
			int count = 0;
			rL.setMaxProgress(species2overlap.length());
			try (PreparedStatement insertStmd = c.prepareStatement("INSERT INTO species2overlap VALUES (?, ?, ?);");
					BufferedReader buf = new BufferedReader(new FileReader(species2overlap));) {
				String line;
				while ((line = buf.readLine()) != null) {
					final String[] tokens = line.split("\t");
					final int source = Integer.parseInt(tokens[0]);
					final int target = Integer.parseInt(tokens[1]);
					final int maxOverlap = Integer.parseInt(tokens[2]);
					insertStmd.setInt(1, source);
					insertStmd.setInt(2, target);
					insertStmd.setInt(3, maxOverlap);
					insertStmd.execute();
					count++;
					rL.reportProgress(line.length() + 1);
				}
			}
			c.commit();
			c.setAutoCommit(true);
			System.out.println(String.format("Table species2overlap: added %,d items", count));
			rL.reportFinish();
			rL.reportRuntime();

			System.out.println(">Indexing SQL Table species2overlap");
			rL.setTime();
			stmt.execute("DROP INDEX IF EXISTS species2overlapIndex");
			stmt.execute("CREATE INDEX species2overlapIndex ON species2overlap (source, target)");
			rL.reportFinish();
			rL.reportRuntime();

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void createFactorsTable(Map<Integer, File> weightFiles) {
		try {

			for (Entry<Integer, File> e : weightFiles.entrySet()) {

				int n = e.getKey();
				File weightFile = e.getValue();

				String table = "factors_" + n;
				String index = "factors_" + n + "_index";

				System.out.println(">Creating SQL Table " + table);
				rL.setTime();
				stmt.execute("CREATE TABLE IF NOT EXISTS " + table + " (acc TEXT, factor REAL)");
				stmt.execute("DELETE FROM " + table);
				c.setAutoCommit(false);
				int count = 0;
				rL.setMaxProgress(weightFile.length());
				try (PreparedStatement insertStmd = c.prepareStatement("INSERT INTO " + table + " VALUES (?, ?);");
						BufferedReader buf = new BufferedReader(new FileReader(weightFile));) {
					String line;
					while ((line = buf.readLine()) != null) {
						final String[] tokens = line.split("\t");
						final String accession = tokens[0];
						final double weight = Double.parseDouble(tokens[1]);
						insertStmd.setString(1, accession);
						insertStmd.setDouble(2, weight);
						insertStmd.execute();
						count++;
						rL.reportProgress(line.length() + 1);
					}
				}
				c.commit();
				c.setAutoCommit(true);
				System.out.println(String.format("Table " + table + ": added %,d items", count));
				rL.reportFinish();
				rL.reportRuntime();

				System.out.println(">Indexing SQL Table " + table);
				rL.setTime();
				stmt.execute("DROP INDEX IF EXISTS " + index);
				stmt.execute("CREATE INDEX " + index + " ON " + table + " (acc)");
				rL.reportFinish();
				rL.reportRuntime();

			}

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void createAcc2DominatorsTable(File genusDominationFile) {
		try {

			String table = "acc2dominator";
			String index = "acc2dominator_index";

			System.out.println(">Creating SQL Table " + table);
			rL.setTime();
			stmt.execute("CREATE TABLE IF NOT EXISTS " + table
					+ " (acc TEXT, dominator TEXT, btop TEXT, qstart INTEGER, sstart INTEGER, slen INTEGER)");
			stmt.execute("DELETE FROM " + table);
			c.setAutoCommit(false);
			int count = 0;
			rL.setMaxProgress(genusDominationFile.length());
			try (PreparedStatement insertStmd = c
					.prepareStatement("INSERT INTO " + table + " VALUES (?, ?, ?, ?, ?, ?)");
					BufferedReader buf = new BufferedReader(new FileReader(genusDominationFile));) {
				String line;
				while ((line = buf.readLine()) != null) {
					final String[] tokens = line.split("\t");
					final String accession = tokens[0];
					final String dominator = tokens[1];
					final String btop = tokens[2];
					final Integer qstart = Integer.parseInt(tokens[3]);
					final Integer sstart = Integer.parseInt(tokens[4]);
					final Integer slen = Integer.parseInt(tokens[5]);
					insertStmd.setString(1, accession);
					insertStmd.setString(2, dominator);
					insertStmd.setString(3, btop);
					insertStmd.setInt(4, qstart);
					insertStmd.setInt(5, sstart);
					insertStmd.setInt(6, slen);
					insertStmd.execute();
					count++;
					rL.reportProgress(line.length() + 1);
				}
			}
			c.commit();
			c.setAutoCommit(true);
			System.out.println(String.format("Table " + table + ": added %,d items", count));
			rL.reportFinish();
			rL.reportRuntime();

			System.out.println(">Indexing SQL Table " + table);
			rL.setTime();
			stmt.execute("DROP INDEX IF EXISTS " + index);
			stmt.execute("CREATE INDEX " + index + " ON " + table + " (dominator)");
			rL.reportFinish();
			rL.reportRuntime();

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void createTree2Newick(String type, String newick) {
		try {

			String table = "tree2newick";
			String index = "tree2newick_index";

			System.out.println(">Creating SQL Table " + table);
			stmt.execute("CREATE TABLE IF NOT EXISTS " + table + " (type TEXT, newick TEXT)");
			stmt.execute("DELETE FROM " + table);
			c.setAutoCommit(false);
			int count = 0;
			try (PreparedStatement insertStmd = c.prepareStatement("INSERT INTO " + table + " VALUES (?, ?);")) {
				insertStmd.setString(1, type);
				insertStmd.setString(2, newick);
				insertStmd.execute();
			}
			c.commit();
			c.setAutoCommit(true);
			System.out.println(String.format("Table " + table + ": added %,d items", count));
			rL.reportFinish();
			rL.reportRuntime();

			System.out.println(">Indexing SQL Table " + table);
			rL.setTime();
			stmt.execute("DROP INDEX IF EXISTS " + index);
			stmt.execute("CREATE INDEX " + index + " ON " + table + " (type)");
			rL.reportFinish();
			rL.reportRuntime();

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
