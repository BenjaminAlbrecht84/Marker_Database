package mairaDatabase.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Statement;
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
	
	public void createFactorsTable(File weightFile, int MAX_PROTEINS_PER_GCF) {
		try {
			
			String table = "factors_"+MAX_PROTEINS_PER_GCF;
			String index = "factors_"+MAX_PROTEINS_PER_GCF+"_index";
			
			System.out.println(">Creating SQL Table factors");
			rL.setTime();
			stmt.execute("CREATE TABLE IF NOT EXISTS "+table+" (acc TEXT, factor REAL)");
			stmt.execute("DELETE FROM "+table);
			c.setAutoCommit(false);
			int count = 0;
			rL.setMaxProgress(weightFile.length());
			try (PreparedStatement insertStmd = c.prepareStatement("INSERT INTO "+table+" VALUES (?, ?);");
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
			System.out.println(String.format("Table "+table+": added %,d items", count));
			rL.reportFinish();
			rL.reportRuntime();

			System.out.println(">Indexing SQL Table "+table);
			rL.setTime();
			stmt.execute("DROP INDEX IF EXISTS "+index);
			stmt.execute("CREATE INDEX "+index+" ON "+table+" (acc)");
			rL.reportFinish();
			rL.reportRuntime();

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}
