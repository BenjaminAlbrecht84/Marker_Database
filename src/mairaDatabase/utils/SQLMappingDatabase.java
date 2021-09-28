package mairaDatabase.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.sql.*;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

public class SQLMappingDatabase {

	private ResourceLoader rL = new ResourceLoader();

	private String databaseFile;
	private File tmpDir;
	private Connection c;
	private Statement stmt;

	public SQLMappingDatabase(String databaseFile, File tmpDir) {
		init(databaseFile, tmpDir);
	}

	public SQLMappingDatabase(SQLMappingDatabase mappingDatabase) {
		init(mappingDatabase.getDatabaseFile(), mappingDatabase.getTmpDir());
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

	public List<Integer> getTaxIdByAcc(String acc) {
		String sql = "SELECT taxid FROM acc2gcf WHERE acc='" + acc + "'";
		try {
			ResultSet rs = stmt.executeQuery(sql);
			return toIntArrayList(rs);
		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "Accession: " + acc);
		}
		return null;
	}

	public List<String> getGCFByAcc(String acc) {
		String sql = "SELECT gcf FROM acc2gcf WHERE acc='" + acc + "'";
		try {
			ResultSet rs = stmt.executeQuery(sql);
			return toStringArrayList(rs);
		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "Accession: " + acc);
		}
		return null;
	}

	public Integer getTaxIDByGCF(String gcf) {
		String sql = "SELECT taxid FROM gcf2taxid WHERE gcf='" + gcf + "'";
		try {
			ResultSet rs = stmt.executeQuery(sql);
			return rs.getInt(1);
		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "GCF: " + gcf);
		}
		return null;
	}
	
	public Integer getAvgSizeByTaxid(int taxid) {
		String sql = "SELECT avg_size FROM species2size WHERE taxid='" + taxid + "'";
		try {
			ResultSet rs = stmt.executeQuery(sql);
			return rs.getInt(1);
		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "Taxid: " + taxid);
		}
		return null;
	}
	
	public Integer getSizeByGCF(String gcf) {
		String sql = "SELECT size FROM gcf2taxid WHERE gcf='" + gcf + "'";
		try {
			ResultSet rs = stmt.executeQuery(sql);
			return rs.getInt(1);
		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "GCF: " + gcf);
		}
		return null;
	}

	private List<String> toStringArrayList(ResultSet rs) throws SQLException {
		List<String> result = new ArrayList<>();
		while (rs.next())
			result.add(rs.getString(1));
		return result;
	}

	private List<Integer> toIntArrayList(ResultSet rs) throws SQLException {
		List<Integer> result = new ArrayList<>();
		while (rs.next())
			result.add(rs.getInt(1));
		return result;
	}
	
	public void createSpecies2sizeTable(File speciesSizeFile) {
		try {

			System.out.println(">Creating SQL Table species2size");
			rL.setTime();
			stmt.execute("CREATE TABLE IF NOT EXISTS species2size (taxid INTEGER, avg_size INTEGER)");
			stmt.execute("DELETE FROM species2size");
			c.setAutoCommit(false);
			int count = 0;
			rL.setMaxProgress(speciesSizeFile.length());
			try (PreparedStatement insertStmd = c.prepareStatement("INSERT INTO species2size VALUES (?, ?);");
					BufferedReader buf = new BufferedReader(new FileReader(speciesSizeFile));) {
				String line;
				while ((line = buf.readLine()) != null) {
					final String[] tokens = line.split("\t");
					final int taxid = Integer.parseInt(tokens[0]);
					final int size = Integer.parseInt(tokens[1]);
					insertStmd.setInt(1, taxid);
					insertStmd.setInt(2, size);
					insertStmd.execute();
					count++;
					rL.reportProgress(line.length() + 1);
				}
			}
			c.commit();
			c.setAutoCommit(true);
			System.out.println(String.format("Table species2size: added %,d items", count));
			rL.reportFinish();
			rL.reportRuntime();

			System.out.println(">Indexing SQL Table species2size");
			rL.setTime();
			stmt.execute("DROP INDEX IF EXISTS species2sizeIndex");
			stmt.execute("CREATE INDEX species2sizeIndex ON species2size (taxid)");
			rL.reportFinish();
			rL.reportRuntime();

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void createGcf2size2taxidTable(File gcf2size2taxidFile) {
		try {

			System.out.println(">Creating SQL Table gcf2taxid");
			rL.setTime();
			stmt.execute("CREATE TABLE IF NOT EXISTS gcf2taxid (gcf TEXT, size INTEGER, taxid INTEGER)");
			stmt.execute("DELETE FROM gcf2taxid");
			c.setAutoCommit(false);
			int count = 0;
			rL.setMaxProgress(gcf2size2taxidFile.length());
			try (PreparedStatement insertStmd = c.prepareStatement("INSERT INTO gcf2taxid VALUES (?, ?, ?);");
					BufferedReader buf = new BufferedReader(new FileReader(gcf2size2taxidFile));) {
				String line;
				while ((line = buf.readLine()) != null) {
					final String[] tokens = line.split("\t");
					final String gcf = tokens[0];
					final int size = Integer.parseInt(tokens[1]);
					final int taxid = Integer.parseInt(tokens[2]);
					insertStmd.setString(1, gcf);
					insertStmd.setInt(2, size);
					insertStmd.setInt(3, taxid);
					insertStmd.execute();
					count++;
					rL.reportProgress(line.length() + 1);
				}
			}
			c.commit();
			c.setAutoCommit(true);
			System.out.println(String.format("Table gcf2taxid: added %,d items", count));
			rL.reportFinish();
			rL.reportRuntime();

			System.out.println(">Indexing SQL Table gcf2taxid");
			rL.setTime();
			stmt.execute("DROP INDEX IF EXISTS gcf2taxidIndex");
			stmt.execute("CREATE INDEX gcf2taxidIndex ON gcf2taxid (gcf)");
			rL.reportFinish();
			rL.reportRuntime();

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void createAcc2gcf2taxidTable(File acc2gcf2taxidFile) {
		try {

			System.out.println(">Creating SQL Table acc2gcf");
			rL.setTime();
			stmt.execute("CREATE TABLE IF NOT EXISTS acc2gcf (acc TEXT, gcf TEXT, taxid integer)");
			stmt.execute("DELETE FROM acc2gcf");
			c.setAutoCommit(false);
			int count = 0;
			rL.setMaxProgress(acc2gcf2taxidFile.length());
			try (PreparedStatement insertStmd = c.prepareStatement("INSERT INTO acc2gcf VALUES (?, ?, ?);");
					BufferedReader buf = new BufferedReader(new FileReader(acc2gcf2taxidFile));) {
				String line;
				while ((line = buf.readLine()) != null) {
					final String[] tokens = line.split("\t");
					final String accession = tokens[0];
					final String gcf = tokens[1];
					final int taxid = Integer.parseInt(tokens[2]);
					insertStmd.setString(1, accession);
					insertStmd.setString(2, gcf);
					insertStmd.setInt(3, taxid);
					insertStmd.execute();
					count++;
					rL.reportProgress(line.length() + 1);
				}
			}
			c.commit();
			c.setAutoCommit(true);
			System.out.println(String.format("Table acc2gcf2taxid: added %,d items", count));
			rL.reportFinish();
			rL.reportRuntime();

			System.out.println(">Indexing SQL Table acc2gcf2taxid");
			rL.setTime();
			stmt.execute("DROP INDEX IF EXISTS acc2gcfIndex");
			stmt.execute("CREATE INDEX acc2gcfIndex ON acc2gcf (acc,gcf)");
			rL.reportFinish();
			rL.reportRuntime();

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public String getDatabaseFile() {
		return databaseFile;
	}

	public File getTmpDir() {
		return tmpDir;
	}

	public void close() {
		try {
			c.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
	}

}
