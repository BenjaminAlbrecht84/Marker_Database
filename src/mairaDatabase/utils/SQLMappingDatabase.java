package mairaDatabase.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.sql.*;
import java.util.ArrayList;
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

	public ArrayList<Integer> getTaxIdByAcc(String acc) {
		String sql = "SELECT taxid FROM acc2gcf WHERE acc='" + acc + "'";
		try {
			ResultSet rs = stmt.executeQuery(sql);
			return toIntArrayList(rs);
		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "Accession: " + acc);
		}
		return null;
	}

	public ArrayList<String> getGCFByAcc(String acc) {
		String sql = "SELECT gcf FROM acc2gcf WHERE acc='" + acc + "'";
		try {
			ResultSet rs = stmt.executeQuery(sql);
			return toStringArrayList(rs);
		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "Accession: " + acc);
		}
		return null;
	}

	public ArrayList<Integer> getTaxIDByGCF(String gcf) {
		String sql = "SELECT taxid FROM acc2gcf WHERE gcf='" + gcf + "'";
		try {
			ResultSet rs = stmt.executeQuery(sql);
			return toIntArrayList(rs);
		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "Accession: " + gcf);
		}
		return null;
	}

	public ArrayList<String> getACCByGCF(String gcf) {
		String sql = "SELECT acc FROM acc2gcf WHERE gcf='" + gcf + "'";
		try {
			ResultSet rs = stmt.executeQuery(sql);
			return toStringArrayList(rs);
		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "Accession: " + gcf);
		}
		return null;
	}

	private ArrayList<String> toStringArrayList(ResultSet rs) throws SQLException {
		ArrayList<String> result = new ArrayList<>();
		while (rs.next())
			result.add(rs.getString(1));
		return result;
	}

	private ArrayList<Integer> toIntArrayList(ResultSet rs) throws SQLException {
		ArrayList<Integer> result = new ArrayList<>();
		while (rs.next())
			result.add(rs.getInt(1));
		return result;
	}

	public void createAcc2gcf2taxidTable(File acc2gcf2taxid) {
		try {

			System.out.println(">Creating SQL Table acc2gcf2taxid");
			rL.setTime();
			stmt.execute("CREATE TABLE IF NOT EXISTS acc2gcf (acc TEXT, gcf TEXT, taxid integer)");
			stmt.execute("DELETE FROM acc2gcf");
			c.setAutoCommit(false);
			int count = 0;
			rL.setMaxProgress(acc2gcf2taxid.length());
			try (PreparedStatement insertStmd = c.prepareStatement("INSERT INTO acc2gcf VALUES (?, ?, ?);");
					BufferedReader buf = new BufferedReader(new FileReader(acc2gcf2taxid));) {
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
