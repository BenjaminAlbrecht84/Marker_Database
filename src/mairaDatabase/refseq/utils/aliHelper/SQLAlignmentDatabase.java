package mairaDatabase.refseq.utils.aliHelper;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.sql.*;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import mairaDatabase.utils.FastaReader;
import mairaDatabase.utils.SQLMappingDatabase;
import mairaDatabase.utils.FastaReader.FastaEntry;

public class SQLAlignmentDatabase {

	private String databaseFile;
	private Connection c;
	private Statement stmt;

	public SQLAlignmentDatabase(String databaseFolder, String genus, File tmpDir)
			throws ClassNotFoundException, SQLException {
		this.databaseFile = databaseFolder + File.separator + genus + ".db";
		Class.forName("org.sqlite.JDBC");
		c = DriverManager.getConnection("jdbc:sqlite:" + this.databaseFile);
		stmt = c.createStatement();
		stmt.execute("PRAGMA temp_store_directory = '" + tmpDir.getAbsolutePath() + "'");
		String clusterTable = genus + "_clusterTable";
		stmt.execute("CREATE TABLE IF NOT EXISTS " + clusterTable
				+ " (qacc TEXT, racc TEXT, qstart INTEGER, qend INTEGER, qlen INTEGER, sstart INTEGER, send INTEGER, slen INTEGER, pident DOUBLE, btop TEXT)");
		String markerTable = genus + "_markerTable";
		stmt.execute("CREATE TABLE IF NOT EXISTS " + markerTable
				+ " (qacc TEXT, racc TEXT, qstart INTEGER, qend INTEGER, qlen INTEGER, sstart INTEGER, send INTEGER, slen INTEGER, pident DOUBLE, btop TEXT)");
	}

	public void close() throws SQLException {
		c.close();
	}

	public boolean containsAcc(String acc, String table) {
		try {
			stmt.execute("CREATE INDEX IF NOT EXISTS " + table + "Index ON " + table + " (qacc)");
			String sql = "SELECT EXISTS(SELECT 1 FROM " + table + " WHERE qacc='" + acc + "' LIMIT 1)";
			ResultSet rs = stmt.executeQuery(sql);
			if (rs.getInt(1) == 1)
				return true;
		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "Genus: " + table, ex);
		}
		return false;
	}

	public void addAlignmentInfo(String table, AlignmentReceiver receiver) {
		String sql = "SELECT * FROM " + table;
		try {
			ResultSet rs = stmt.executeQuery(sql);
			while (rs.next())
				receiver.addAlignment(new AlignmentInfo(rs));
		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "Genus: " + table, ex);
		}
	}

	public ResultSet getAllDistinctAlignments(String table) {
		String sql = "SELECT DISTINCT * FROM " + table;
		try {
			ResultSet rs = stmt.executeQuery(sql);
			return rs;
		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "Genus: " + table, ex);
		}
		return null;
	}

	public List<AlignmentInfo> getAlignments(String acc, String table) {
		List<AlignmentInfo> alis = new ArrayList<>();
		String sql = "SELECT * FROM " + table + " WHERE qacc='" + acc + "'";
		try {
			ResultSet rs = stmt.executeQuery(sql);
			while (rs.next())
				alis.add(new AlignmentInfo(rs));
		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "Genus: " + table, ex);
		}
		return alis;
	}

	public int getAlignmentCount(String acc, String table) {
		String sql = "SELECT COUNT(*) FROM " + table + " WHERE qacc='" + acc + "'";
		try {
			ResultSet rs = stmt.executeQuery(sql);
			return rs.getInt(1);
		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "Genus: " + table, ex);
		}
		return 0;
	}

	public class AlignmentInfo {

		private String query, ref, btop;
		private double identity;
		private double qstart, qend, qlen;
		private double sstart, send, slen;

		public AlignmentInfo(ResultSet rs) throws SQLException {
			this.query = rs.getString(1);
			this.ref = rs.getString(2);
			this.qstart = rs.getInt(3);
			this.qend = rs.getInt(4);
			this.qlen = rs.getInt(5);
			this.sstart = rs.getInt(6);
			this.send = rs.getInt(7);
			this.slen = rs.getInt(8);
			this.identity = rs.getDouble(9);
			this.btop = rs.getString(10);
		}

		public int getQueryStart() {
			return (int) qstart;
		}

		public int getSubjectStart() {
			return (int) sstart;
		}

		public int getQueryLen() {
			return (int) qlen;
		}

		public int getSubjectLen() {
			return (int) slen;
		}

		public String getQuery() {
			return query;
		}

		public String getRef() {
			return ref;
		}

		public double getIdentity() {
			return identity;
		}

		public double getQueryCoverage() {
			return (Math.abs(qstart - qend) / qlen) * 100.;
		}

		public double getRefCoverage() {
			return (Math.abs(sstart - send) / slen) * 100.;
		}

		public String getBtop() {
			return btop;
		}

		public String toString() {
			return query + " " + ref + " " + getQueryCoverage() + " " + getRefCoverage() + " " + identity;
		}

	}

	public synchronized void addAlignmentTable(String tableName, File src, File tab, boolean dropTable) {
		if (!tab.exists())
			return;
		try {

			long time = System.currentTimeMillis();
			if (dropTable) {
				stmt.execute("DROP TABLE IF EXISTS " + tableName);
				stmt.execute("DROP INDEX IF EXISTS " + tableName + "Index");
			}
			stmt.execute("CREATE TABLE IF NOT EXISTS " + tableName
					+ " (qacc TEXT, racc TEXT, qstart INTEGER, qend INTEGER, qlen INTEGER, sstart INTEGER, send INTEGER, slen INTEGER, pident DOUBLE, btop TEXT)");

			long runtime = (System.currentTimeMillis() - time) / 1000;
			if (src != null) {
				int count = addSelfAlignments(tableName, src);
				System.err.println(
						String.format("SQL>Table " + tableName + ": added %,d items", count) + " (" + runtime + "s)");
			}
			int count = addAlignments(tableName, tab);
			System.err.println(
					String.format("SQL>Table " + tableName + ": added %,d items", count) + " (" + runtime + "s)");
			if (dropTable) {
				time = System.currentTimeMillis();
				stmt.execute("CREATE INDEX " + tableName + "Index ON " + tableName + " (qacc)");
				runtime = (System.currentTimeMillis() - time) / 1000;
				System.err.println("SQL>Table " + tableName + " indexed (" + runtime + "s)");
			}

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private int addSelfAlignments(String tableName, File src) {
		int count = 0;
		try {
			c.setAutoCommit(false);
			try (PreparedStatement insertStmd = c
					.prepareStatement("INSERT INTO " + tableName + " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?);")) {
				System.err.println("SQL> " + tableName + " " + src.getAbsolutePath());
				ArrayList<FastaEntry> tokens = FastaReader.read(src);
				System.err.println("SQL> " + tokens.size() + " proteins adding to " + tableName);
				for (FastaEntry token : tokens) {
					final String acc = token.getName();
					final String seq = token.getSequence();
					final String qacc = acc;
					final String racc = acc;
					final int qstart = 0;
					final int qend = seq.length();
					final int qlen = seq.length();
					final int sstart = 0;
					final int send = seq.length();
					final int slen = seq.length();
					final double pident = 100;
					final String btop = seq.length() + "M";
					insertStmd.setString(1, qacc);
					insertStmd.setString(2, racc);
					insertStmd.setInt(3, qstart);
					insertStmd.setInt(4, qend);
					insertStmd.setInt(5, qlen);
					insertStmd.setInt(6, sstart);
					insertStmd.setInt(7, send);
					insertStmd.setInt(8, slen);
					insertStmd.setDouble(9, pident);
					insertStmd.setString(10, btop);
					insertStmd.execute();
					count++;
				}
			}
			c.commit();
			c.setAutoCommit(true);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return count;
	}

	private int addAlignments(String tableName, File tab) {
		int count = 0;
		try {
			c.setAutoCommit(false);
			try (PreparedStatement insertStmd = c
					.prepareStatement("INSERT INTO " + tableName + " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?);")) {
				try (BufferedReader buf = new BufferedReader(new FileReader(tab))) {
					String line;
					while ((line = buf.readLine()) != null) {
						final String[] tokens = line.split("\t");
						final String qacc = tokens[0];
						final String racc = tokens[1];
						final int qstart = Integer.parseInt(tokens[2]);
						final int qend = Integer.parseInt(tokens[3]);
						final int qlen = Integer.parseInt(tokens[4]);
						final int sstart = Integer.parseInt(tokens[5]);
						final int send = Integer.parseInt(tokens[6]);
						final int slen = Integer.parseInt(tokens[7]);
						final double pident = Double.parseDouble(tokens[8]);
						final String btop = tokens.length < 10 || pident < 99 ? "" : tokens[9];
						insertStmd.setString(1, qacc);
						insertStmd.setString(2, racc);
						insertStmd.setInt(3, qstart);
						insertStmd.setInt(4, qend);
						insertStmd.setInt(5, qlen);
						insertStmd.setInt(6, sstart);
						insertStmd.setInt(7, send);
						insertStmd.setInt(8, slen);
						insertStmd.setDouble(9, pident);
						insertStmd.setString(10, btop);
						insertStmd.execute();
						count++;
					}
				}
			}
			c.commit();
			c.setAutoCommit(true);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return count;
	}

	public Statement getStmt() {
		return stmt;
	}

}
