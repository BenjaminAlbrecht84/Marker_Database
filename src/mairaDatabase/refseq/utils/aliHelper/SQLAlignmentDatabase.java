package mairaDatabase.refseq.utils.aliHelper;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import mairaDatabase.utils.FastaReader;
import mairaDatabase.utils.FastaReader.FastaEntry;
import mairaDatabase.utils.SQLMappingDatabase;

public class SQLAlignmentDatabase {

	private File tmpDir;
	private String databaseFile;
	private Connection c;
	private Statement stmt;
	private String genus;
	private String clusterTable, markerTable, accTable, disjoinTable, gcfTable;

	public SQLAlignmentDatabase(SQLAlignmentDatabase aligmentDatabase) throws ClassNotFoundException, SQLException {
		init(aligmentDatabase.getDatabaseFile(), aligmentDatabase.getGenus(), aligmentDatabase.getTmpDir());
	}
	
	public SQLAlignmentDatabase(SQLAlignmentDatabase aligmentDatabase, SQLMappingDatabase mappingDatabase) throws ClassNotFoundException, SQLException {
		init(aligmentDatabase.getDatabaseFile(), aligmentDatabase.getGenus(), aligmentDatabase.getTmpDir());
		String stmt_attach = "ATTACH '" + mappingDatabase.getDatabaseFile() + "' AS mapping";
		Statement attach_stmt = c.createStatement();
		attach_stmt.executeUpdate(stmt_attach);
	}

	public SQLAlignmentDatabase(String databaseFolder, String genus, File tmpDir)
			throws ClassNotFoundException, SQLException {
		String databaseFile = databaseFolder + File.separator + genus + ".db";
		init(databaseFile, genus, tmpDir);
	}

	private void init(String databaseFile, String genus, File tmpDir) throws ClassNotFoundException, SQLException {
		this.databaseFile = databaseFile;
		this.tmpDir = tmpDir;
		this.genus = genus;
		Class.forName("org.sqlite.JDBC");
		c = DriverManager.getConnection("jdbc:sqlite:" + this.databaseFile);
		stmt = c.createStatement();
//		stmt.execute("PRAGMA temp_store_directory = '" + tmpDir.getAbsolutePath() + "'");

		clusterTable = genus + "_clusterTable";
		stmt.execute("CREATE TABLE IF NOT EXISTS " + clusterTable
				+ " (qacc_id INTEGER, racc_id INTEGER, qstart INTEGER, qend INTEGER, qlen INTEGER, sstart INTEGER, send INTEGER, slen INTEGER, pident DOUBLE, btop TEXT)");
		stmt.execute("CREATE INDEX IF NOT EXISTS " + clusterTable + "_qaccIndex ON " + clusterTable + " (qacc_id)");
		stmt.execute("CREATE INDEX IF NOT EXISTS " + clusterTable + "_raccIndex ON " + clusterTable + " (racc_id)");
		stmt.execute("CREATE INDEX IF NOT EXISTS " + clusterTable + "_pidentIndex ON " + clusterTable + " (pident)");

		markerTable = genus + "_markerTable";
		stmt.execute("CREATE TABLE IF NOT EXISTS " + markerTable
				+ " (qacc_id INTEGER, racc_id INTEGER, qstart INTEGER, qend INTEGER, qlen INTEGER, sstart INTEGER, send INTEGER, slen INTEGER, pident DOUBLE, btop TEXT)");
		stmt.execute("CREATE INDEX IF NOT EXISTS " + markerTable + "_qaccIndex ON " + markerTable + " (qacc_id)");
		stmt.execute("CREATE INDEX IF NOT EXISTS " + markerTable + "_raccIndex ON " + markerTable + " (racc_id)");
		stmt.execute("CREATE INDEX IF NOT EXISTS " + markerTable + "_pidentIndex ON " + markerTable + " (pident)");

		disjoinTable = genus + "_disjoinTable";
		stmt.execute("CREATE TABLE IF NOT EXISTS " + disjoinTable
				+ " (qgcf_id INTEGER, rgcf_id INTEGER, qspecies INTEGER, rspecies INTEGER, disjoin INTEGER)");
		stmt.execute(
				"CREATE INDEX IF NOT EXISTS " + disjoinTable + "_gcfIndex ON " + disjoinTable + " (qgcf_id, rgcf_id)");
		stmt.execute("CREATE INDEX IF NOT EXISTS " + disjoinTable + "_speciesIndex ON " + disjoinTable
				+ " (qspecies, rspecies)");

		accTable = genus + "_accTable";
		stmt.execute("CREATE TABLE IF NOT EXISTS " + accTable + " (id INTEGER PRIMARY KEY, acc TEXT)");
		stmt.execute("CREATE INDEX IF NOT EXISTS " + accTable + "_accIndex ON " + accTable + " (acc)");

		gcfTable = genus + "_gcfTable";
		stmt.execute("CREATE TABLE IF NOT EXISTS " + gcfTable + " (id INTEGER PRIMARY KEY, gcf TEXT)");
		stmt.execute("CREATE INDEX IF NOT EXISTS " + gcfTable + "_gcfIndex ON " + gcfTable + " (gcf)");

		String clusterTmpTable = genus + "_tmp_clusterTable";
		stmt.execute("DROP TABLE IF EXISTS " + clusterTmpTable);
		String markerTmpTable = genus + "_tmp_markerTable";
		stmt.execute("DROP TABLE IF EXISTS " + markerTmpTable);
	}

	public void close() {
		try {
			c.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
	}

	public boolean containsAcc(String acc) {
		try {
			String sql = "SELECT EXISTS(SELECT 1 FROM " + accTable + " WHERE acc ='" + acc + "')";
			ResultSet rs = stmt.executeQuery(sql);
			if (rs.getInt(1) == 1)
				return true;
		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "Genus: " + accTable, ex);
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
		if (containsAcc(acc))
			return getAlignments(getAccessionId(acc), table);
		return Collections.emptyList();
	}
	
	public List<AlignmentInfo> getGenomeProteinAlignments(int gcfId, String table) {
		List<AlignmentInfo> alis = new ArrayList<>();
		String sql = "SELECT a.*, qat.acc, rat.acc "
				+ " FROM gcf2taxid"
				+ " JOIN acc2gcf AS at USING(gcf_id)"
				+ " JOIN " + accTable + " AS qat USING(acc)"
				+ " JOIN " + table + " AS a ON qat.id = a.qacc_id"
				+ " JOIN " + accTable + " AS rat ON rat.id = a.racc_id"
				+ " WHERE gcf_id =" + gcfId + "";
		try {
			ResultSet rs = stmt.executeQuery(sql);
			while (rs.next()) {
				AlignmentInfo aliInfo = new AlignmentInfo(rs);
				if (aliInfo.isOkay(table))
					alis.add(aliInfo);
			}
		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "Genus: " + table, ex);
		}
		return alis;
	}
	
	public List<AlignmentInfo> getAlignments(int accId, String table) {
		List<AlignmentInfo> alis = new ArrayList<>();
		String sql = "SELECT a.*, qat.acc, rat.acc "
				+ " FROM " + table + " AS a "
				+ " JOIN " + accTable + " AS qat ON qat.id = a.qacc_id"
				+ " JOIN " + accTable + " AS rat ON rat.id = a.racc_id"
				+ " WHERE qacc_id ='" + accId + "'";
		try {
			ResultSet rs = stmt.executeQuery(sql);
			while (rs.next()) {
				AlignmentInfo aliInfo = new AlignmentInfo(rs);
				if (aliInfo.isOkay(table))
					alis.add(aliInfo);
			}
		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "Genus: " + table, ex);
		}
		return alis;
	}

	public Integer getAccessionId(String acc) {
		try {
			String accTable = genus + "_accTable";
			ResultSet rs = stmt.executeQuery("SELECT id FROM " + accTable + " WHERE acc = '" + acc + "'");
			if (rs.next())
				return rs.getInt(1);
		} catch (SQLException e) {
			e.printStackTrace();
		}
		return null;
	}

	public String getAccession(int id) {
		try {
			ResultSet rs = stmt.executeQuery("SELECT acc FROM " + accTable + " WHERE id = " + id);
			if (rs.next())
				return rs.getString(1);
		} catch (SQLException e) {
			e.printStackTrace();
		}
		return null;
	}

	public boolean containsGcf(String gcf) {
		try {
			String sql = "SELECT EXISTS(SELECT 1 FROM " + gcfTable + " WHERE gcf ='" + gcf + "')";
			ResultSet rs = stmt.executeQuery(sql);
			if (rs.getInt(1) == 1)
				return true;
		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "Genus: " + accTable, ex);
		}
		return false;
	}

	public Integer getGenomeId(String gcf) {
		try {
			ResultSet rs = stmt.executeQuery("SELECT id FROM " + gcfTable + " WHERE gcf = '" + gcf + "'");
			if (rs.next())
				return rs.getInt(1);
			throw new SQLException("ERROR: Genome unknown: " + gcf);
		} catch (SQLException e) {
			e.printStackTrace();
		}
		return null;
	}

	public String getGenome(int id) {
		try {
			ResultSet rs = stmt.executeQuery("SELECT gcf FROM " + gcfTable + " WHERE id = " + id);
			if (rs.next())
				return rs.getString(1);
		} catch (SQLException e) {
			e.printStackTrace();
		}
		return null;
	}

	public boolean containsDisjoin(int id1, int id2) {
		try {
			String sql = "SELECT EXISTS(SELECT 1 FROM " + disjoinTable + " WHERE qgcf_id =" + id1 + " AND rgcf_id = "
					+ id2 + ")";
			ResultSet rs = stmt.executeQuery(sql);
			if (rs.getInt(1) == 1)
				return true;
		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "Genus: " + accTable, ex);
		}
		return false;
	}
	
	public Integer getAccessionCount() {
		try {
			ResultSet rs = stmt.executeQuery("SELECT COUNT(*) FROM " + accTable);
			if (rs.next())
				return rs.getInt(1);
			throw new SQLException("ERROR: Unknown accession count");
		} catch (SQLException e) {
			e.printStackTrace();
		}
		return null;
	}

	public List<int[]> getProteinGenomeInfo(Integer accId) {

		try {

			//@formatter:off
			String stmt_info = "SELECT "
							+ "    c.species_id, "
							+ "    c.gcf_id "
							+ " FROM " + accTable + " AS a"
							+ " JOIN mapping.acc2gcf AS b"
							+ "   on a.acc = b.acc "
							+ " JOIN mapping.gcf2taxid AS c "
							+ "   ON b.gcf_id = c.gcf_id "
							+ " WHERE a.id = " + accId;
			//@formatter:on 
			PreparedStatement infoStmt = c.prepareStatement(stmt_info);
			ResultSet rs = infoStmt.executeQuery();
			
			List<int[]> proteinGenomeInfo = new ArrayList<>();
			while (rs.next()) {
				int[] info = { rs.getInt(1), rs.getInt(2) };
				proteinGenomeInfo.add(info);
			}
			return proteinGenomeInfo;

		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "Protein: " + accId, ex);
		}
		return null;
	}

	public int getAlignmentCount(String acc, String table) {
		if (containsAcc(acc)) {
			String sql = "SELECT COUNT(*) FROM " + table + " WHERE qacc_id ='" + getAccessionId(acc) + "'";
			try {
				ResultSet rs = stmt.executeQuery(sql);
				return rs.getInt(1);
			} catch (SQLException ex) {
				Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "Genus: " + table, ex);
			}
		}
		return 0;
	}

	public class AlignmentInfo {

		private int queryId, refId;
		private String query, ref, btop;
		private double identity;
		private double qstart, qend, qlen;
		private double sstart, send, slen;

		public AlignmentInfo(ResultSet rs) throws SQLException {
			this.queryId = rs.getInt(1);
			this.refId = rs.getInt(2);
			this.qstart = rs.getInt(3);
			this.qend = rs.getInt(4);
			this.qlen = rs.getInt(5);
			this.sstart = rs.getInt(6);
			this.send = rs.getInt(7);
			this.slen = rs.getInt(8);
			this.identity = rs.getDouble(9);
			this.btop = rs.getString(10);
			this.query = rs.getString(11);
			this.ref = rs.getString(12);
		}

		public boolean isOkay(String table) {
			if (query == null)
				System.err.println("ERROR: Unknown query_id " + queryId + " for table " + table);
			if (ref == null)
				System.err.println("ERROR: Unknown ref_id " + refId + " for table " + table);
			return query != null && ref != null;
		}

		public int getQueryId() {
			return queryId;
		}

		public int getRefId() {
			return refId;
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
			int count = 0;
			if (src != null)
				count += addSelfAlignments(tableName, src);
			count += addAlignments(tableName, tab);
			long runtime = (System.currentTimeMillis() - time) / 1000;
			System.err.println(
					String.format("SQL>Table " + tableName + ": added %,d items", count) + " (" + runtime + "s)");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private int addSelfAlignments(String tableName, File src) {

		int aliCounter = 0;
		try {

			long time = System.currentTimeMillis();

			c.setAutoCommit(false);

			// @formatter:off
			String insert_acc_stmt = "INSERT INTO " + accTable + "(acc) VALUES (?);";
			// @formatter:on

			int accCounter = 0;
			try (PreparedStatement accStmt = c.prepareStatement(insert_acc_stmt)) {
				System.err.println("SQL> " + accTable + " " + src.getAbsolutePath());
				ArrayList<FastaEntry> tokens = FastaReader.read(src);
				System.err.println("SQL> " + tokens.size() + " accessions adding to " + accTable);
				for (FastaEntry token : tokens) {
					final String acc = token.getName();
					if (!containsAcc(acc)) {
						accStmt.setString(1, acc);
						accStmt.execute();
					}
				}
			}

			long runtime = (System.currentTimeMillis() - time) / 1000;
			System.err.println(String.format("SQL>Table " + tableName + ": %,d self-accessions added", accCounter)
					+ " (" + runtime + "s)");

			c.commit();

			// @formatter:off
			String insert_ali_stmt = "INSERT INTO " + tableName + "("
					+ " qacc_id, "
					+ " racc_id, "
					+ " qstart, "
					+ " qend, "
					+ " qlen, "
					+ " sstart, "
					+ " send, "
					+ " slen, "
					+ " pident, "
					+ " btop) "
					+ " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);";
			// @formatter:on

			time = System.currentTimeMillis();
			try (PreparedStatement aliStmt = c.prepareStatement(insert_ali_stmt)) {

				System.err.println("SQL> " + tableName + " " + src.getAbsolutePath());
				ArrayList<FastaEntry> tokens = FastaReader.read(src);
				System.err.println("SQL> " + tokens.size() + " proteins adding to " + tableName);
				for (FastaEntry token : tokens) {
					final String acc = token.getName();
					final String seq = token.getSequence();
					final int qaccId = getAccessionId(acc);
					final int raccId = getAccessionId(acc);
					final int qstart = 0;
					final int qend = seq.length();
					final int qlen = seq.length();
					final int sstart = 0;
					final int send = seq.length();
					final int slen = seq.length();
					final double pident = 100;
					final String btop = seq.length() + "M";
					aliStmt.setInt(1, qaccId);
					aliStmt.setInt(2, raccId);
					aliStmt.setInt(3, qstart);
					aliStmt.setInt(4, qend);
					aliStmt.setInt(5, qlen);
					aliStmt.setInt(6, sstart);
					aliStmt.setInt(7, send);
					aliStmt.setInt(8, slen);
					aliStmt.setDouble(9, pident);
					aliStmt.setString(10, btop);
					aliStmt.execute();
					aliCounter++;
				}
			}

			c.commit();

			runtime = (System.currentTimeMillis() - time) / 1000;
			System.err.println(String.format("SQL>Table " + tableName + ": %,d self-alignments added", aliCounter)
					+ " (" + runtime + "s)");

			c.setAutoCommit(true);

		} catch (Exception e) {
			e.printStackTrace();
		}
		return aliCounter;
	}

	private int addAlignments(String tableName, File tab) {
		int aliCounter = 0;
		try {

			c.setAutoCommit(false);

			// @formatter:off
			String insert_acc_stmt = "INSERT INTO " + accTable + "(acc) VALUES (?);";
			// @formatter:on

			long time = System.currentTimeMillis();
			int accCounter = 0;
			try (PreparedStatement accStmt = c.prepareStatement(insert_acc_stmt)) {
				try (BufferedReader buf = new BufferedReader(new FileReader(tab))) {
					String line;
					while ((line = buf.readLine()) != null) {
						final String[] tokens = line.split("\t");
						final String qacc = tokens[0];
						if (!containsAcc(qacc)) {
							accStmt.setString(1, qacc);
							accStmt.execute();
							accCounter++;
						}
						final String racc = tokens[1];
						if (!containsAcc(racc)) {
							accStmt.setString(1, racc);
							accStmt.execute();
							accCounter++;
						}
					}
				}
			}

			c.commit();

			long runtime = (System.currentTimeMillis() - time) / 1000;
			System.err.println(String.format("SQL>Table " + tableName + ": %,d accessions", accCounter) + " added ("
					+ runtime + "s)");

			// @formatter:off
			String insert_ali_stmt = "INSERT INTO " + tableName + "("
					+ " qacc_id, "
					+ " racc_id, "
					+ " qstart, "
					+ " qend, "
					+ " qlen, "
					+ " sstart, "
					+ " send, "
					+ " slen, "
					+ " pident, "
					+ " btop) "
					+ " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?);";
			// @formatter:on

			time = System.currentTimeMillis();
			try (PreparedStatement aliStmt = c.prepareStatement(insert_ali_stmt)) {
				try (BufferedReader buf = new BufferedReader(new FileReader(tab))) {
					String line;
					while ((line = buf.readLine()) != null) {
						final String[] tokens = line.split("\t");
						final int qaccId = getAccessionId(tokens[0]);
						final int raccId = getAccessionId(tokens[1]);
						final int qstart = Integer.parseInt(tokens[2]);
						final int qend = Integer.parseInt(tokens[3]);
						final int qlen = Integer.parseInt(tokens[4]);
						final int sstart = Integer.parseInt(tokens[5]);
						final int send = Integer.parseInt(tokens[6]);
						final int slen = Integer.parseInt(tokens[7]);
						final double pident = Double.parseDouble(tokens[8]);
						final String btop = tokens.length < 10 || pident < 99 ? "" : tokens[9];
						aliStmt.setInt(1, qaccId);
						aliStmt.setInt(2, raccId);
						aliStmt.setInt(3, qstart);
						aliStmt.setInt(4, qend);
						aliStmt.setInt(5, qlen);
						aliStmt.setInt(6, sstart);
						aliStmt.setInt(7, send);
						aliStmt.setInt(8, slen);
						aliStmt.setDouble(9, pident);
						aliStmt.setString(10, btop);
						aliStmt.execute();
						aliCounter++;
					}
				}
			}
			c.commit();

			runtime = (System.currentTimeMillis() - time) / 1000;
			System.err.println(String.format("SQL>Table " + tableName + ": %,d alignments added", aliCounter) + " ("
					+ runtime + "s)");

			c.setAutoCommit(true);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return aliCounter;
	}

	public Connection getConnection() {
		return c;
	}

	public Statement getStmt() {
		return stmt;
	}

	public File getTmpDir() {
		return tmpDir;
	}

	public String getDatabaseFile() {
		return databaseFile;
	}

	public String getGenus() {
		return genus;
	}

}
