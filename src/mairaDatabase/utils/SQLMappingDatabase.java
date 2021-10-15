package mairaDatabase.utils;

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
//			stmt.execute("PRAGMA page_size = 10485760");
//			stmt.execute("PRAGMA temp_store_directory = '" + tmpDir.getAbsolutePath() + "'");
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

	public List<String> getAvgGCFsByAcc(String acc) {
		String sql = "SELECT gcf a FROM acc2gcf" 
				+ " JOIN gcf2taxid b USING(gcf_id)"
				+ " JOIN species2size c ON (c.taxid=b.species_id)" 
				+ " WHERE acc='" + acc + "'"
				+ " AND b.size >= c.avg_size";
		try {
			ResultSet rs = stmt.executeQuery(sql);
			List<String> arr = toStringArrayList(rs);
			return arr;
		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "Accession: " + acc);
		}
		return null;
	}

	public List<String> getGCFByAcc(String acc) {
		String sql = "SELECT gcf FROM acc2gcf"
				+ " JOIN gcf2taxid USING(gcf_id)"
				+ " WHERE acc='" + acc + "'";
		try {
			ResultSet rs = stmt.executeQuery(sql);
			return toStringArrayList(rs);
		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "Accession: " + acc);
		}
		return null;
	}

	public List<String> getAccByGCF(String gcf) {
		String sql = "SELECT acc FROM acc2gcf "
				+ " JOIN gcf2taxid USING(gcf_id)"
				+ " WHERE gcf='" + gcf + "'";
		try {
			ResultSet rs = stmt.executeQuery(sql);
			return toStringArrayList(rs);
		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "Genome: " + gcf);
		}
		return null;
	}

	public List<String> getGCFByGenus(int genusId) {
		String sql = "SELECT gcf FROM gcf2taxid WHERE genus_id=" + genusId;
		try {
			ResultSet rs = stmt.executeQuery(sql);
			return toStringArrayList(rs);
		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "GenusId: " + genusId);
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
	
	public Integer getGcfId(String gcf) {
		String sql = "SELECT gcf_id FROM gcf2taxid WHERE gcf='" + gcf + "'";
		try {
			ResultSet rs = stmt.executeQuery(sql);
			return rs.getInt(1);
		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "GCF: " + gcf);
		}
		return null;
	}
	
	public String getGcf(int gcfId) {
		String sql = "SELECT gcf FROM gcf2taxid WHERE gcf_id='" + gcfId + "'";
		try {
			ResultSet rs = stmt.executeQuery(sql);
			if(rs.next())
				return rs.getString(1);
		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "GCF: " + gcfId);
		}
		return null;
	}
	
	public Integer getSpeciesIdByGCF(int gcfId) {
		String sql = "SELECT species_id FROM gcf2taxid WHERE gcf_id='" + gcfId + "'";
		try {
			ResultSet rs = stmt.executeQuery(sql);
			if(rs.next())
				return rs.getInt(1);
		} catch (SQLException ex) {
			ex.printStackTrace();
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "GCF: " + gcfId);
		}
		return null;
	}
	
	public Integer maxSpeciesCount() {
		String sql = "SELECT COUNT(DISTINCT species_id) AS count FROM gcf2taxid GROUP BY genus_id ORDER BY count DESC LIMIT 1";
		try {
			ResultSet rs = stmt.executeQuery(sql);
			if(rs.next())
				return rs.getInt(1);
		} catch (SQLException ex) {
			ex.printStackTrace();
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "");
		}
		return null;
	}
	
	public List<Integer> getSpeciesIdByGenusId(int genusId) {
		String sql = "SELECT DISTINCT species_id FROM gcf2taxid WHERE genus_id='" + genusId + "'";
		try {
			ResultSet rs = stmt.executeQuery(sql);
			return toIntArrayList(rs);
		} catch (SQLException ex) {
			ex.printStackTrace();
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "GCF: " + genusId);
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
	
	public Integer getSizeByGCF(int gcfId) {
		String sql = "SELECT size FROM gcf2taxid WHERE gcf_id='" + gcfId + "'";
		try {
			ResultSet rs = stmt.executeQuery(sql);
			return rs.getInt(1);
		} catch (SQLException ex) {
			Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, "GCF_id: " + gcfId);
		}
		return null;
	}

	private List<String> toStringArrayList(ResultSet rs) throws SQLException {
		List<String> result = new ArrayList<>();
		while (rs.next())
			result.add(rs.getString(1));
		return result;
	}

	private List<Object[]> toGenomeInfoList(ResultSet rs) throws SQLException {
		List<Object[]> result = new ArrayList<>();
		while (rs.next()) {
			Object[] o = { rs.getString(1), rs.getInt(2), rs.getInt(3) };
			result.add(o);
		}
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
			stmt.execute("CREATE TABLE IF NOT EXISTS species2size ("
					+ "		taxid INTEGER PRIMARY KEY, "
					+ "		avg_size INTEGER)");
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
			stmt.execute(
					"CREATE TABLE IF NOT EXISTS gcf2taxid ("
					+ "gcf_id INTEGER PRIMARY KEY AUTOINCREMENT, "
					+ "gcf TEXT, "
					+ "size INTEGER, "
					+ "taxid INTEGER, "
					+ "species_id INTEGER, "
					+ "genus_id INTEGER)");
			stmt.execute("DELETE FROM gcf2taxid");
			c.setAutoCommit(false);
			int count = 0;
			rL.setMaxProgress(gcf2size2taxidFile.length());
			try (PreparedStatement insertStmd = c.prepareStatement("INSERT INTO gcf2taxid VALUES (?, ?, ?, ?, ?, ?);");
					BufferedReader buf = new BufferedReader(new FileReader(gcf2size2taxidFile));) {
				String line;
				while ((line = buf.readLine()) != null) {
					final String[] tokens = line.split("\t");
					final String gcf = tokens[0];
					final int size = Integer.parseInt(tokens[1]);
					final int taxid = Integer.parseInt(tokens[2]);
					final int speciesId = Integer.parseInt(tokens[3]);
					final int genusId = Integer.parseInt(tokens[4]);
					insertStmd.setObject(1, null);
					insertStmd.setString(2, gcf);
					insertStmd.setInt(3, size);
					insertStmd.setInt(4, taxid);
					insertStmd.setInt(5, speciesId);
					insertStmd.setInt(6, genusId);
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
			stmt.executeUpdate("DROP INDEX IF EXISTS gcf2taxid_gcfIndex");
			stmt.executeUpdate("CREATE INDEX gcf2taxid_gcfIndex ON gcf2taxid (gcf)");
			stmt.executeUpdate("DROP INDEX IF EXISTS gcf2taxid_speciesIdIndex");
			stmt.executeUpdate("CREATE INDEX gcf2taxid_taxidIndex ON gcf2taxid (species_id)");
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
			stmt.execute("CREATE TABLE IF NOT EXISTS acc2gcf ("
					+ "acc_id INTEGER PRIMARY KEY AUTOINCREMENT, "
					+ "acc TEXT, "
					+ "gcf_id INTEGER, "
					+ "taxid integer)");
			stmt.execute("DELETE FROM acc2gcf");
			c.setAutoCommit(false);
			int count = 0;
			rL.setMaxProgress(acc2gcf2taxidFile.length());
			try (PreparedStatement insertStmd = c.prepareStatement("INSERT INTO acc2gcf VALUES (?, ?, ?, ?);");
					BufferedReader buf = new BufferedReader(new FileReader(acc2gcf2taxidFile));) {
				String line;
				while ((line = buf.readLine()) != null) {
					final String[] tokens = line.split("\t");
					final String accession = tokens[0];
					final int gcf_id = getGcfId(tokens[1]);
					final int taxid = Integer.parseInt(tokens[2]);
					insertStmd.setObject(1, null);
					insertStmd.setString(2, accession);
					insertStmd.setInt(3, gcf_id);
					insertStmd.setInt(4, taxid);
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
			stmt.execute("DROP INDEX IF EXISTS acc2gcf_accIndex");
			stmt.execute("CREATE INDEX acc2gcf_accIndex ON acc2gcf (acc)");
			stmt.execute("DROP INDEX IF EXISTS acc2gcf_gcfIndex");
			stmt.execute("CREATE INDEX acc2gcf_gcfIndex ON acc2gcf (gcf_id)");
			stmt.execute("DROP INDEX IF EXISTS acc2gcf_taxidIndex");
			stmt.execute("CREATE INDEX acc2gcf_taxidIndex ON acc2gcf (taxid)");
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
