package refseq.utils.aliHelper;

import utils.FastaReader;
import utils.SQLMappingDatabase;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.sql.*;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

public class SQLAlignmentDatabase {

    private String databaseFile;
    private Connection c;
    private Statement stmt;

    public SQLAlignmentDatabase(String databaseFolder, String genus, File tmpDir) {
        try {
            this.databaseFile = databaseFolder + File.separator + genus + ".db";
            Class.forName("org.sqlite.JDBC");
            c = DriverManager.getConnection("jdbc:sqlite:" + this.databaseFile);
            stmt = c.createStatement();
            stmt.execute("PRAGMA temp_store_directory = '" + tmpDir.getAbsolutePath() + "'");
            String clusterTable = genus + "_clusterTable";
            stmt.execute("CREATE TABLE IF NOT EXISTS " + clusterTable + " (qacc TEXT, racc TEXT, qstart integer, qend integer, qlen integer, sstart integer, send integer, slen integer, pident double, length integer)");
            String markerTable = genus + "_markerTable";
            stmt.execute("CREATE TABLE IF NOT EXISTS " + markerTable + " (qacc TEXT, racc TEXT, qstart integer, qend integer, qlen integer, sstart integer, send integer, slen integer, pident double, length integer)");
        } catch (ClassNotFoundException ex) {
            Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, null, ex);
        } catch (SQLException ex) {
            Logger.getLogger(SQLMappingDatabase.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void close() {
        try {
            c.close();
        } catch (SQLException e) {
            e.printStackTrace();
        }
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

    public ArrayList<AlignmentInfo> getAlignments(String acc, String table) {
        ArrayList<AlignmentInfo> alis = new ArrayList<>();
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

        private String query, ref;
        private double identity;
        private double qstart, qend, qlen;
        private double sstart, send, slen;

        public AlignmentInfo(ResultSet rs) throws SQLException {
            this.query = rs.getString(1);
            this.ref = rs.getString(2);
            this.qstart = rs.getDouble(3);
            this.qend = rs.getDouble(4);
            this.qlen = rs.getDouble(5);
            this.sstart = rs.getDouble(6);
            this.send = rs.getDouble(7);
            this.slen = rs.getDouble(8);
            this.identity = rs.getDouble(9);
        }

        public double getQlen() {
            return qlen;
        }

        public double getSlen() {
            return slen;
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
            stmt.execute("CREATE TABLE IF NOT EXISTS " + tableName + " (qacc TEXT, racc TEXT, qstart integer, qend integer, qlen integer, sstart integer, send integer, slen integer, pident double, length integer)");

            long runtime = (System.currentTimeMillis() - time) / 1000;
            if (src != null) {
                int count = addSelfAlignments(tableName, src);
                System.err.println(String.format("SQL>Table " + tableName + ": added %,d items", count) + " (" + runtime + "s)");
            }
            int count = addAlignments(tableName, tab);
            System.err.println(String.format("SQL>Table " + tableName + ": added %,d items", count) + " (" + runtime + "s)");
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
            try (PreparedStatement insertStmd = c.prepareStatement("INSERT INTO " + tableName + " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?);")) {
                System.err.println("SQL> " + tableName + " " + src.getAbsolutePath());
                ArrayList<Object[]> tokens = FastaReader.read(src);
                System.err.println("SQL> " + tokens.size() + " proteins adding to " + tableName);
                for (Object[] token : tokens) {
                    final String acc = token[0].toString();
                    final String seq = token[1].toString();
                    final String qacc = acc;
                    final String racc = acc;
                    final int qstart = 0;
                    final int qend = seq.length();
                    final int qlen = seq.length();
                    final int sstart = 0;
                    final int send = seq.length();
                    final int slen = seq.length();
                    final double pident = 100;
                    final int length = seq.length();
                    insertStmd.setString(1, qacc);
                    insertStmd.setString(2, racc);
                    insertStmd.setInt(3, qstart);
                    insertStmd.setInt(4, qend);
                    insertStmd.setInt(5, qlen);
                    insertStmd.setInt(6, sstart);
                    insertStmd.setInt(7, send);
                    insertStmd.setInt(8, slen);
                    insertStmd.setDouble(9, pident);
                    insertStmd.setInt(10, length);
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
            try (PreparedStatement insertStmd = c.prepareStatement("INSERT INTO " + tableName + " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?);")) {
                BufferedReader buf = new BufferedReader(new FileReader(tab));
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
                    final int length = Integer.parseInt(tokens[9]);
                    insertStmd.setString(1, qacc);
                    insertStmd.setString(2, racc);
                    insertStmd.setInt(3, qstart);
                    insertStmd.setInt(4, qend);
                    insertStmd.setInt(5, qlen);
                    insertStmd.setInt(6, sstart);
                    insertStmd.setInt(7, send);
                    insertStmd.setInt(8, slen);
                    insertStmd.setDouble(9, pident);
                    insertStmd.setInt(10, length);
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

}
