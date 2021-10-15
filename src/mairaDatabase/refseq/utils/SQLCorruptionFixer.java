package mairaDatabase.refseq.utils;

import java.io.File;
import java.sql.Connection;
import java.sql.ResultSet;

import mairaDatabase.refseq.utils.aliHelper.SQLAlignmentDatabase;

public class SQLCorruptionFixer {

	public static void run(String genus, File aliDir) {

		SQLAlignmentDatabase db;
		try {

			db = new SQLAlignmentDatabase(aliDir.getAbsolutePath(), genus, aliDir);
			db.getStmt().setFetchSize(1000000);
			Connection c = db.getConnection();
			c.setAutoCommit(false);

			ResultSet rs = db.getStmt().executeQuery("PRAGMA integrity_check");
			if (rs.next() && !rs.getString(1).equals("ok")) {
				ProcessBuilder builder = new ProcessBuilder();
				String dbFile = aliDir.getAbsolutePath() + File.separator + genus + ".db";
				String tmpFile = aliDir.getAbsolutePath() + File.separator + genus + "_tmp.db";
				new File(tmpFile).delete();
				builder.command("sqlite3 " + dbFile + " \".recover\" | sqlite3 " + tmpFile + " | mv -f " + tmpFile + " "
						+ dbFile);
				Process process = builder.start();
				process.waitFor();
			}

		} catch (Exception e) {
			e.printStackTrace();
		}

	}

}
