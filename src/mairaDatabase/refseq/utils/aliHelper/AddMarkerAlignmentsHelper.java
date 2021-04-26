package mairaDatabase.refseq.utils.aliHelper;

import java.io.File;
import java.sql.SQLException;

public class AddMarkerAlignmentsHelper {

	private static File markerOutputFolder;

	public static void main(String[] args) {

		if (args.length != 4) {
			System.out.println("src|alis|tmp|cores");
			System.exit(0);
		}

		try {

			String src = args[0];
			String aliPath = args[1];
			String tmpPath = args[2];

			File tmpDir = new File(tmpPath);
			File aliFolder = new File(aliPath);
			markerOutputFolder = new File(src + File.separatorChar + "genus_marker_proteins");

			System.out.println(">Processing DIAMOND results");
			File diamondTabFile = new File(markerOutputFolder.getAbsolutePath() + File.separator + "allProteins.tab");
			SQLAlignmentDatabase aliDatabase = createAlignmentDatabase(aliFolder.getAbsolutePath(), "Genus", tmpDir);
			aliDatabase.addAlignmentTable("Genus" + "_markerTable", null, diamondTabFile, false);
			
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	private static synchronized SQLAlignmentDatabase createAlignmentDatabase(String aliFolder, String genus,
			File tmpFile) throws ClassNotFoundException, SQLException {
		return new SQLAlignmentDatabase(aliFolder, genus, tmpFile);
	}

}
