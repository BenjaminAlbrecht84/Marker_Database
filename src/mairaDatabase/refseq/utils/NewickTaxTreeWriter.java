package mairaDatabase.refseq.utils;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.zip.GZIPOutputStream;

import mairaDatabase.utils.SQLMairaDatabase;
import mairaDatabase.utils.taxTree.TaxDump;
import mairaDatabase.utils.taxTree.TaxTree;

public class NewickTaxTreeWriter {

	public static void run(File src, SQLMairaDatabase mairaDatabase) {

		File treeFile = new File(src.getAbsolutePath() + File.separator + "ncbi_taxonomy.tre");
		treeFile.deleteOnExit();
		File gzTreeFile = new File(treeFile.getAbsolutePath() + ".gz");

		TaxTree tree = new TaxDump().parse(src.getAbsolutePath());
		String newick = tree.toNewick();
		
		// writing newick string to file
		try {
			FileWriter fW = new FileWriter(treeFile);
			fW.write(newick);
			fW.close();
			compressGzipFile(treeFile.getAbsolutePath(), gzTreeFile.getAbsolutePath());
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		// writing newick string into database
		mairaDatabase.createTree2Newick("tax", newick);

	}

	private static void compressGzipFile(String file, String gzipFile) throws IOException {
		FileInputStream fis = new FileInputStream(file);
		FileOutputStream fos = new FileOutputStream(gzipFile);
		GZIPOutputStream gzipOS = new GZIPOutputStream(fos);
		byte[] buffer = new byte[1024];
		int len;
		while ((len = fis.read(buffer)) != -1) {
			gzipOS.write(buffer, 0, len);
		}
		// close resources
		gzipOS.close();
		fos.close();
		fis.close();
	}

}
