package mairaDatabase.refseq.utils;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import mairaDatabase.utils.FileUtils;

public class Cleaner {

	public static void apply(File src, File database) {

		List<File> files = Stream.of(src.listFiles()).filter(f -> f.isFile()).collect(Collectors.toList());

		// deleting all tab and txt files
		for (File f : files) {
			if (f.getName().endsWith(".tab"))
				f.delete();
			if (f.getName().endsWith(".txt"))
				f.delete();
		}
		database.delete();
		
		// deleting folders
		FileUtils.deleteDirectory(src.getAbsolutePath() + File.separator + "taxdump");
		
	}

}
