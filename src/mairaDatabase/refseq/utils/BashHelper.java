package mairaDatabase.refseq.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class BashHelper {

	public static void apply(String src, File markerDir, File genusDir) {

		try {

			File scriptDir = new File(src + File.separator + "scripts");
			scriptDir.mkdir();

			File diamondBash = new File(scriptDir.getAbsolutePath() + File.separator + "diamond_makedb.sh");
			File lastBash = new File(scriptDir.getAbsolutePath() + File.separator + "last_makedb.sh");
			File ellaBash = new File(scriptDir.getAbsolutePath() + File.separator + "ella_makedb.sh");
			File[] shFiles = { diamondBash, lastBash, ellaBash };
			for (File f : shFiles) {
				f.delete();
				f.createNewFile();
			}

			BufferedWriter diamondWriter = new BufferedWriter(new FileWriter(diamondBash));
			BufferedWriter lastWriter = new BufferedWriter(new FileWriter(lastBash));
			BufferedWriter ellaWriter = new BufferedWriter(new FileWriter(ellaBash));
			BufferedWriter[] writers = { diamondWriter, lastWriter, ellaWriter };

			try {

				// writing header
				for (BufferedWriter w : writers) {
					w.write("#!/bin/sh\n");
					w.write("exe=$1\n");
					w.write("cpus=$2\n");
				}
				lastWriter.write("size=$3\n");
				ellaWriter.write("size=$3\n");

				// writing marker db commands
				for (File f : markerDir.listFiles((dir, name) -> name.toLowerCase().endsWith(".faa"))) {

					String diamondDb = replaceSuffix(f, "dmnd").replace("genus", "diamond");
					diamondWriter.write("$exe makedb --in " + f.getAbsolutePath() + " -d " + diamondDb + " -p $cpus\n");

					String lastdDb = removeSuffix(f).replace("genus", "last");
					lastWriter.write("mkdir " + lastdDb + "\n");
					lastWriter.write("$exe -P $cpus -s $size -p " + lastdDb + "/last_db " + f.getAbsolutePath() + "\n");

					String ellaDb = replaceSuffix(f, "ella").replace("genus", "ella");
					ellaWriter
							.write("$exe makedb -i " + f.getAbsolutePath() + " -d " + ellaDb + " -p $cpus -s $size\n");

				}

				// writing genus db commands
				for (File dir : Stream.of(genusDir.listFiles()).filter(d -> d.isDirectory())
						.collect(Collectors.toList())) {

					for (File f : Stream.of(dir.listFiles()).filter(f -> f.isFile() && f.getName().endsWith(".faa"))
							.collect(Collectors.toList())) {

						String diamondDb = replaceSuffix(f, "dmnd");
						diamondWriter
								.write("$exe makedb --in " + f.getAbsolutePath() + " -d " + diamondDb + " -p $cpus\n");

						String lastdDb = removeSuffix(f) + "_last_db";
						lastWriter.write("mkdir " + lastdDb + "\n");
						lastWriter.write("$exe -P $cpus -s $size -p " + lastdDb + "/last_db " + f.getAbsolutePath() + "\n");

						String ellaDb = replaceSuffix(f, "ella");
						ellaWriter.write("$exe makedb -i " + f.getAbsolutePath() + " -d " + ellaDb + " -p $cpus\n");

					}

				}

			} finally {
				for (BufferedWriter w : writers)
					w.close();
			}

		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	private static String replaceSuffix(File f, String newExtension) {
		String path = f.getAbsolutePath();
		int i = path.lastIndexOf('.');
		if (i > 0 && i < path.length() - 1)
			path = path.substring(0, i + 1);
		return path + newExtension;
	}

	private static String removeSuffix(File f) {
		String path = f.getAbsolutePath();
		int i = path.lastIndexOf('.');
		if (i > 0 && i < path.length() - 1)
			path = path.substring(0, i);
		return path;
	}

}
