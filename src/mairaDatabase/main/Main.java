package mairaDatabase.main;

import java.io.File;

import mairaDatabase.refseq.RefseqManager;

public class Main {

	public static void main(String[] args) {

		if (args.length != 5 && args.length != 6) {
			System.out.println("SRC|TMP|DB|CORES|MEMORY|Taxa");
			return;
		}

		File src = new File(args[0]);
		File tmp = new File(args[1]);
		String aliFolder = args[2];
		int cores = Integer.parseInt(args[3]);
		int memory = Integer.parseInt(args[4]);
		String genera = args.length > 5 ? args[5] : null;
		new RefseqManager().run(src, tmp, aliFolder, cores, memory, genera.trim().split(","));

	}

}
