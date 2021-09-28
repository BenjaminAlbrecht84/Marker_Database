package mairaDatabase.refseq.utils;

import java.io.*;

public class DiamondRunner {

	public static File view(File daa, int cores, String diamondBin) {
		String out = daa.getAbsolutePath().replaceAll("\\.daa", ".tab");
		String cmd = diamondBin + " view -a " + daa.getAbsolutePath()
				+ " -f 6 qseqid sseqid qstart qend qlen sstart send slen pident length -o " + out + " -k 0 -p " + cores;
		executingCommand(cmd);
		return new File(out);
	}

	public static File makedb(File faaFile, int cores, String diamondBin) {
		String in = faaFile.getAbsolutePath();
		String db = faaFile.getAbsolutePath().replaceAll("\\.faa", ".dmnd");
		String cmd = diamondBin + " makedb --in " + in + " -d " + db + " -p " + cores;
		executingCommand(cmd);
		return new File(db);
	}

	public static File blastp(File dbFile, File faaFile, File tmpFile, int identity, double blockSize, int cores,
			String diamondBin) {
		double b = blockSize;
		String out = faaFile.getAbsolutePath().replaceAll("\\.faa", ".tab");
		String db = dbFile.getAbsolutePath(), in = faaFile.getAbsolutePath();
		String cmd = diamondBin + " blastp -d " + db + " -q " + in + " --shape-mask 111101011101111 -b " + b
				+ " -c 1 --id " + identity
				+ " -k 0 -f 6 qseqid sseqid qstart qend qlen sstart send slen pident btop -o " + out + " -p " + cores;
		if (tmpFile != null)
			cmd += " -t " + tmpFile.getAbsolutePath();
		executingCommand(cmd);
		return new File(out);
	}

	private static int executingCommand(String command) {
		try {

			System.err.println("DIAMOND>Executing " + command);

			Process proc = Runtime.getRuntime().exec(command);
			StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");
			StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

			errorGobbler.start();
			outputGobbler.start();
			int exitVal = proc.waitFor();

			return exitVal;

		} catch (Exception e) {
			e.printStackTrace();
		}

		return 1;
	}

	private static class StreamGobbler extends Thread {

		private final InputStream is;
		private final String type;

		StreamGobbler(InputStream is, String type) {
			this.is = is;
			this.type = type;
		}

		public void run() {
			try {
				BufferedReader br = new BufferedReader(new InputStreamReader(is));
				String line;
				if (type.equals("ERROR")) {
					while ((line = br.readLine()) != null)
						System.err.println("DIAMOND_ERR>" + line);
				}
				if (type.equals("OUTPUT")) {
					while ((line = br.readLine()) != null)
						System.err.println("DIAMOND_OUT>" + line);
				}
			} catch (IOException ioe) {
				ioe.printStackTrace();
			}
		}

	}

}
