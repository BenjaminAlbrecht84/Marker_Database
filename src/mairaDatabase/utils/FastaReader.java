package mairaDatabase.utils;

import java.io.*;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipException;

public class FastaReader {

	public static ArrayList<FastaEntry> read(File fastaFile) {

		ArrayList<FastaEntry> readInfo = new ArrayList<FastaEntry>();
		HashSet<SparseString> readNames = new HashSet<SparseString>();
		try {

			BufferedReader buf;
			try {
				buf = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fastaFile))));
			} catch (ZipException | EOFException e) {
				buf = new BufferedReader(new FileReader(fastaFile));
			}

			String line, name = "";
			boolean readSequence = false;
			StringBuilder seqBuilder = new StringBuilder("");
			while ((line = buf.readLine()) != null) {
				if (line.startsWith(">")) {
					if (seqBuilder.length() != 0 && !name.isEmpty()) {
						FastaEntry entry = new FastaEntry(name, seqBuilder.toString());
						if (checkFASTFile(entry, readNames))
							readInfo.add(entry);
					}
					seqBuilder = new StringBuilder("");
					name = line.substring(1).split(" ")[0];
					readSequence = true;
				} else if (line.startsWith("+")) {
					readSequence = false;
				} else if (readSequence) {
					seqBuilder.append(line);
				}

			}
			if (seqBuilder.length() != 0 && !name.isEmpty()) {
				FastaEntry entry = new FastaEntry(name, seqBuilder.toString());
				if (checkFASTFile(entry, readNames))
					readInfo.add(entry);
			}
			buf.close();

		} catch (Exception e) {
			//e.printStackTrace();
			return null;
		}
		readNames = null;

		return readInfo;

	}
	
	public static class FastaEntry{
		
		private SparseString name, sequence;
		
		public FastaEntry(String name, String sequence) {
			super();
			this.name = new SparseString(name);
			this.sequence = new SparseString(sequence);
		}

		public String getName() {
			return name.toString();
		}
		
		public SparseString getSparseName() {
			return name;
		}

		public String getSequence() {
			return sequence.toString();
		}

		public int getSequenceLength() {
			return sequence.toString().length();
		}
		
	}

	private static boolean checkFASTFile(FastaEntry e, HashSet<SparseString> readNames) {
		if (readNames.contains(e.getName())) {
			System.err.println("WARNING: read " + e.getName() + " occurs multiple times in FASTA file.");
			return false;
		}
		readNames.add(new SparseString(e.getName()));
		return true;
	}

}
