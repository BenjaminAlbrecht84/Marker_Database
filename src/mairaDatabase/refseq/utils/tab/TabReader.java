package mairaDatabase.refseq.utils.tab;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;

public class TabReader {

	public static ArrayList<TabHit> run(File tabFile, GraphManipulator receiver) throws IOException {

		ArrayList<TabHit> hits = new ArrayList<TabHit>();
		StringBuffer word = new StringBuffer();
		String[] words = new String[10];

		InputStream is = new BufferedInputStream(new FileInputStream(tabFile));
		byte[] buffer = new byte[1024 * 1024];
		int readCharacters;
		int pos = 0;
		while ((readCharacters = is.read(buffer)) != -1) {
			for (int i = 0; i < readCharacters; i++) {
				char c = (char) buffer[i];
				if (c == '\n') {
					words[pos++] = word.toString();
					word = new StringBuffer();
					String[] a = words;
					TabHit h = new TabHit(a[0], a[1], Double.parseDouble(a[2]), Double.parseDouble(a[3]), Double.parseDouble(a[4]),
							Double.parseDouble(a[5]), Double.parseDouble(a[6]), Double.parseDouble(a[7]), Double.parseDouble(a[8]));
					hits.add(h);
					if (receiver != null && hits.size() % 1000000 == 0) {
						receiver.addGraphNodes(hits);
						hits.clear();
					}
					pos = 0;
				} else if (c == '\t') {
					words[pos++] = word.toString();
					word = new StringBuffer();
				} else {
					word.append(c);
				}
			}
		}
		if (receiver != null) {
			receiver.addGraphNodes(hits);
			hits.clear();
		}
		is.close();
		return hits;

	}

	public static class TabHit {

		private String query, ref;
		private double identity;
		private double qstart, qend, qlen;
		private double sstart, send, slen;

		public TabHit(String query, String ref, double qstart, double qend, double qlen, double sstart, double send, double slen, double identity) {
			this.query = query;
			this.ref = ref;
			this.qstart = qstart;
			this.qend = qend;
			this.qlen = qlen;
			this.sstart = sstart;
			this.send = send;
			this.slen = slen;
			this.identity = identity;
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

}
