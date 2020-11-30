package mairaDatabase.refseq.step1_clustering;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import mairaDatabase.refseq.RefseqManager;
import mairaDatabase.refseq.utils.aliHelper.SQLAlignmentDatabase;
import mairaDatabase.refseq.utils.aliHelper.SQLAlignmentDatabase.AlignmentInfo;
import mairaDatabase.utils.FastaReader;
import mairaDatabase.utils.SparseString;
import mairaDatabase.utils.FastaReader.FastaEntry;

public class Clustering {

	public enum ClusteringMode {
		MARKER_DB, GENUS_DB
	};

	private int ID_THRESHOLD, COV_THRESHOLD;
	private final static int MIN_PROTEINS_CLUSTER = 1000;

	public void run(String genus, SQLAlignmentDatabase alignmentDatabase, File faaFile, File proteinOutFile,
			BufferedWriter dominationWriter, int MIN_ID, ClusteringMode mode) {

		String table = genus + "_clusterTable";
		COV_THRESHOLD = MIN_ID;
		ID_THRESHOLD = MIN_ID;
		List<FastaEntry> genusProteins = FastaReader.read(faaFile);
		long time = System.currentTimeMillis();

		List<ClusterNode> clusterNodes = new ArrayList<>(genusProteins.size());
		for (FastaEntry protein : genusProteins) {
			SparseString acc = new SparseString(protein.getName());
			int count = alignmentDatabase.getAlignmentCount(acc.toString(), table);
			clusterNodes.add(new ClusterNode(acc, count));
		}
		Map<String, ClusterNode> acc2node = new HashMap<>(clusterNodes.size());
		clusterNodes.stream().forEach(v -> acc2node.put(v.getAcc(), v));
		Collections.sort(clusterNodes);

		int selectedNodes = clusterNodes.size();
		for (ClusterNode v : clusterNodes) {
			if (!v.isDominated()) {
				List<AlignmentInfo> alis = alignmentDatabase.getAlignments(v.getAcc(), table);
				for (AlignmentInfo ali : alis) {
					double qLen = ali.getQlen(), slen = ali.getSlen();
					ClusterNode w = acc2node.get(ali.getRef());
					if (w == null || w.isDominated() || qLen < RefseqManager.MIN_LENGTH
							|| slen < RefseqManager.MIN_LENGTH)
						continue;
					boolean isSelfHit = ali.getQuery().equals(ali.getRef());
					if (!isSelfHit && ali.getIdentity() > ID_THRESHOLD && ali.getQueryCoverage() > COV_THRESHOLD) {
						w.setDominatedBy(v);
						selectedNodes--;
					}
				}
			}
			if (selectedNodes < MIN_PROTEINS_CLUSTER)
				break;
		}

		long written = 0;
		try {
			BufferedWriter proteinsWriter = new BufferedWriter(new FileWriter(proteinOutFile));
			try {
				for (FastaEntry protein : genusProteins) {
					String acc = protein.getName();
					String seq = protein.getSequence();
					ClusterNode v = acc2node.get(acc);
					if (!v.isDominated()) {
						proteinsWriter.write(">" + acc + "\n" + seq + "\n");
						written++;
					} else if (mode == ClusteringMode.GENUS_DB) {
						for (ClusterNode w : v.getDominatedBy())
							dominationWriter.write(v.getAcc() + "\t" + w.getAcc() + "\n");
					}
				}
			} finally {
				proteinsWriter.close();
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		long runtime = (System.currentTimeMillis() - time) / 1000;
		System.err.println(genus + ": " + written + " proteins reported (" + runtime + "s)");

	}

	public class ClusterNode implements Comparable<ClusterNode> {

		private int outDegree;
		private SparseString acc;
		private Set<ClusterNode> dominatedBy = new HashSet<>();

		public ClusterNode(SparseString acc, int outDegree) {
			this.acc = acc;
			this.outDegree = outDegree;
		}

		@Override
		public int compareTo(ClusterNode v) {
			return Integer.compare(v.outDegree, outDegree);
		}

		public String getAcc() {
			return acc.toString();
		}

		public void setDominatedBy(ClusterNode v) {
			dominatedBy.add(v);
		}

		public boolean isDominated() {
			return !dominatedBy.isEmpty();
		}

		public Set<ClusterNode> getDominatedBy() {
			return dominatedBy;
		}

	}

}
