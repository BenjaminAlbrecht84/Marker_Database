package mairaDatabase.refseq.step1_clustering;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;

import jloda.util.Pair;
import mairaDatabase.refseq.RefseqManager;
import mairaDatabase.refseq.utils.aliHelper.SQLAlignmentDatabase;
import mairaDatabase.refseq.utils.aliHelper.SQLAlignmentDatabase.AlignmentInfo;
import mairaDatabase.utils.FastaReader;
import mairaDatabase.utils.FastaReader.FastaEntry;
import mairaDatabase.utils.SQLMappingDatabase;
import mairaDatabase.utils.SparseString;
import mairaDatabase.utils.taxTree.TaxNode;
import mairaDatabase.utils.taxTree.TaxTree;

public class Clustering {

	public enum ClusteringMode {
		MARKER_DB, GENUS_DB
	};

	private int ID_THRESHOLD, COV_THRESHOLD;
	private final static int MIN_PROTEINS_CLUSTER = 1000;

	public void run(String genus, SQLAlignmentDatabase alignmentDatabase, SQLMappingDatabase mappingDatabase,
			TaxTree taxTree, File faaFile, File proteinOutFile, BufferedWriter dominationWriter, int MIN_ID,
			ClusteringMode mode) throws Exception {

		long time = System.currentTimeMillis();
		String table = genus + "_clusterTable";
		COV_THRESHOLD = MIN_ID;
		ID_THRESHOLD = MIN_ID;
		List<FastaEntry> genusProteins = FastaReader.read(faaFile);

		if (mode == ClusteringMode.MARKER_DB)
			genusProteins = genusProteins.stream().filter(e -> hasUniqueGenus(e.getName(), mappingDatabase, taxTree))
					.collect(Collectors.toList());
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
					double qLen = ali.getQueryLen(), slen = ali.getSubjectLen();
					ClusterNode w = acc2node.get(ali.getRef());
					if (w == null)
						continue;
					if (w.isDominated())
						continue;
					boolean isSelfHit = ali.getQuery().equals(ali.getRef());
					if (!isSelfHit && ali.getIdentity() > ID_THRESHOLD && ali.getQueryCoverage() > COV_THRESHOLD) {
						w.setDominatedBy(v, ali);
						selectedNodes--;
					}
				}
			}
			if (selectedNodes < MIN_PROTEINS_CLUSTER)
				break;
		}

		long written = 0;

		try (BufferedWriter proteinsWriter = new BufferedWriter(new FileWriter(proteinOutFile))) {
			for (FastaEntry protein : genusProteins) {
				String acc = protein.getName();
				String seq = protein.getSequence();
				ClusterNode v = acc2node.get(acc);
				if (!v.isDominated()) {
					if (seq.length() > RefseqManager.MIN_LENGTH) {
						proteinsWriter.write(">" + acc + "\n" + seq + "\n");
						written++;
					}
				} else if (mode == ClusteringMode.GENUS_DB) {
					Pair<ClusterNode, AlignmentInfo> pair = v.getDominatedBy();
					ClusterNode dominator = pair.getFirst();
					AlignmentInfo aliInfo = pair.getSecond();
					dominationWriter.write(v.getAcc() + "\t" + dominator.getAcc() + "\t" + aliInfo.getBtop() + "\t"
							+ aliInfo.getQueryStart() + "\t" + aliInfo.getSubjectStart() + "\t"
							+ aliInfo.getSubjectLen() + "\n");
				}
			}
		}

		long runtime = (System.currentTimeMillis() - time) / 1000;
		System.err.println(genus + ": " + written + " proteins reported (" + runtime + "s)");

	}

	private boolean hasUniqueGenus(String acc, SQLMappingDatabase mappingDatabase, TaxTree taxTree) {
		List<Integer> taxids = mappingDatabase.getTaxIdByAcc(acc);
		Set<TaxNode> nodes = taxids.stream().map(id -> taxTree.getNode(id)).filter(Objects::nonNull)
				.collect(Collectors.toSet());
		Set<String> genera = nodes.stream().map(v -> v.getAncestorAtRank("genus")).filter(Objects::nonNull)
				.map(v -> v.getName()).collect(Collectors.toSet());
		return genera.size() == 1;
	}

	public class ClusterNode implements Comparable<ClusterNode> {

		private int outDegree;
		private SparseString acc;
		private Pair<ClusterNode, AlignmentInfo> dominatedBy;

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

		public void setDominatedBy(ClusterNode v, AlignmentInfo aliInfo) {
			dominatedBy = new Pair<>(v, aliInfo);
		}

		public boolean isDominated() {
			return dominatedBy != null;
		}

		public Pair<ClusterNode, AlignmentInfo> getDominatedBy() {
			return dominatedBy;
		}

	}

}
