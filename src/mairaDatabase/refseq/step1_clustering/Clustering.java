package mairaDatabase.refseq.step1_clustering;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import jloda.util.Pair;
import mairaDatabase.refseq.utils.aliHelper.SQLAlignmentDatabase;
import mairaDatabase.refseq.utils.aliHelper.SQLAlignmentDatabase.AlignmentInfo;
import mairaDatabase.utils.FastaReader;
import mairaDatabase.utils.FastaReader.FastaEntry;
import mairaDatabase.utils.SQLMappingDatabase;
import mairaDatabase.utils.SparseString;
import mairaDatabase.utils.Statistics;
import mairaDatabase.utils.taxTree.TaxTree;

public class Clustering {

	public enum ClusteringMode {
		MARKER_DB, GENUS_DB
	};

	private int ID_THRESHOLD, COV_THRESHOLD;
	private final static int MIN_PROTEINS_CLUSTER = 1000;

	public void run(String genus, SQLAlignmentDatabase alignmentDatabase, SQLMappingDatabase mappingDatabase,
			TaxTree taxTree, File faaFile, File proteinOutFile, BufferedWriter overlapWriter,
			BufferedWriter dominationWriter, int MIN_ID, ClusteringMode mode) throws Exception {

		long time = System.currentTimeMillis();
		String table = genus + "_clusterTable";
		COV_THRESHOLD = MIN_ID;
		ID_THRESHOLD = MIN_ID;
		List<FastaEntry> genusProteins = FastaReader.read(faaFile);

		List<ClusterNode> clusterNodes = new ArrayList<>(genusProteins.size());
		for (FastaEntry protein : genusProteins) {
			SparseString acc = new SparseString(protein.getName());
			int count = alignmentDatabase.getAlignmentCount(acc.toString(), table);
			clusterNodes.add(new ClusterNode(acc, count));
		}
		Map<String, ClusterNode> acc2node = new HashMap<>(clusterNodes.size());
		clusterNodes.stream().forEach(v -> acc2node.put(v.getAcc(), v));
		Collections.sort(clusterNodes);

		Map<String, GCFOverlap> gcf2overlap = new HashMap<>();
		int selectedNodes = clusterNodes.size();
		for (ClusterNode v : clusterNodes) {
			if (!v.isDominated()) {
				List<AlignmentInfo> alis = alignmentDatabase.getAlignments(v.getAcc(), table);
				for (AlignmentInfo ali : alis) {
					ClusterNode w = acc2node.get(ali.getRef());
					if (w == null || w.isDominated() || w.isDominator())
						continue;
					boolean isSelfHit = ali.getQuery().equals(ali.getRef());
					if (!isSelfHit && ali.getIdentity() > ID_THRESHOLD && ali.getRefCoverage() > COV_THRESHOLD) {
						w.setDominatedBy(v, ali);
						selectedNodes--;
						if (mode == ClusteringMode.MARKER_DB) {
							List<String> queryGenomes = mappingDatabase.getGCFByAcc(ali.getQuery());
							List<String> refGenomes = mappingDatabase.getGCFByAcc(ali.getRef());
							queryGenomes.forEach(gcf -> gcf2overlap.computeIfAbsent(gcf, key -> new GCFOverlap(gcf))
									.addOverlap(refGenomes));
						}
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
					proteinsWriter.write(">" + acc + "\n" + seq + "\n");
					written++;
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

		if (mode == ClusteringMode.MARKER_DB) {
			Map<Pair<Integer, Integer>, List<Integer>> speciesOverlapList = new HashMap<>();
			gcf2overlap.values().stream().forEach(o -> {
				for (Entry<Pair<Integer, Integer>, Integer> e : o.getMaxSpeciesOverlaps(taxTree, mappingDatabase)
						.entrySet())
					speciesOverlapList.computeIfAbsent(e.getKey(), key -> new ArrayList<Integer>()).add(e.getValue());
			});
			speciesOverlapList.keySet().stream().forEach(p -> {
				try {
					overlapWriter.write(p.getFirst() + "\t" + p.getSecond() + "\t"
							+ Statistics.getMaxFromIntegers(speciesOverlapList.get(p)) + "\n");
				} catch (IOException e) {
					e.printStackTrace();
				}
			});
		}

		long runtime = (System.currentTimeMillis() - time) / 1000;
		System.err.println(genus + ": " + written + " proteins reported (" + runtime + "s)");

	}

	private Integer getSpeciesTaxidByGCF(String gcf, TaxTree taxTree, SQLMappingDatabase mappingDatabase) {
		return taxTree.getNode(mappingDatabase.getTaxIDByGCF(gcf)).getAncestorAtRank("species").getTaxid();
	}

	private Double getGenomeCompleteness(String gcf, TaxTree taxTree, SQLMappingDatabase mappingDatabase) {
		double size = mappingDatabase.getSizeByGCF(gcf);
		double avgSize = mappingDatabase.getAvgSizeByTaxid(getSpeciesTaxidByGCF(gcf, taxTree, mappingDatabase));
		return size / avgSize;
	}

	public class GCFOverlap {

		private String gcf;
		private Map<String, Integer> overlaps = new HashMap<String, Integer>();

		public GCFOverlap(String gcf) {
			this.gcf = gcf;
		}

		public void addOverlap(List<String> genomes) {
			genomes.stream().forEach(g -> overlaps.put(g, overlaps.computeIfAbsent(g, key -> 0) + 1));
		}

		public Map<Pair<Integer, Integer>, Integer> getMaxSpeciesOverlaps(TaxTree taxTree, SQLMappingDatabase mapping) {
			Map<Pair<Integer, Integer>, Integer> speciesOverlaps = new HashMap<>();
			int total = mapping.getSizeByGCF(gcf);
			int s1 = getSpeciesTaxidByGCF(gcf, taxTree, mapping);
			double c1 = getGenomeCompleteness(gcf, taxTree, mapping);
			for (Entry<String, Integer> e : overlaps.entrySet()) {
				int overlap = (int) Math.round(((double) e.getValue() / (double) total) * 100.);
				int s2 = getSpeciesTaxidByGCF(e.getKey(), taxTree, mapping);
				double c2 = getGenomeCompleteness(e.getKey(), taxTree, mapping);
				if (c1 > 0.9 && c2 > 0.9 && s1 != s2 && overlap > 0) {
					int maxOverlap = speciesOverlaps.computeIfAbsent(new Pair<Integer, Integer>(s1, s2), key -> 0);
					if (overlap > maxOverlap)
						speciesOverlaps.put(new Pair<Integer, Integer>(s1, s2), overlap);
				}
			}
			return speciesOverlaps;
		}

		public String getGcf() {
			return gcf;
		}

	}

	public class ClusterNode implements Comparable<ClusterNode> {

		private int outDegree;
		private SparseString acc;
		private Pair<ClusterNode, AlignmentInfo> dominatedBy;
		private boolean isDominator = false;

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
			v.setDominator(true);
			dominatedBy = new Pair<>(v, aliInfo);
		}

		public boolean isDominated() {
			return dominatedBy != null;
		}

		public Pair<ClusterNode, AlignmentInfo> getDominatedBy() {
			return dominatedBy;
		}

		public boolean isDominator() {
			return isDominator;
		}

		public void setDominator(boolean isDominator) {
			this.isDominator = isDominator;
		}

	}

}
