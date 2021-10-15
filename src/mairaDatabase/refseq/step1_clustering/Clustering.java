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
import java.util.Objects;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;

import jloda.util.Pair;
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
		Map<SparseString, ClusterNode> acc2node = new TreeMap<>();
		clusterNodes.stream().forEach(v -> acc2node.put(new SparseString(v.getAcc()), v));
		Collections.sort(clusterNodes);

		int selectedNodes = clusterNodes.size();
		Map<SparseString, int[]> genomeToInfo = new TreeMap<>();
		for (ClusterNode v : clusterNodes) {
			List<AlignmentInfo> alis = alignmentDatabase.getAlignments(v.getAcc(), table);
			for (AlignmentInfo ali : alis) {
				ClusterNode w = acc2node.get(new SparseString(ali.getRef()));
				boolean isSelfHit = ali.getQuery().equals(ali.getRef());
				if (w != null && !isSelfHit && ali.getIdentity() > ID_THRESHOLD
						&& ali.getRefCoverage() > COV_THRESHOLD) {
					if (!v.isDominated() && !w.isDominated() && !w.isDominator()) {
						w.setDominatedBy(v, ali);
						selectedNodes--;
					}
					if (mode == ClusteringMode.MARKER_DB) {
						Set<Integer> vSpecies = getSpeciesTaxidByACC(v.getAcc(), taxTree, mappingDatabase);
						Set<Integer> wSpecies = getSpeciesTaxidByACC(w.getAcc(), taxTree, mappingDatabase);
						if (vSpecies.stream().anyMatch(id -> !wSpecies.contains(id))) {
							v.setGcfInfo(getGenomeInfo(v.getAcc(), taxTree, mappingDatabase, genomeToInfo));
							w.setGcfInfo(getGenomeInfo(w.getAcc(), taxTree, mappingDatabase, genomeToInfo));
							v.addCrossSpeciesAlignedNode(w);
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
				ClusterNode v = acc2node.get(new SparseString(acc));
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

			Map<Pair<Integer, Integer>, Integer> maxSpeciesOverlaps = new HashMap<>();
			for (Entry<SparseString, int[]> gcfInfo : genomeToInfo.entrySet()) {

				String gcf = gcfInfo.getKey().toString();
				int genomeId = gcfInfo.getValue()[0];
				int s1 = gcfInfo.getValue()[1];

				Map<Pair<Integer, Integer>, Integer> genomeOverlaps = new HashMap<>();
				List<ClusterNode> gcfNodes = clusterNodes.stream()
						.filter(v -> v.getGcfInfo().stream().anyMatch(p -> p[0] == genomeId))
						.collect(Collectors.toList());

				List<ClusterNode> crossSpeciesAlignedNodes = gcfNodes.stream().map(v -> v.getCrossSpeciesAlignedNodes())
						.flatMap(List::stream).distinct().collect(Collectors.toList());
				for (ClusterNode w : crossSpeciesAlignedNodes) {
					w.getGcfInfo().stream().filter(wInfo -> {
						int s2 = wInfo[1];
						return s1 != s2;
					}).forEach(wInfo -> genomeOverlaps.put(new Pair<>(s1, wInfo[1]),
							genomeOverlaps.computeIfAbsent(new Pair<>(s1, wInfo[1]), key -> 0) + 1));
				}

				int total = mappingDatabase.getSizeByGCF(gcf);
				for (Entry<Pair<Integer, Integer>, Integer> e : genomeOverlaps.entrySet()) {
					int overlap = (int) Math.round(((double) e.getValue() / (double) total) * 100.);
					int maxOverlap = maxSpeciesOverlaps.computeIfAbsent(e.getKey(), key -> 0);
					if (overlap > maxOverlap)
						maxSpeciesOverlaps.put(e.getKey(), overlap);
				}

			}

			maxSpeciesOverlaps.entrySet().stream().forEach(e -> {
				try {
					overlapWriter.write(e.getKey().getFirst() + "\t" + e.getKey().getSecond() + "\t" + e.getValue() + "\n");
				} catch (IOException ex) {
					ex.printStackTrace();
				}
			});

		}

		long runtime = (System.currentTimeMillis() - time) / 1000;
		System.err.println(genus + ": " + written + " proteins reported (" + runtime + "s)");

	}

	private List<int[]> getGenomeInfo(String acc, TaxTree taxTree, SQLMappingDatabase mappingDatabase,
			Map<SparseString, int[]> genomeToInfo) {
		List<int[]> genomeInfo = new ArrayList<>();
		for (String genome : getCompleteGenomes(acc, 90, taxTree, mappingDatabase, genomeToInfo)) {
			int[] gcfInfo = genomeToInfo.computeIfAbsent(new SparseString(genome), key -> {
				int[] info = new int[3];
				info[0] = genomeToInfo.size();
				info[1] = getSpeciesTaxidByGCF(genome, taxTree, mappingDatabase);
				info[2] = getGenomeCompleteness(genome, taxTree, mappingDatabase);
				return info;
			});
			int speciesId = gcfInfo[1];
			if (speciesId != -1) {
				int genomeId = gcfInfo[0];
				int[] info = { genomeId, speciesId };
				genomeInfo.add(info);
			}
		}
		return genomeInfo;
	}

	private List<String> getCompleteGenomes(String gcf, int minCompleteness, TaxTree taxTree,
			SQLMappingDatabase mappingDatabase, Map<SparseString, int[]> genomeToInfo) {
		return mappingDatabase.getGCFByAcc(gcf).stream().filter(g -> {
			int c = genomeToInfo.containsKey(new SparseString(gcf)) ? genomeToInfo.get(new SparseString(gcf))[2]
					: getGenomeCompleteness(g, taxTree, mappingDatabase);
			return c > minCompleteness;
		}).collect(Collectors.toList());
	}

	private int getSpeciesTaxidByGCF(String gcf, TaxTree taxTree, SQLMappingDatabase mappingDatabase) {
		TaxNode v = taxTree.getNode(mappingDatabase.getTaxIDByGCF(gcf));
		TaxNode w;
		if (v != null && (w = v.getAncestorAtRank("species")) != null)
			return w.getTaxid();
		return -1;
	}

	private Set<Integer> getSpeciesTaxidByACC(String acc, TaxTree taxTree, SQLMappingDatabase mappingDatabase) {
		return mappingDatabase.getTaxIdByAcc(acc).stream().filter(Objects::nonNull).map(id -> taxTree.getNode(id))
				.filter(Objects::nonNull).map(v -> v.getAncestorAtRank("species")).filter(Objects::nonNull)
				.map(TaxNode::getTaxid).collect(Collectors.toSet());
	}

	private int getGenomeCompleteness(String gcf, TaxTree taxTree, SQLMappingDatabase mappingDatabase) {
		Integer speciesId = getSpeciesTaxidByGCF(gcf, taxTree, mappingDatabase);
		Integer size = mappingDatabase.getSizeByGCF(gcf);
		if (speciesId != -1 && size != null) {
			Integer avgSize = mappingDatabase.getAvgSizeByTaxid(speciesId);
			return (int) (Math.round((double) size / (double) avgSize) * 100);
		}
		return -1;
	}

	public class ClusterNode implements Comparable<ClusterNode> {

		private int outDegree;
		private SparseString acc;
		private Pair<ClusterNode, AlignmentInfo> dominatedBy;
		private boolean isDominator = false;
		private List<int[]> gcfInfo;
		private List<ClusterNode> crossSpeciesAlignedNodes;

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

		public void addCrossSpeciesAlignedNode(ClusterNode w) {
			if (crossSpeciesAlignedNodes == null)
				crossSpeciesAlignedNodes = new ArrayList<>();
			crossSpeciesAlignedNodes.add(w);
		}

		public List<ClusterNode> getCrossSpeciesAlignedNodes() {
			if (crossSpeciesAlignedNodes == null)
				return Collections.emptyList();
			return crossSpeciesAlignedNodes;
		}

		public List<Integer> getGenomeSpecies() {
			if (gcfInfo == null)
				return Collections.emptyList();
			return gcfInfo.stream().map(info -> info[1]).collect(Collectors.toList());
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

		public void setGcfInfo(List<int[]> gcfInfo) {
			this.gcfInfo = gcfInfo;
		}

		public List<int[]> getGcfInfo() {
			if (gcfInfo == null)
				return Collections.emptyList();
			return gcfInfo;
		}

	}

}
