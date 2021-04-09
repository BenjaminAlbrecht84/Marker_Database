package mairaDatabase.refseq.step2_filtering;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

import mairaDatabase.refseq.RefseqManager;
import mairaDatabase.refseq.utils.aliHelper.SQLAlignmentDatabase;
import mairaDatabase.utils.FastaReader;
import mairaDatabase.utils.SQLMappingDatabase;
import mairaDatabase.utils.SparseString;
import mairaDatabase.utils.Statistics;
import mairaDatabase.utils.FastaReader.FastaEntry;
import mairaDatabase.utils.taxTree.TaxTree;

public class Filtering {

	private int ID_THRESHOLD, COV_THRESHOLD;
	private SQLMappingDatabase mappingDatabase;
	private SQLAlignmentDatabase alignmentDatabase;
	private String table;

	public void run(File faaFile, String genus, BufferedWriter factorWriter, BufferedWriter markerWriter,
			TaxTree taxTree, SQLMappingDatabase mappingDatabase, SQLAlignmentDatabase alignmentDatabase,
			int MAX_PROTEINS_PER_GCF, int MIN_ID) {

		this.table = genus + "_clusterTable";
		this.COV_THRESHOLD = MIN_ID;
		this.ID_THRESHOLD = MIN_ID;
		this.mappingDatabase = mappingDatabase;
		this.alignmentDatabase = alignmentDatabase;
		long time = System.currentTimeMillis();

		List<FastaEntry> markerProteins = FastaReader.read(faaFile).stream()
				.filter(p -> p.getSequenceLength() > RefseqManager.MIN_LENGTH).collect(Collectors.toList());
		List<MarkerNode> markerNodes = new ArrayList<>(markerProteins.size());
		Map<String, MarkerNode> acc2node = new HashMap<>(markerProteins.size());
		for (FastaEntry protein : markerProteins) {
			SparseString acc = protein.getSparseName();
			int count = alignmentDatabase.getAlignmentCount(acc.toString(), table);
			MarkerNode v = new MarkerNode(acc, count);
			markerNodes.add(v);
			acc2node.put(acc.toString(), v);
		}
		Collections.sort(markerNodes);

		int selectedNodes = markerNodes.size();
		Map<String, Integer> gcf2Counts = new HashMap<>();
		try {
			for (MarkerNode v : markerNodes) {
				Map<String, Integer> localCounts = new HashMap<>();
				for (String gcf : getCoveredGenomes(v, acc2node)) {
					gcf2Counts.putIfAbsent(gcf, 0);
					localCounts.putIfAbsent(gcf, 0);
					localCounts.put(gcf, localCounts.get(gcf) + 1);
					if (gcf2Counts.get(gcf) < MAX_PROTEINS_PER_GCF)
						v.setSelected(true);
				}
				if (v.isSelected()) {
					for (String gcf : localCounts.keySet()) {
						gcf2Counts.putIfAbsent(gcf, 0);
						gcf2Counts.put(gcf, gcf2Counts.get(gcf) + localCounts.get(gcf));
					}
				} else
					selectedNodes--;
			}
		} catch (Exception e) {
			e.printStackTrace();
		}

		// writing out marker proteins
		for (FastaEntry protein : markerProteins) {
			String acc = protein.getName();
			if (acc2node.get(acc).isSelected()) {
				try {
					markerWriter.write(">" + acc + "\n" + protein.getSequence() + "\n");
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}

		// writing out marker protein weights
		for (MarkerNode v : acc2node.values()) {
			if (v.isSelected()) {
				List<Double> gcfFactors = new ArrayList<>();
				for (String gcf : getCoveredGenomes(v, acc2node))
					gcfFactors.add(1. / (double) gcf2Counts.get(gcf));
				try {
					factorWriter.write(v.getAcc() + "\t" + Statistics.getMean(gcfFactors) + "\n");
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		}

		long runtime = (System.currentTimeMillis() - time) / 1000;
		System.err.println(genus + ": " + selectedNodes + "/" + markerProteins.size() + " marker proteins selected ("
				+ runtime + "s)");

	}

	private Set<String> getCoveredGenomes(MarkerNode v, Map<String, MarkerNode> acc2node) {
		List<SQLAlignmentDatabase.AlignmentInfo> alis = alignmentDatabase.getAlignments(v.getAcc(), table);
		Set<String> coveredGenomes = new HashSet<>();
		for (SQLAlignmentDatabase.AlignmentInfo ali : alis) {
			if (ali.getIdentity() > ID_THRESHOLD && ali.getQueryCoverage() > COV_THRESHOLD) {
				for (String gcf : mappingDatabase.getGCFByAcc(ali.getRef()))
					coveredGenomes.add(gcf);
			}
		}
		return coveredGenomes;
	}

	public class MarkerNode implements Comparable<MarkerNode> {

		private int outDegree;
		private SparseString acc;
		private boolean selected = false;

		public MarkerNode(SparseString acc, int outDegree) {
			this.acc = acc;
			this.outDegree = outDegree;
		}

		@Override
		public int compareTo(MarkerNode v) {
			return Integer.compare(v.outDegree, outDegree);
		}

		public String getAcc() {
			return acc.toString();
		}

		public void setSelected(boolean selected) {
			this.selected = selected;
		}

		public boolean isSelected() {
			return selected;
		}

	}

}
