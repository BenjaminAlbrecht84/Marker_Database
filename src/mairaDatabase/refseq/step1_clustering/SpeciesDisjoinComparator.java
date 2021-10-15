package mairaDatabase.refseq.step1_clustering;

import java.io.BufferedWriter;
import java.io.File;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Objects;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;

import jloda.util.FastAFileIterator;
import jloda.util.Pair;
import mairaDatabase.refseq.utils.aliHelper.SQLAlignmentDatabase;
import mairaDatabase.refseq.utils.aliHelper.SQLAlignmentDatabase.AlignmentInfo;
import mairaDatabase.utils.ResourceLoader;
import mairaDatabase.utils.SQLMappingDatabase;

public class SpeciesDisjoinComparator {

	private final Thread.UncaughtExceptionHandler uncaughtExceptionHandler = new Thread.UncaughtExceptionHandler() {
		@Override
		public void uncaughtException(Thread th, Throwable ex) {
			System.out.println("Uncaught exception: " + ex);
		}
	};

	private Map<IntPair, Double> minSpeciesDisjoint = new TreeMap<>();
	private Map<IntPair, Integer> genomeOverlaps = new TreeMap<>();
	private ResourceLoader rL = new ResourceLoader();
	private FastAFileIterator proteinIterator;
	private Iterator<Entry<IntPair, Integer>> genomeOverlapsIterator;
	
	private Set<Integer> sourceSpeciesIds, targetSpeciesIds;

	public void run(File faaFile, String genus, int genusId, SQLAlignmentDatabase alignmentDatabase,
			SQLMappingDatabase mappingDatabase, BufferedWriter speciesDisjoinWriter, int MIN_ID, int MIN_COV,
			int cores) {

		System.out.println(">Assessing species disjoints for " + genus);

		try {

			long time = System.currentTimeMillis();
			String table = genus + "_clusterTable";

			List<Integer> speciesIds = mappingDatabase.getSpeciesIdByGenusId(genusId);
			List<Set<Integer>> speciesBatches = new ArrayList<>();
			int n = Math.floorDiv(mappingDatabase.maxSpeciesCount(), 2);
			for (int i = 0; i < speciesIds.size(); i += n)
				speciesBatches.add(new HashSet<>(speciesIds.subList(i, Math.min(i + n, speciesIds.size()))));
			
			int nBatches = speciesBatches.size() * speciesBatches.size();
			int batchCounter = 1;
			for (Set<Integer> b1 : speciesBatches) {
				for (Set<Integer> b2 : speciesBatches) {
					
					System.out.println("Batch " + (batchCounter++) +"/"+nBatches);
					
					sourceSpeciesIds = b1;
					targetSpeciesIds = b2;
					
					proteinIterator = new FastAFileIterator(faaFile.getAbsolutePath());
					List<Runnable> genomeOverlapThreads = new ArrayList<>();
					for (int i = 0; i < cores; i++) {
						Thread t = new GenomeOverlapThread(table, mappingDatabase, alignmentDatabase, MIN_ID, MIN_COV);
						t.setUncaughtExceptionHandler(uncaughtExceptionHandler);
						genomeOverlapThreads.add(t);
					}
					rL.runThreads(cores, genomeOverlapThreads, alignmentDatabase.getAccessionCount());

					genomeOverlapsIterator = genomeOverlaps.entrySet().stream().filter(e -> e.getValue() > 0).iterator();
					List<Runnable> speciesDisjointThreads = new ArrayList<>();
					for (int i = 0; i < cores; i++) {
						Thread t = new SpeciesDisjointThread(mappingDatabase);
						t.setUncaughtExceptionHandler(uncaughtExceptionHandler);
						speciesDisjointThreads.add(t);
					}
					rL.runThreads(cores, speciesDisjointThreads, genomeOverlaps.size());
					genomeOverlaps.clear();;

					for (Entry<IntPair, Double> e : minSpeciesDisjoint.entrySet())
						speciesDisjoinWriter
								.write(e.getKey().getFirst() + "\t" + e.getKey().getSecond() + "\t" + e.getValue() + "\n");

					long runtime = (System.currentTimeMillis() - time) / 1000;
					System.err.println(
							genus + ": " + minSpeciesDisjoint.size() + " species disjoints reported (" + runtime + "s)");
					minSpeciesDisjoint.clear();;
					
				}
			}

		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	private synchronized List<Entry<IntPair, Integer>> nextGenomeOverlapEntry() {
		List<Entry<IntPair, Integer>> overlaps = new ArrayList<>();
		while (genomeOverlapsIterator.hasNext() && overlaps.size() < 1000) {
			overlaps.add(genomeOverlapsIterator.next());
		}
		return overlaps;
	}

	private synchronized void submitSpeciesDisjoints(TreeMap<IntPair, Double> localMinSpeciesDisjoint) {
		for (Entry<IntPair, Double> e : localMinSpeciesDisjoint.entrySet()) {
			IntPair speciesPair = e.getKey();
			double curDisjointPerc = minSpeciesDisjoint.computeIfAbsent(speciesPair, key -> 100.);
			double disjointPerc1 = e.getValue();
			if (disjointPerc1 < curDisjointPerc)
				minSpeciesDisjoint.put(speciesPair, disjointPerc1);
		}
	}

	private class SpeciesDisjointThread extends Thread {

		private SQLMappingDatabase mappingDatabase;

		public SpeciesDisjointThread(SQLMappingDatabase mappingDatabase) {
			this.mappingDatabase = new SQLMappingDatabase(mappingDatabase);
		}

		@Override
		public void run() {
			List<Entry<IntPair, Integer>> overlaps = null;
			while (!(overlaps = nextGenomeOverlapEntry()).isEmpty()) {
				TreeMap<IntPair, Double> localMinSpeciesDisjoint = new TreeMap<>();
				for (Entry<IntPair, Integer> e : overlaps) {
					IntPair genomePair = e.getKey();
					int overlap1 = e.getValue();
					double size1 = mappingDatabase.getSizeByGCF(genomePair.getFirst());
					double disjoint1 = size1 - overlap1;
					double disjointPerc1 = (disjoint1 / size1) * 100.;
					Integer species1 = mappingDatabase.getSpeciesIdByGCF(genomePair.getFirst());
					Integer species2 = mappingDatabase.getSpeciesIdByGCF(genomePair.getSecond());
					if (species1 != null && species2 != null) {
						IntPair speciesPair = new IntPair(species1, species2);
						double curDisjointPerc = localMinSpeciesDisjoint.computeIfAbsent(speciesPair, key -> 100.);
						if (disjointPerc1 < curDisjointPerc)
							localMinSpeciesDisjoint.put(speciesPair, disjointPerc1);
					}
					rL.reportProgress(1);
				}
				submitSpeciesDisjoints(localMinSpeciesDisjoint);
			}
			mappingDatabase.close();
			rL.countDown();
		}

	}

	private synchronized List<String> nextGenusProtein() {
		List<String> proteins = new ArrayList<>();
		while (proteinIterator.hasNext() && proteins.size() < 1000) {
			Pair<String, String> entry = proteinIterator.next();
			String acc = entry.getFirst().substring(1);
			proteins.add(acc);
		}
		return proteins;
	}

	private synchronized void submitGenomeOverlaps(TreeMap<IntPair, Integer> localGenomeOverlaps) {
		for (IntPair genomePair : localGenomeOverlaps.keySet()) {
			int curOverlap = genomeOverlaps.computeIfAbsent(genomePair, key -> 0);
			genomeOverlaps.put(genomePair, curOverlap + localGenomeOverlaps.get(genomePair));
		}
	}

	private class GenomeOverlapThread extends Thread {

		private String table;
		private SQLAlignmentDatabase alignmentDatabase;
		private int MIN_ID, MIN_COV;

		public GenomeOverlapThread(String table, SQLMappingDatabase mappingDatabase,
				SQLAlignmentDatabase alignmentDatabase, int MIN_ID, int MIN_COV)
				throws ClassNotFoundException, SQLException {
			this.table = table;
			this.alignmentDatabase = new SQLAlignmentDatabase(alignmentDatabase, mappingDatabase);
			this.MIN_ID = MIN_ID;
			this.MIN_COV = MIN_COV;
		}

		@Override
		public void run() {
			List<String> proteins = null;
			while (!(proteins = nextGenusProtein()).isEmpty()) {
				TreeMap<IntPair, Integer> localGenomeOverlaps = new TreeMap<>();
				for (String acc : proteins) {
					Integer accId = alignmentDatabase.getAccessionId(acc);
					if (accId == null)
						continue;
					if (!hasSourceSpeciesIds(accId))
						continue;
					AlignedProtein p = new AlignedProtein(accId);
					List<AlignmentInfo> alis = alignmentDatabase.getAlignments(accId, table).stream()
							.filter(a -> a.getRefCoverage() > MIN_COV && a.getIdentity() > MIN_ID)
							.collect(Collectors.toList());
					for (AlignmentInfo ali : alis)
						p.addDominatingNode(ali.getRefId(), alignmentDatabase);
					for (IntPair genomePair : p.getGenomeOverlapInfo(alignmentDatabase)) {
						int curOverlap = localGenomeOverlaps.computeIfAbsent(genomePair, key -> 0);
						localGenomeOverlaps.put(genomePair, curOverlap + 1);
					}
					p.freeMemory();
					rL.reportProgress(1);
				}
				submitGenomeOverlaps(localGenomeOverlaps);
			}
			alignmentDatabase.close();
			rL.countDown();
		}
		
		private boolean hasSourceSpeciesIds(int accId) {
			return alignmentDatabase.getProteinGenomeInfo(accId).stream().mapToInt(i -> i[0])
					.anyMatch(id -> sourceSpeciesIds.contains(id));
		}

	}

	public class IntPair implements Comparable<IntPair> {

		private int first, second;

		public IntPair(int first, int second) {
			this.first = first;
			this.second = second;
		}

		public int getFirst() {
			return first;
		}

		public int getSecond() {
			return second;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + getEnclosingInstance().hashCode();
			result = prime * result + Objects.hash(first, second);
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj) {
				return true;
			}
			if (obj == null) {
				return false;
			}
			if (getClass() != obj.getClass()) {
				return false;
			}
			IntPair other = (IntPair) obj;
			if (!getEnclosingInstance().equals(other.getEnclosingInstance())) {
				return false;
			}
			return first == other.first && second == other.second;
		}

		private SpeciesDisjoinComparator getEnclosingInstance() {
			return SpeciesDisjoinComparator.this;
		}

		@Override
		public int compareTo(IntPair p) {
			if (first != p.first)
				return Integer.compare(first, p.first);
			return Integer.compare(second, p.second);
		}

	}

	public class AlignedProtein {

		private int accId;
		private Set<Integer> dominatingNodeIds = new HashSet<>();

		public AlignedProtein(int accId) {
			this.accId = accId;
		}

		public Set<IntPair> getGenomeOverlapInfo(SQLAlignmentDatabase alignmentDatabase) {
			List<int[]> vInfo = alignmentDatabase.getProteinGenomeInfo(accId);
			Set<IntPair> genomeOverlapInfo = getGenomeOverlaps(vInfo, accId, alignmentDatabase);
			for (int domId : dominatingNodeIds)
				genomeOverlapInfo.addAll(getGenomeOverlaps(vInfo, domId, alignmentDatabase));
			return genomeOverlapInfo;
		}

		public void addDominatingNode(int refId, SQLAlignmentDatabase alignmentDatabase) {
			List<int[]> vInfo = alignmentDatabase.getProteinGenomeInfo(accId);
			if (!getGenomeOverlaps(vInfo, refId, alignmentDatabase).isEmpty())
				dominatingNodeIds.add(refId);
		}

		private Set<IntPair> getGenomeOverlaps(List<int[]> vInfo, int wId, SQLAlignmentDatabase alignmentDatabase) {
			Set<IntPair> genomeOverlaps = new HashSet<>();
			List<int[]> wInfo = alignmentDatabase.getProteinGenomeInfo(wId);
			int[] vSpecies = vInfo.stream().mapToInt(i -> i[0]).toArray();
			int[] vGenomes = vInfo.stream().mapToInt(i -> i[1]).toArray();
			int[] wSpecies = wInfo.stream().mapToInt(i -> i[0]).toArray();
			int[] wGenomes = wInfo.stream().mapToInt(i -> i[1]).toArray();
			for (int i = 0; i < vGenomes.length; i++) {
				for (int j = 0; j < wGenomes.length; j++) {
					int g1 = vGenomes[i];
					int g2 = wGenomes[j];
					int s1 = vSpecies[i];
					int s2 = wSpecies[j];
					IntPair info = new IntPair(g1, g2);
					if (g1 != g2 && s1 != s2 && sourceSpeciesIds.contains(s1) && targetSpeciesIds.contains(s2))
						genomeOverlaps.add(info);
				}
			}
			return genomeOverlaps;
		}

		public int getAccId() {
			return accId;
		}

		public void freeMemory() {
			dominatingNodeIds = null;
		}

	}

}
