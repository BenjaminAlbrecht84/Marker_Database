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
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

import mairaDatabase.refseq.utils.aliHelper.SQLAlignmentDatabase;
import mairaDatabase.refseq.utils.aliHelper.SQLAlignmentDatabase.AlignmentInfo;
import mairaDatabase.utils.ResourceLoader;
import mairaDatabase.utils.SQLMappingDatabase;

public class SpeciesDisjoinComparator2 {

	private final Thread.UncaughtExceptionHandler uncaughtExceptionHandler = new Thread.UncaughtExceptionHandler() {
		@Override
		public void uncaughtException(Thread th, Throwable ex) {
			System.out.println("Uncaught exception: " + ex);
		}
	};

	private Map<IntPair, Double> minSpeciesDisjoint = new ConcurrentHashMap<>();
	private ResourceLoader rL = new ResourceLoader();
	private Iterator<String> genomesIterator;

	public void run(File faaFile, String genus, int genusId, SQLAlignmentDatabase alignmentDatabase,
			SQLMappingDatabase mappingDatabase, BufferedWriter speciesDisjoinWriter, int MIN_ID, int MIN_COV,
			int cores) {

		System.out.println(">Assessing species disjoints for " + genus);

		try {

			long time = System.currentTimeMillis();
			String table = genus + "_clusterTable";

			List<String> genusGenomes = mappingDatabase.getGCFByGenus(genusId);
			genomesIterator = genusGenomes.iterator();
			List<Runnable> speciesDisjointThreads = new ArrayList<>();
			for (int i = 0; i < cores; i++) {
				Thread t = new SpeciesDisjointThread(table, mappingDatabase, alignmentDatabase, MIN_ID, MIN_COV);
				t.setUncaughtExceptionHandler(uncaughtExceptionHandler);
				speciesDisjointThreads.add(t);
			}
			rL.runThreads(cores, speciesDisjointThreads, genusGenomes.size());

			for (Entry<IntPair, Double> e : minSpeciesDisjoint.entrySet())
				speciesDisjoinWriter
						.write(e.getKey().getFirst() + "\t" + e.getKey().getSecond() + "\t" + e.getValue() + "\n");

			long runtime = (System.currentTimeMillis() - time) / 1000;
			System.err.println(
					genus + ": " + minSpeciesDisjoint.size() + " species disjoints reported (" + runtime + "s)");
			minSpeciesDisjoint.clear();

		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	private synchronized List<String> nextGenusGenomes() {
		List<String> genomes = new ArrayList<>();
		while (genomesIterator.hasNext() && genomes.size() < 10) {
			genomes.add(genomesIterator.next());
		}
		return genomes;
	}

	private class SpeciesDisjointThread extends Thread {

		private String table;
		private SQLMappingDatabase mappingDatabase;
		private SQLAlignmentDatabase alignmentDatabase;
		private int MIN_ID, MIN_COV;

		public SpeciesDisjointThread(String table, SQLMappingDatabase mappingDatabase,
				SQLAlignmentDatabase alignmentDatabase, int MIN_ID, int MIN_COV)
				throws ClassNotFoundException, SQLException {
			this.table = table;
			this.mappingDatabase = mappingDatabase;
			this.alignmentDatabase = alignmentDatabase;
			this.MIN_ID = MIN_ID;
			this.MIN_COV = MIN_COV;
		}

		@Override
		public void run() {
			List<String> genomeBatch = null;
			while (!(genomeBatch = nextGenusGenomes()).isEmpty()) {

				for (String gcf : genomeBatch) {

					try {

						mappingDatabase = new SQLMappingDatabase(mappingDatabase);
						alignmentDatabase = new SQLAlignmentDatabase(alignmentDatabase, mappingDatabase);

						Integer gcfId = mappingDatabase.getGcfId(gcf);
						if (gcfId == null)
							continue;
						Integer speciesId = mappingDatabase.getSpeciesIdByGCF(gcfId);
						if (speciesId == null)
							continue;

						TreeMap<IntPair, Integer> genomeOverlaps = new TreeMap<>();
						Map<Integer, List<AlignmentInfo>> proteinAlignments = alignmentDatabase
								.getGenomeProteinAlignments(gcfId, table).stream()
								.collect(Collectors.groupingBy(AlignmentInfo::getQueryId));

						for (Entry<Integer, List<AlignmentInfo>> e : proteinAlignments.entrySet()) {
							if (e.getKey() == null)
								continue;
							int accId = e.getKey();
							List<AlignmentInfo> alis = e.getValue().stream()
									.filter(a -> a.getRefCoverage() > MIN_COV && a.getIdentity() > MIN_ID)
									.collect(Collectors.toList());
							AlignedProtein p = new AlignedProtein(accId, gcfId, speciesId);
							for (AlignmentInfo ali : alis)
								p.addDominatingNode(ali.getRefId(), alignmentDatabase);
							for (IntPair genomePair : p.getGenomeOverlapInfo(alignmentDatabase)) {
								int curOverlap = genomeOverlaps.computeIfAbsent(genomePair, key -> 0);
								genomeOverlaps.put(genomePair, curOverlap + 1);
							}
							p.freeMemory();

						}

						for (Entry<IntPair, Integer> e : genomeOverlaps.entrySet()) {
							IntPair genomePair = e.getKey();
							int overlap1 = e.getValue();
							double size1 = mappingDatabase.getSizeByGCF(genomePair.getFirst());
							double disjoint1 = size1 - overlap1;
							double disjointPerc1 = (disjoint1 / size1) * 100.;
							Integer species1 = mappingDatabase.getSpeciesIdByGCF(genomePair.getFirst());
							Integer species2 = mappingDatabase.getSpeciesIdByGCF(genomePair.getSecond());
							if (species1 != null && species2 != null) {
								IntPair speciesPair = new IntPair(species1, species2);
								double curDisjointPerc = minSpeciesDisjoint.computeIfAbsent(speciesPair, key -> 100.);
								if (disjointPerc1 < curDisjointPerc)
									minSpeciesDisjoint.put(speciesPair, disjointPerc1);
							}
						}

						rL.reportProgress(1);

					} catch (ClassNotFoundException | SQLException e1) {
						e1.printStackTrace();
					}

					alignmentDatabase.close();
					mappingDatabase.close();

				}

			}

			rL.countDown();

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

		private SpeciesDisjoinComparator2 getEnclosingInstance() {
			return SpeciesDisjoinComparator2.this;
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
		private int gcfId;
		private int speciesId;
		private Set<Integer> dominatingNodeIds = new HashSet<>();

		public AlignedProtein(int accId, int gcfId, int speciesId) {
			this.accId = accId;
			this.gcfId = gcfId;
			this.speciesId = speciesId;
		}

		public Set<IntPair> getGenomeOverlapInfo(SQLAlignmentDatabase alignmentDatabase) {
			Set<IntPair> genomeOverlapInfo = getGenomeOverlaps(accId, alignmentDatabase);
			for (int domId : dominatingNodeIds)
				genomeOverlapInfo.addAll(getGenomeOverlaps(domId, alignmentDatabase));
			return genomeOverlapInfo;
		}

		public void addDominatingNode(int refId, SQLAlignmentDatabase alignmentDatabase) {
			if (!getGenomeOverlaps(refId, alignmentDatabase).isEmpty())
				dominatingNodeIds.add(refId);
		}

		private Set<IntPair> getGenomeOverlaps(int wId, SQLAlignmentDatabase alignmentDatabase) {
			Set<IntPair> genomeOverlaps = new HashSet<>();
			List<int[]> wInfo = alignmentDatabase.getProteinGenomeInfo(wId);
			int[] wSpecies = wInfo.stream().mapToInt(i -> i[0]).toArray();
			int[] wGenomes = wInfo.stream().mapToInt(i -> i[1]).toArray();
			for (int j = 0; j < wGenomes.length; j++) {
				int g1 = gcfId;
				int g2 = wGenomes[j];
				int s1 = speciesId;
				int s2 = wSpecies[j];
				IntPair info = new IntPair(g1, g2);
				if (g1 != g2 && s1 != s2)
					genomeOverlaps.add(info);
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
