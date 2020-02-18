package refseq.step1_filtering;

import refseq.utils.DiamondRunner;
import refseq.utils.aliHelper.SQLAlignmentDatabase;
import utils.FastaReader;
import utils.SQLMappingDatabase;
import utils.taxTree.TaxNode;
import utils.taxTree.TaxTree;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicLong;

public class SuppressorManager {

    private CountDownLatch latch;
    private ExecutorService executor;

    private long maxProgress = 0, lastProgress = 0;
    private AtomicLong progress = new AtomicLong();
    private long time;

    public void runClustering(String rank, String srcPath, String aliFolder, TaxTree taxTree, SQLMappingDatabase mappingDatabase, int cores, int memory, int identity, File tmpFile) {

        File outFolder = new File(srcPath + File.separator + rank + "_marker_proteins_clustered");
        outFolder.mkdir();
        File[] initFaaFiles = new File(srcPath + File.separator + "genus_database_proteins").listFiles((dir, name) -> name.endsWith(".faa") && !name.endsWith("_new.faa"));
        ArrayList<File> faaFiles = new ArrayList<>(Arrays.asList(initFaaFiles));
        Collections.sort(faaFiles, Comparator.comparingLong(File::length).reversed());
        faaFiles.stream().forEach(f -> maxProgress += f.length());
        initFaaFiles = faaFiles.toArray(new File[faaFiles.size()]);

        System.out.println(">Assessing new proteins for " + initFaaFiles.length + " protein sets");
        ArrayList<Runnable> clusterNewFileThreads = new ArrayList<>();
        for (File faaFile : initFaaFiles)
            clusterNewFileThreads.add(new ClusterNewFileThread(faaFile, tmpFile, aliFolder));
        executor = Executors.newFixedThreadPool(cores);
        time = System.currentTimeMillis();
        runInParallel(clusterNewFileThreads);
        executor.shutdown();
        reportFinish();
        System.out.println("Runtime: " + getUptime());

        System.out.println(">Running DIAMOND on " + initFaaFiles.length + " protein sets");
        ArrayList<Runnable> clusterThreadsPart1 = new ArrayList<>();
        for (File faaFile : initFaaFiles)
            clusterThreadsPart1.add(new ClusterThreadPart1(faaFile, tmpFile, outFolder, aliFolder, cores, memory, identity, taxTree, mappingDatabase));
        executor = Executors.newFixedThreadPool(1);
        time = System.currentTimeMillis();
        runInParallel(clusterThreadsPart1);
        executor.shutdown();
        reportFinish();
        System.out.println("Runtime: " + getUptime());

        System.out.println(">Clustering " + initFaaFiles.length + " protein sets");
        ArrayList<Runnable> clusterThreadsPart2 = new ArrayList<>();
        for (File faaFile : initFaaFiles)
            clusterThreadsPart2.add(new ClusterThreadPart2(faaFile, tmpFile, outFolder, aliFolder, identity, taxTree, mappingDatabase));
        executor = Executors.newFixedThreadPool(cores);
        time = System.currentTimeMillis();
        runInParallel(clusterThreadsPart2);
        executor.shutdown();
        reportFinish();
        System.out.println("Runtime: " + getUptime());

    }

    private class ClusterNewFileThread implements Runnable {

        private File faaFile, tmpFile;
        private String aliFolder;

        public ClusterNewFileThread(File faaFile, File tmpFile, String aliFolder) {
            this.faaFile = faaFile;
            this.tmpFile = tmpFile;
            this.aliFolder = aliFolder;
        }

        @Override
        public void run() {
            try {
                String genus = faaFile.getName().replaceAll("\\.faa", "").replaceAll("_", "");
                SQLAlignmentDatabase alignmentDatabase = createAlignmentDatabase(aliFolder, genus, tmpFile);
                File newFile = new File(faaFile.getAbsolutePath().replaceAll("\\.faa", "_new.faa"));
                newFile.delete();
                newFile.createNewFile();
                BufferedWriter writer = new BufferedWriter(new FileWriter(newFile));
                ArrayList<Object[]> tokens = FastaReader.read(faaFile);
                for (Object[] token : tokens) {
                    String acc = token[0].toString();
                    if (!alignmentDatabase.containsAcc(acc, genus + "_clusterTable"))
                        writer.write(">" + acc + "\n" + token[1].toString() + "\n");
                }
                writer.close();
                alignmentDatabase.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
            reportProgress(faaFile.length());
            latch.countDown();
        }

    }

    private class ClusterThreadPart1 implements Runnable {

        private File faaFile, tmpFile, outFolder;
        private String aliFolder;
        private int cores, memory, identity;
        private TaxTree taxTree;
        private SQLMappingDatabase mappingDatabase;

        public ClusterThreadPart1(File faaFile, File tmpFile, File outFolder, String aliFolder, int cores, int memory, int identity, TaxTree taxTree, SQLMappingDatabase mappingDatabase) {
            this.faaFile = faaFile;
            this.tmpFile = tmpFile;
            this.outFolder = outFolder;
            this.aliFolder = aliFolder;
            this.cores = cores;
            this.memory = memory;
            this.identity = identity;
            this.taxTree = taxTree;
            this.mappingDatabase = createMappingDatabase(mappingDatabase);
        }

        @Override
        public void run() {

            String genus = faaFile.getName().replaceAll("\\.faa", "").replaceAll("_", "");
            SQLAlignmentDatabase sqlAliDatabase = createAlignmentDatabase(aliFolder, genus, tmpFile);

            try {

                File newFile = new File(faaFile.getAbsolutePath().replaceAll("\\.faa", "_new.faa"));
                newFile.createNewFile();
                if (newFile.exists() && newFile.length() > 0) {
                    File dbFile1 = DiamondRunner.makedb(faaFile, cores);
                    File tabFile1 = DiamondRunner.blastp(dbFile1, newFile, tmpFile, identity, memory, cores);
                    sqlAliDatabase.addAlignmentTable(genus + "_clusterTable", null, tabFile1, false);
                    dbFile1.delete();
                    tabFile1.delete();

                    File dbFile2 = DiamondRunner.makedb(newFile, cores);
                    File tabFile2 = DiamondRunner.blastp(dbFile2, faaFile, tmpFile, identity, memory, cores);
                    sqlAliDatabase.addAlignmentTable(genus + "_clusterTable", null, tabFile2, false);
                    dbFile2.delete();
                    tabFile2.delete();
                }
                newFile.delete();

            } catch (Exception e) {
                e.printStackTrace();
            }

            mappingDatabase.close();
            sqlAliDatabase.close();
            reportProgress(faaFile.length());
            latch.countDown();

        }
    }

    private class ClusterThreadPart2 implements Runnable {

        private String aliFolder, outFolder, genus;
        private File faaFile, tmpFile;
        private TaxTree taxTree;
        private int identity;
        private SQLMappingDatabase mappingDatabase;

        public ClusterThreadPart2(File faaFile, File tmpFile, File outFolder, String aliFolder, int identity, TaxTree taxTree, SQLMappingDatabase mappingDatabase) {
            this.faaFile = faaFile;
            this.tmpFile = tmpFile;
            this.outFolder = outFolder.getAbsolutePath();
            this.identity = identity;
            this.aliFolder = aliFolder;
            this.taxTree = taxTree;
            this.mappingDatabase = createMappingDatabase(mappingDatabase);
            this.genus = faaFile.getName().replaceAll("\\.faa", "").replaceAll("_", "");
        }

        @Override
        public void run() {
            SQLAlignmentDatabase aliDatabase = createAlignmentDatabase(aliFolder, genus, tmpFile);
            File outFile = new File(outFolder + File.separator + faaFile.getName().replaceAll("\\.faa", "_clustered.faa"));
            outFile.delete();
            new Clustering().run(genus, aliDatabase, faaFile, outFile, identity, identity);
            mappingDatabase.close();
            aliDatabase.close();
            reportProgress(faaFile.length());
            latch.countDown();
        }

    }

    public void runMarker(String rank, String srcPath, String aliFolder, TaxTree taxTree, SQLMappingDatabase mappingDatabase, int cores, int memory, int identity, File tmpFile) {

        try {

            File clusterFolder = new File(srcPath + File.separator + rank + "_marker_proteins_clustered");
            File outFolder = new File(srcPath + File.separator + rank + "_marker_proteins");
            outFolder.mkdir();
            ArrayList<File> faaFiles = new ArrayList<>(Arrays.asList(clusterFolder.listFiles((dir, name) -> name.endsWith(".faa") && !name.endsWith("_new.faa"))));
            Collections.sort(faaFiles, Comparator.comparingLong(File::length).reversed());
            maxProgress = 0;
            faaFiles.stream().forEach(f -> maxProgress += f.length());

            System.out.println(">Assessing newly computed representatives for " + faaFiles.size() + " protein sets");
            time = System.currentTimeMillis();
            ArrayList<Runnable> newFileThreads = new ArrayList<>();
            for (File faaFile : faaFiles)
                newFileThreads.add(new MarkerNewFileThread(faaFile, tmpFile, outFolder.getAbsolutePath(), aliFolder));
            executor = Executors.newFixedThreadPool(cores);
            runInParallel(newFileThreads);
            executor.shutdown();
            reportFinish();
            System.out.println("Runtime: " + getUptime());

            System.out.println(">Aligning new vs all proteins using DIAMOND");
            time = System.currentTimeMillis();
            File clusteredProteins = new File(outFolder + File.separator + "cluster.faa");
            clusteredProteins.createNewFile();
            clusteredProteins.deleteOnExit();
            for (File f : clusterFolder.listFiles((dir, name) -> name.endsWith("_clustered.faa")))
                appendToFile(f, clusteredProteins);
            File clusterDb = DiamondRunner.makedb(clusteredProteins, cores);
            clusteredProteins.delete();
            clusterDb.deleteOnExit();
            File newProteins = new File(outFolder + File.separator + "newProteins.faa");
            newProteins.createNewFile();
            newProteins.deleteOnExit();
            for (File f : outFolder.listFiles((dir, name) -> name.endsWith("_new.faa")))
                appendToFile(newProteins, f);
            File tabFile1 = DiamondRunner.blastp(clusterDb, newProteins, tmpFile, identity, memory, cores);
            newProteins.delete();
            clusterDb.delete();

            System.out.println(">Aligning all vs new proteins using DIAMOND");
            File newClusteredProteins = new File(outFolder + File.separator + "newCluster.faa");
            newClusteredProteins.createNewFile();
            newClusteredProteins.deleteOnExit();
            for (File f : outFolder.listFiles((dir, name) -> name.endsWith("_new.faa")))
                appendToFile(f, newClusteredProteins);
            File newClusterDb = newClusteredProteins.exists() ? DiamondRunner.makedb(newClusteredProteins, cores) : null;
            newClusterDb.deleteOnExit();
            newClusteredProteins.delete();
            File allProteins = new File(outFolder + File.separator + "allProteins.faa");
            allProteins.createNewFile();
            allProteins.deleteOnExit();
            for (File f : faaFiles)
                appendToFile(allProteins, f);
            File tabFile2 = DiamondRunner.blastp(newClusterDb, allProteins, tmpFile, identity, memory, cores);
            allProteins.delete();
            newClusterDb.delete();
            System.out.println("Runtime: " + getUptime());

            System.out.println(">Processing DIAMOND result for " + faaFiles.size() + " genera");
            time = System.currentTimeMillis();
            ArrayList<Runnable> markerThreadsPart1 = new ArrayList<>();
            File[] tabFiles = {tabFile1, tabFile2};
            for (File faaFile : faaFiles)
                markerThreadsPart1.add(new MarkerThreadPart1(faaFile, tabFiles, aliFolder, tmpFile, taxTree, mappingDatabase, outFolder));
            executor = Executors.newFixedThreadPool(1);
            runInParallel(markerThreadsPart1);
            executor.shutdown();
            reportFinish();
            System.out.println("Runtime: " + getUptime());

            System.out.println(">Computing marker proteins for " + faaFiles.size() + " protein sets");
            time = System.currentTimeMillis();
            ArrayList<Runnable> markerThreadsPart2 = new ArrayList<>();
            for (File faaFile : faaFiles)
                markerThreadsPart2.add(new MarkerThreadPart2(faaFile, tmpFile, identity, aliFolder, outFolder, taxTree, mappingDatabase));
            executor = Executors.newFixedThreadPool(cores);
            runInParallel(markerThreadsPart2);
            executor.shutdown();
            reportFinish();
            System.out.println("Runtime: " + getUptime());

            // removing new files
            for (Runnable t : newFileThreads)
                ((MarkerNewFileThread) t).newFile.delete();

        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    private class MarkerThreadPart2 implements Runnable {

        private String outFolder, aliFolder;
        private File faaFile, tmpFile;
        private int identity;
        private TaxTree taxTree;
        private SQLMappingDatabase mappingDatabase;

        public MarkerThreadPart2(File faaFile, File tmpFile, int identity, String aliFolder, File outFolder, TaxTree taxTree, SQLMappingDatabase mappingDatabase) {
            this.faaFile = faaFile;
            this.tmpFile = tmpFile;
            this.identity = identity;
            this.aliFolder = aliFolder;
            this.outFolder = outFolder.getAbsolutePath();
            this.taxTree = taxTree;
            this.mappingDatabase = createMappingDatabase(mappingDatabase);
        }

        @Override
        public void run() {
            String genus = faaFile.getName().replaceAll("_clustered\\.faa", "").replaceAll("_", "");
            SQLAlignmentDatabase aliDatabase = createAlignmentDatabase(aliFolder, genus, tmpFile);
            File outFile = new File(outFolder + File.separator + genus + "_marker.faa");
            new Selecting().run(genus, taxTree, mappingDatabase, aliDatabase, faaFile, outFile, identity, identity);
            aliDatabase.close();
            mappingDatabase.close();
            reportProgress(faaFile.length());
            latch.countDown();
        }

    }

    private class MarkerThreadPart1 implements Runnable {

        private File[] tabFiles;
        private File faaFile, tmpFile, outFolder;
        private String aliFolder;
        private TaxTree taxTree;
        private SQLMappingDatabase mappingDatabase;

        public MarkerThreadPart1(File faaFile, File[] tabFiles, String aliFolder, File tmpFile, TaxTree taxTree, SQLMappingDatabase mappingDatabase, File outFolder) {
            this.faaFile = faaFile;
            this.tabFiles = tabFiles;
            this.aliFolder = aliFolder;
            this.tmpFile = tmpFile;
            this.taxTree = taxTree;
            this.mappingDatabase = createMappingDatabase(mappingDatabase);
            this.outFolder = outFolder;
        }

        @Override
        public void run() {

            String genus = faaFile.getName().replaceAll("_clustered\\.faa", "").replaceAll("_", "");
            File tabFile = new File(outFolder + File.separator + genus + "_marker.tab");
            tabFile.deleteOnExit();
            try {
                BufferedWriter writer = new BufferedWriter(new FileWriter(tabFile));
                for (File tab : tabFiles) {
                    if (tab == null || !tab.exists() || tab.length() == 0)
                        continue;
                    BufferedReader buf = new BufferedReader(new FileReader(tab));
                    String line;
                    while ((line = buf.readLine()) != null) {
                        final String[] tokens = line.split("\t");
                        final String qacc = tokens[0];
                        final String g = getRank(mappingDatabase.getTaxIdByAcc(qacc), "genus");
                        if (g != null && g.equals(genus))
                            writer.write(line + "\n");
                    }
                    buf.close();
                }
                writer.close();
            } catch (Exception e) {
                e.printStackTrace();
            }

            SQLAlignmentDatabase aliDatabase = createAlignmentDatabase(aliFolder, genus, tmpFile);
            aliDatabase.addAlignmentTable(genus + "_markerTable", null, tabFile, false);
            tabFile.delete();
            aliDatabase.close();
            mappingDatabase.close();
            reportProgress(faaFile.length());
            latch.countDown();

        }

        private String getRank(ArrayList<Integer> taxids, String rank) {
            String name = null;
            for (int taxid : taxids) {
                TaxNode w = taxTree.getNode(taxid);
                while (w != null) {
                    if (w.getRank().equals(rank)) {
                        if (name == null)
                            name = w.getName();
                        else if (!name.equals(w.getName()))
                            return null;
                        break;
                    }
                    w = w.getParent();
                }
            }
            if (name != null)
                return name.replaceAll("_", "");
            return null;
        }
    }

    private class MarkerNewFileThread implements Runnable {

        private File faaFile, tmpFile;
        private String outFolder, aliFolder;
        private File newFile;

        public MarkerNewFileThread(File faaFile, File tmpFile, String outFolder, String aliFolder) {
            this.faaFile = faaFile;
            this.tmpFile = tmpFile;
            this.outFolder = outFolder;
            this.aliFolder = aliFolder;
        }

        @Override
        public void run() {

            String genus = faaFile.getName().replaceAll("_clustered\\.faa", "").replaceAll("_", "");
            SQLAlignmentDatabase sqlAliDatabase = createAlignmentDatabase(aliFolder, genus, tmpFile);

            try {

                newFile = new File(outFolder + File.separator + genus + "_new.faa");
                newFile.createNewFile();
                BufferedWriter writer = new BufferedWriter(new FileWriter(newFile));
                ArrayList<Object[]> tokens = FastaReader.read(faaFile);
                for (Object[] o : tokens) {
                    String acc = o[0].toString();
                    if (!sqlAliDatabase.containsAcc(acc, genus + "_markerTable"))
                        writer.write(">" + acc + "\n" + o[1] + "\n");
                }
                writer.close();
            } catch (Exception e) {
                e.printStackTrace();
            }

            sqlAliDatabase.close();
            reportProgress(faaFile.length());
            latch.countDown();
        }
    }

    private void runInParallel(ArrayList<Runnable> threads) {
        latch = new CountDownLatch(threads.size());
        for (Runnable t : threads)
            executor.submit(t);
        try {
            latch.await();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    private synchronized SQLMappingDatabase createMappingDatabase(SQLMappingDatabase mappingDatabase) {
        return new SQLMappingDatabase(mappingDatabase);
    }

    private synchronized SQLAlignmentDatabase createAlignmentDatabase(String aliFolder, String genus, File tmpFile) {
        return new SQLAlignmentDatabase(aliFolder, genus, tmpFile);
    }

    private void appendToFile(File source, File target) {
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(target, true));
            for (Object[] o : FastaReader.read(source))
                writer.write(">" + o[0].toString() + "\n" + o[1].toString() + "\n");
            writer.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private synchronized void reportProgress(long delta) {
        progress.getAndAdd(delta);
        int p = ((int) ((((double) progress.get() / (double) maxProgress)) * 100) / 5) * 5;
        if (p > lastProgress && p < 100) {
            lastProgress = p;
            System.out.print(p + "% (" + getUptime() + ") ");
        }
    }

    private void reportFinish() {
        progress.set(0);
        lastProgress = 0;
        System.out.print(100 + "%\n");
    }

    private String getUptime() {
        long runtime = (System.currentTimeMillis() - time) / 1000;
        return runtime + "s";
    }

}
