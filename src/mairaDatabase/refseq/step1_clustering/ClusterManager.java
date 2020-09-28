package mairaDatabase.refseq.step1_clustering;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import mairaDatabase.refseq.utils.DiamondRunner;
import mairaDatabase.refseq.utils.Formatter;
import mairaDatabase.refseq.utils.aliHelper.SQLAlignmentDatabase;
import mairaDatabase.utils.FastaReader;
import mairaDatabase.utils.ResourceLoader;
import mairaDatabase.utils.SQLMappingDatabase;
import mairaDatabase.utils.FastaReader.FastaEntry;
import mairaDatabase.utils.taxTree.TaxTree;

public class ClusterManager {
	
	private ResourceLoader rL = new ResourceLoader();

    public void runClustering(String rank, String srcPath, String aliFolder, TaxTree taxTree, SQLMappingDatabase mappingDatabase, int cores, int memory, int identity, File tmpFile) {

        File outFolder = new File(srcPath + File.separator + rank + "_marker_proteins_clustered");
        outFolder.mkdir();
        File[] initFaaFiles = new File(srcPath + File.separator + "genus_database_proteins").listFiles((dir, name) -> name.endsWith(".faa") && !name.endsWith("_new.faa"));
        ArrayList<File> faaFiles = new ArrayList<>(Arrays.asList(initFaaFiles));
        Collections.sort(faaFiles, Comparator.comparingLong(File::length).reversed());
        long totalFileLength = 0L;
        for(File f : faaFiles)
        	totalFileLength += f.length();
        initFaaFiles = faaFiles.toArray(new File[faaFiles.size()]);

        System.out.println(">Assessing new proteins for " + initFaaFiles.length + " protein sets");
        List<Runnable> collectNewProteinsThreads = new ArrayList<>();
        for (File faaFile : initFaaFiles)
            collectNewProteinsThreads.add(new CollectNewProteinsThread(faaFile, tmpFile, aliFolder));
        rL.runThreads(cores, collectNewProteinsThreads, totalFileLength);

        System.out.println(">Running DIAMOND on " + initFaaFiles.length + " protein sets");
        List<Runnable> alignProteinsThreads = new ArrayList<>();
        for (File faaFile : initFaaFiles)
            alignProteinsThreads.add(new AlignProteinsThread(faaFile, tmpFile, aliFolder, cores, memory, identity, mappingDatabase));
        rL.runThreads(1, alignProteinsThreads, totalFileLength);

        System.out.println(">Clustering " + initFaaFiles.length + " protein sets");
        List<Runnable> clusterProteinsThread = new ArrayList<>();
        for (File faaFile : initFaaFiles)
            clusterProteinsThread.add(new ClusterProteinsThread(faaFile, tmpFile, outFolder, aliFolder, identity, mappingDatabase));
        rL.runThreads(cores, clusterProteinsThread, totalFileLength);

    }

    private class CollectNewProteinsThread implements Runnable {

        private File faaFile, tmpFile;
        private String aliFolder;

        public CollectNewProteinsThread(File faaFile, File tmpFile, String aliFolder) {
            this.faaFile = faaFile;
            this.tmpFile = tmpFile;
            this.aliFolder = aliFolder;
        }

        @Override
        public void run() {
            try {
                String genus = Formatter.removeNonAlphanumerics(faaFile.getName().replaceAll("\\.faa", ""));
                SQLAlignmentDatabase alignmentDatabase = createAlignmentDatabase(aliFolder, genus, tmpFile);
                File newFile = new File(faaFile.getAbsolutePath().replaceAll("\\.faa", "_new.faa"));
                newFile.delete();
                newFile.createNewFile();
                BufferedWriter writer = new BufferedWriter(new FileWriter(newFile));
                ArrayList<FastaEntry> tokens = FastaReader.read(faaFile);
                for (FastaEntry token : tokens) {
                    String acc = token.getName();
                    if (!alignmentDatabase.containsAcc(acc, genus + "_clusterTable"))
                        writer.write(">" + acc + "\n" + token.getSequence() + "\n");
                }
                writer.close();
                alignmentDatabase.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
            rL.reportProgress(faaFile.length());
            rL.countDown();
        }

    }

    private class AlignProteinsThread implements Runnable {

        private File faaFile, tmpFile;
        private String aliFolder;
        private int cores, memory, identity;
        private SQLMappingDatabase mappingDatabase;

        public AlignProteinsThread(File faaFile, File tmpFile, String aliFolder, int cores, int memory, int identity, SQLMappingDatabase mappingDatabase) {
            this.faaFile = faaFile;
            this.tmpFile = tmpFile;
            this.aliFolder = aliFolder;
            this.cores = cores;
            this.memory = memory;
            this.identity = identity;
            this.mappingDatabase = createMappingDatabase(mappingDatabase);
        }

        @Override
        public void run() {

            String genus = Formatter.removeNonAlphanumerics(faaFile.getName().replaceAll("\\.faa", ""));
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
            rL.reportProgress(faaFile.length());
            rL.countDown();

        }
    }

    private class ClusterProteinsThread implements Runnable {

        private String aliFolder, outFolder, genus;
        private File faaFile, tmpFile;
        private int identity;
        private SQLMappingDatabase mappingDatabase;

        public ClusterProteinsThread(File faaFile, File tmpFile, File outFolder, String aliFolder, int identity, SQLMappingDatabase mappingDatabase) {
            this.faaFile = faaFile;
            this.tmpFile = tmpFile;
            this.outFolder = outFolder.getAbsolutePath();
            this.identity = identity;
            this.aliFolder = aliFolder;
            this.mappingDatabase = createMappingDatabase(mappingDatabase);
            this.genus = Formatter.removeNonAlphanumerics(faaFile.getName().replaceAll("\\.faa", ""));
        }

        @Override
        public void run() {
            SQLAlignmentDatabase aliDatabase = createAlignmentDatabase(aliFolder, genus, tmpFile);
            File outFile = new File(outFolder + File.separator + faaFile.getName().replaceAll("\\.faa", "_clustered.faa"));
            outFile.delete();
            new Clustering().run(genus, aliDatabase, faaFile, outFile, identity, identity);
            mappingDatabase.close();
            aliDatabase.close();
            rL.reportProgress(faaFile.length());
            rL.countDown();
        }

    }

    private synchronized SQLMappingDatabase createMappingDatabase(SQLMappingDatabase mappingDatabase) {
        return new SQLMappingDatabase(mappingDatabase);
    }

    private synchronized SQLAlignmentDatabase createAlignmentDatabase(String aliFolder, String genus, File tmpFile) {
        return new SQLAlignmentDatabase(aliFolder, genus, tmpFile);
    }

}
