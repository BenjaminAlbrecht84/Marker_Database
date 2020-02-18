package refseq.step2_selecting;

import refseq.utils.aliHelper.SQLAlignmentDatabase;
import utils.SQLMappingDatabase;
import utils.taxTree.TaxTree;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

public class FilterManager {

    private CountDownLatch latch;
    private ExecutorService executor;

    private int maxProgress, lastProgress = 0;
    private AtomicInteger progress = new AtomicInteger();
    private long time;

    private BufferedWriter factorWriter, markerWriter;
    private TaxTree taxTree;
    private int NUM_OF_PROTEINS, MIN_ID;
    private String aliFolder;
    private File tmpDir;

    public void run(String rank, String srcPath, String aliFolder, File tmpDir, TaxTree taxTree, SQLMappingDatabase mappingDatabase, int NUM_OF_PROTEINS, int MIN_ID, int cores) {

        this.taxTree = taxTree;
        this.NUM_OF_PROTEINS = NUM_OF_PROTEINS;
        this.MIN_ID = MIN_ID;
        this.aliFolder = aliFolder;
        this.tmpDir = tmpDir;
        this.time = System.currentTimeMillis();

        File markerDir = new File(srcPath + File.separator + rank + "_marker_proteins");

        try {

            File markerFile = new File(srcPath + File.separator + "genus_marker_db_" + NUM_OF_PROTEINS + ".faa");
            markerFile.delete();
            markerWriter = new BufferedWriter(new FileWriter(markerFile, true));
            File out = new File(srcPath + File.separator + "genus_marker_db_" + NUM_OF_PROTEINS + "_weights.tab");
            out.delete();
            factorWriter = new BufferedWriter(new FileWriter(out, true));

            System.out.println(">Filtering marker proteins");
            File[] faaFiles = markerDir.listFiles();
            ArrayList<Runnable> threads = new ArrayList<>();
            for (File faaFile : faaFiles)
                threads.add(new SelectorThread(faaFile, mappingDatabase));
            maxProgress = threads.size();
            executor = Executors.newFixedThreadPool(cores);
            runInParallel(threads);
            executor.shutdown();
            markerWriter.close();
            factorWriter.close();

            reportFinish();
            System.out.println("Runtime: " + getUptime());


        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    private class SelectorThread implements Runnable {

        private File faaFile;
        private SQLMappingDatabase mappingDatabase;
        private SQLAlignmentDatabase alignmentDatabase;
        private String genus;

        public SelectorThread(File faaFile, SQLMappingDatabase mappingDatabase) {
            this.faaFile = faaFile;
            this.mappingDatabase = createMappingDatabase(mappingDatabase);
            this.genus = faaFile.getName().replaceAll("_marker\\.faa", "").replaceAll("_", "");
            this.alignmentDatabase = createAlignmentDatabase(genus);
        }

        @Override
        public void run() {
            new Filtering().run(faaFile, genus, factorWriter, markerWriter, taxTree, mappingDatabase, alignmentDatabase, NUM_OF_PROTEINS, MIN_ID);
            mappingDatabase.close();
            alignmentDatabase.close();
            reportProgress(1);
            latch.countDown();
        }
    }

    private synchronized SQLMappingDatabase createMappingDatabase(SQLMappingDatabase mappingDatabase) {
        return new SQLMappingDatabase(mappingDatabase);
    }

    private synchronized SQLAlignmentDatabase createAlignmentDatabase(String genus) {
        return new SQLAlignmentDatabase(aliFolder, genus, tmpDir);
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

    private synchronized void reportProgress(int delta) {
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
