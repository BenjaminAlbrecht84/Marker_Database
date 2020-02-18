package refseq.utils.aliHelper;

import refseq.utils.DiamondRunner;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

public class AddAlignmentsHelper {

    private static File[] daaFiles;
    private static int daaPointer = 0;
    private static File tmpDir;

    private static long time = System.currentTimeMillis();
    private static int maxProgress, lastProgress = 0;
    private static AtomicInteger progress = new AtomicInteger();
    private static CountDownLatch latch;
    private static ExecutorService executor;

    public static void main(String[] args) {

        File clusterAlignmentFolder = new File(args[0]);
        File markerAlignmentFolder = new File(args[1]);
        File databaseFolder = new File(args[2]);
        tmpDir = new File(args[3]);
        int cores = Integer.parseInt(args[4]);
        executor = Executors.newFixedThreadPool(cores);

        daaFiles = clusterAlignmentFolder.listFiles((dir, name) -> name.endsWith(".daa"));
        maxProgress = daaFiles.length;
        ArrayList<Runnable> threads = new ArrayList<>();
        daaPointer = 0;
//        for (int i = 0; i < cores; i++)
//            threads.add(new AlignmentsAdder(databaseFolder.getAbsolutePath(), null, "clusterTable", 1));
//        runInParallel(threads);
//        reportFinish();

        daaFiles = markerAlignmentFolder.listFiles((dir, name) -> name.endsWith(".daa") && !name.endsWith("2015.daa"));
        maxProgress = daaFiles.length;
        threads = new ArrayList<>();
        daaPointer = 0;
        for (int i = 0; i < cores; i++)
            threads.add(new AlignmentsAdder(databaseFolder.getAbsolutePath(), clusterAlignmentFolder.getAbsolutePath(), "markerTable", 1));
        runInParallel(threads);
        reportFinish();

        executor.shutdown();

    }

    private static synchronized File nextDaaFile() {
        if (daaPointer < daaFiles.length)
            return daaFiles[daaPointer++];
        return null;
    }

    private static synchronized SQLAlignmentDatabase createDatabase(String databaseFile, String genus, File tmpDir) {
        return new SQLAlignmentDatabase(databaseFile, genus, tmpDir);
    }

    private static class AlignmentsAdder implements Runnable {

        private String dbFolder, srcFolder, dbType;
        private int cores;

        public AlignmentsAdder(String dbFolder, String srcFolder, String dbType, int cores) {
            this.cores = cores;
            this.dbType = dbType;
            this.dbFolder = dbFolder;
            this.srcFolder = srcFolder;
        }

        @Override
        public void run() {
            File daa;
            while ((daa = nextDaaFile()) != null) {
                String genus = daa.getName().replaceAll("_clustered", "");
                genus = genus.replaceAll("\\.daa", "");
                genus = genus.replaceAll("_", "");
                SQLAlignmentDatabase aliDb = createDatabase(dbFolder, genus, tmpDir);
                File tab = DiamondRunner.view(daa, cores);
                File src = srcFolder != null ? new File(srcFolder + File.separator + daa.getName().replaceAll("\\.daa", ".faa")) : null;
                aliDb.addAlignmentTable(genus + "_" + dbType, src, tab, true);
                tab.delete();
                reportProgress(1);
            }
            latch.countDown();
        }
    }

    private static void runInParallel(ArrayList<Runnable> threads) {
        latch = new CountDownLatch(threads.size());
        for (Runnable t : threads)
            executor.submit(t);
        try {
            latch.await();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    private static synchronized void reportProgress(int delta) {
        progress.getAndAdd(delta);
        int p = ((int) ((((double) progress.get() / (double) maxProgress)) * 100) / 5) * 5;
        if (p > lastProgress && p < 100) {
            lastProgress = p;
            System.out.print(p + "% (" + getUptime() + ") ");
        }
    }

    private static void reportFinish() {
        progress.set(0);
        lastProgress = 0;
        System.out.print(100 + "%\n");
    }

    private static String getUptime() {
        long runtime = (System.currentTimeMillis() - time) / 1000;
        return runtime + "s";
    }


}
