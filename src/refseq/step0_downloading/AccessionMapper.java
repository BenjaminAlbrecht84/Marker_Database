package refseq.step0_downloading;

import refseq.utils.AssemblyParser;
import utils.FastaReader;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Objects;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

public class AccessionMapper {

    private int maxProgress, lastProgress = 0;
    private AtomicInteger progress = new AtomicInteger();
    private CountDownLatch latch;
    private ExecutorService executor;
    private long time = System.currentTimeMillis();

    private ArrayList<File> faaFiles;
    private int filePointer = 0;

    private HashMap<String, Integer> gcfToTaxID;
    private BufferedWriter fW;

    public void run(String src, File downloadFolder, File outFile, int cores) throws IOException {

        gcfToTaxID = new AssemblyParser().getGcf2Taxid(src);
        if (outFile.exists())
            outFile.delete();
        fW = new BufferedWriter(new FileWriter(outFile, true));

        try {
            executor = Executors.newFixedThreadPool(cores);
            System.out.println(">Creating mapping file");
            ArrayList<Runnable> mappers = new ArrayList<>();
            for (int i = 0; i < cores; i++)
                mappers.add(new Mapper());
            faaFiles = new ArrayList<>();
            for (File dir : Objects.requireNonNull(downloadFolder.listFiles(pathname -> pathname.isDirectory()))) {
                for (File f : Objects.requireNonNull(dir.listFiles((dir1, name) -> name.endsWith(".faa.gz"))))
                    faaFiles.add(f);
            }
            maxProgress = faaFiles.size();
            runInParallel(mappers);
            reportFinish();
            System.out.println("Runtime: " + getRuntime());

        } finally {
            fW.close();
            executor.shutdown();
        }

    }

    private synchronized File nextFile() {
        if (filePointer < faaFiles.size())
            return faaFiles.get(filePointer++);
        return null;
    }

    private class Mapper implements Runnable {

        @Override
        public void run() {
            File faaFile;
            while ((faaFile = nextFile()) != null) {
                String searchID = "GCF_" + faaFile.getName().split("_")[1];
                Integer taxID = gcfToTaxID.get(searchID);
                if (taxID != null) {
                    for (Object[] o : FastaReader.read(faaFile)) {
                        String acc = o[0].toString();
                        try {
                            fW.write(acc + "\t" + searchID + "\t" + taxID + "\n");
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    }
                }
                reportProgress(1);
            }
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
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    private void reportProgress(int delta) {
        progress.addAndGet(delta);
        int p = ((int) ((((double) progress.get() / (double) maxProgress)) * 100) / 5) * 5;
        if (p > lastProgress && p < 100) {
            lastProgress = p;
            System.out.print(p + "% (" + getRuntime() + ") ");
        }
    }

    private void reportFinish() {
        progress.set(0);
        lastProgress = 0;
        System.out.print(100 + "%\n");
    }

    private String getRuntime() {
        long runtime = (System.currentTimeMillis() - time) / 1000;
        return runtime + "s";
    }

}
