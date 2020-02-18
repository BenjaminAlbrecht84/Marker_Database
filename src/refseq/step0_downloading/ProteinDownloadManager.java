package refseq.step0_downloading;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

import refseq.utils.AssemblyParser;
import refseq.utils.Downloader;
import utils.FastaReader;
import utils.SQLMappingDatabase;
import utils.taxTree.TaxNode;
import utils.taxTree.TaxTree;

public class ProteinDownloadManager {

    private String[] majorRanks = {"species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom", "root"};
    private ArrayList<String> rankList = new ArrayList<>(Arrays.asList(majorRanks));

    private TaxTree taxTree;

    private int DOWNLOAD_THREADS = 8;
    private int maxProgress, lastProgress = 0;
    private AtomicInteger progress = new AtomicInteger();

    private CountDownLatch latch;
    private ExecutorService executor;
    private ArrayList<Runnable> threads;

    private ArrayList<TaxNode> ftpNodes;
    private int ftpPointer = 0;
    private File[] proteinDirs;
    private int dirPointer = 0;

    private long time = System.currentTimeMillis();

    public void run(String srcPath, SQLMappingDatabase mappingDatabase, TaxTree taxTree, int cores) {

        this.DOWNLOAD_THREADS = 2 * cores;
        this.taxTree = taxTree;

        try {

            File refseqProteins = new File(srcPath + File.separator + "refseq_proteins");

            // adding FTP address to leaves
            System.out.println(">Parsing assembly summary table ");
            time = System.currentTimeMillis();
            HashMap<Integer, ArrayList<String>> taxidToFTP = new AssemblyParser().getTaxidToFTP(srcPath);
            ftpNodes = new ArrayList<>();
            int totalFtpLinks = 0;
            for (TaxNode v : taxTree.getNodes()) {
                if (taxidToFTP.containsKey(v.getTaxid()) && isBacteriaAcc(v)) {
                    v.setInfo(taxidToFTP.get(v.getTaxid()));
                    ftpNodes.add(v);
                    totalFtpLinks += taxidToFTP.get(v.getTaxid()).size();
                    v.reportFTPLeaf();
                }
            }
            System.out.println("Runtime: " + getUptime());

            // downloading protein data
            executor = Executors.newFixedThreadPool(DOWNLOAD_THREADS);
            System.out.println(">Downloading protein information for " + ftpNodes.size() + " taxonomic node(s)");
            time = System.currentTimeMillis();
            refseqProteins.mkdir();
            ConcurrentHashMap<String, File> proteinToFile = new ConcurrentHashMap<String, File>();
            listFilesRec(refseqProteins, proteinToFile);
            threads = new ArrayList<>();
            for (int i = 0; i < DOWNLOAD_THREADS; i++)
                threads.add(new ProteinDownloadingThread(refseqProteins, proteinToFile));
            maxProgress = totalFtpLinks;
            runInParallel(threads);
            reportFinish();
            executor.shutdown();
            System.out.println("Runtime: " + getUptime());

            // update mapping database
            File acc2gcf2taxid = new File(srcPath + File.separator + "acc2gcfs2taxid.tab");
            new AccessionMapper().run(srcPath, refseqProteins, acc2gcf2taxid, cores);
            mappingDatabase.createAcc2gcf2taxidTable(acc2gcf2taxid);

            // collecting protein data
            executor = Executors.newFixedThreadPool(cores);
            System.out.println(">Creating " + refseqProteins.listFiles(File::isDirectory).length + " genus protein sets");
            time = System.currentTimeMillis();
            File proteinFolder = new File(srcPath + File.separator + "genus_database_proteins");
            proteinFolder.mkdir();
            threads = new ArrayList<>();
            proteinDirs = refseqProteins.listFiles(File::isDirectory);
            dirPointer = 0;
            for (int i = 0; i < cores; i++)
                threads.add(new ProteinCollector(proteinFolder, mappingDatabase));
            maxProgress = 0;
            for (File dir : refseqProteins.listFiles(File::isDirectory))
                maxProgress += dir.listFiles().length;
            runInParallel(threads);
            reportFinish();
            System.out.println("Runtime: " + getUptime());
            executor.shutdown();

        } catch (Exception e) {
            e.printStackTrace();
        }


    }

    private synchronized File nextProteinDir() {
        if (dirPointer < proteinDirs.length)
            return proteinDirs[dirPointer++];
        return null;
    }

    public class ProteinCollector implements Runnable {

        private File proteinFolder;
        private SQLMappingDatabase mapping;

        public ProteinCollector(File proteinFolder, SQLMappingDatabase mairaDatabase) {
            this.proteinFolder = proteinFolder;
            this.mapping = new SQLMappingDatabase(mairaDatabase.getDatabaseFile(), mairaDatabase.getTmpDir());
        }

        @Override
        public void run() {

            File dir;
            while ((dir = nextProteinDir()) != null) {

                HashMap<String, HashMap<Integer, HashSet<String>>> seenAccessionMap = new HashMap<>();
                File out = new File(proteinFolder + "/" + dir.getName().replace("\\s+", "_") + ".faa");
                if (out.exists()) {
                    ArrayList<Object[]> tokens = FastaReader.read(out);
                    for (Object[] o : tokens) {
                        if (!o[1].toString().isEmpty())
                            addSeenAccession(seenAccessionMap, o[0].toString());
                    }
                }

                try {
                    BufferedWriter fW = new BufferedWriter(new FileWriter(out, true));
                    File[] files = dir.listFiles();
                    try {
                        for (File f : files) {
                            if (f.isFile()) {
                                ArrayList<Object[]> info = FastaReader.read(f);
                                for (Object[] o : info) {
                                    String acc = o[0].toString();
                                    String seq = o[1].toString();
                                    if (!isSeenAccession(seenAccessionMap, acc)) {
                                        addSeenAccession(seenAccessionMap, acc);
                                        fW.write(">" + acc + "\n" + seq + "\n");
                                    }
                                }
                            }
                            reportProgress(1);
                        }
                    } finally {
                        fW.close();
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                }

            }
            latch.countDown();
        }

        private void addSeenAccession(HashMap<String, HashMap<Integer, HashSet<String>>> seenAccessionMap, String acc) {
            String prefix = acc.split("_")[0];
            int id = Integer.parseInt(acc.split("_")[1].replaceAll("\\.", ""));
            int bucket = id % 10;
            seenAccessionMap.putIfAbsent(prefix, new HashMap<>());
            seenAccessionMap.get(prefix).putIfAbsent(bucket, new HashSet<>());
            seenAccessionMap.get(prefix).get(bucket).add(acc);
        }

        private boolean isSeenAccession(HashMap<String, HashMap<Integer, HashSet<String>>> seenAccessionMap, String acc) {
            String prefix = acc.split("_")[0];
            int id = Integer.parseInt(acc.split("_")[1].replaceAll("\\.", ""));
            int bucket = id % 10;
            if (!seenAccessionMap.containsKey(prefix) || !seenAccessionMap.get(prefix).containsKey(bucket) || !seenAccessionMap.get(prefix).get(bucket).contains(acc))
                return false;
            return true;
        }

    }

    private void listFilesRec(File file, ConcurrentHashMap<String, File> proteinToFile) {
        if (file.isDirectory()) {
            for (File f : file.listFiles())
                listFilesRec(f, proteinToFile);
        } else
            proteinToFile.put(file.getName(), file);
    }

    private boolean isBacteriaAcc(TaxNode v) {
        TaxNode p = v;
        while (p != null) {
            if (p.getName().equals("Bacteria"))
                return true;
            p = p.getParent();
        }
        return false;
    }

    private synchronized TaxNode nextFTPNode() {
        if (ftpPointer < ftpNodes.size())
            return ftpNodes.get(ftpPointer++);
        return null;
    }

    public class ProteinDownloadingThread implements Runnable {

        private ConcurrentHashMap<String, File> proteinToFile;
        private File proteinFolder;

        public ProteinDownloadingThread(File proteinFolder,
                                        ConcurrentHashMap<String, File> proteinToFile) {
            this.proteinFolder = proteinFolder;
            this.proteinToFile = proteinToFile;
        }

        @Override
        public void run() {

            try {

                TaxNode v;
                while ((v = nextFTPNode()) != null) {

                    ArrayList<String> ftps = (ArrayList<String>) v.getInfo();
                    if (ftps == null)
                        continue;
                    String genus = getRank(v, "genus");
                    if (genus == null)
                        continue;

                    for (String ftp : ftps) {

                        boolean isReadable = false;
                        while (!isReadable) {

                            String[] split = ftp.split("/");
                            String outFile = split[split.length - 1] + "_protein.ftp";
                            String remoteFaaFile = split[split.length - 1] + "_protein.faa.gz";
                            File localFolder = new File(proteinFolder + File.separator + genus.replaceAll("\\s+", "_"));
                            if (!proteinToFile.containsKey(remoteFaaFile)) {
                                localFolder.mkdir();
                                Downloader.getFtpFile(ftp, remoteFaaFile, localFolder, outFile);
                            }

                            File localFaaFile = new File(localFolder + File.separator + remoteFaaFile);
                            if (FastaReader.read(localFaaFile) == null) {
                                proteinToFile.remove(remoteFaaFile);
                                localFaaFile.delete();
                            } else
                                isReadable = true;

                        }

                    }

                    reportProgress(ftps.size());

                }

            } catch (Exception e) {
                latch.countDown();
                e.printStackTrace();
            }

            latch.countDown();

        }

    }

    private String getRank(TaxNode v, String rank) {
        TaxNode w = v;
        while (w != null) {
            if (w.getRank().equals(rank))
                return w.getName();
            w = w.getParent();
        }
        return null;
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
