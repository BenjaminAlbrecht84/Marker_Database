package refseq.step0_downloading;

import refseq.utils.Downloader;

import java.io.File;

public class NCBIDownloader {

    public void run(String srcPath) {

        System.out.println(">Downloading NCBI taxdump files");
        String remoteFile = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz";
        String out = "new_taxdump.out";
        Downloader.getFile(remoteFile, new File(srcPath), out);

        String tarGzPathLocation = srcPath + File.separator + "new_taxdump.tar.gz";
        ProcessBuilder builder = new ProcessBuilder();
        File taxdump = new File(srcPath + File.separator + "taxdump");
        taxdump.mkdir();
        builder.command("sh", "-c", String.format("tar xfz %s -C %s", tarGzPathLocation, taxdump.getAbsolutePath()));
        builder.directory(new File(srcPath));
        try {
            Process process = builder.start();
            process.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }

        new File(srcPath + File.separator + out).delete();
        new File(tarGzPathLocation).delete();
        for (File f : taxdump.listFiles()) {
            if (f.getName().equals("names.dmp") || f.getName().equals("nodes.dmp"))
                continue;
            f.delete();
        }

        System.out.println(">Downloading assembly summary table");
        new File(srcPath + File.separator + "assembly_summary_refseq.txt").delete();
        remoteFile = "ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt";
        out = "assembly_summary_refseq.out";
        Downloader.getFile(remoteFile, new File(srcPath), out);
        new File(srcPath + File.separator + out).delete();


    }

}
