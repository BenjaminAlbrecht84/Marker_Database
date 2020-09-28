package mairaDatabase.refseq.step0_downloading;

import java.io.File;

import mairaDatabase.refseq.utils.Downloader;

public class NCBIDownloader {

    private static final String NCBI_TAXDUMP_FILE = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz";
    private static final String NCBI_ASSEMBLY_FILE = "ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt";
    
    public void run(String srcPath) {
    	
    	// downloading zipped taxdump files
        System.out.println(">Downloading NCBI taxdump files");
        String taxdumpOut = "new_taxdump.out";
        Downloader.getFile(NCBI_TAXDUMP_FILE, new File(srcPath), taxdumpOut);
        new File(srcPath + File.separator + taxdumpOut).delete();
        
        // extracting zipped taxdump files
        String taxdumoZippedPath = srcPath + File.separator + "new_taxdump.tar.gz";
        File taxdump = new File(srcPath + File.separator + "taxdump");
        taxdump.mkdir();
        ProcessBuilder builder = new ProcessBuilder();
        builder.command("sh", "-c", String.format("tar xfz %s -C %s", taxdumoZippedPath, taxdump.getAbsolutePath()));
        builder.directory(new File(srcPath));
        try {
            Process process = builder.start();
            process.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
        
        // only keep the files names.dmp and nodes.dmp
        new File(taxdumoZippedPath).delete();
        for (File f : taxdump.listFiles()) {
            if (f.getName().equals("names.dmp") || f.getName().equals("nodes.dmp"))
                continue;
            f.delete();
        }
        
        // downloading assembly summary table
        System.out.println(">Downloading assembly summary table");
        new File(srcPath + File.separator + "assembly_summary_refseq.txt").delete();
        String assemblyOut = "assembly_summary_refseq.out";
        Downloader.getFile(NCBI_ASSEMBLY_FILE, new File(srcPath), assemblyOut);
        new File(srcPath + File.separator + assemblyOut).delete();


    }

}
