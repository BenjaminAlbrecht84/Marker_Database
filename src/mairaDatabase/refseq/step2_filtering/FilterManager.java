package mairaDatabase.refseq.step2_filtering;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;

import mairaDatabase.refseq.utils.Formatter;
import mairaDatabase.refseq.utils.aliHelper.SQLAlignmentDatabase;
import mairaDatabase.utils.FileUtils;
import mairaDatabase.utils.ResourceLoader;
import mairaDatabase.utils.SQLMappingDatabase;
import mairaDatabase.utils.taxTree.TaxTree;

public class FilterManager {
	
	private ResourceLoader rL = new ResourceLoader();
	
    private BufferedWriter factorWriter, markerWriter;
    private TaxTree taxTree;
    private int NUM_OF_PROTEINS, MIN_ID;
    private String aliFolder;
    private File tmpDir, weightFile;

    public void run(String rank, String srcPath, String aliFolder, File tmpDir, File markerDir, TaxTree taxTree, SQLMappingDatabase mappingDatabase, int NUM_OF_PROTEINS, int MIN_ID, int cores) {

        this.taxTree = taxTree;
        this.NUM_OF_PROTEINS = NUM_OF_PROTEINS;
        this.MIN_ID = MIN_ID;
        this.aliFolder = aliFolder;
        this.tmpDir = tmpDir;
        rL.setTime();

        try {
        	
        	File markerDatabase = new File(srcPath+File.separator+"marker_db");
            File markerFile = new File(markerDatabase.getAbsolutePath() + File.separator + "genus_marker_db_" + NUM_OF_PROTEINS + ".faa");
            markerFile.delete();
            markerWriter = new BufferedWriter(new FileWriter(markerFile, true));
           
            weightFile = new File(srcPath + File.separator + "genus_marker_db_" + NUM_OF_PROTEINS + "_weights.tab");
            weightFile.delete();
            factorWriter = new BufferedWriter(new FileWriter(weightFile, true));

            System.out.println(">Filtering marker proteins");
            File[] faaFiles = markerDir.listFiles();
            ArrayList<Runnable> threads = new ArrayList<>();
            for (File faaFile : faaFiles)
                threads.add(new FilterThread(faaFile, mappingDatabase));
            rL.runThreads(cores, threads, threads.size());
  
            markerWriter.close();
            factorWriter.close();

        } catch (Exception e) {
            e.printStackTrace();
        }
        
        FileUtils.deleteDirectory(markerDir.getAbsolutePath());

    }

    private class FilterThread implements Runnable {

        private File faaFile;
        private SQLMappingDatabase mappingDatabase;
        private SQLAlignmentDatabase alignmentDatabase;
        private String genus;

        public FilterThread(File faaFile, SQLMappingDatabase mappingDatabase) {
            this.faaFile = faaFile;
            this.mappingDatabase = createMappingDatabase(mappingDatabase);
            this.genus = Formatter.removeNonAlphanumerics(faaFile.getName().replaceAll("_marker\\.faa", ""));
            this.alignmentDatabase = createAlignmentDatabase(genus);
        }

        @Override
        public void run() {
            new Filtering().run(faaFile, genus, factorWriter, markerWriter, taxTree, mappingDatabase, alignmentDatabase, NUM_OF_PROTEINS, MIN_ID);
            mappingDatabase.close();
            alignmentDatabase.close();
            rL.reportProgress(1);
            rL.countDown();
        }
    }

    private synchronized SQLMappingDatabase createMappingDatabase(SQLMappingDatabase mappingDatabase) {
        return new SQLMappingDatabase(mappingDatabase);
    }

    private synchronized SQLAlignmentDatabase createAlignmentDatabase(String genus) {
        return new SQLAlignmentDatabase(aliFolder, genus, tmpDir);
    }

	public File getWeightFile() {
		return weightFile;
	}

}
