package mairaDatabase.main;

import java.io.File;

import jloda.fx.util.ArgsOptions;
import jloda.fx.util.ResourceManagerFX;
import jloda.util.Basic;
import jloda.util.CanceledException;
import jloda.util.PeakMemoryUsageMonitor;
import jloda.util.ProgramProperties;
import jloda.util.UsageException;
import mairaDatabase.refseq.RefseqManager;

public class MairaDatabase {

	public static void main(String[] args) {
		try {
			ResourceManagerFX.addResourceRoot(MairaDatabase.class, "mairaDatabase.resources");
			ProgramProperties.setProgramName("MairaDatabase");
            ProgramProperties.setProgramVersion(Version.SHORT_DESCRIPTION);
            Basic.setDebugMode(true);
            
            PeakMemoryUsageMonitor.start();
            (new MairaDatabase()).run(args);
            System.err.println("Total time:  " + PeakMemoryUsageMonitor.getSecondsSinceStartString());
            System.err.println("Peak memory: " + PeakMemoryUsageMonitor.getPeakUsageString());
            System.exit(0);
		} catch (Exception ex) {
			Basic.caught(ex);
			System.exit(1);
		}
	}

	private void run(String[] args) throws CanceledException, UsageException {
		
        final ArgsOptions options = new ArgsOptions(args, this, "MAIRA Database - automatically computes a complete current version of the MAIRA database.");
        options.setVersion(ProgramProperties.getProgramVersion());
        options.setLicense("Copyright (C) 2020 Algorithms in Bioinformatics, University of Tuebingen. This program comes with ABSOLUTELY NO WARRANTY.");
        options.setAuthors("Benjamin Albrecht");
        options.setCommandMandatory(true);
        
        String srcPath;
		String tmpPath;
		String aliPath;
		int cores;
		int memory;
		String generaInput;
		String diamondBin;
		
		srcPath = options.getOptionMandatory("-f", "srcFolder", "Folder for the computed database files", "");
		tmpPath = options.getOptionMandatory("-t", "tmpFolder", "Temporary folder used for speeding-up the computation (eg. /dev/shm)", "");
		aliPath = options.getOptionMandatory("-a", "aliFolder", "Alignment folder capturing alignment information", "");
		cores = options.getOption("-p", "threads", "Number of threads (default all)", Runtime.getRuntime().availableProcessors());
		memory = options.getOption("-m", "memory", "Available working memory in GB (default max)", (int) Math.round(Runtime.getRuntime().maxMemory() / Math.pow(10,9)));
		generaInput = options.getOption("-g", "genera", "Genera to be considered in the database (default all)", "");
		diamondBin = options.getOption("-d", "diamond", "Path of DIAMOND binary", "diamond");
		
		File src = new File(srcPath);
		File tmp = new File(tmpPath);
		String[] genera = generaInput.isEmpty() ? null : generaInput.trim().split(",");
		new RefseqManager().run(src, tmp, aliPath, cores, memory, genera, diamondBin);
		 
	}

}
