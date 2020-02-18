package startUp;

import refseq.RefseqManager;

import java.io.File;

public class Main {

    public static void main(String[] args) {

        if (args.length != 5) {
            System.out.println("SRC|TMP|DB|CORES|MEMORY");
            return;
        }

        File src = new File(args[0]);
        File tmp = new File(args[1]);
        String aliFolder = args[2];
        int cores = Integer.parseInt(args[3]);
        int memory = Integer.parseInt(args[4]);
        new RefseqManager().run(src, tmp, aliFolder, cores, memory);

    }

}
