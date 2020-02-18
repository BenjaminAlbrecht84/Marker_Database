package utils;

import java.io.*;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipException;

public class FastaReader {

    public static ArrayList<Object[]> read(File fastaFile) {

        ArrayList<Object[]> readInfo = new ArrayList<Object[]>();
        HashSet<SparseString> readNames = new HashSet<SparseString>();
        try {

            BufferedReader buf;
            try {
                buf = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fastaFile))));
            } catch (ZipException | EOFException e) {
                buf = new BufferedReader(new FileReader(fastaFile));
            }

            String line, id = "";
            boolean readSequence = false;
            StringBuilder seq = new StringBuilder("");
            while ((line = buf.readLine()) != null) {
                if (line.startsWith("@") || line.startsWith(">")) {
                    if (seq.length() != 0 && !id.isEmpty()) {
                        Object[] o = {new SparseString(id), new SparseString(seq.toString()), seq.length()};
                        if (checkFASTFile(o, readNames))
                            readInfo.add(o);
                    }
                    seq = new StringBuilder("");
                    id = line.substring(1).split(" ")[0];
                    readSequence = true;
                } else if (line.startsWith("+")) {
                    readSequence = false;
                } else if (readSequence) {
                    seq.append(line);
                }

            }
            if (seq.length() != 0 && !id.isEmpty()) {
                Object[] o = {new SparseString(id), new SparseString(seq.toString()), seq.length()};
                if (checkFASTFile(o, readNames))
                    readInfo.add(o);
            }
            buf.close();

        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
        readNames = null;

        return readInfo;

    }

    private static boolean checkFASTFile(Object[] o, HashSet<SparseString> readNames) {
        if (readNames.contains(o[0])) {
            System.err.println("WARNING: read " + o[0].toString() + " occurs multiple times in FASTA file.");
            return false;
        }
        readNames.add((SparseString) o[0]);
        return true;
    }

}
