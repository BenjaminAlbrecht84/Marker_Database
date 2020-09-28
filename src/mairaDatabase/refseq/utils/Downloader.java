package mairaDatabase.refseq.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

public class Downloader {

    public static File getFile(String remoteFile, File localDir, String outFile) {
        File file = null;
        try {

            int returnValue = 1;
            while (returnValue != 0) {

                file = new File(localDir + File.separator + remoteFile);
                file.delete();
                File out = new File(localDir + File.separator + outFile);
                out.delete();
                 String cmd = "wget --timeout 10 " + remoteFile + " -P " + localDir;
//                String cmd = "/usr/bin/wget --timeout 10 " + remoteFile + " -P " + localDir;

                returnValue = executingCommand(cmd, out);

            }

        } catch (Exception e) {
            e.printStackTrace();
        }
        return file;

    }

    public static File getFtpFile(String ftp, String remoteFile, File localDir, String outFile) {
        File file = null;
        try {

            int returnValue = 1;
            while (returnValue != 0) {

                file = new File(localDir + File.separator + remoteFile);
                file.delete();
                File out = new File(localDir + File.separator + outFile);
                out.delete();
                String cmd = "wget " + ftp + "/" + remoteFile + " -P " + localDir;
//                String cmd = "/usr/bin/wget --timeout 10 " + ftp + "/" + remoteFile + " -P " + localDir;

                returnValue = executingCommand(cmd, out);

            }

        } catch (Exception e) {
            e.printStackTrace();
        }
        return file;

    }

    private static int executingCommand(String command, File out) {
        try {

            Runtime rt = Runtime.getRuntime();
            Process proc = rt.exec(command);

            // checking error messages
            StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), out);
            errorGobbler.start();

            int exitVal = proc.waitFor();

            FileWriter fW = new FileWriter(out);
            fW.write(errorGobbler.getMessages());
            fW.close();

            return exitVal;

        } catch (Exception e) {
            e.printStackTrace();
        }

        return 1;
    }

    private static class StreamGobbler extends Thread {
        InputStream is;
        StringBuilder build = new StringBuilder();

        StreamGobbler(InputStream is, File out) {
            this.is = is;
        }

        public void run() {
            try {
                InputStreamReader isr = new InputStreamReader(is);
                BufferedReader br = new BufferedReader(isr);
                String line = null;
                while ((line = br.readLine()) != null)
                    build.append(line + "\n");
            } catch (IOException ioe) {
                ioe.printStackTrace();
            }
        }

        public String getMessages() {
            return build.toString();
        }

    }

    private static File gunzip(File infile) {

        byte[] buffer = new byte[1024];

        try {

            GZIPInputStream gzis = new GZIPInputStream(new FileInputStream(infile));
            File outfile = new File(infile.getAbsolutePath().replaceAll(".gz", ""));
            FileOutputStream out = new FileOutputStream(outfile);

            int len;
            while ((len = gzis.read(buffer)) > 0) {
                out.write(buffer, 0, len);
            }

            gzis.close();
            out.close();

            infile.delete();
            return outfile;

        } catch (IOException ex) {
            ex.printStackTrace();
        }

        return null;
    }

}
