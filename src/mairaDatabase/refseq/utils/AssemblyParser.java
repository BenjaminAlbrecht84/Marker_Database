package mairaDatabase.refseq.utils;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class AssemblyParser {

	public Map<Integer, ArrayList<String>> getTaxidToFTP(String src) {

		HashMap<Integer, ArrayList<String>> taxidToFTP = new HashMap<>();
		try {
			BufferedReader buf = new BufferedReader(
					new FileReader(new File(src + File.separator + "assembly_summary_refseq.txt")));
			String l;
			while ((l = buf.readLine()) != null) {
				if (!l.startsWith("#")) {
					String[] split = l.split("\t");
					int taxID = Integer.parseInt(split[5]);
					String ftp = split[19];
					if (ftp.startsWith("ftp")) {
						taxidToFTP.putIfAbsent(taxID, new ArrayList<String>());
						taxidToFTP.get(taxID).add(ftp);
					}
				}
			}
			buf.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return taxidToFTP;
	}

	public Map<String, Integer> getGcf2Taxid(String src) {

		HashMap<String, Integer> gcfToTaxID = new HashMap<>();
		try {
			BufferedReader buf = new BufferedReader(
					new FileReader(new File(src + File.separator + "assembly_summary_refseq.txt")));
			String l;
			while ((l = buf.readLine()) != null) {
				if (!l.startsWith("#")) {
					String[] split = l.split("\t");
					int taxID = Integer.parseInt(split[5]);
					String gcfAcc = split[0];
					gcfToTaxID.put(gcfAcc, taxID);
				}
			}
			buf.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return gcfToTaxID;
	}

	public Map<String, String> getGcfToYear(String src) {

		HashMap<String, String> gcfToYear = new HashMap<>();
		try {
			BufferedReader buf = new BufferedReader(
					new FileReader(new File(src + File.separator + "assembly_summary_refseq.txt")));
			String l;
			while ((l = buf.readLine()) != null) {
				if (!l.startsWith("#")) {
					String[] split = l.split("\t");
					String gcf = split[0];
					String year = split[14].split("/")[0];
					gcfToYear.put(gcf, year);
				}
			}
			buf.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return gcfToYear;
	}
	
	public Map<String, String> getGcfToAssembyLevel(String src) {

		HashMap<String, String> gcfToAssemblyLevel = new HashMap<>();
		try {
			BufferedReader buf = new BufferedReader(
					new FileReader(new File(src + File.separator + "assembly_summary_refseq.txt")));
			String l;
			while ((l = buf.readLine()) != null) {
				if (!l.startsWith("#")) {
					String[] split = l.split("\t");
					String gcf = split[0];
					String assemblyLevel = split[11];
					gcfToAssemblyLevel.put(gcf, assemblyLevel);
				}
			}
			buf.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return gcfToAssemblyLevel;
	}

}
