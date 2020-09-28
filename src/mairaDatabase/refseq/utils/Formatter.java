package mairaDatabase.refseq.utils;

public class Formatter {
	
	public static String removeNonAlphanumerics(String genus) {
		return genus.replaceAll("[^A-Za-z0-9]", "");
	}
	
}
