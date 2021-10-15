module MAIRA_database {
	
	exports mairaDatabase.refseq.step2_filtering;
	exports mairaDatabase.refseq.utils.aliHelper;
	exports mairaDatabase.utils.taxTree;
	exports mairaDatabase.main;
	exports mairaDatabase.utils;
	exports mairaDatabase.refseq.step0_downloading;
	exports mairaDatabase.refseq.utils.tab;
	exports mairaDatabase.refseq;
	exports mairaDatabase.refseq.utils;
	exports mairaDatabase.refseq.step1_clustering;

	requires java.logging;
	requires transitive java.sql;
	requires transitive jloda;
	
}