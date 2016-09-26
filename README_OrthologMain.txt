README FILE
################################################################################################
Get_Ortholog_By_RBH (Get Orthologous proteins By Reciprocal Best Hits)

Script authors: 	Teeratas Kijpornyongpan and Lucio Navarro
			Purdue University
			7 February 2015

# For the impatient, this is the simplest command line to run the tool:

	perl Ortholog_main.pl -r <reference.fasta> <sp1.fasta> <sp2.fasta> . . . <spN.fasta>

However, in the "How to run the program" section you will find more details about the
program options and parameters.

################################################################################################

### Description ###

This tool was created to predict orthologous proteins/genes across two or more genomes.
Get_Ortholog_By_RBH may use a well annotated protein set as a reference genome to find 
possible ortholog sequences among multiple species. The tool can be also used for simply 
predicting orthologs between just two species. Additionally, duplicated or paralogous 
sequences may be also predicted during the analysis.


### The algorithm ###

The tool relies on BLAST-based methods for its purpose, specifically on Reciprocal Best 
Hits (RBH) from blastp searches between every genome and the reference (for orthologs) and 
within each genome (for possible paralogs). 

The analysis involve three basic steps:

	1- BLAST searches: reciprocal blastp are run between each genome of interest and the 
	reference genome, and within each genome. The between-genome BLAST only reports the best 
	target sequence hit in the output (option -max_target_seqs is set to 1 in the blastp),
	while the within-BLAST reports up to 200 targets (-max_target_seqs is set to 200).
	 
	2- BLAST output filtering: blastp hits with alignment identity >= 30% and protein 
	coverage >= 50% for both the query and subject are selected and used in the next step.
	The minimum alignment identity and minimum coverage parameters can be modified by options 
	-i and -c respectively. See below in "How to run the program". 
	
	3- Reciprocal best hit identification: After filtering, blast searches whose query and 
	subjects are reciprocally found as best hits between the genome of interest and the 
	reference are selected as ortholog sequences. In parallel, using the within-genome BLAST
	searches the self best hit proteins are filtered out and any hit with score >=65% of
	the original self best hit score are collected as possible paralogs. The score for 
	selecting paralogs can be modified by option -s. See below in "How to run the program".


### The inputs ###

The input files must be protein sequences in fasta format. One file for each genome is 
placed in the same folder containing the main tool and the packages. The file names must 
follow the format "speciesname.fasta". We recommend to use short species denotation in the 
file names in order to facilitate protein set identification during the analysis and the 
final outputs. A simple example would be: Ecoli.fasta, Dmel.fasta, Athaliana.fasta.


### The outputs ###

After completing the analysis, a folder called "results" will be created. Here you can find
a "Summary_Table.txt" showing each predicted ortholog set (using reference protein IDs as 
identifiers) and the ortholog and possible paralog protein IDs found for each species. Another
file called "Ortholog_sets_sequences.fasta" will contain the protein sequences in fasta for
members of each ortholog set (paralog sequences are not included in this fasta file). The
protein IDs of each ortholog set in the fasta file are identified using the common reference
protein ID followed by the species protein IDs.

Additionally, if you are using a SQL database, tablenames "Orthologsets" and "Duplicatesets"
will be created as well as tables for each species or genomes with protein ID, additional
documentation (if any), and sequences. SQL database can be enabled with option "-q 1", which
connects the mySQL server at "thr.genomics.purdue.edu". Since a username and password are 
required to connect the mySQL server you will need to provide this information with options
-u <username> and -p <password>. 


### How to run the program ###

	perl Ortholog_main.pl -r <reference.fasta> [-i %cutoffidentity -c %cutoffcoverage 
	-s %cutoffscore -q 1 -u <user> -p <password>] <sp1.fasta> <sp2.fasta> ... <spN.fasta>
	
	Program options:
	
	-r		reference species genome set.
	
	-i		%identity cutoff value (Default %identity >= 30) during BLAST output
			filtering.
				
	-c		%coverage cutoff value (Default %coverage > 50) during BLAST output 
			filtering. Higher values may reduce the possibilities of finding more 
			orthologs.
	
	-s		%score cutoff value (Default %score > 65) for selecting paralogs. Higher
			values may select for more similar or closer sequences (recent paralogs?).

	-q 1		use "-q 1" to enable the connection to the mySQL sever (thr.genomics.purdue.edu)
			By default the SQL connection is disabled.

	-u		username for connecting mySQL server. Use this option only with "-q 1". 

	-p		password for connecting mySQL server. Use this option only with "-q 1".
	

### External tools required ###

- BLAST+ (can be downloaded from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- Perl

# Optional tools:
- DBI Perl package
- SQL database (user and password to connect the server "thr.genomics.purdue.edu").

################################################################################################



