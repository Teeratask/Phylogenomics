#!/usr/bin/perl -w
use strict;
use Blastresult;
use Localblast;
use Data::Dumper;
use Getopt::Std;

our ( $opt_b, $opt_o, $opt_r, $opt_i, $opt_c, $opt_s, $opt_q, $opt_u, $opt_p );
getopt ( 'b:o:r:i:c:s:u:p:q' );
#------------------------------------------------------------------------------#
#The main program#
# Teeratas Kijpornyongpan and Lucio Navarro
# 30 April 2014
#
# Input: Sequences in .fasta file of each species 
#									(reference species and comparing spp.)
# How to use program: perl program.pl -r referencesp._file
#					   	-i %cutoffidentity -c %cutoffcoverage -s %cutoffscore
#						[-q 1 -u username -p password] sp.1_file sp.2_file ...
#	Note: [] --> If want to create a result in mySQL database 
# Program options: -b --> enable Localblast (1 for enable localBlast)
#				   -o --> enable finding ortholog seqs
#				   -r --> reference sp. genome
#				   -i --> %identity cutoff value in getResult
#						  Default %identity > 30%
#				   -c --> %coverage cutoff value in getResult+getDuplicateset
#						  Default %coverage > 50%
#				   -s --> %score cutoff value in getDuplicateset 
#						  Default %score > 65%
#				   -q --> option if want to store the result in mySQL database
#						  at the server 'thr.genomics.purdue.edu'
#						   (use 1 for yes, do not fill it for no)
#				-u, -p -->put username and password to connect the SQL server
#						  (if tying -q 1 to create a result in the database)
#------------------------------------------------------------------------------#

# To check whether options are received from command line:
#print "r = $opt_r\ni = $opt_i\nc = $opt_c\ns = $opt_s\nq = $opt_q\nu = 
#$opt_u\np = $opt_p\n";

my ( @inputfile_list, @species_list );

# store reference species inputfile and species name as the first position [0] in the arrays
push @inputfile_list, $opt_r;
my ( $name, $junk ) = split ".fast", $opt_r, 2;
push @species_list, $name;

foreach my $filename ( @ARGV ) {
	push @inputfile_list, $filename;
	( $name, $junk ) = split ".fast", $filename, 2;
	push @species_list, $name;
}

# To check the content:: 	
#	print "Input list:\n";
#	foreach my $i ( @inputfile_list ) {
#		print $i, "\n";
#	}
#	print "Species list:\n";
#	foreach my $i ( @species_list ) {
#		print $i, "\n";
#	}

if ( $opt_b ) {
	LocalBlast ( \@inputfile_list, \@species_list );
}

# Make a loop of find ortholog set for eacg pair (reference --> each_sp.)
# $resultpage = result from Reference --> Subject
# $subjectfile = fasta seq file of each subject species
# $r_resultpage = result from Subject --> Reference
# $r_subjectfile = fasta seq file of reference species

my $orthologset = {};
my $duplicateset = [];

if ( $opt_o ) {
	foreach my $i ( 1 .. $#species_list ) {
		my $resultpage = "blastp\.".$species_list[0]."\.".$species_list[$i]."\.out";
		my $r_resultpage = "blastp\.".$species_list[$i]."\.".$species_list[0]."\.out";
		my $subjectfile = $inputfile_list[$i];
		my $r_subjectfile = $inputfile_list[0];
		$orthologset = getOrthologSet ( $orthologset, $resultpage,
						$subjectfile, $r_resultpage, $r_subjectfile, $opt_i, $opt_c );
	}
	# To check $orthologset content:: print Dumper ( $orthologset );
	
	foreach my $i ( 0 .. $#species_list ) {
		my $resultpage = "blastp\.".$species_list[$i]."\.".$species_list[$i]."\.out";
		my $subjectfile = $inputfile_list[$i];
		$duplicateset = getDuplicateSet ( $duplicateset, $resultpage, $subjectfile, 
										$opt_i, $opt_c, $opt_s );
	}
	print Dumper ( $duplicateset );
	# To check $duplicateset content:: print Dumper ( $duplicateset );


	# Then, we will create output files in the new directory called results
	# The result files are: Ortholog set sequences file, and the summary table
	
	my $dirhandle = "./results";
	mkdir $dirhandle;
	
	writeOrthologSeq( $orthologset, $duplicateset, \@species_list, $dirhandle );
	writeSummaryTable( $orthologset, $duplicateset, \@species_list, $dirhandle );
}	
# If you key options to create the result in mySQL database
if ( $opt_q ) {
	createOrthologTable( $orthologset, \@species_list, $opt_u, $opt_p );
	createDuplicateTable( $duplicateset, $orthologset, \@species_list, $opt_u, $opt_p );
}

exit 0;
