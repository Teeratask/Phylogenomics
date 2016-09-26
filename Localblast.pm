#!/usr/bin/perl

package Localblast;
require Exporter;
@ISA = "Exporter";
@EXPORT = qw( LocalBlast );

use strict ;

# Format the databases for blast search
#sub LocalBlast
# Make a local BLAST (blastp program) in the PC between sequences
# To use: LocalBlast ( $inputfile_list, $species_list )
sub LocalBlast {
	my ( $inputfile_list, $species_list ) = @_;
	my $file_number = @$inputfile_list;
	foreach my $i ( @$inputfile_list ) {
		system("makeblastdb -in $i -dbtype prot");
	}

	foreach my $i ( 1 .. ($file_number - 1) ) {
		system("blastp -query $$inputfile_list[$i] -db $$inputfile_list[0] -evalue 1e-6 -out blastp.$$species_list[$i].$$species_list[0].out -use_sw_tback -max_target_seqs 1 -seg no -soft_masking true -outfmt 6");
		system("blastp -query $$inputfile_list[0] -db $$inputfile_list[$i] -evalue 1e-6 -out blastp.$$species_list[0].$$species_list[$i].out -use_sw_tback -max_target_seqs 1 -seg no -soft_masking true -outfmt 6");
		system("blastp -query $$inputfile_list[$i] -db $$inputfile_list[$i] -evalue 1e-20 -out blastp.$$species_list[$i].$$species_list[$i].out -use_sw_tback -max_target_seqs 200 -seg no -soft_masking true -outfmt 6");
	}
		system("blastp -query $$inputfile_list[0] -db $$inputfile_list[0] -evalue 1e-20 -out blastp.$$species_list[0].$$species_list[0].out -use_sw_tback -max_target_seqs 200 -seg no -soft_masking true -outfmt 6");
}
1;