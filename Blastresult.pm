#!/usr/bin/perl -w

# Blastresult.pm
# Run in parallel with the Ortholog_main.pl
# Last edited on Feb 7, 2015
# Bug found on 13 Jan --> problems about hash ref when found 2 or more identical genes in genome
# Bug found on Feb 5 --> Cannot print duplicateset of other species (but correctly find duplicateset)
# All bugs are fixed, Feb 7, 2015

package Blastresult;
require Exporter;
@ISA = "Exporter";
@EXPORT = qw( 
			  getResult getOrthologSet 
			  getDuplicateSet createOrthologTable createDuplicateTable
			  writeOrthologSeq writeSummaryTable
			);

use strict;
use Data::Dumper;
use DBI;

# Subroutine: getResult
# Scrap the result from the BLAST result page and go to the subject sequences file
# to find the hit gene and store additional info (e.g. documentation and FASTA)
# Filtering of the possible hit: %identity and %coverage cutoff value (obtain from
#	command line), otherwise the default values is 30% identity and 50% coverage
# To use: $best_hit = getResult( $resultpage, $subjectfile, $identity, $coverage );
# $resultpage --> BLAST output result page (might be a file or a string variable)
# $subjectfile = name of the subject file for finding the hit gene
sub getResult {
	my ( $resultline, $subjectfile, $identity, $coverage ) = @_;
	my $hit = {};
	my $rest_blastresult;
	
	# To check :: print "$resultpage\n";
	# Extract result from Blast output, get only first line of each unique query because
	# it's the best hit, also check for the %identity cutoff value
	# If it's not passed the cutoff --> throw away and don't further consider them
	
	chomp $resultline;
	( $$hit{query}, $$hit{subject}, $$hit{identity},
	       $$hit{alignlength}, $rest_blastresult ) = split "\t", $resultline, 5;
	my @rest_blastresult = split "\t", $rest_blastresult;
	$$hit{score} = pop @rest_blastresult;
	$$hit{score} =~ s/\s+//;
	
	# If there's no assigned %identity cutoff value, use 30.
	unless ( $identity ) {
		$identity = 30;
	}
	
	if ( $$hit{identity} < $identity ) {
		return 0;
	}
	
	# Next, we will find the subject gene in the subject sequence file called $subjectfile
	# and check whether the alignment length reach over than %coverage cutoff value
	# (50 as a default value)--> If not, throw it away
	# If pass this criteria: store the sequence
	
	$hit = getSequence ( $hit, $subjectfile );
	unless ( $coverage ) {
		$coverage = 50;
	}	
	unless ( $$hit{subject_fheader} ) {
		die "ERROR!! CANNOT FOUND BEST HIT SEQUENCE\n";
	}
	elsif ( $$hit{alignlength}/length( $$hit{sequence})*100 >= $coverage ) {
		return $hit;
	}
	else {
		return 0;
	}
}

# Subroutine: getSequence
# get protein or nucleotide from subjectfile using seq ID
# To use: \%hit (or \%best_hit) = getSequence ( \%hit, $subjectfile ) 
sub getSequence {
	my ( $hit, $subjectfile ) = @_;
	open ( my $in, "<", $subjectfile )  || die "$subjectfile cannot be opened\n";
	my ( $header, $doc ) = ( undef, undef );
	while ( my $line = <$in> ) {
		chomp $line;
		my $hit_subject = $$hit{subject};
		$hit_subject =~ s/\|//g;
		if ( $line =~ /^>/ ) {
			( $header, $doc ) = split " ", $line, 2;
			$header =~ s/>|\|//g;
		}
		if ( $hit_subject eq $header ) {
			( $$hit{subject_fheader}, $$hit{doc} ) = split " ", $line, 2;
			unless ( !( $$hit{doc} ) ) {
				chomp $$hit{doc};
			}
			while ( my $line = <$in> ) {
				if ( $line =~/^>/ ) {
					chomp $$hit{sequence};
					last;
				}
				else {
					chomp $line;
					$$hit{sequence} .= $line;
				}
			}
		last;
		}
	}
	close $in;
	return $hit;
}



# Subroutine: getOrthologSet
# Store orthologset information in a complex struc called 'orthologset' 
# $orthologset = { refgene1 => [{id=>.., doc=>.., sequence=>..},{ },{ }, ..] , gene2 =>[ ], .. }
#                       orthologset gene1                                       gene2   ..
#                                  gene1 hit w/ sp.A            spB spC ..
# To use: $orthologset = getOrthologSet ( $orthologset, $resultpage
#							$subjectfile, $r_resultpage, $r_subjectfile, $opt_i, $opt_c );
{
	my $idx = 1;	# index will be used for positioning of each species
					# in $orthologset; position = 0 is used for reference sp.
					
	sub getOrthologSet {
		
		my ( $orthologset, $resultpage, $subjectfile, $r_resultpage, 
			 $r_subjectfile, $identity, $coverage ) = @_;
		my $besthit_list = [];
		my $r_besthit_list = [];
		
		# First, check whether there's any assigned cutoff value
		unless ( $identity ) {
			$identity = undef;
		}
		unless ( $coverage ) {
			$coverage = undef;
		}
		
		# First, we will get best_hit list from forward BLAST result
		# (Reference --> Subject)
		# Store result in an array of hashes called $besthit_list
		open ( my $in, "<", $resultpage ) || die "$resultpage cannot be opened\n";
		while ( my $line = <$in> ) {
			my $best_hit = getResult( $line, $subjectfile, $identity, $coverage );
			unless ( $best_hit == 0 ) {
				push @$besthit_list, $best_hit;
			}
		}
		close $in;
		
		# Then, we do the same thing for reverse BLAST result
		# (Subject --> Reference)
		# Store result in an array of hashes called $r_besthit_list
		open ( $in, "<", $r_resultpage ) || die "$r_resultpage cannot be opened\n";
		while ( my $line = <$in> ) {
			my $r_best_hit = getResult( $line, $r_subjectfile, $identity, $coverage );
			unless ( $r_best_hit == 0 ) {
				push @$r_besthit_list, $r_best_hit;
			}
		}
		close $in;
			
		# To test Data structure: print "BEST HIT\n", Dumper ( $besthit_list );
		# 						  print "R_BEST HIT\n", Dumper ( $r_besthit_list );	
	
	
		# After that, we will find the reciprocal hits and store it in the
		# $orthologset, using reference gene name as a hash key and the value
		# is an array of hash ( store id, doc and sequence of orthog genes for each species
		# for position 0 in an array is spared for reference genome
		
		foreach my $hit ( @$besthit_list ) {
			foreach my $r_hit ( @$r_besthit_list ) {
				if ( $hit->{query} eq $r_hit->{subject} &&
					 $hit->{subject} eq $r_hit->{query} ) {
					unless ( $orthologset->{$hit->{query}} ) {
						$orthologset->{$hit->{query}}->[0]->{id} = $r_hit->{subject};
						$orthologset->{$hit->{query}}->[0]->{doc} = $r_hit->{doc};
						$orthologset->{$hit->{query}}->[0]->{sequence} = $r_hit->{sequence};
					}
					$orthologset->{$hit->{query}}->[$idx]->{id} = $hit->{subject};
					$orthologset->{$hit->{query}}->[$idx]->{doc} = $hit->{doc};
					$orthologset->{$hit->{query}}->[$idx]->{sequence} = $hit->{sequence};
					last;
				}
			}
		}
		
		# Next, check for each orthog set, for the species that does not have 
		# the reciprocal hits with reference genome, we will leave the position 
		# of that species as N/A
		foreach my $each_refgene ( keys %$orthologset ) {
			foreach my $i ( 1 .. $idx ) {
				unless ( $orthologset->{$each_refgene}->[$i]->{id} ) {
					$orthologset->{$each_refgene}->[$i]->{id} = 'N/A';
				}
			}
		}
		
		$idx++;
		#To check the data structure:: print Dumper ( $orthologset );
		return $orthologset;
	}
}



# Subroutine: getDuplicateSet
# Retreive all possible duplicated genes using %identity %coverage %score cutoff value
# 	(30%, 50% and 65% as the default values, respectively
# from %Blast_result hashes store in a complex struc called '$duplicateset' (array of array of array of hashes)
# $paralogset = [ [ [ {head=>, doc=>, sequence=>}, {head=>...}, ....], [ {}, {} ], .. ],  [ [ {},{}], [{},{}] ], .. ]
#								             species1											   sp.2
#						            Dup_geneset1						Dup_geneset2
#							Dup_geneset1.1			 gene1.2
#						 hashes head doc seq
# To use: $duplicateset = getDuplicateSet( $duplicateset, $resultpage, $subjectfile
#											, $opt_i, $opt_c, $opt_s )
sub getDuplicateSet {
	my ( $duplicateset, $resultpage, $subjectfile, 
			$identity, $coverage, $score ) = @_;
	my ( $best_hit, $hit ) = ( 0, 0 );
	my $dupset_idx = 0;
	my $dupset = [];
	# $position will be used as an index for each gene inside paralog	
	my $position = 1;
	
	# First, check whether there's any defined cutoff value
	unless ( $identity ) {
		$identity = undef;
	}
	unless ( $coverage ) {
		$coverage = 50;
	}
	unless ( $score ) {
		$score = 65;
	}
	
	# Then parsing hit information from the result file
	open ( my $in, "<", $resultpage ) || die "$resultpage cannot be opened\n";
	my $first_hit = <$in>;
	$hit = getResult ( $first_hit, $subjectfile, $identity, $coverage );
	$best_hit = $hit;
	while ( my $line = <$in> ) {
		$hit = getResult ( $line, $subjectfile, $identity, $coverage );
		unless ( $hit ) {
			next;
		}
		# To check whether it is the self hit (same gene
		#	, which will have the highest BLAST score) --> to calculate %score
		if ( ( $$hit{query} eq $$hit{subject} ) && $$hit{identity} == 100.00 ) {
			$best_hit = $hit;
			if ( $position != 1 ) {
				$position = 1;
				$dupset_idx++;
			}
			next;
		}
		if ( $best_hit == 0 ) {
			next;
		}
		if ( ( $$hit{alignlength}/$$best_hit{alignlength}*100 >= $coverage )
				&& ( $$hit{score}/$$best_hit{score}*100  >= $score )
				&& ( $$hit{query} eq $$best_hit{query} ) ) {
			
			unless ( $dupset->[$dupset_idx]->[0]->{id} ) {
			# Check whether this new gene is a unique duplicateset						
				$dupset->[$dupset_idx]->[0]->{id} = $$best_hit{subject};
				$dupset->[$dupset_idx]->[0]->{doc} = $$best_hit{doc};
				$dupset->[$dupset_idx]->[0]->{sequence} = $$best_hit{sequence};		
			}
			
			$dupset->[$dupset_idx]->[$position]->{id} = $$hit{subject};
			$dupset->[$dupset_idx]->[$position]->{doc} = $$hit{doc};
			$dupset->[$dupset_idx]->[$position]->{sequence} = $$hit{sequence};
			$position++;
			next;
			#To check:	print Dumper ($best_hit);
			#			print Dumper ($hit);
		}
		$position = 1;
		# To check whether there's any content before moving to next dupset_idx
		if ( $dupset->[$dupset_idx]->[0]->{id} ) {
				$dupset_idx++;
		}
	}
	close $in;
	pop @$dupset;
	#To check the content: print Dumper ($dupset) all duplicates in single sp.;
	push @$duplicateset, $dupset;
	return $duplicateset;
}


# Subroutine: createOrthologTable
# Pick the orthologous set result and crate tables in mySQL database 
# site: (thr.genomics.purdue.edu, mysql->user)
# Tables: one table for orthologous set of each species
# and other tables storing all ortholog genes of each species (ID, doc and seq)
# Input: $orthologset	Output: database in mysql
# To use: createOrthologTable( $orthologset, $species_list, $user, $password )
sub createOrthologTable {
	my ( $orthologset, $species_list, $user, $password ) = @_;
	
	# First, we will print the table of orthologous sets in tabular format
	# Each Col: each species, Each row: each orthologous set (ortholog of each sp.)
	my @ortholog_table;
	foreach my $each_species ( @$species_list ) {
		$ortholog_table[0] .= $each_species."\t";
	}
	chop $ortholog_table[0];
	# To check the header row (Col. header) print $ortholog_table[0], "\n";
	
	my $i = 1;
	my $sp_num = @$species_list;
	foreach my $keys ( keys %$orthologset ) {
		$ortholog_table[$i] = "";
		foreach my $idx ( 0 .. ( $sp_num - 1 ) ) {
			$ortholog_table[$i] .= $orthologset->{$keys}->[$idx]->{id}."\t";
		}
		chop $ortholog_table[$i];
		$i++;
	}
	
	# To check the Tabular format: 	foreach my $line ( @ortholog_table ) {
	#									print "$line\n";
	#								}
	
	# Next, we will use this tabular result --> transform into the database
	# One table is for orthologous set table
	# First, we will connect the server: thr.genomics.purdue.edu
	# User and password come from the input from command line (-u and -p, see README)
	
	my $dsn = "DBI:mysql:$user:thr.genomics.purdue.edu:3306";
	my $dbh = DBI->connect( $dsn, $user, $password ) or dir $DBI::errstr;
	my $sth;

	# Then, we will create a tablename called 'Orthologsets' + initial names of spp.
	
	my $tablename = "Orthologsets";
	foreach my $each_species ( @$species_list ) {
		if ( $each_species =~ /.{2}/ ) {
			$tablename .= "_$&";
		}
	}
	
	# To check : print $tablename, "\n";
	
	my $newtable = "CREATE TABLE $tablename ( Number SMALLINT,";
	foreach my $each_species ( @$species_list ) {
		$newtable .= " $each_species TINYTEXT,";
	}
	
	$newtable .= " Input_time TIMESTAMP)";

	# To check correctness of the statement: print $newtable, "\n";
	$sth = $dbh->prepare( $newtable );
	$sth->execute or die "PROCESS \"CREATE NEW TABLE\" cannot be done\n";
	
	# To check whether the table's already created in SQL Database::
	#	$sth = $dbh->prepare( "show tables" );
	#	$sth->execute;
	#	while ( my @result = $sth->fetchrow_array ) {
	#		print "@result\n";
	#	}
	
	# After that, we will add content of each row in Orthologsets table
	my $newrow;
	$newrow = qq{INSERT INTO $tablename ( Number,};
	foreach my $each_species ( @$species_list ) {
		$newrow .= " $each_species,";
	}
	chop $newrow;
	$newrow .= ") VALUES ( ?,";
	foreach ( @$species_list ) {
		$newrow .= " ?,";
	}
	chop $newrow;
	$newrow .= ")";
	# To check statement: print $newrow, "\n";
	
	my $num = 1;
	foreach my $i ( 1 .. $#ortholog_table ) {
		my @each_row = split "\t", $ortholog_table[$i];
		my $new_eachrow = $newrow;
		$new_eachrow =~ s/\?/$num/;
		foreach my $i ( @each_row ) {
			$new_eachrow =~ s/\?/\"$i\"/;
		}
	# To check each statement:	print $new_eachrow, "\n";
		$sth = $dbh->prepare( $new_eachrow );
		$sth->execute or die "PROCESS \"INSERT NEW ROW\" cannot be done\n";
		$num++;
	}
	
	# To check the content in the table:
	#my $showtabcontent = "SELECT * FROM $tablename";
	#$sth = $dbh->prepare( $showtabcontent );
	#$sth->execute or die "PROCESS \"SELECT * FROM $tablename\" cannot be done\n";;
	#while ( my @result = $sth->fetchrow_array ) {
	#	print "@result\n";
	#}
	# End: Check the content in the table
	
	
	# The next step: we will create SQL tables of all sequences of each species
	my $species_number = @$species_list;
	foreach my $species_num ( 0 .. ( $species_number - 1 ) ) {
		# First, we will create each table
		$tablename = $$species_list[$species_num];
		my $newtable = "CREATE TABLE $tablename ( ID TINYTEXT, 
						Documentation TINYTEXT, Sequence TEXT, Input_time TIMESTAMP )";
		$sth = $dbh->prepare( $newtable );
		$sth->execute or die "PROCESS \"CREATE NEW TABLE\" cannot be done\n";
		
		# After that, we will add content of each row in Orthologsets table
		$newrow = qq{INSERT INTO $tablename ( ID, Documentation, Sequence ) VALUES ( ?, ?, ?) };
		$sth = $dbh->prepare( $newrow );
		foreach my $eachorth ( keys %$orthologset ) {
			if ( $$orthologset{$eachorth}->[$species_num]->{id} eq 'N/A' ) {
				next;
			}
			else {
				$sth->execute( $$orthologset{$eachorth}->[$species_num]->{id},
						   $$orthologset{$eachorth}->[$species_num]->{doc},
						   $$orthologset{$eachorth}->[$species_num]->{sequence})
						   or die "PROCESS \"INSERT NEW ROW\" cannot be done\n";
			}
		}
	}
	
	# To check whether the table's already added in SQL Database::
	#	$sth = $dbh->prepare( "show tables" );
	#	$sth->execute;
	#	while ( my @result = $sth->fetchrow_array ) {
	#		print "@result\n";
	#	}
	
	$sth->finish;
	$dbh->disconnect;
	return 1;
}


# Subroutine: createDuplicateTable
# Pick the possible duplicate gene sets and crate tables in mySQL database 
# site: (thr.genomics.purdue.edu, mysql->user)
# Tables: one table for duplicate gene set of each species
# and add the duplicate gene and info in each species table
# (additional sequences other than sequences from orthologset)
# Input: $duplicateset	Output: database in mysql
# To use: createDuplicateTable( $duplicateset, $orthologset, $species_list, $user, $password )
sub createDuplicateTable {
	my ( $duplicateset, $orthologset, $species_list, $user, $password ) = @_;
	
	# First, we will print the table of duplicate sets in tabular format
	# Each Col: each species, Each row: each orthologous set (ortholog of each sp.)
	my @duplicate_table;
	# @duplicate_table will have headers of each column as shown below
	# Species	Set_no	Query	Duplicate_genes
	# To check the header row (Col. header) print $ortholog_table[0], "\n";
	
	# Create a loop to format each row in the table using three indices
	#	$i is a sp. num index, $j is a duplicate set num index for each sp.
	#	$k is an index for each duplicate gene in the same duplicateset
	my $sp_num = @$species_list;
	foreach my $i ( 0 .. ( $sp_num - 1 ) ) {	
		my $dupset_num = @{$duplicateset->[$i]};
		foreach my $j ( 0 .. ( $dupset_num - 1 ) ) {
			my $dupgene_num = @{$duplicateset->[$i]->[$j]};
			my $eachrow;
			my $query = $duplicateset->[$i]->[$j]->[0]->{id};
			my $dupgene_list = "";
			foreach my $k ( 1 .. ( $dupgene_num - 1 ) ) {
				$dupgene_list .= $duplicateset->[$i]->[$j]->[$k]->{id}.",";
			}
			chop $dupgene_list;
			$eachrow = $$species_list[$i]."\t".($j+1)."\t".$query."\t".$dupgene_list;
			push @duplicate_table, $eachrow;
		}
	}	
	
	# To check the Tabular format: 	foreach my $line ( @duplicate_table ) {
	#									print "$line\n";
	#								}
	
	
	# Next, we will use this tabular result --> transform into the database
	# r orthologous set table
	# First, we will connect the server: thr.genomics.purdue.edu
	# User and password come from the input from command line (-u and -p, see README)
	
	my $dsn = "DBI:mysql:$user:thr.genomics.purdue.edu:3306";
	my $dbh = DBI->connect( $dsn, $user, $password ) or dir $DBI::errstr;
	my $sth;

	my $tablename = "Duplicatesets";
	
	my $newtable = qq{CREATE TABLE $tablename ( Species TINYTEXT, Set_no INT, 
					Query TINYTEXT, Duplicates TEXT, Input_time TIMESTAMP)};
	
	# To check correctness of the statement: print $newtable, "\n";
	$sth = $dbh->prepare( $newtable );
	$sth->execute or die "PROCESS \"CREATE NEW TABLE\" cannot be done\n";
	
	my $newrow = qq{INSERT INTO $tablename ( Species, Set_no, Query, Duplicates ) 
					VALUES ( ?, ?, ?, ?) };
	$sth = $dbh->prepare( $newrow );
	
	foreach my $eachrow ( @duplicate_table ) {
		my ( $species, $set_no, $query, $duplicates ) = split "\t", $eachrow;
			$sth->execute( $species, $set_no, $query, $duplicates )
				or die "PROCESS \"INSERT NEW ROW\" cannot be done\n";
	}
	
		
	# Next we will add the extra duplicate gene info into a species name table
	# So we have to check whether each gene in $duplicateset already exists in
	# $orthologset, if it is, just skipped it, otherwise it will be added into the
	# table.
	# $i is an index for species
		
	foreach my $i ( 0 .. ( $sp_num - 1 ) ) {	
		$tablename = $$species_list[$i];
		$newrow = qq{INSERT INTO $tablename ( ID, Documentation, Sequence ) VALUES ( ?, ?, ?) };
		$sth = $dbh->prepare( $newrow );
		my @newseq_list = ( 'id' );
		foreach my $each_dupset ( @{$duplicateset->[$i]} ) {
			foreach my $each_dupgene ( @$each_dupset ) {
				my $repeat;			# $repeat used to check whether the seq already exists
									# in the sp. table/ or already added into the table
				foreach my $eachorth ( keys %$orthologset ) {
					foreach my $newseq ( @newseq_list ) {
						if ( $each_dupgene->{id} eq $newseq ) {
							$repeat = 1;
							last;
						}
					}
					if ( $repeat ) {
						last;
					}
					elsif ( $each_dupgene->{id} eq $orthologset->{$eachorth}->[$i]->{id} ) {
						$repeat = 1;
						last;
					}
				}
				if ( $repeat ) {
					next;
				}
				$sth->execute( $each_dupgene->{id}, $each_dupgene->{doc},
						$each_dupgene->{sequence})
						or die "PROCESS \"INSERT NEW ROW\" cannot be done\n";
				push @newseq_list, $each_dupgene->{id};
			}
		}
	}
		
	$dbh->disconnect;
	return 1;
}



# Subroutine: writeOrthologSeq
# Write orthologous sequences for each orthologous set
# genes that are in the same ortholog sets will have the same prefix tag as ref gene of that
# orthologset, and also write possible duplicate genes if it has
# Write the result in the .fasta file output
# To use: writeOrthologSeq( $orthologset, $duplicateset, $species_list, $dirhandle ) 
sub writeOrthologSeq {
	
	my ( $orthologset, $duplicateset, $species_list, $dirhandle ) = @_;
	open (my $SEQFASTA, '>', "$dirhandle/Ortholog_sets_sequences.fasta"); 
	foreach my $key (keys %$orthologset) {
	
		foreach my $idx ( 0 .. ( @$species_list - 1) ) {
			unless ($orthologset->{$key}->[$idx]->{id} eq "N/A") {
				my $header = ">"."$key\_$orthologset->{$key}->[$idx]->{id}";
				my $sequence = 	"$orthologset->{$key}->[$idx]->{sequence}";
				my $length = "80";
		
				print $SEQFASTA "$header\n";
	
				for ( my $pos = 0 ; $pos < length($sequence) ; $pos += $length ) {
					print $SEQFASTA substr($sequence, $pos, $length), "\n";
        		}
				
				# Then, print all possible duplicate genes
				my $dupset_num = @{$duplicateset->[$idx]};
				foreach my $j ( 0 .. ( $dupset_num - 1)) {
					my $dupgene_num = @{$duplicateset->[$idx]->[$j]};
					if ( $orthologset->{$key}->[$idx]->{id} eq $duplicateset->[$idx]->[$j]->[0]->{id} ) {
						foreach my $l ( 1 .. ($dupgene_num - 1 )) {
							unless ( $orthologset->{$key}->[$idx]->{id} eq $duplicateset->[$idx]->[$j]->[$l]->{id} ) {
								$header = ">"."$key\_$duplicateset->[$idx]->[$j]->[$l]->{id}";
								$sequence = $duplicateset->[$idx]->[$j]->[$l]->{sequence};
								
								print $SEQFASTA "$header\n";
	
								for ( my $pos = 0 ; $pos < length($sequence) ; $pos += $length ) {
									print $SEQFASTA substr($sequence, $pos, $length), "\n";
								}	
							}
						}
						last;	
					}
				}	
				
			}
		}
		print $SEQFASTA "\n";
	
	}

	close ($SEQFASTA);
	return 1;
}


# Subroutine: writeSummaryTable
# Write Summary table: orthologous sequences for each species and all possible
# duplicated genes in each species
# Write the result in the .txt file output
# To use: writeSummaryTable( $orthologset, $duplicateset, $species_list, $dirhandle ); 
sub writeSummaryTable {
	my ( $orthologset, $duplicateset, $species_list, $dirhandle ) = @_;
	open (my $TABLETXT, '>', "$dirhandle/Summary_Table.txt"); 
	print $TABLETXT "----------ORTHOLOGOUS GROUPS---------\n";
	
	my $spp_columns = "";
	foreach my $idx ( 0 .. ( @$species_list - 1 ) ) {
		$spp_columns .= "$$species_list[$idx]\t$$species_list[$idx]\_duplicates\t";
	}
	print $TABLETXT "$spp_columns\n";

	foreach my $key (keys %$orthologset) {
		foreach my $idx ( 0 .. ( @$species_list - 1 ) ) {
			my $duplicate;
			print $TABLETXT "$orthologset->{$key}->[$idx]->{id}\t";
			my $dupset_num = @{$duplicateset->[$idx]};
			foreach my $j ( 0 .. ( $dupset_num - 1)) {
				my $dupgene_num = @{$duplicateset->[$idx]->[$j]};
				if ( $orthologset->{$key}->[$idx]->{id} eq $duplicateset->[$idx]->[$j]->[0]->{id} ) {
					my $dupgenes = "";
					foreach my $l ( 1 .. ($dupgene_num - 1 )) {
						unless ( $orthologset->{$key}->[$idx]->{id} eq $duplicateset->[$idx]->[$j]->[$l]->{id} ) {
							$dupgenes .= $duplicateset->[$idx]->[$j]->[$l]->{id}.",";
						}
					}
					chop $dupgenes;
					print $TABLETXT "$dupgenes\t";
					$duplicate = 1;
					last;
				}
			}
			unless ( $duplicate ) {
				print $TABLETXT "\t";
			}
		} 
		print $TABLETXT "\n";
	}
}

1;
