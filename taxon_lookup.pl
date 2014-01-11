#!/usr/bin/perl
#############################################################################################################################################
# taxon_lookup.pl
# Retrieves the taxonomy of any taxon at any taxonomic level. Works in both a single organism mode (-s) and batch mode with multiple
# taxa listed in a file, one organism per line (-i). Names must be identical to those in NCBI's database.
# This script requires NCBI's taxdump db which can be downloaded from here: ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
# Change the paths below to the correct path of the files.
#
# Written by Itai Sharon, itai.sharon@gmail.com (08/25/2013)
#############################################################################################################################################
use strict;

my $nodes_dmp = '/work/bio_db/taxdb/taxdump/nodes.dmp';
my $names_dmp = '/work/bio_db/taxdb/taxdump/names.dmp';

(-e $nodes_dmp) || die "\nCannot find nodes dump file under $nodes_dmp. You may need to update the path for variable \$nodes_dmp\n\n";
(-e $names_dmp) || die "\nCannot find names dump file under $names_dmp. You may need to update the path for variable \$names_dmp\n\n";

my %id2info = ();
my %lookup = ();

($#ARGV == 1) || die "\nUsage: $0 {-s <sequence-name> | -i <sequence-list-file>}" . 
			"\nOutput is in the format:\n<Taxon name>\t:\tRoot > Superkingdom > Phylum > Class > Order > Family > Genus > species\n\n";
if($ARGV[0] eq '-s') {
	$lookup{$ARGV[1]} = undef;
}
elsif($ARGV[0] eq '-i') {
	open(IN, $ARGV[1]) || die "\nCannot read $ARGV[1]\n\n";
	while(<IN>) {
		chomp;
		$lookup{$_} = undef;
	}
}
close(IN);

open(IN, $names_dmp) || die "\nCannot read $names_dmp\n\n";
while(<IN>) {
	# 1218	|	Prochlorococcus	|		|	scientific name	|
	my ($id, $name) = split(/\s+\|\s+/);
	$lookup{$name} = $id if(exists($lookup{$name}));
	$id2info{$id}{'name'} = $name if($_ =~ /scientific name/);
}
close(IN);

open(IN, $nodes_dmp) || die "\nCannot read $names_dmp\n\n";
while(<IN>) {
	# 1219	|	1218	|	species	|	PM	| ..
	my ($id, $parent_id, $taxonomy) = split(/\s+\|\s+/);
	$id2info{$id}{'parent'} = $parent_id;
	$id2info{$id}{'taxonomy'} = $taxonomy;
}
close(IN);

foreach my $org (keys %lookup) {
	my %hier = ();
	my $id = $lookup{$org};
	if(!defined($id)) {
		print "$org\t:\tNO RESULTS\n";
		next;
	}
	$hier{$id2info{$id}{'taxonomy'}} = $id;
	while(($id = $id2info{$id}{'parent'}) != 1) {
		$hier{$id2info{$id}{'taxonomy'}} = $id;
	}
	print "$org\t:\tRoot";

	foreach my $tax(reverse('species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom')) {
		if(exists($hier{$tax})) {
			print " > ", $id2info{$hier{$tax}}{'name'};
		}
		else {
			print " > N/A";
		}
	}
	print "\n";
}
