#!/usr/bin/perl
use strict;
use Bio::SearchIO;

($#ARGV == 0) || die "\nUsage: $0 <blast-file>\n\n\tThis script will create a tab delimited summary file of a blast file with the following information:\n\tQuery\tQuery length\thit\thit length\t% identity\te-value\talignment length\tquery start\tquery end\thit start\thit end\thit description\n\tNote that hit start and hit end are reversed if alignment is on the opposite strand.\n\n";

my $in = new Bio::SearchIO(-file => $ARGV[0]);

while(my $result = $in->next_result) {
	my $hit = $result->next_hit;
	if(!defined($hit)) {
		print $result->query_name, "\t", $result->query_length, "\t*** No hits found ***\n";
		next;
	}
	do {
		while(my $hsp = $hit->next_hsp) {
			my ($hit_start, $hit_end) = ($hsp->hit->strand == 1)? ($hsp->hit->start, $hsp->hit->end) : ($hsp->hit->end, $hsp->hit->start);
			print $result->query_name, "\t", $result->query_length, "\t", $hit->name, "\t", $hit->length, "\t", int(100*$hsp->percent_identity)/100, "\t", $hsp->evalue, "\t",
				$hsp->length, "\t", $hsp->query->start, "\t", $hsp->query->end, "\t$hit_start\t$hit_end\t", $hit->description, "\n";
		}
	} while($hit = $result->next_hit);
}
