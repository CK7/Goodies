#!/usr/bin/perl
#############################################################################################################
# This script lists the positions of all Ns in the input fasta file
#
# Written by Itai Sharon (itai.sharon@gmail.con)
#############################################################################################################

use strict;

my $IN_ref = \*STDIN;
if($#ARGV == 0) {
	open(IN, $ARGV[0]) || die "\nCannot read $ARGV[0]\n\n";
	$IN_ref = \*IN;
}

my ($header, $seq) = (undef, undef);

while(<$IN_ref>) {
	if($_ =~ /^>(\S+)/) {
		my $nheader = $1;
		if(defined($header)) {
			print "\n* $header (", length($seq), " bps) *\n";
			$seq =~ tr/acgtn/ACGTN/;
			while($seq =~ /(N+)/g) {
				my $subseq = $1;
				my $type = ($subseq =~ /^N+$/)? "N-segment" : "Non-DNA";
				my ($start, $end) = ((pos($seq)-length($subseq)+1), pos($seq));
				print "$start-$end\t$type, ", ($end-$start+1), " bps\n";
			}
		}
		($header, $seq) = ($nheader, '');
	}
	else {
		chomp;
		$seq .= $_;
	}
}

if(defined($header)) {
	print "\n* $header (", length($seq), " bps) *\n";
	$seq =~ tr/acgtn/ACGTN/;
	while($seq =~ /(N+)/g) {
		my $subseq = $1;
		my $type = ($subseq =~ /^N+$/)? "N-segment" : "Non-DNA";
		my ($start, $end) = ((pos($seq)-length($subseq)+1), pos($seq));
		print "$start-$end\t$type, ", ($end-$start+1), " bps\n";
	}
}


close($IN_ref) if($ARGV == 0);
