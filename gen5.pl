#!/usr/bin/perl
use strict;
use Bio::SeqIO;
use Math::Random::OO::Normal;

my $seed = 42;

(@ARGV >= 6) || die "\n$0 is a very simple NGS PE-read generator with no error models etc.\n\nUsage: $0 <fasta-file> <out-prefix> <coverage> <insert-size> <insert-size_sd> <read-length> [<seed>]\n\n";

my ($in_file, $out_prefix, $coverage, $insert_size, $insert_size_sd, $read_length) = @ARGV;
$seed = $ARGV[6] if(@ARGV == 7);

open(LOG, ">$out_prefix.log") || die "\nCannot write to $out_prefix.log\n\n";
print LOG "Input file    \t$in_file\n";
print LOG "Out prefix    \t$out_prefix\n";
print LOG "Coverage      \t$coverage\n";
print LOG "Insert size   \t$insert_size\n";
print LOG "Insert size_sd\t$insert_size_sd\n";
print LOG "Read length   \t$read_length\n\n";
print LOG "Seed          \t$seed\n\n";

srand($seed);

open(R1, ">$out_prefix" . "_R1.fq") || die "\nCannot write to $out_prefix" . "_R1.fq\n\n";
open(R2, ">$out_prefix" . "_R2.fq") || die "\nCannot write to $out_prefix" . "_R2.fq\n\n";

my $in = new Bio::SeqIO(-file => $in_file);

my $index = 1;
while(my $seq_obj = $in->next_seq) {
	my $length = $seq_obj->length;
	my $npairs = int($length*$coverage/($read_length*2));

	if($length < $insert_size ) {
		print LOG $seq_obj->display_id, "\t$length\tShorter than insert-size ($insert_size), skipping\n";
		next;
	}
	print LOG $seq_obj->display_id, "\t$length\t$npairs\n";	

	my @r1_starts = split((0 x $npairs), //);
	my $last_start = $length-$insert_size;
	foreach my $x (0 .. ($npairs-1)) {
		$r1_starts[$x] = int(rand($last_start)+1);
	}
	@r1_starts = sort {$a <=> $b} @r1_starts;
	my $isize_generator = Math::Random::OO::Normal->new($insert_size, $insert_size_sd);
	$isize_generator->seed($seed);
	foreach my $r1_start (@r1_starts) {
		my $sense = (rand() > 0.5)? 1 : 0;
		my $isize = int($isize_generator->next);
		$isize = $length-$r1_start+1 if($r1_start+$isize-1 > $length);
		$isize = $read_length if($isize < $read_length);

		my $r2_start = $r1_start+$isize-$read_length;
#print "$r1_start\t$r2_start\t$isize\t$insert_size, $insert_size_sd\n";
#next;

		if(($r1_start <= 0) || ($r1_start+$read_length-1 > $length)) {
			print LOG "Skipping: impossible values for read 1: ($r1_start, ", ($r1_start+$read_length-1), "), sequence ", $seq_obj->display_id, ", length=$length\n";
			next;
		}
		if(($r2_start <= 0) || ($r2_start+$read_length-1 > $length)) {
			print LOG "Skipping: impossible values for read 2: ($r2_start, ", ($r2_start+$read_length-1), "), sequence ", $seq_obj->display_id, ", length=$length\n";
			next;
		}
		my $r1 = $seq_obj->subseq($r1_start, $r1_start+$read_length-1);
		my $r2 = $seq_obj->subseq($r2_start, $r2_start+$read_length-1);
		$r2 = reverse($r2);
		$r2 =~ tr/ACGT/TGCA/;

		if(rand() > 0.5) {
			print R1 "\@$index", "_1 seq_id=", $seq_obj->display_id, " start=$r1_start insert_size=$isize Sense\n$r1\n+\n", ('I' x $read_length), "\n";
			print R2 "\@$index", "_2 seq_id=", $seq_obj->display_id, " start=$r2_start insert_size=$isize Anti-sense\n$r2\n+\n", ('I' x $read_length), "\n";
		}
		else {
			print R2 "\@$index", "_2 seq_id=", $seq_obj->display_id, " start=$r1_start insert_size=$isize Sense\n$r1\n+\n", ('I' x $read_length), "\n";
			print R1 "\@$index", "_1 seq_id=", $seq_obj->display_id, " start=$r2_start insert_size=$isize Anti-sense\n$r2\n+\n", ('I' x $read_length), "\n";
		}
		$index++;
	}
}

close(R1);
close(R2);
close(LOG);
