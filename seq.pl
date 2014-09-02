#!/usr/bin/perl
use strict;

###############################################################################################################
# seq.pl
# This program retrieves one sequence from a FASTA file or the standard input and allows operations om it,
# including getting its reverse complement or just a segment of the sequence.
#
# Written by Itai Sharon, itai.sharon@gmail.com (09/02/2014)
###############################################################################################################

my $VERSION = "v1.00, 09/02/2014";

my ($start, $end) = (undef, undef);
my ($seq_id, $input_file, $output_file) = (undef, undef, undef);
my $line_length = 80;
my $reverse = 1;
my $IN_ref = \*STDIN;

(($ARGV[0] ne '-h') && ($ARGV[0] ne '--help')) || usage();
#die "\nUsage: $0 [-s <query-seq>] [-d <input-file>] [-L <start-coord>,<end-coord>] [-S <1|2>] [-l <line-length>] [-o <output-file>]\n\n";

while(@ARGV > 0) {
	my $flag = shift(@ARGV);
	my $arg = shift(@ARGV);
	if($flag eq '-S') {
		$reverse = $arg;
	}
	elsif($flag eq '-s') {
		$seq_id = $arg;
	}
	elsif($flag eq '-d') {
		$input_file = $arg;
	}
	elsif($flag eq '-L') {
		($start, $end) = split(/,/, $arg);
		$start = undef if($start == 0);
		$end = undef if($end == 0);
		die "Error: end < start (-L)\n\n" if(defined($start) && defined($end) && ($end < $start));
	}
	elsif($flag eq '-e') {
		$end = $arg;
	}
	elsif($flag eq '-l') {
		$line_length = $arg;
	}
	elsif($flag eq '-o') {
		$output_file = $arg;
	}
	else {
		die "\nError: unknown flag $flag\n\n";
	}
	defined($arg) || die "\nError: flag $flag requires an argument\n\n";
}

if(defined($input_file)) {
	open(IN, $input_file) || die "\nError: could not read file $input_file\n\n";
	$IN_ref = \*IN;
}

# Look for the sequence. If input stream contains more than one sequence and seq_id was not specified then the first sequence will be pickes
my ($header, $seq) = (undef, undef);
while(<$IN_ref>) {
	chomp;
	if($_ =~ />(\S+)/) {
		last if(defined($header));
		if(($1 eq $seq_id) || !defined($seq_id)) {
			$header = $_; 
			$seq_id = $1;
		}
	}
	elsif(defined($header)) {
		$seq .= $_;
	}
}

$end = length($seq) if(defined($end) && ($end > length($seq)));

if(defined($start) || defined($end)) {
	if(defined($start) && defined($end)) {
		$seq_id .= ":$start-$end";
		$seq = substr($seq, $start-1, $end-$start+1);
	}
	elsif(defined($start)) {
		$seq_id .= ":$start-" . length($seq);
		$seq = substr($seq, $start-1);
	}
	else {
		$seq_id .= ":1-$end";
		$seq = substr($seq, 0, $end);
	}
	$header = ">$seq_id";
}
if($reverse == 2) {
	$seq =~ tr/ACGTacgt/TGCAtgca/;
	$seq = reverse($seq);
}

$seq =~ s/(.{$line_length})/$1\n/g;
if(defined($output_file)) {
	open(OUT, ">$output_file") || die "\nCannot write to $output_file\n\n";
	print OUT "$header\n$seq\n";
	close(OUT);
}
else {
	print "$header\n$seq\n";
}

##########################################################################################################################################
sub usage {
	print STDERR "\nseq.pl $VERSION\n\n";
	print STDERR "This program is similar to NCBI's fastacmd command but it works on FASTA files instead of BLAST databases.\n";
	print STDERR "Input can be taken either from a file or from the standard input.\n";
	print STDERR "\nArguments:\n\n";	
	print STDERR "  -d\tInput file (defualt: input will be taken from the standard input)\n";	
	print STDERR "  -s\tSequence to retrieve (default: first sequence in input stream will be chosen)\n";	
	print STDERR "  -l\tLine length (default: 80 characters)\n";	
	print STDERR "  -o\tOutput file (default: output will be written to the standard output)\n";	
	print STDERR "  -L\t1-offset range of sequence to extract in the format start,end.\n";
	print STDERR "    \t0 in start refers to beginning of sequence, 0 in end refers to the end of the sequence. Defualt: 0,0\n";
	print STDERR "  -S\tStrand. 1 referse to forward, 2 referse to reverse (default: 1)\n";
	die "\n";
}
