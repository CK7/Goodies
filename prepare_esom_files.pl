#!/usr/bin/perl
use strict;
use File::Copy;
use Bio::SeqIO;

######################################################################################################################################################
# Version 1.00, 06/25/2012, written by Itai Sharon (itaish@berkeley.edu)
# Version 1.01, 10/23/2012, Fixed a tiny bug in step 3 (search for "v1.01") 
# Version 1.02, 12/28/2012, (search for "V1.02") 
# Version 1.03, 09/03/2013, Fixed a bug that made lower-case sequences to be ignored (look for "V1.03"). Also fixed time/date printing and added 
#                           command line printout
######################################################################################################################################################
# Copyright (C) 2012 Itai Sharon
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), 
# to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
# DEALINGS IN THE SOFTWARE.
######################################################################################################################################################


######################################################################################################################################################
# User-defined parameters
######################################################################################################################################################
my $VERSION = "prepare_esom_files.pl V1.03";
my $kmer_size = 4;			# 0 means no DNA signature
my $normalize_abundance_pattern = 1;	# 2: log-transform, 1: normalize abundance pattern rows, 0: don't normalize
my $output_coverage = 0;		# 3: Add log-transformed coverage, 2: Add normalized coverage, 1: add coverage as a separate column, 0: don't add coverage
#my $output_normalized_coverage = 0;	# Same as $output_coverage but if (1) the coverage is normalized by the highest coverage 
#my $output_log_transformed_coverage = 0;# Same as $output_coverage but if (1) the coverage is normalized by the highest coverage 
my $assembly_file = undef;		# file that contains all sequences
my $annotation_file = undef;		# Annotation file for the assembly file. Optional, may contain only some of the contig names.
my $min_size = 3000;			# Minimum size for considering scaffolds
my $window_size = 3000;			
my $read_size = 0;
my $out_directory = undef;		# Output directory, mandatory 
my %column2files = ();
my %scaf2class = ();
my %class2annotation = ();

######################################################################################################################################################
sub usage {
	my $prog = $0;
	$prog =~ s/.+\///;
	print STDERR "\nUsage: $prog [-k <k-mer size>] [--raw_abundances | --log_transform] [--coverage | --normalized_coverage | --log_transform_coverage]\n";
	print STDERR "                             [-w <window-size>] [-a <annotation-file>] [-m <min-segment-size>]\n";
	print STDERR "                             [-sg <sam-files-glob>] [-sa <id> <sam-file-1> ... <sam-file-m>]\n";
	print STDERR "                             <out-directory> <assembly-file>\n\n";
	print STDERR "Option                    \t| Meaning\n";
	print STDERR "-------------------------------------------------------------------------------------------------------------------------------------\n";
	print STDERR "-k                        \t| Specify K-mer size (must be >=0, 0 means no DNA signature statistics to be used. Default=4)\n";
	print STDERR "-sg                       \t| Specify SAM file glob path for the computation of abundance pattern. Each file will be assigned\n";
	print STDERR "                          \t| a different dimensions. This parameter can be specified multiple times for different glob paths\n";
	print STDERR "                          \t| and is not mutually exclusive with the -si parameters.\n";
	print STDERR "-sa                       \t| Specify SAM files for specific abundance pattern column. id is the identifier of the column\n";
	print STDERR "                          \t| and all arguments after and until the next flag (or the last 2 parameters) are files that\n";
	print STDERR "                          \t| will be used for that column.\n";
	print STDERR "--log_transform           \t| log-transform each abundance pattern value (default: normalize rows)\n";
	print STDERR "--raw_abundances          \t| Do not normalize abundance pattern over all samples (default: normalize rows)\n";
	print STDERR "--coverage                \t| Add the sum of abundance pattern columns (before normalization) as an additional column (default: no)\n";
	print STDERR "--normalized_coverage     \t| Same as --coverage but values are normalized by the highest coverage in the coverage column (default: no)\n";
	print STDERR "--log_transformed_coverage\t| Same as --coverage but log10(values) reported on coverage column (default: no)\n";
	print STDERR "-w                        \t| Window size (default: 3,000 bp)\n";
	print STDERR "-m                        \t| Minimum segment size to be considered (default: 3,000 bp)\n";
	print STDERR "-a                        \t| Annotation file in the format of <scaf>\t<annotation>\n";
	print STDERR "\nLast two arguments must be the output directory and the assembly fasta file.\n";
	die "\n";
}

######################################################################################################################################################
# 1. Read command line
######################################################################################################################################################
my @timeDate = localtime(time);
print STDERR "\n$VERSION\n";
my $timestr = sprintf("%02d:%02d:%02d %02d/%02d/%02d", $timeDate[2], $timeDate[1], $timeDate[0], ($timeDate[4]+1), $timeDate[3], (1900+$timeDate[5]));	# V1.03
print STDERR "$timestr\n";	# V1.03

($#ARGV >= 1) || usage;

print STDERR "Command line: $0 ", join(" ", @ARGV), "\n";	# V1.03

$assembly_file = pop(@ARGV);
$out_directory = pop(@ARGV);

my %files = ();
while($ARGV[0] =~ /^\-/) {
	my $flag = shift @ARGV;
	if($flag eq '-k') {
		$kmer_size = shift(@ARGV);
	}
	elsif($flag eq '-w') {
		$window_size = shift(@ARGV);
	}
	elsif($flag eq '-m') {
		$min_size = shift(@ARGV);
	}
	elsif($flag eq '-a') {
		$annotation_file = shift(@ARGV);
	}
	elsif($flag eq '-sa') {
		my $id = shift(@ARGV);
		my $sam_file = shift(@ARGV);

		while(defined($sam_file) && ($sam_file !~ /^\-/)) {
			((-e $sam_file) && (-f $sam_file)) || die "\nError: file $sam_file specified with -sa option either does not exist or is not a file\n\n";
			!exists($files{$sam_file}) || die "\nError: file $sam_file specified more than once with either -sa or -sg options\n\n";
			$column2files{$id}{$sam_file} = 1;
			$files{$sam_file} = 1;
			$sam_file = shift(@ARGV);
		}
		if($sam_file =~ /^\-/) {
			unshift(@ARGV, $sam_file);
		}
	}
	elsif($flag eq '-sg') {
		my $sam_files_glob = shift(@ARGV);
		foreach my $sam_file (glob($sam_files_glob)) {
			if(!(-f $sam_file)) {
				print STDERR "Warning: $sam_file specified with the -sg option is not a file, skipping\n\n";
				next;
			}
			!exists($files{$sam_file}) || die "\nError: file $sam_file specified more than once with either -sa or -sg options\n\n";
			$column2files{$sam_file}{$sam_file} = 1;
			$files{$sam_file} = 1;
		}
	}
	elsif($flag eq '--raw_abundances') {
		$normalize_abundance_pattern = 0;
	}
	elsif($flag eq '--log_transform') {
		$normalize_abundance_pattern = 2;
	}
	elsif($flag eq '--coverage') {
		$output_coverage = 1;
#		$output_normalized_coverage = 0;
#		$output_log_transformed_coverage = 0;
	}
	elsif($flag eq '--normalized_coverage') {
		$output_coverage = 2;
#		$output_normalized_coverage = 1;
#		$output_log_transformed_coverage = 0;
	}
	elsif($flag eq '--log_transformed_coverage') {
		$output_coverage = 3;
#		$output_normalized_coverage = 0;
#		$output_log_transformed_coverage = 1;
	}
	else {
		die "\nError: unrecognized option $flag\n\n";
	}
}

(-e $assembly_file) || die "\nCannot read $assembly_file\n\n";
mkdir($out_directory);
((-e $out_directory) && (-d $out_directory)) || die "\nCannot create directory $out_directory\n\n";
(!$output_coverage || ((keys %column2files) > 0)) || die "\nError: cannot use --coverage or --normalized_coverage without specifying SAM files\n\n";

######################################################################################################################################################
# 1.5 Read annotation file
######################################################################################################################################################
if(defined($annotation_file)) {
	print STDERR "Reading annotation file $annotation_file ... ";
	open(IN, $annotation_file) || die "\nCannot read $annotation_file (specified with the -a option)\n\n";
	my %annotation2class = ();
	while(<IN>) {
		chomp;
		my ($scaf, $annotation) = split(/\t/);
# die "($scaf, $annotation)\n";
		if(!exists($annotation2class{$annotation})) {
			$annotation2class{$annotation} = ((keys %annotation2class)+1);
			$class2annotation{$annotation2class{$annotation}} = $annotation;
		}
		$scaf2class{$scaf} = $annotation2class{$annotation};
	}
	close(IN);
	print STDERR "ok\n";
}

######################################################################################################################################################
# 2. Read SAM files
######################################################################################################################################################
my @sample_order = ();
my %read_starts = ();
my %normalization_factor = ();
my %scaf2length = ();

my $first_total = undef;
foreach my $column_id (sort keys %column2files) {
	my $total = 0;
	push(@sample_order, $column_id);
	foreach my $sam_file (keys %{$column2files{$column_id}}) {
		print STDERR "Reading $sam_file (column $column_id) ... ";
		my $file_total = 0;
		open(IN, $sam_file) || die "\nCannot read $sam_file\n\n";
		while(<IN>) {
			if($_ =~ /^\@SQ/) {
				# @SQ	SN:NODE_840_length_578_cov_29.223183.36.514	LN:479
				($_ =~ /\@SQ\s+SN:(\S+)\s+LN:(\d+)/) || die "Line with unknown format: $_\n";
				if(!exists($scaf2length{$1})) {
					$scaf2length{$1} = $2;
				}
				# Sanity check, just make sure that the info in all SAM files is consistent
				elsif($scaf2length{$1} != $2) {
					die "\nFatal error: length for $1 is different in this file and one of the previous\n\n";
				}
			}
			elsif($_ !~ /^\@/) {
				chomp;
				$file_total++;	
				# HWI-ST330_0096:1:28:7341:66385#GATCAG/1	16	NODE_126_length_55165_cov_80.104744	2	255	100M	*	0	0	GAGAGTTTATAAAAACTACTTGGGAAGGTATTAAGACTTTAATTTCAACAGTTCTTGATGCAATAAAGGTAAAAGTTGAGACTATTTGGAATGGACTAAA	FEEEGEHHHHEFFHCHHFHHHHDGHBEEEECGEFFCEEDCHHCEHHHHHHHFDEHFHHHCHHFHHHHHHHHGHHHHHHHHEHHHHHHHHHHHHHHHHHHH	XA:i:0	MD:Z:100	NM:i:0
				my @fs = split(/\t/); 
				my ($scaf, $pos, $l) = ($fs[2], $fs[3], length($fs[9]));
				$read_size = $l if($l > $read_size);
				next if($scaf eq '*');
				$read_starts{$column_id}{$scaf}{$pos}++;
			}
		}
		close(IN);
		print STDERR "ok, $file_total mapped reads\n";
		$total += $file_total;
	}
	$first_total = $total if(!defined($first_total));
	$normalization_factor{$column_id} = $first_total/$total;
	print STDERR "Read $total reads overall for column $column_id, normalization factor is ", $normalization_factor{$column_id}, "\n\n";
}

######################################################################################################################################################
# 3. Read assembly file and determine segments
######################################################################################################################################################
print STDERR "Determining segments ... ";
my %segments = ();
my $nsegments = 0;

my $in = new Bio::SeqIO(-file => $assembly_file);

while(my $scaf = $in->next_seq) {
	next if($scaf->length < $min_size);
	my $s = 'X' . $scaf->seq;	# V1.02
	$s =~ tr/acgtn/ACGTN/;
	# If the sequence is shorter than a window size then we will treat it as just one piece, even with N's. Sequence will be excluded 
	# if the number of non-N's is smaller than $min_size;
	if($scaf->length < $window_size) {
		my @subseqs = split(/N+/, $s);
		my $nbps = ($s =~ tr/ACGT/ACGT/);
		($nbps > $min_size) || next;
		$segments{$scaf->display_id}{1} = [$scaf->length-1, scalar(@subseqs), $nbps, $scaf->display_id . '_1'];	# V1.02
		$nsegments++;
		next;
	}

	$s .= 'NX';	# This will ensure that the loop below will cover all segments, including the last one 
	my $start = 1;	# V1.02
	my $index = 1;
	while($s =~ m/(N+)/g) {
		my $end = pos($s)-length($1)-1;
		# v1.01: We will use the segment if it's bigger than either $min_size or $window_size 
		if(($end-$start+1 < $window_size) && ($end-$start+1 < $min_size)) {
			$start = pos($s);
			next;
		}
		my $ncsegments = int(($end-$start+1)/$window_size);
		$ncsegments = 1 if($ncsegments == 0);
		my $cwindow_size = int(($end-$start+1)/$ncsegments);
		foreach my $i (1 .. ($ncsegments-1)) {
			my $subseq = substr($s, $start, $cwindow_size);
			my $nbps = ($subseq =~ tr/ACGT/ACGT/);
			$segments{$scaf->display_id}{$start} = [$start+$cwindow_size-1, 0, $nbps, $scaf->display_id . "_$index"];	# V1.02
			$start += $cwindow_size;
			$index++;
		}
		my $subseq = substr($s, $start, $end-$start+1);	# V1.02
		my $nbps = ($subseq =~ tr/ACGT/ACGT/);
		$segments{$scaf->display_id}{$start} = [$end, 1, $nbps, $scaf->display_id . "_$index"];
		$start = pos($s);
		$index++;
		$nsegments += $ncsegments;
	}
}
print STDERR "ok, $nsegments segments assigned\n";

######################################################################################################################################################
# 4. Compute abundance pattern
######################################################################################################################################################
print STDERR "Computing abundance patterns ... " if((keys %column2files) > 0);
my %abundance_pattern_columns = ();
my %coverage_column = ();

foreach my $scaf (keys %segments) {
	foreach my $start (keys %{$segments{$scaf}}) {
		my $end = $segments{$scaf}{$start}[0];
		my $possible_starting_points = $segments{$scaf}{$start}[2] - $segments{$scaf}{$start}[1]*($read_size-1);
		foreach my $column_id (keys %read_starts) {
			my $coverage = 0;
			foreach my $i ($start .. $end) {
				$coverage += $read_starts{$column_id}{$scaf}{$i};
			}
			$abundance_pattern_columns{$scaf}{$start}{$column_id} = $coverage * $normalization_factor{$column_id} / $possible_starting_points;
			# We will either use this or not
			$coverage_column{$scaf}{$start} += $abundance_pattern_columns{$scaf}{$start}{$column_id};
		}
	}
} 

if($normalize_abundance_pattern != 0) {
	foreach my $scaf (keys %abundance_pattern_columns) {
		foreach my $start (keys %{$abundance_pattern_columns{$scaf}}) {
			foreach my $column_id (keys %{$abundance_pattern_columns{$scaf}{$start}}) {
				if($normalize_abundance_pattern == 1) {
					$abundance_pattern_columns{$scaf}{$start}{$column_id} /= $coverage_column{$scaf}{$start} if($coverage_column{$scaf}{$start} > 0);
				}
				elsif($normalize_abundance_pattern == 2) {
					$abundance_pattern_columns{$scaf}{$start}{$column_id} = log($abundance_pattern_columns{$scaf}{$start}{$column_id}+1)/log(10);
				}
			}
		}
	}
}

######################################################################################################################################################
# 5. If we need to compute also a normalized coverage column then do it here 
######################################################################################################################################################
if($output_coverage == 3) {
	foreach my $scaf (keys %coverage_column) {
		foreach my $start (keys %{$coverage_column{$scaf}}) {
			$coverage_column{$scaf}{$start} = log($coverage_column{$scaf}{$start})/log(10);	
		}
	}
}
elsif($output_coverage == 2) {
	my $normalizer = undef;
	foreach my $scaf (keys %coverage_column) {
		foreach my $start (keys %{$coverage_column{$scaf}}) {
			$normalizer = $coverage_column{$scaf}{$start} if($coverage_column{$scaf}{$start} > $normalizer);	
		}
	}
	foreach my $scaf (keys %coverage_column) {
		foreach my $start (keys %{$coverage_column{$scaf}}) {
			$coverage_column{$scaf}{$start} /= $normalizer;	
		}
	}
}
print STDERR "ok\n" if((keys %column2files) > 0);

######################################################################################################################################################
# 6. Now go back to the assembly file and compute the DNA signature
######################################################################################################################################################
# Begin by creating the list of k-mers. Reverse complementory kmers will be unified under the k-mer that comes earlier in an 
# alphabetical order
my @mer_order = ();
my %DNA_signature_columns = ();
if($kmer_size > 0) {
	print STDERR "Computing DNA signature ... ";
	my %mer_dictionary = ();
	make_list_of_possible_tetramers('', $kmer_size, \%mer_dictionary, \@mer_order);

	my $in = new Bio::SeqIO(-file => $assembly_file);

	while(my $seq = $in->next_seq) {
		next if($seq->length < $min_size);
		foreach my $start (keys %{$segments{$seq->display_id}}) {
			my $end = $segments{$seq->display_id}{$start}[0];
			my $s = $seq->subseq($start, $end);	# V1.02
			$s =~ tr/acgtn/ACGTN/;	# V1.03
			my @mer_dist = ();
			calc_kmer_dist($s, \@mer_dist, \%mer_dictionary, \@mer_order, $kmer_size);
			$DNA_signature_columns{$seq->display_id}{$start} = \@mer_dist;	
		}
	}
	print STDERR "ok\n";
}

######################################################################################################################################################
# 7. And now, write the output files
######################################################################################################################################################
print STDERR "Writing output to $out_directory ... ";
open(LRN, ">$out_directory/esom.lrn") || die "\nCannot write to $out_directory/esom.lrn\n\n";
open(CLS, ">$out_directory/esom.cls") || die "\nCannot write to $out_directory/esom.cls\n\n";
open(NAMES, ">$out_directory/esom.names") || die "\nCannot write to $out_directory/esom.names\n\n";

my $ndimensions = @mer_order;
$ndimensions += @sample_order;
$ndimensions++ if($output_coverage > 0);

print CLS "\% $nsegments\n";
my $step = int(255 / (1+int((keys %class2annotation) ** (1/3))));
my $cls = 1;
my ($r, $g, $bl) = (255, 255, 255);
while($cls <= (keys %class2annotation)) {
	print CLS "\%$cls\t", $class2annotation{$cls}, "\t$r\t$g\t$bl\n";
	$cls++;
	if(($bl-$step) > 0) {
		$bl -= $step;
	}
	else {
		$bl = 255;
		if(($g-$step) > 0) {
			$g -= $step;
		}
		else {
			$g = 255;
			$r -= $step;
		}
	}
}

print NAMES "\% $nsegments\n";
print LRN "\% $nsegments\n";
print LRN "\% ", ($ndimensions+1), "\n";
print LRN "\% 9", ("\t1" x $ndimensions), "\n";
print LRN "\% key";
foreach my $column_id (@sample_order) {
	print LRN "\t$column_id";
}
print LRN "\tCoverage" if($output_coverage == 1);
print LRN "\tNormalized coverage" if($output_coverage == 2);
print LRN "\tLog-transformed coverage" if($output_coverage == 3);
foreach my $mer (@mer_order) {
	print LRN "\t$mer";
}
print LRN "\n";

my $dp = 0;
foreach my $scaf (keys %segments) {
	foreach my $start (sort {$a <=> $b} keys %{$segments{$scaf}}) {
		my $end = $segments{$scaf}{$start}[0];
		$dp++;
		print NAMES "$dp\t", $segments{$scaf}{$start}[3], "\t$scaf:($start,$end), ", $segments{$scaf}{$start}[1], " segment(s), ",  $segments{$scaf}{$start}[2], 
				'/', ($end-$start+1), " non-N bps\n";
		my $cls = exists($scaf2class{$scaf})? $scaf2class{$scaf} : 0;
		print CLS "$dp\t$cls\n";
		print LRN "$dp";
		foreach my $column_id (@sample_order) {
			print LRN "\t", $abundance_pattern_columns{$scaf}{$start}{$column_id};
		}
		print LRN "\t", $coverage_column{$scaf}{$start} if($output_coverage > 0);
		print LRN "\t", join("\t", @{$DNA_signature_columns{$scaf}{$start}}) if(@mer_order>0);
		print LRN "\n";
		
	}
}

print STDERR "ok\n";
print STDERR "\nFinished successfully\n";

#=======================================================================================================================================#
sub make_list_of_possible_tetramers {
	my ($mer, $k, $mer_dictionary_ref, $mer_order_ref) = @_;

	if($k == 0) {
                my $rc_mer = reverse($mer);
		$rc_mer =~ tr/ACGT/TGCA/;
		push(@{$mer_order_ref}, $mer) if($rc_mer ge $mer);
		$mer_dictionary_ref->{$mer} = ($rc_mer ge $mer)? $mer : $rc_mer;
                return;
        }

	my @bases = ("A", "T", "C", "G");
	foreach my $na (@bases) {
		make_list_of_possible_tetramers("$mer$na", $k-1, $mer_dictionary_ref, $mer_order_ref);
	}
}

#=======================================================================================================================================#
sub calc_kmer_dist {
	# Note: sequences may contain N's (see above in the segments computation part)
	my ($seq, $mer_dist_ref, $mer_dictionary_ref, $mer_order_ref, $k) = @_;
	my @subseqs = split(/N+/, $seq);
	my $total = length($seq)-@subseqs*($k+1);
	my %dist = ();
	foreach my $s (@subseqs) {
		my $n = length($s)-$k;
		foreach my $i (0 .. $n) {
			$dist{$mer_dictionary_ref->{substr($s, $i, $k)}}++;
		}
	}
	foreach my $mer (@{$mer_order_ref}) {
		push(@{$mer_dist_ref}, $dist{$mer}/$total);
	}
}

