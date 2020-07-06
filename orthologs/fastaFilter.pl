#!/usr/bin/perl

## filter poor sequences

## ------------------------------------------------------
## code modified from orthoMCL 'orthomclFilterFasta' tool
## from the OrthoMCL software
## https://orthomcl.org/common/downloads/software/v2.0/
## Li et al. (2003) Genome Research 13: 2178-2189
## doi: 10.1101/gr.1224503
## ------------------------------------------------------

use strict;

&usage() unless scalar(@ARGV) == 3 || scalar(@ARGV) == 5;

my $input = shift(@ARGV);
my $abbrev = $input;
$abbrev =~ s/(\w+)(\.fa$|\.fasta$|\.faa$)/$1/ || die "File $input does not have a common fasta file name ending\n";

my $minLength = shift(@ARGV);
my $maxStopPercent = shift(@ARGV);

my $goodFile = $input;
$goodFile =~ s/$abbrev/$abbrev\_good/;
my $poorFile = $input;
$poorFile =~ s/$abbrev/$abbrev\_poor/;
if(scalar(@ARGV)==2){
    $goodFile = shift(@ARGV);
    $poorFile = shift(@ARGV);
}


open(IN, $input) || die "Can't open input file $input\n";

my $rejectRates = [];
open(GOOD, ">$goodFile") || die "Can't open output file $goodFile\n";
open(BAD, ">$poorFile") || die "Can't open output file $poorFile\n";

my $seqCount;
my $rejectSeqCount;
my $currentSeq;
my $currentLen;
my $currentStopCnt;

# process lines of one file
while (<IN>) {
    chomp;
    # handle prev seq
    if (/\>/) {
	(/\>([^|]+)|/ && $1 eq $abbrev) ||
	    die "The ID on def line '$_' is missing the prefix '$abbrev|' '$1'\n" unless $1 eq $abbrev;
	if ($currentSeq) {
	    die "Error: zero length sequence in file $input.  Look near line '$_'\n" if $currentLen == 0;
	    $seqCount++;
	    $rejectSeqCount += &handleSeq($currentSeq, $currentLen, $currentStopCnt);
	    $currentSeq = "";
	    $currentLen = 0;
	    $currentStopCnt = 0;
	}
    } else {
	$currentLen += length($_);
	$currentStopCnt += tr/[^A-Za-z]//; # this removes the stop codon from $_
    }
    $currentSeq .= "$_\n";
}
$rejectSeqCount += &handleSeq($currentSeq, $currentLen, $currentStopCnt);
$seqCount++;

# add file stats to reject count if it qualifies
if ($rejectSeqCount) {
    my $pct = $rejectSeqCount/$seqCount * 100;
    if ($pct > 10) {
	push(@$rejectRates, [$input, $pct]);
    }
}

close(IN);
close(GOOD);
close(BAD);


if (scalar(@$rejectRates)) {
  print "\nProteomes with > 10% poor proteins:\n";
  my @sortedRR = sort {$b->[1] <=> $a->[1]} @$rejectRates;
  foreach my $reject (@sortedRR) {
    my $intPct = int($reject->[1]);
    print "  $reject->[0]\t$intPct%\n";
  }
}

sub handleSeq {
  my ($seq, $len, $stopCnt) = @_;
  my $isBad = 0;
  my $stopPercent = (($len - $stopCnt)/$len)* 100;

  if ($len < $minLength || $stopPercent > $maxStopPercent) {
    print BAD $seq;
    $isBad = 1;
  } else {
    print GOOD $seq;
  }
  return $isBad;
}

sub usage {
  print STDERR "
Filter protein sequences in fasta file by minimum length and max percent stop codons.

Usage:
  fastaFilter.pl input_file min_length max_percent_stops [good_proteins_file poor_proteins_file]

where:
  input_file:              fasta input file
  min_length:              minimum allowed length of proteins.  (suggested: 10)
  max_percent_stop:        maximum percent stop codons.  (suggested 0)
  good_proteins_file:      optional. Output file good proteins
  poor_proteins_file:      optional. Output file poor proteins

EXAMPLE: fastaFilter.pl SPECIES.fasta 10 0

";
  exit(1);
}
