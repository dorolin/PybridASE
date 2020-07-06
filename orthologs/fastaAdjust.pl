#!/usr/bin/perl

## clean up fasta file

## ------------------------------------------------------
## code modified from orthoMCL 'orthomclAdjustFasta' tool
## from the OrthoMCL software
## https://orthomcl.org/common/downloads/software/v2.0/
## Li et al. (2003) Genome Research 13: 2178-2189
## doi: 10.1101/gr.1224503
## ------------------------------------------------------

use strict;
use warnings;

&usage() unless scalar(@ARGV) >= 4;

my $input = shift(@ARGV);
my $taxon = shift(@ARGV);
my $delim = shift(@ARGV);
my $field = shift(@ARGV);
my @xtrstr = @ARGV;

open(IN, $input) || die "Cannot open input file '$input'\n";
open(OUT, ">$taxon.fasta") || die "Cannot open output file '$taxon.fasta'\n";


my %ids;
while(<IN>) {
    chomp $_;
    if(/\>/){
	s/^\>\s*//;
	my @a = split(/$delim/);
	my $id = $a[$field-1];
	for(my $i = 0; $i < scalar(@xtrstr); $i++){
	    $id =~ s/$xtrstr[$i]//;
	}
	die "Fasta file '$input' contains a duplicate id: $id\n" if $ids{$id};
	$ids{$id} = 1;
	print OUT ">$taxon\|$id\n";
    }else{
	s/\*+$//; ## remove asterisk(s) at end of line
	##s/\*/X/g; ## replace asterisk(s) in middle of sequence by X
	##          ## those might be stop-codons -> don't do
	unless($_ =~ /^\s*$/){
	    print OUT "$_\n";
	}
    }
}


sub usage {
print STDERR "
Adjust definition line fields in .fasta file and remove '*'s at end of
sequences, if any.

Usage: fastaAdjust.pl inputfile taxon fielddelimiter field [strings to remove]

Example: fastaAdjust.pl file.fa SP1 \" \" 5 \"transcript:\"

";
exit(1);
}
