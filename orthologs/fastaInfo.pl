#!/usr/bin/perl

## ----------------------------------------------------

## Extract info from fasta file
##   (gene info and length of transcript)

## Example: ./fastaInfo.pl MEL_good.fasta 2

## ----------------------------------------------------

use warnings;
use strict;

die "Usage: fastaInfo.pl <infile.fa> <definition line fields to keep>\n" unless @ARGV == 2;


my $fa = shift(@ARGV);
my $nfields = shift(@ARGV);

my $out = $fa;
$out =~ s/(\.fa$|\.fasta$|\.faa$)//;
$out = $out."_info.txt";

## open outfile
open (OUT, ">", "$out") or die "Couldn't open $out file to write to\n";

## open infile
open (FA, $fa) or die "Couldn't read $fa file\n";

## first entry
my $line = <FA>;
chomp $line;
if ($line !~ /^>/){
    die "first line in $fa file is not a definition line\n";
}
$line =~ s/^>//;
my @fields = split (/[\s\|]+/,$line); ## split by whitespace or '|'
my $m = 0;
my $field = '';
while ($m < $nfields){
    $field = shift(@fields);
    print OUT $field." ";
    $m++;
}

## remaining lines
my $count = 0;
my $base = '';
while ($line = <FA>){
    chomp $line;
    if ($line =~ /^>/){ ## new gene
	print OUT $count."\n";
	$count = 0;
	$line =~ s/^>//;
	@fields = split (/[\s\|]+/,$line); ## split by whitespace or '|'
	$m = 0;
	$field = '';
	while ($m < $nfields){
	    $field = shift(@fields);
	    print OUT $field." ";
	    $m++;
	}
    }
    else {
	$base = $line;
	$base =~ s/[A-Za-z]//g;
	if(length($base)>0){
	    die "unexpected symbol $base\n" ;
	}
	$count += length($line);
    }
}
## print last count
print OUT $count."\n";

close (FA);
close (OUT);


exit;
