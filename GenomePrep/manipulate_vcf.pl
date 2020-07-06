#!/usr/bin/env perl

## manipulate vcf file that only contains alt/alt homozygous variants
## for a single individual and re-code them as ref/alt heterozygotes,
## and only keep GT field (also adjust AC, AF, MLEAC, MLEAF)
## (assuming that the following fields are present: GT:AD:DP:GQ:PL)

use strict;
use warnings;
use Getopt::Long;
use File::Basename;


my ($infile, $outfile);
GetOptions(
    'i|I=s'  => \$infile,
    'o|O=s'  => \$outfile,
    'h|help' => \&usage
)
or (print "\nERROR: Argument invalid\n" and &usage);

&usage if (!defined($infile));


if (!defined($outfile)){
	$outfile=$infile;
	$outfile=~ s/\.vcf$/\_het.vcf/;
}


open (IN, $infile) or die ("\nCan't open input file $infile\n\n");
open (OUT, ">$outfile") or die ("\nCan't write to output file: $outfile\n\n");

while (<IN>){
    ## remove description lines
    next if (/^\#\#FORMAT\=<ID\=AD,/);
    next if (/^\#\#FORMAT\=<ID\=DP,/);
    next if (/^\#\#FORMAT\=<ID\=GQ,/);
    next if (/^\#\#FORMAT\=<ID\=PL,/);

    chomp;

    if (/^\#/){ ## print header (description or column names)
	print OUT $_."\n";
    }
    else{ ## process SNP info
	my @line=split(/\s/,$_);
	my @info=split(/;/,$line[7]);

	foreach my $i (0..$#info){
	    if($info[$i] =~ m/^AC\=/){
		$info[$i] = "AC=1";
	    }
	    elsif($info[$i] =~ m/^AF\=/){
		$info[$i] = "AF=0.5";
	    }
	    elsif($info[$i] =~ m/^MLEAC\=/){
		$info[$i] = "MLEAC=1";
	    }
	    elsif($info[$i] =~ m/^MLEAF\=/){
		$info[$i] = "MLEAF=0.5";
	    }
	}

	foreach my $i (0..6){
	    print OUT $line[$i]."\t";
	}
	foreach my $i (0..$#info){
	    print OUT $info[$i].";";
	}
	print OUT "\t"."GT"."\t"."0\/1"."\n";
    }
}

close (IN);
close (OUT);

print "\nvcf file manipulation finished. File saved as: $outfile\n\n";


# Show usage
# ==============================================================================
sub usage{
	print "\n";
    print "  Usage:\n";
    print "    ".basename($0)."\n";
    print "      -i <SNPs input file (vcf format for single individual)>\n";
    print "      -o <output file> (optional)\n";
    print "\n";
	print "  Example:\n";
	print "      ".basename($0)." -i test.vcf\n";
    print "\n\n";
    exit;
}
# ==============================================================================
