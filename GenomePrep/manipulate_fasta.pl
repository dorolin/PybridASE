#!/usr/bin/env perl

## manipulate fasta output file from GATK FastaAlternateReferenceMaker,
## which changes original scaffold names from e.g.
## ">Peax302Chr1" to ">1 Peax302Chr1:1-189397759"

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
	$outfile=$infile.".edited";
}


open (IN, $infile) or die ("\nCan't open input file $infile\n\n");
open (OUT, ">$outfile") or die ("\nCan't write to output file: $outfile\n\n");

while (<IN>){

    chomp;

    if (/^>/){ ## header

	my @line=split(/\s/,$_);
	my @scaf=split(/:/,$line[1]);

	print OUT ">".$scaf[0]."\n";

    }
    else{ ## sequence
	print OUT $_."\n";
    }
}

close (IN);
close (OUT);

print "\nfasta file manipulation finished. File saved as: $outfile\n\n";


# Show usage
# ==============================================================================
sub usage{
	print "\n";
    print "  Usage:\n";
    print "    ".basename($0)."\n";
    print "      -i <fasta file from GATK FastaAlternateReferenceMaker>\n";
    print "      -o <output file> (optional)\n";
    print "\n";
	print "  Example:\n";
	print "      ".basename($0)." -i test.fa\n";
    print "\n\n";
    exit;
}
# ==============================================================================
