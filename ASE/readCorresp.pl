#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;


## Example:
## ## SLURM_ARRAY_TASK_ID=35
## ## samples=$HOME/Style/step02_AxEx/files_AxEx.txt
## ## filebase=$(cat $samples | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')
## ## FILE1=$HOME/Style/step02_AxEx/Pax_genome/secondPass/${filebase}/Aligned.sortedByCoord.out.bam
## ## FILE2=$HOME/Style/step02_AxEx/Pex_genome/secondPass/${filebase}/Aligned.sortedByCoord.out.bam
## $HOME/Style/scripts/readCorresp.pl -f $FILE1 -s $FILE2 -o ${filebase}_aligninfo.txt


## Output:
## QNAME FILE1(FLAG RNAME POS MAPQ NH:i HI:i AS:i) FILE2(FLAG RNAME POS MAPQ NH:i HI:i AS:i)
## (where entries for multimappers are concatenated with ',')

&usage if (scalar(@ARGV)!=6);

my ($firstsam, $secondsam, $outfile);
GetOptions(
    'f|F=s'     => \$firstsam,
    's|S=s'     => \$secondsam,
    'o|O=s'     => \$outfile,
    'h|help'    => \&usage
)
or &usage("\nERROR: Argument invalid (see above)\n\n");


## external programs
my $samtools='samtools';

## three additional sam flags for output
my $f1=11; ## NH:i: (number of alignments: 1=unique)
my $f2=12; ## HI:i: (query hit index: enumerates multiple alignments)
my $f3=13; ## AS:i: (alignment score)

my %allreads=();
my @line=();
my @oldinfo=();
my @newinfo=();
my $concatinfo='';
my $read='';

## Read first infile
## ----------------------------------------------------------
if($firstsam =~ /\.bam$/){
    open (FILE, "$samtools view $firstsam | ")
	or die ("\nCan't open input file $firstsam\n\n");
}
else{
    open (FILE, "$firstsam")
	or die ("\nCan't open input file $firstsam\n\n");
}

while (<FILE>){
    @line = split(/\s+/,$_);

    $read=$line[0];
    if (!exists($allreads{$read})){
	## initialize hash of arrays
	$allreads{$read} = [ ". . . . . . .", ". . . . . . ."];
	$allreads{$read}->[0] = "$line[1] $line[2] $line[3] $line[4] $line[$f1] $line[$f2] $line[$f3]";
    }
    else{ ## append info
	@oldinfo = split(/ /,$allreads{$read}->[0]);
	@newinfo = ($line[1],$line[2],$line[3],$line[4],$line[$f1],$line[$f2],$line[$f3]);
	$concatinfo = "$oldinfo[0],$newinfo[0]";
	foreach my $i (1 .. $#oldinfo) {
	    $concatinfo .= " $oldinfo[$i],$newinfo[$i]";
	}
	$allreads{$read}->[0] = $concatinfo;
    }
}

close (FILE);


## Read second infile
## ----------------------------------------------------------
if($secondsam =~ /\.bam$/){
    open (FILE, "$samtools view $secondsam | ")
	or die ("\nCan't open input file $secondsam\n\n");
}
else{
    open (FILE, "$secondsam")
	or die ("\nCan't open input file $secondsam\n\n");
}

while (<FILE>){
    @line = split(/\s+/,$_);

    $read=$line[0];
    if (!exists($allreads{$read})){
	## initialize hash of arrays
	$allreads{$read} = [ ". . . . . . .", ". . . . . . ."];
	$allreads{$read}->[1] = "$line[1] $line[2] $line[3] $line[4] $line[$f1] $line[$f2] $line[$f3]";
    }
    else{
	if ($allreads{$read}->[1] =~ /^\. \. \. \. \. \. \.$/){ ## no info
	    $allreads{$read}->[1] = "$line[1] $line[2] $line[3] $line[4] $line[$f1] $line[$f2] $line[$f3]";
	}
	else{ ## append info
	    @oldinfo = split(/ /,$allreads{$read}->[1]);
	    @newinfo = ($line[1],$line[2],$line[3],$line[4],$line[$f1],$line[$f2],$line[$f3]);
	    $concatinfo = "$oldinfo[0],$newinfo[0]";
	    foreach my $i (1 .. $#oldinfo) {
		$concatinfo .= " $oldinfo[$i],$newinfo[$i]";
	    }
	    $allreads{$read}->[1] = $concatinfo;
	}
    }
}

close (FILE);


## Write outfile (this is not sorted anymore)
## ----------------------------------------------------------
open (OUT, ">$outfile")
    or die ("\nCan't write to file $outfile\n\n");

for $read (keys %allreads) {
    print OUT "$read\t$allreads{$read}->[0]\t$allreads{$read}->[1]\n";
}

close (OUT);




## =================================================================
sub usage{
    my $txt='';
    $txt=shift if (scalar(@_)>0);
    print "\n";
    print "  Usage:\n";
    print "    ".basename($0)."\n";
    print "      -f    <first input file (sam format)>\n";
    print "      -s    <second input file (sam format)>\n";
    print "      -o    <output file>\n";
    print "\n\n";
    print $txt;
    exit;
}
## =================================================================


## then sort the output with bash commands, see 'readCorresp.sh'
