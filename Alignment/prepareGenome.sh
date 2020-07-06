#!/bin/bash

#SBATCH --mail-user=me@work
#SBATCH --mail-type=end,fail
#SBATCH --job-name="prepareGenome"
#SBATCH --output=%x-%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=50G
#SBATCH --partition=all


module load vital-it/7
module load UHTS/Aligner/STAR/2.6.0c


#$HOME/bin/gffread/gffread Peax304_final.gff -E -T -o Peax304_final.gtf

## mkdir Pax_304_generated

STAR  --runMode genomeGenerate --runThreadN 4 --genomeDir Pax_304_generated --genomeFastaFiles Petunia_axillaris_304.faa --sjdbGTFfile Peax304_final.gtf --limitGenomeGenerateRAM 100000000000




## STAR

## basic options:

# --runThreadN NumberOfThreads
# --runMode genomeGenerate
# --genomeDir /path/to/genomeDir
# --genomeFastaFiles /path/to/genome/fasta1 /path/to/genome/fasta2 ...
# --sjdbGTFfile /path/to/annotations.gtf
# --sjdbOverhang ReadLength-1

## to get an estimate for sjdbOverhang, e.g. run
## awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' AxExMS2_2_cleaned.fastq
## -> read length is 101 for many, so default of 100 works fine

## advanced options:

## annotations in GFF format:
# --sjdbGTFfeatureExon
# --sjdbGTFtagExonParentTranscript (=transcript id by default)
# --sjdbGTFtagExonParentGene


## from the STAR manual 2.6.0a:

# In addition to the aforementioned options, for GFF3 formatted annotations you need to use --sjdbGTFtagExonParentTranscript Parent. In general, for --sjdbGTFfile files STAR only processes lines which have --sjdbGTFfeatureExon (=exon by default) in the 3rd field (col- umn). The exons are assigned to the transcripts using parent-child relationship defined by the --sjdbGTFtagExonParentTranscript (=transcript id by default) GTF/GFF attribute.

# [although it might be easier to convert gff to gtf...]

## genome with a large number of references:
# If you are using a genome with a large (>5,000) number of references (chro- somes/scaffolds), you may need to reduce the --genomeChrBinNbits to reduce RAM consumption. The following scaling is recommended: --genomeChrBinNbits = min(18,log2[max(GenomeLength/NumberOfReferences,ReadLength)]). 3 gigaBase genome with 100,000 chromosomes/scaffolds, this is equal to 15.

