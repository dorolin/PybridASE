#!/bin/bash

#SBATCH --mail-user=me@work
#SBATCH --mail-type=end,fail
#SBATCH --job-name="snps_preprocess"
#SBATCH --output=%x-%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=71:50:00 ## 31:50:00 should be enough
#SBATCH --mem-per-cpu=20G ## 10G should be enough
#SBATCH --partition=all


module load vital-it/7
module load UHTS/Analysis/picard-tools/2.18.11
module load UHTS/Analysis/GenomeAnalysisTK/4.1.0.0
module load UHTS/Analysis/samtools/1.8


mkdir -p tmp


## may want to sort bam first

picard-tools SortSam \
	     I=$HOME/Style/CoGe_data/Exreads_on_Pax304.unsorted.bam \
	     O=$HOME/Style/CoGe_data/Exreads_on_Pax304.sorted.bam \
	     SORT_ORDER=coordinate \
	     VALIDATION_STRINGENCY=SILENT \
	     TMP_DIR=tmp


## first, add read groups info to bam file, if not exist
## (check with samtools view -H sample.bam | grep '@RG')

picard-tools AddOrReplaceReadGroups \
	     I=$HOME/Style/CoGe_data/Exreads_on_Pax304.sorted.bam \
	     O=$HOME/Style/CoGe_data/Exreads_on_Pax304_RG.bam \
	     RGID=1 \
	     RGLB=lib1 \
	     RGPL=illumina \
	     RGPU=unit1 \
	     RGSM=pex1 \
	     VALIDATION_STRINGENCY=SILENT \
	     TMP_DIR=tmp

## giving mostly dummy names here as I don't have the original info at hand
## see http://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups

picard-tools MarkDuplicates \
	     I=$HOME/Style/CoGe_data/Exreads_on_Pax304_RG.bam \
	     O=Exreads_on_Pax304_marked_duplicates.bam \
	     M=Exreads_on_Pax304_dup_metrics.txt \
	     DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES \
	     VALIDATION_STRINGENCY=SILENT \
	     CLEAR_DT=false \
	     MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
	     TMP_DIR=tmp


## problem is that this file is not correctly ordered
picard-tools ReorderSam \
	     I=Exreads_on_Pax304_marked_duplicates.bam \
	     O=Exreads_on_Pax304_marked_duplicates_reordered.bam \
	     R=$HOME/Style/Pax_genome/Petunia_axillaris_304.fa \
	     CREATE_INDEX=TRUE \
	     VALIDATION_STRINGENCY=SILENT \
	     TMP_DIR=tmp


## remove intermediate files
rm -f Exreads_on_Pax304.sorted.bam
rm -f Exreads_on_Pax304_RG.bam
rm -f Exreads_on_Pax304_marked_duplicates.bam*


## if runtime problems, consider changing MAX_OPTICAL_DUPLICATE_SET_SIZE
## if memory problems, consider SORTING_COLLECTION_SIZE_RATIO
## if memory problems, consider MAX_RECORDS_IN_RAM
