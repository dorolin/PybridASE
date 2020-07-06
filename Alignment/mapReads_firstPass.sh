#!/bin/bash

#SBATCH --mail-user=me@work
#SBATCH --mail-type=end,fail
#SBATCH --job-name="mapReads_firstPass"
#SBATCH --output=%x-%A-%a.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1:50:00
#SBATCH --mem-per-cpu=12G
#SBATCH --array=1-39
#SBATCH --partition=all


## run in $HOME/Style/step02_AxEx/Pax_genome/firstPass/

## added here
## --outFilterScoreMinOverLread 0.3
## --outFilterMatchNminOverLread 0.3
## as recommended in
## https://github.com/alexdobin/STAR/issues/169
## when many reads were not aligned because they are "too short"
## (not sure if this will have conflicts with
##  --outFilterMismatchNoverReadLmax)

## added --outFilterMatchNmin 20

## changed --outFilterMismatchNoverReadLmax from 0.06 to 0.04


module load vital-it/7
module load UHTS/Aligner/STAR/2.6.0c


# get fastq file
infiles=$HOME/Style/step02_AxEx/files_AxEx.txt
filebase=$(cat $infiles | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')
FILE=$(basename ${filebase}"_cleaned.fastq")
ID="ID:"$filebase

mkdir -p "./"$filebase

STAR --genomeDir $HOME/Style/Pax_genome/Pax_304_generated \
     --runThreadN 4 \
     --readFilesPrefix $HOME/Style/step01/ \
     --readFilesIn $FILE \
     --outFileNamePrefix "./"$filebase"/" \
     --clip5pNbases 6 \
     --alignIntronMin 20 \
     --alignIntronMax 1500000 \
     --outFilterMismatchNmax 999 \
     --outFilterMismatchNoverReadLmax 0.04 \
     --outFilterScoreMinOverLread 0.3 \
     --outFilterMatchNminOverLread 0.3 \
     --outFilterMatchNmin 20 \
     --outFilterMultimapNmax 20 \
     --outFilterMultimapScoreRange 2 \
     --outFilterType BySJout \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattrRGline $ID \
     --outSAMprimaryFlag AllBestScore \
     --outSAMattributes NH HI AS nM NM MD vA vG vW \
     --waspOutputMode SAMtag \
     --varVCFfile $HOME/Style/SNPs/Pex_on_Pax/Pax_304/Exreads_on_Pax304_snps_filtered_2_sel.vcf



## STAR

## basic options:

# --runThreadN NumberOfThreads
# --genomeDir /path/to/genomeDir
# --readFilesIn /path/to/read1 [/path/to/read2]

## (some) advanced options:

# --readFilesPrefix prefix for read files names
# --sjdb* options from genome prep step can alternatively be specified
#                 at mapping step
# --outFilterType BySJout (reduces the number of ”spurious” junctions)
# --outFilterMultimapNmax int (max number of multiple alignments allowed
# 			       for a read)
# --outSAMtype BAM SortedByCoordinate (similar to samtools sort command)
# --twopassMode Basic (include novel junctions)
# --quantMode GeneCounts (as with htseq-count)
# --outReadsUnmapped Fastx (default is None)
# --runRNGseed int (random number generator seed)


## various options to control SAM/BAM output with --outSAMattr*

## specify multiple read group (RG) tags to be added in BAM output:
# --outSAMattrRGline ID:xxx , ID:zzz
