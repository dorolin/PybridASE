#!/bin/bash

#SBATCH --mail-user=me@work
#SBATCH --mail-type=fail
#SBATCH --job-name="qualityControl"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=03:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --array=1-72%20
#SBATCH --partition=all


module load vital-it/7
module load UHTS/Quality_control/fastqc/0.11.7
module load UHTS/Analysis/trimmomatic/0.36


#ln -s $HOME/Style/data/*.fastq ./
#ln -s $HOME/Style/data/TruSeq3-SE ./

#mkdir raw
#mkdir filtered

#ls -1 *.fastq | sed 's/\.fastq@//g' >files.txt

param_store=files.txt

file=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')

## check quality with fastqc
fastqc -t 2 -o raw $file".fastq"


## trim
trimmomatic SE -phred33 -threads 2 $file".fastq" $file"_cleaned.fastq" ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:15 TRAILING:15 SLIDINGWINDOW:4:20 MINLEN:60 AVGQUAL:25


## re-check quality
fastqc -t 2 -o filtered $file"_cleaned.fastq"


# trimmomatic
# -----------

## usage:
# trimmomatic SE [-threads <threads>] [-phred33 | -phred64] [-trimlog <logFile>] <input> <output> <step 1> ...

## options:
# ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>
# SLIDINGWINDOW:<windowSize>:<requiredQuality>
# MAXINFO:<targetLength>:<strictness>
# LEADING:<quality>
# TRAILING:<quality>
# CROP:<length>
# HEADCROP:<length>
# MINLEN:<length>
# AVGQUAL

## processing steps occur in the order they are specified

## example:
# trimmomatic SE -phred33 input.fq.gz output.fq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

## info:
# http://www.usadellab.org/cms/?page=trimmomatic
# https://github.com/Jeanielmj/bioinformatics-workshop/wiki/Adapter-Removal,-Trimming,-and-Filtering

# -phred33 because '/' and '>' characters don't exist in phred66

