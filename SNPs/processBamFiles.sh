#!/bin/bash

#SBATCH --mail-user=me@work
#SBATCH --mail-type=end,fail
#SBATCH --job-name="processBamFiles"
#SBATCH --output=%x-%A-%a.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:25:00
#SBATCH --mem-per-cpu=3G
#SBATCH --array=1-39
#SBATCH --partition=all


module load vital-it/7
module load UHTS/Analysis/picard-tools/2.18.11
#module load UHTS/Analysis/GenomeAnalysisTK/4.1.3.0



## run in $HOME/Style/step02_AxEx/Pax_genome/secondPass


infiles=$HOME/Style/step02_AxEx/files_AxEx.txt
filebase=$(cat $infiles | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')
TMPDIR=$filebase"/tmp"

mkdir -p $TMPDIR


picard-tools MarkDuplicates \
	     I=$filebase"/Aligned.sortedByCoord.out.bam" \
	     O=$filebase"/Aligned.sortedByCoord.md.bam" \
	     M=$filebase"/Aligned.sortedByCoord.metrics" \
	     CREATE_INDEX=true \
	     CLEAR_DT=false \
	     VALIDATION_STRINGENCY=SILENT \
	     DUPLICATE_SCORING_STRATEGY=RANDOM \
	     READ_NAME_REGEX=null \
	     MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
	     TMP_DIR=$TMPDIR
