#!/bin/bash

#SBATCH --mail-user=me@work
#SBATCH --mail-type=end,fail
#SBATCH --job-name="mapReads_merge"
#SBATCH --output=%x-%A-%a.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:50:00
#SBATCH --mem-per-cpu=12G
#SBATCH --partition=all


module load vital-it/7
##module load UHTS/Analysis/GenomeAnalysisTK/4.1.3.0
module load UHTS/Analysis/picard-tools/2.18.11
##module load UHTS/Analysis/samtools/1.8


## run in $HOME/Style/step02_AxEx/Pax_genome/secondPass/


cohort=$HOME/Style/step02_AxEx/Pax_genome/secondPass/samples.list

infiles=$HOME/Style/step02_AxEx/files_AxEx.txt

#cat $infiles | sed 's/^/$HOME\/Style\/step02_AxEx\/Pax_genome\/secondPass\//g' | sed 's/$/\/Aligned.sortedByCoord.out.bam/g' > $cohort

bamlist=$(cat ${cohort} | sed 's/^/I=/g' | tr '\n' ' ')

TMPDIR="tmp"
mkdir -p $TMPDIR

picard-tools MergeSamFiles \
             ${bamlist} \
             O=Aligned.sortedByCoord.out.all.bam \
	     CREATE_INDEX=true \
             VALIDATION_STRINGENCY=SILENT \
	     TMP_DIR=$TMPDIR






