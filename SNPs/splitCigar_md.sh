#!/bin/bash

#SBATCH --mail-user=me@work
#SBATCH --mail-type=end,fail
#SBATCH --job-name="splitCigar"
#SBATCH --output=%x-%A-%a.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:50:00
#SBATCH --mem-per-cpu=12G
#SBATCH --array=1-39
#SBATCH --partition=all



module load vital-it/7
module load UHTS/Analysis/GenomeAnalysisTK/4.1.3.0


## run in $HOME/Style/step02_AxEx/Pax_genome/splitCigar_md


infiles=$HOME/Style/step02_AxEx/files_AxEx.txt
filebase=$(cat $infiles | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')
TMPDIR=${filebase}/tmp

mkdir -p $TMPDIR


## input bam need to be indexed
## (automatically when running 'processBamFiles.sh')

GenomeAnalysisTK SplitNCigarReads \
		 -I $HOME/Style/step02_AxEx/Pax_genome/secondPass/${filebase}/Aligned.sortedByCoord.md.bam \
		 -O ${filebase}/Aligned.sortedByCoord.md.sp.bam \
		 -R $HOME/Style/Pax_genome/Petunia_axillaris_304_PexIUPAC.fa \
		 --max-mismatches-in-overhang 2 \
		 --tmp-dir $TMPDIR

