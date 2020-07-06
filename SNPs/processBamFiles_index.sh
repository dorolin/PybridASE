#!/bin/bash

#SBATCH --mail-user=me@work
#SBATCH --mail-type=end,fail
#SBATCH --job-name="index"
#SBATCH --output=%x-%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:30:00
#SBATCH --mem-per-cpu=5000M
#SBATCH --array=1-39
#SBATCH --partition=all


## the specified resources are a huge overkill,
## runs in 30s and needs 1M


module load vital-it/7
module load UHTS/Analysis/samtools/1.8


## run in $HOME/Style/step02_AxEx/Pax_genome/secondPass


infiles=$HOME/Style/step02_AxEx/files_AxEx.txt
filebase=$(cat $infiles | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')

cd $filebase

samtools index -@ 3 Aligned.sortedByCoord.out.bam

