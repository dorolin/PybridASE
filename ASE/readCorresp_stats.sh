#!/bin/bash

#SBATCH --mail-user=me@work
#SBATCH --mail-type=end,fail
#SBATCH --job-name="readCorresp_stats"
#SBATCH --output=%x-%A-%a.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=23:50:00
#SBATCH --mem-per-cpu=30G
#SBATCH --array=34,35,36,37,38,39,31,32,33
#SBATCH --partition=all


## get correct array task ids for samples of interest
## grep -n "AxExMS1" $samples ## 34,35,36
## grep -n "AxExMS2" $samples ## 37,38,39
## grep -n "AxExLS" $samples ## 31,32,33

module load vital-it/7
##module load R/latest
module load R/3.6.1

## run in '$HOME/Style/step04_AxEx/Pax_genome/ASE/'
## ## run in '$HOME/Style/step04_AxEx/Pex_genome/ASE/'

## annotated variants from ANNOVAR
refgenefile=$HOME/Style/Pax_genome/annovar/Peax304db/Peax304_refGene.txt
## refgenefile=$HOME/Style/Pex_genome/annovar/Peex304db/Peex304_refGene.txt

## GeneiASE in/out dir
indir=$HOME/Style/step04_AxEx/Pax_genome/ASE/
## indir=$HOME/Style/step04_AxEx/Pex_genome/ASE/

## outdir
outdir=${indir}

samples=$HOME/Style/step02_AxEx/files_AxEx.txt
## sample
filebase=$(cat $samples | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')

## outfile prefix
pfx="/" ## "/" is 'none'

## read correspondence
correspfile=$HOME/Style/step02_AxEx/corresp/${filebase}_aligninfo_Pax.txt
## correspfile=$HOME/Style/step02_AxEx/corresp/${filebase}_aligninfo_Pex.txt


hostname; date
start_time="$(date -u +%s)"
echo "processing $filebase"


Rscript --vanilla $HOME/Style/scripts/readCorresp_stats.R \
	${refgenefile} ${correspfile} ${indir}/${filebase} \
	${outdir}/${filebase} ${pfx}
## note: 'readCorresp_stats.R' can be optimized for speed


end_time="$(date -u +%s)"
elaps=$(dc -e "$end_time $start_time - p")
elapstime=`printf '%dh:%dm:%ds\n' $(($elaps/3600)) $(($elaps%3600/60)) $(($elaps%60))`
echo "---------------------------"
echo "readCorresp_stats took $elapstime"
echo "---------------------------"

