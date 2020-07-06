#!/bin/bash

#SBATCH --mail-user=me@work
#SBATCH --mail-type=end,fail
#SBATCH --job-name="geneiase"
#SBATCH --output=%x-%A-%a.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=8:50:00
#SBATCH --mem-per-cpu=2G
#SBATCH --array=34,35,36,37,38,39,31,32,33
#SBATCH --partition=all


## get correct array task ids for samples of interest
## grep -n "AxExMS1" $samples ## 34,35,36
## grep -n "AxExMS2" $samples ## 37,38,39
## grep -n "AxExLS" $samples ## 31,32,33

module load vital-it/7
##module load R/latest
module load R/3.5.1 ## require library("binom")


mygeneiase=$HOME/bin/geneiase-1.0.1/bin/geneiase

## run in '$HOME/Style/step04_AxEx/Pax_genome/ASE'


## bam files require read groups: these have already been added when
##  calling SNPs in step03

samples=$HOME/Style/step02_AxEx/files_AxEx.txt
filebase=$(cat $samples | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')
FILE=input.tab
BINO=`cat ${filebase}/bias.txt`

hostname; date
start_time="$(date -u +%s)"
echo "processing $filebase"


$mygeneiase -t static -i ${filebase}/${FILE} -p $BINO -m 1 \
	    -o ${filebase}/geneiase.out


end_time="$(date -u +%s)"
elaps=$(dc -e "$end_time $start_time - p")
elapstime=`printf '%dh:%dm:%ds\n' $(($elaps/3600)) $(($elaps%3600/60)) $(($elaps%60))`
echo "------------------------"
echo "GeneiASE took $elapstime"
echo "------------------------"

