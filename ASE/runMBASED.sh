#!/bin/bash

#SBATCH --mail-user=me@work
#SBATCH --mail-type=end,fail
#SBATCH --job-name="mbased"
#SBATCH --output=%x-%A-%a.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=23:50:00
#SBATCH --mem-per-cpu=12G
#SBATCH --array=34,35,36,37,38,39,31,32,33
#SBATCH --partition=all

## Note: time and memory requirements strongly depend on
##       number of simulations in 'runMBASED.R' to
##       compute p-values: <2h and <1.5G for nsims=10^5


## get correct array task ids for samples of interest
## grep -n "AxExMS1" $samples ## 34,35,36
## grep -n "AxExMS2" $samples ## 37,38,39
## grep -n "AxExLS" $samples ## 31,32,33

module load vital-it/7
##module load R/latest
module load R/3.6.1 ## require library("MBASED")

## run in '$HOME/Style/step04_AxEx/Pax_genome/ASE/'
## ## run in '$HOME/Style/step04_AxEx/Pex_genome/ASE/'

## set number of simulations
nsims=100000

## set rho, the dispersion parameter
## default is 0: uses binomial rather than betabinomial
## 0.002 in tutorial and 0.004 in paper
## note: geneiase default is -r 0.012
##       (the greater the value the more conservative the test)
rho=0.012

## annotated variants from ANNOVAR
annosnpsfile=$HOME/Style/Pax_genome/annovar/snps_filtered_sel_anno.variant_function
##annosnpsfile=$HOME/Style/Pex_genome/annovar/snps_filtered_sel_anno.variant_function

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


hostname; date
start_time="$(date -u +%s)"
echo "processing $filebase"


Rscript --vanilla $HOME/Style/scripts/runMBASED.R \
	${annosnpsfile} ${indir}/${filebase} \
	${outdir}/${filebase} ${pfx} \
	${nsims} ${rho}
## a couple of MBASED settings are hard-coded in 'runMBASED.R'


end_time="$(date -u +%s)"
elaps=$(dc -e "$end_time $start_time - p")
elapstime=`printf '%dh:%dm:%ds\n' $(($elaps/3600)) $(($elaps%3600/60)) $(($elaps%60))`
echo "---------------------------"
echo "runMBASED took $elapstime"
echo "---------------------------"

