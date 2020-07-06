#!/bin/bash

#SBATCH --mail-user=me@work
#SBATCH --mail-type=end,fail
#SBATCH --job-name="stats_multiQC"
#SBATCH --output=%x-%A-%a.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --partition=all


module load vital-it/7
module load UHTS/Analysis/MultiQC/1.8

## run in $HOME/Style/step02_AxEx/Pax_genome/secondPass/

## compute stats on all collected QC results

multiqc --dirs --filename multiQC --outdir stats .


## multiqc
## Usage: multiqc [OPTIONS] <analysis directory>



