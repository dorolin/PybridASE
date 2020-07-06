#!/bin/bash

#SBATCH --mail-user=me@work
#SBATCH --mail-type=end,fail
#SBATCH --job-name="count_reads"
#SBATCH --output=%x-%A-%a.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:50:00
#SBATCH --mem-per-cpu=10G
#SBATCH --array=1-39
#SBATCH --partition=all


module load vital-it/7
module load UHTS/Analysis/GenomeAnalysisTK/4.1.3.0


## run in '$HOME/Style/step04_AxEx/Pax_genome/ASE'


## bam files require read groups: these have already been added when
##  calling SNPs in step03

samples=$HOME/Style/step02_AxEx/files_AxEx.txt
bamdir=$HOME/Style/step03_AxEx/Pax_genome
filebase=$(cat $samples | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')
##FILE=${bamdir}/${filebase}/Aligned.sortedByCoord.sp_RG.bam
FILE=${bamdir}/${filebase}/Aligned.sortedByCoord.md.sp_RG.bam


mkdir -p ${filebase}/tmp



hostname; date
start_time="$(date -u +%s)"
echo "processing $filebase"


## not removing duplicates, as recommended
GenomeAnalysisTK ASEReadCounter \
		 -R $HOME/Style/Pax_genome/Petunia_axillaris_304.fa \
		 -I $FILE \
		 -V $HOME/Style/step04_AxEx/Pax_genome/SNPs/cohort_md_snps_filtered_sel.vcf.gz \
		 -O ${filebase}/counts.txt \
		 --disable-read-filter NotDuplicateReadFilter \
		 --min-mapping-quality 60 \
		 --min-base-quality 25 \
		 --tmp-dir ${filebase}/tmp


## also making test run with actually removing duplicates
GenomeAnalysisTK ASEReadCounter \
		 -R $HOME/Style/Pax_genome/Petunia_axillaris_304.fa \
		 -I $FILE \
		 -V $HOME/Style/step04_AxEx/Pax_genome/SNPs/cohort_md_snps_filtered_sel.vcf.gz \
		 -O ${filebase}/counts_md.txt \
		 --min-mapping-quality 60 \
		 --min-base-quality 25 \
		 --tmp-dir ${filebase}/tmp




end_time="$(date -u +%s)"
elaps=$(dc -e "$end_time $start_time - p")
elapstime=`printf '%dh:%dm:%ds\n' $(($elaps/3600)) $(($elaps%3600/60)) $(($elaps%60))`
echo "------------------------------"
echo "ASEReadcounter took $elapstime"
echo "------------------------------"

