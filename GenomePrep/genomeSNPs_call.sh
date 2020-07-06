#!/bin/bash

#SBATCH --mail-user=me@work
#SBATCH --mail-type=end,fail
#SBATCH --job-name="snps_call"
#SBATCH --output=%x-%A-%a.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=71:50:00
#SBATCH --mem-per-cpu=10G
#SBATCH --array=0-7
#SBATCH --partition=all


module load vital-it/7
module load UHTS/Analysis/picard-tools/2.18.11
module load UHTS/Analysis/GenomeAnalysisTK/4.1.3.0
module load UHTS/Analysis/samtools/1.8


## run in '$HOME/Style/SNPs/Pex_on_Pax/Pax_304'

# ## [check version of genome used for mapping]
# samtools view $HOME/Style/CoGe_data/Exreads_on_Pax304_marked_duplicates.bam | cut -f3 | sort | uniq > Exreads_on_Pax_scaffolds.txt

# ## prepare reference genome (if not done yet)
# cd $HOME/Style/Pax_genome/
# ## change extension of genome reference file
# scp -p Petunia_axillaris_304.faa Petunia_axillaris_304.fa
# ## create index
# samtools faidx Petunia_axillaris_304.fa
# ## create dictionary
# picard-tools CreateSequenceDictionary \
# 	     R=Petunia_axillaris_304.fa


cd $HOME/Style/SNPs/Pex_on_Pax/Pax_304


intervaldir=$HOME/Style/Pax_genome/intervals_chr/intervalLists

interval=$(ls ${intervaldir}/*scattered.interval_list | sed -n "$(( ${SLURM_ARRAY_TASK_ID} + 1))p")

OUTDIR=./interval_${SLURM_ARRAY_TASK_ID}

mkdir -p ${OUTDIR}/tmp

cd $OUTDIR



hostname; date
echo "array task $SLURM_ARRAY_TASK_ID for $(basename $interval)"
echo "output directory ${OUTDIR}"
echo ""
start_time="$(date -u +%s)"


## HaplotypeCaller in single sample mode
## -------------------------------------
GenomeAnalysisTK HaplotypeCaller \
		 -I $HOME/Style/CoGe_data/Exreads_on_Pax304_marked_duplicates_reordered.bam \
 		 -O Exreads_on_Pax304.vcf.gz \
		 -R $HOME/Style/Pax_genome/Petunia_axillaris_304.fa \
		 -L $interval \
		 -ERC NONE \
		 --native-pair-hmm-threads 8 \
		 --tmp-dir tmp



# consider options
#     --tmp-dir
#     --emit-ref-confidence / -ERC NONE (the default)
#          [ for -ERC NONE, name output file *.vcf.gz ]
#          [ for -ERC GVCF, name output file *.g.vcf.gz ]


