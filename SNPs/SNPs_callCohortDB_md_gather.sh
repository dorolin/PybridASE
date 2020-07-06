#!/bin/bash

#SBATCH --mail-user=me@work
#SBATCH --mail-type=end,fail
#SBATCH --job-name="snps_gatherVcfs"
#SBATCH --output=%x-%A-%a.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=3:50:00 ## only takes a minute
#SBATCH --mem-per-cpu=50G ## 20G should be enough
#SBATCH --partition=all



module load vital-it/7
module load UHTS/Analysis/picard-tools/2.18.11
module load UHTS/Analysis/samtools/1.8
module load UHTS/Analysis/GenomeAnalysisTK/4.1.3.0


## run in '$HOME/Style/step03_AxEx/Pax_genome/cohortDB_md/'


# ## Pax genome:
# > vcfInput.list

# for i in {0..587}; do
#     echo "./interval_${i}/cohort_md.vcf.gz" >> vcfInput.list
# done


# ## Pex genome:
# > vcfInput.list

# for i in {0..514}; do
#     echo "./interval_${i}/cohort_md.vcf.gz" >> vcfInput.list
# done



# picard-tools GatherVcfs \
# 	     I=vcfInput.list \
# 	     O=cohort_md_all.vcf.gz \
# 	     VALIDATION_STRINGENCY=SILENT

## => actually no indexes are created (WARNING)
## -> better use MergeVcfs instead


picard-tools MergeVcfs \
	     I=vcfInput.list \
	     O=cohort_md_all.vcf.gz \
	     VALIDATION_STRINGENCY=SILENT

## if not full list of contigs in vcf header, add
## --SEQUENCE_DICTIONARY=$HOME/Style/Pex_genome/Petunia_exserta_303_CoGe.dict \
