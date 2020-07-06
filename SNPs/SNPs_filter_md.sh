#!/bin/bash

#SBATCH --mail-user=me@work
#SBATCH --mail-type=end,fail
#SBATCH --job-name="snps_filter_md"
#SBATCH --output=%x-%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:50:00 ## 1:00:00 should be enough
#SBATCH --mem-per-cpu=10G
#SBATCH --partition=all


module load vital-it/7
module load UHTS/Analysis/picard-tools/2.18.11
module load UHTS/Analysis/GenomeAnalysisTK/4.1.3.0
module load UHTS/Analysis/samtools/1.8

## NOTE:
## had unsolvable issues with UHTS/Analysis/GenomeAnalysisTK/4.1.0.0
## for VariantFiltration, which seems to run just fine with local
## $HOME/bin/gatk-4.1.3.0/gatk


## run in '$HOME/Style/step04_AxEx/Pax_genome/SNPs'


rawvcf="$HOME/Style/step03_AxEx/Pax_genome/cohortDB_md/cohort_md_all.vcf.gz"

reference="$HOME/Style/Pax_genome/Petunia_axillaris_304.fa"


## subset to biallelic SNPs-only

GenomeAnalysisTK SelectVariants \
		 -V $rawvcf \
		 --select-type-to-include SNP \
		 --restrict-alleles-to BIALLELIC \
		 -O cohort_md_snps.vcf.gz



## hard-filter SNPs

# # # ## as recommended in https://gatkforums.broadinstitute.org/gatk/discussion/3891/calling-variants-in-rnaseq
# # # $HOME/bin/gatk-4.1.3.0/gatk VariantFiltration \
# # # 		 -R $reference \
# # # 		 -V cohort_md_snps.vcf.gz \
# # # 		 --window 35 --cluster 3 \
# # # 		 --filter-name "FS30" --filter-expression "FS > 30.0" \
# # # 		 --filter-name "QD2" --filter-expression "QD < 2.0" \
# # # 		 -O cohort_md_snps_filtered.vcf.gz

$HOME/bin/gatk-4.1.3.0/gatk VariantFiltration \
		 -R $reference \
		 -V cohort_md_snps.vcf.gz \
		 --window 35 --cluster 4 \
		 --filter-name "FS30" --filter-expression "FS > 30.0" \
		 --filter-name "QD2" --filter-expression "QD < 2.0" \
		 --filter-name "MQ40" --filter-expression "MQ < 40.0" \
		 --filter-name "DP20" --filter-expression "DP < 20" \
		 --filter-name "AF25" --filter-expression "AF < 0.25" \
		 --filter-name "AF75" --filter-expression "AF > 0.75" \
		 -O cohort_md_snps_filtered.vcf.gz


## variants to table for filter stats

GenomeAnalysisTK VariantsToTable \
		 -V cohort_md_snps_filtered.vcf.gz \
		 -F CHROM -F POS -F FILTER -F AC -F AF -F AN -F DP \
		 -F QD -F QUAL -F SOR -F FS -F MQ -F MQRankSum \
		 --show-filtered TRUE \
		 -O cohort_md_snps_filtered.tsv



## select variants

$HOME/bin/gatk-4.1.3.0/gatk SelectVariants \
	       -V cohort_md_snps_filtered.vcf.gz \
               -O cohort_md_snps_filtered_sel.vcf.gz \
	       --exclude-filtered TRUE \
	       -select "DP >= 20" \
	       -select "QD >= 2.0"
## (-select "DP >= 20" to exclude the ones with missing data
##  that received a PASS during VariantFiltration)



## note: a better control might be possible with bcftools annotate -e,
##       see https://samtools.github.io/bcftools/bcftools.html#annotate
##       and https://samtools.github.io/bcftools/bcftools.html#expressions
