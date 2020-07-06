#!/bin/bash

#SBATCH --mail-user=me@work
#SBATCH --mail-type=end,fail
#SBATCH --job-name="snps_filter"
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



## run in '$HOME/Style/SNPs/Pex_on_Pax/Pax_304'

mkdir -p tmp


## subset to SNPs-only

GenomeAnalysisTK SelectVariants \
		 -V Exreads_on_Pax304_all.vcf.gz \
		 -select-type SNP \
		 -O Exreads_on_Pax304_snps.vcf.gz


# # # ## subset to indels-only

# # # GenomeAnalysisTK SelectVariants \
# # # 		 -V Exreads_on_Pax304_all.vcf.gz \
# # # 		 -select-type INDEL \
# # # 		 -O Exreads_on_Pax304_indels.vcf.gz


## hard-filter SNPs

$HOME/bin/gatk-4.1.3.0/gatk VariantFiltration \
		 -V Exreads_on_Pax304_snps.vcf.gz \
		 -O Exreads_on_Pax304_snps_filtered.vcf.gz \
		 --filter-name "QD2" --filter-expression "QD < 2.0" \
		 --filter-name "QUAL30" --filter-expression "QUAL < 30.0" \
		 --filter-name "SOR3" --filter-expression "SOR > 3.0" \
		 --filter-name "FS60" --filter-expression "FS > 60.0" \
		 --filter-name "MQ42" --filter-expression "MQ < 42.0" \
		 --filter-name "MQRankSum12.5" --filter-expression "MQRankSum < -12.5" \
		 --filter-name "ReadPosRankSum8" --filter-expression "ReadPosRankSum < -8.0" \
		 --filter-name "DP30" --filter-expression "DP < 30" \
		 --filter-name "DP200" --filter-expression "DP > 200"



# # # ## hard-filter indels

# # # GenomeAnalysisTK VariantFiltration \
# # # 		 -V Exreads_on_Pax304_indels.vcf.gz \
# # # 		 --filter-name "QD2" --filter-expression "QD < 2.0" \
# # # 		 --filter-name "QUAL30" --filter-expression "QUAL < 30.0" \
# # # 		 --filter-name "FS200" --filter-expression "FS > 200.0" \
# # # 		 --filter-name "ReadPosRankSum20" --filter-expression "ReadPosRankSum < -20.0" \
# # # 		 -O Exreads_on_Pax304_indels_filtered.vcf.gz



# # # ## combine filtered variants

# # # picard-tools SortVcf \
# # # 	     I=Exreads_on_Pax304_snps_filtered.vcf.gz \
# # # 	     I=Exreads_on_Pax304_indels_filtered.vcf.gz \
# # # 	     O=Exreads_on_Pax304_filtered.vcf.gz \
# # # 	     VALIDATION_STRINGENCY=SILENT \
# # # 	     TMP_DIR=tmp


## variants to table for filter stats

$HOME/bin/gatk-4.1.3.0/gatk VariantsToTable \
		 -V Exreads_on_Pax304_snps_filtered.vcf.gz \
		 -F CHROM -F POS -F FILTER -F AC -F AF -F AN -F DP \
		 -F QD -F QUAL -F SOR -F FS -F MQ -F MQRankSum \
		 --show-filtered TRUE \
		 -O Exreads_on_Pax304_snps_filtered.tsv



## select variants

$HOME/bin/gatk-4.1.3.0/gatk SelectVariants \
	       -V Exreads_on_Pax304_snps_filtered.vcf.gz \
               -O Exreads_on_Pax304_snps_filtered_sel.vcf.gz \
	       --exclude-filtered TRUE \
	       -select "DP > 29"
## (-select "DP > 29" to exclude the ones with missing data
##  that received a PASS during VariantFiltration)


## second hard-filter SNPs

$HOME/bin/gatk-4.1.3.0/gatk VariantFiltration \
	       -V Exreads_on_Pax304_snps_filtered_sel.vcf.gz \
	       -O Exreads_on_Pax304_snps_filtered_2.vcf.gz \
	       --filter-name "RefHet" --filter-expression "AC == 2" \
	       --invert-filter-expression


## second select variants

$HOME/bin/gatk-4.1.3.0/gatk SelectVariants \
	       -V Exreads_on_Pax304_snps_filtered_2.vcf.gz \
               -O Exreads_on_Pax304_snps_filtered_2_sel.vcf.gz \
	       --exclude-filtered TRUE


## need uncompressed version for STAR
bcftools view -O v Exreads_on_Pax304_snps_filtered_2_sel.vcf.gz >Exreads_on_Pax304_snps_filtered_2_sel.vcf


## note: a better control might be possible with bcftools annotate -e,
##       see https://samtools.github.io/bcftools/bcftools.html#annotate
##       and https://samtools.github.io/bcftools/bcftools.html#expressions
