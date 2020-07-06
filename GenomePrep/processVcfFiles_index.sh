#!/bin/bash

#SBATCH --mail-user=me@work
#SBATCH --mail-type=end,fail
#SBATCH --job-name="index"
#SBATCH --output=%x-%j.out
#SBATCH --nodes=1
#SBATCH --time=01:30:00
#SBATCH --mem-per-cpu=5000M
#SBATCH --partition=all


module load vital-it/7
module load UHTS/Analysis/samtools/1.8
module load UHTS/Analysis/vcftools/0.1.15

## run in $HOME/Style/SNPs/Pex_on_Pax/Pax_304

## manipulate vcf file used for STAR mapping in WASP-mode
##  to re-code the ALT/ALT homs to REF/ALT hets
$HOME/Style/scripts/manipulate_vcf.pl \
    -i Exreads_on_Pax304_snps_filtered_2_sel.vcf


# ## require compressed vcf
# bgzip -c Exreads_on_Pax304_snps_filtered_2_sel_het.vcf \
#       > Exreads_on_Pax304_snps_filtered_2_sel_het.vcf.gz
bcftools view -O z Exreads_on_Pax304_snps_filtered_2_sel_het.vcf \
	      > Exreads_on_Pax304_snps_filtered_2_sel_het.vcf.gz


## make index
bcftools index -t Exreads_on_Pax304_snps_filtered_2_sel_het.vcf.gz
