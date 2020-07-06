#!/bin/bash

#SBATCH --mail-user=me@work
#SBATCH --mail-type=end,fail
#SBATCH --job-name="snps_call_md"
#SBATCH --output=%x-%A-%a.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=9:50:00
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-39
#SBATCH --partition=all


module load vital-it/7
module load UHTS/Analysis/picard-tools/2.18.11
#module load UHTS/Analysis/samtools/1.8
module load UHTS/Analysis/GenomeAnalysisTK/4.1.3.0


## run in '$HOME/Style/step03_AxEx/Pax_genome/'


infiles=$HOME/Style/step02_AxEx/files_AxEx.txt
filebase=$(cat $infiles | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')
TMPDIR=$filebase"/tmp"

file="$HOME/Style/step02_AxEx/Pax_genome/splitCigar_md/"$filebase"/Aligned.sortedByCoord.md.sp.bam"


mkdir -p $TMPDIR


## require sample name in header
picard-tools AddOrReplaceReadGroups \
             I=$file \
             O=$filebase"/Aligned.sortedByCoord.md.sp_RG.bam" \
             RGSM=$filebase \
	     RGID=$filebase \
	     RGLB=laneX \
             RGPL=illumina \
             RGPU=barcodeX \
	     CREATE_INDEX=true \
             VALIDATION_STRINGENCY=SILENT \
             TMP_DIR=$TMPDIR




## HaplotypeCaller GVCF workflow
## -----------------------------
GenomeAnalysisTK HaplotypeCaller \
		 -I $filebase"/Aligned.sortedByCoord.md.sp_RG.bam" \
		 --sample-name $filebase \
 		 -O $filebase"/Variants_md.g.vcf.gz" \
		 -R $HOME/Style/Pax_genome/Petunia_axillaris_304.fa \
		 -ERC GVCF \
		 --dont-use-soft-clipped-bases true \
		 --native-pair-hmm-threads 8 \
		 --tmp-dir $TMPDIR


## When working with PCR-free data, be sure to set `-pcr_indel_model NONE`

## In GVCF mode, -stand_call_conf is automatically set to zero.
## Call confidence thresholding will then be performed in the subsequent
## GenotypeGVCFs command

## From best practice: adjust the minimum phred-scaled confidence threshold
## for calling variants to 20 (i.e., this will then be done later)


## don't use -R $HOME/Style/Pax_genome/Petunia_axillaris_304_PexIUPAC.fa (only one REF allele allowed in .vcf)



## consider options
##     --tmp-dir
##     --emit-ref-confidence / -ERC NONE (the default)
##          [ for -ERC NONE, name output file *.vcf.gz ]
##          [ for -ERC GVCF, name output file *.g.vcf.gz ]


## Output:
# Either a VCF or GVCF file with raw, unfiltered SNP and indel calls. Regular VCFs must be filtered either by variant recalibration (Best Practice) or hard-filtering before use in downstream analyses. If using the GVCF workflow, the output is a GVCF file that must first be run through GenotypeGVCFs and then filtering before further analysis.
