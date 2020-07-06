#!/bin/bash

#SBATCH --mail-user=me@work
#SBATCH --mail-type=end,fail
#SBATCH --job-name="manipulateRef"
#SBATCH --output=%x-%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:50:00
#SBATCH --mem-per-cpu=12G
#SBATCH --partition=all


module load vital-it/7
module load UHTS/Analysis/GenomeAnalysisTK/4.1.3.0
module load UHTS/Analysis/picard-tools/2.18.11
module load UHTS/Analysis/samtools/1.8


## run in $HOME/Style/Pax_genome/

## name of modified fasta
file="Petunia_axillaris_304_PexIUPAC.fa"


## note that SNPs file needs to be indexed, which can be done with
## bcftools index -t file.vcf.gz ('processVcfFiles_index.sh')

GenomeAnalysisTK FastaAlternateReferenceMaker \
	         -R $HOME/Style/Pax_genome/Petunia_axillaris_304.fa \
		 -O $file \
                 -V $HOME/Style/SNPs/Pex_on_Pax/Pax_304/Exreads_on_Pax304_snps_filtered_2_sel_het.vcf.gz \
                 --use-iupac-sample pex1 \
                 --line-width 80


## FastaAlternateReferenceMaker changes original fasta header lines
## -> have to modify output accordingly

$HOME/Style/scripts/manipulate_fasta.pl -i $file

## makes $file".edited" as default output
mv -f $file".edited" $file


## make new .dict and .fai files
rm -f "$(basename ${file} .fa).dict"
rm -f $file".fai"

picard-tools CreateSequenceDictionary \
	     R=$file \
	     O="$(basename ${file} .fa).dict"

samtools faidx $file

