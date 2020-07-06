#!/bin/bash

#SBATCH --mail-user=me@work
#SBATCH --mail-type=end,fail
#SBATCH --job-name="snps_cohortDB_md"
#SBATCH --output=%x-%A-%a.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=3:50:00
#SBATCH --mem-per-cpu=50G
#SBATCH --array=0-587%20 ## not getting >11 in (Resources)
#SBATCH --partition=all


## most of the runs finish quickly and also do not seem to require
## a lot of real memory, however the virtual memory can be huge

## failed runs may not result in any slurm error message
## i.e., grep "OutOfMemoryError" or "IllegalStateException"


module load vital-it/7
module load UHTS/Analysis/picard-tools/2.18.11
#module load UHTS/Analysis/samtools/1.8
module load UHTS/Analysis/GenomeAnalysisTK/4.1.3.0


## run in '$HOME/Style/step03_AxEx/Pax_genome/cohortDB_md/'



intervaldir=$HOME/Style/Pax_genome/intervals/intervalLists

interval=$(ls ${intervaldir}/*scattered.interval_list | sed -n "$(( ${SLURM_ARRAY_TASK_ID} + 1))p")

cohort=$HOME/Style/step03_AxEx/Pax_genome/cohortDB_md/samples_md.txt

infiles=$HOME/Style/step02_AxEx/files_AxEx.txt

## for --sample-name-map argument with GenomicsDBImport
## (file ending should be .txt)
#paste -d "\t" $infiles <(cat $infiles | sed 's/^/$HOME\/Style\/step03_AxEx\/Pax_genome\//g' | sed 's/$/\/Variants_md.g.vcf.gz/g') > $cohort

OUTDIR=./interval_${SLURM_ARRAY_TASK_ID}
TMPDIR=${OUTDIR}/tmp
DB=${OUTDIR}/DB

mkdir -p $TMPDIR


hostname; date
echo "array task $SLURM_ARRAY_TASK_ID for $(basename $interval)"
echo "output directory ${OUTDIR}"
echo ""
start_time="$(date -u +%s)"

## genotyping with GVCF workflow
# ## -----------------------------
# GenomeAnalysisTK GenomicsDBImport \
# 		 -R $HOME/Style/Pax_genome/Petunia_axillaris_304.fa \
# 		 --sample-name-map $cohort \
# 		 --validate-sample-name-map true \
# 		 -L $interval \
# 		 --merge-input-intervals true \
# 		 --interval-padding 0 \
# 		 --max-num-intervals-to-import-in-parallel 20 \
# 		 --genomicsdb-workspace-path $DB \
# 		 --tmp-dir $TMPDIR



## to run with --java-options use local version of gatk
~/bin/gatk-4.1.3.0/gatk --java-options "-Xmx150G -XX:+UseParallelGC -XX:ParallelGCThreads=4" GenomicsDBImport \
			-R $HOME/Style/Pax_genome/Petunia_axillaris_304.fa \
			--sample-name-map $cohort \
			--validate-sample-name-map true \
			-L $interval \
			--merge-input-intervals true \
			--interval-padding 0 \
			--max-num-intervals-to-import-in-parallel 4 \
			--genomicsdb-workspace-path $DB \
			--tmp-dir $TMPDIR





end_time="$(date -u +%s)"
elaps=$(dc -e "$end_time $start_time - p")
elapstime=`printf '%dh:%dm:%ds\n' $(($elaps/3600)) $(($elaps%3600/60)) $(($elaps%60))`
echo "--------------------------------"
echo "GenomicsDBImport took $elapstime"
echo "--------------------------------"
start_time2="$(date -u +%s)"


GenomeAnalysisTK GenotypeGVCFs \
		 -R $HOME/Style/Pax_genome/Petunia_axillaris_304.fa \
		 -V gendb://$DB \
		 -O ${OUTDIR}/cohort_md.vcf.gz \
		 -L $interval \
		 --merge-input-intervals true \
		 --interval-padding 0 \
		 --standard-min-confidence-threshold-for-calling 20 \
		 --tmp-dir $TMPDIR


end_time="$(date -u +%s)"
elaps=$(dc -e "$end_time $start_time2 - p")
elapstime=`printf '%dh:%dm:%ds\n' $(($elaps/3600)) $(($elaps%3600/60)) $(($elaps%60))`
echo "-----------------------------"
echo "GenotypeGVCFs took $elapstime"
echo "-----------------------------"


date


## In GVCF mode, -stand_call_conf is automatically set to zero.
## Call confidence thresholding will then be performed in the subsequent
## GenotypeGVCFs command

## From best practice: adjust the minimum phred-scaled confidence threshold
## for calling variants to 20 (i.e., this will then be done later)

## don't use -R $HOME/Style/Pax_genome/Petunia_axillaris_304_PexIUPAC.fa (only one REF allele allowed in .vcf)



## intervals GenomicsDBImport:
## --interval-merging-rule ALL: side-by-side but non-overlapping intervals are merged
## --merge-input-intervals false: set to true if using large lists of intervals
## --interval-set-rule UNION: merging approach for intervals
## --interval-padding 0: padding (in bp) to add to each interval
## --max-num-intervals-to-import-in-parallel 1: set to higher value to improve performance, but require more memory


## intervals GenotypeGVCFs:
## --interval-merging-rule ALL: side-by-side but non-overlapping intervals are merged
## --merge-input-intervals false: set to true if using large lists of intervals
## --interval-set-rule UNION: merging approach for intervals
## --interval-padding 0: padding (in bp) to add to each interval


## specify list of intervals:
## https://software.broadinstitute.org/gatk/documentation/article?id=11009
## -L intervals.list (or intervals.interval_list, or intervals.bed) when specifying a text file containing intervals
## GATK-style .list or .intervals: <chr>:<start>-<stop> (only <chr> is strictly required)


