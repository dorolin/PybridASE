#!/bin/bash

#SBATCH --mail-user=me@work
#SBATCH --mail-type=end,fail
#SBATCH --job-name="split_intervals"
#SBATCH --output=%x-%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:50:00
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=all



module load vital-it/7
module load UHTS/Analysis/picard-tools/2.18.11
module load UHTS/Analysis/samtools/1.8
module load UHTS/Analysis/GenomeAnalysisTK/4.1.3.0


## run in [A] '$HOME/Style/Pax_genome/intervals/'
##        -> use --scatter-count 200 (used for RNA-seq SNP calling)
## OR run in [B] '$HOME/Style/Pax_genome/intervals_chr/'
##        -> use --scatter-count 8 (used for DNA-seq SNP calling)

reference=$HOME/Style/Pax_genome/Petunia_axillaris_304.fa
## not using PexIUPAC here to avoid potential problems with not
## having only ACGT or N characters (not tested), this will work as
## long no indels were added to PexIUPAC (both .fa have the same wc)


scaffolds=$HOME/Style/Pax_genome/intervals/output.interval_list ## [A]
##scaffolds=$HOME/Style/Pax_genome/intervals_chr/output.interval_list ## [B]

# ## create dictionary (if not yet existing; not sure if needed here)
# dict="${reference%.fa}.dict"
# picard-tools CreateSequenceDictionary \
#              R=$reference \
# 	       O=$dict


# ## create fasta index (if not yet existing)
# samtools faidx $reference


## [A]

## split reference by N
## --------------------
picard-tools ScatterIntervalsByNs \
	     REFERENCE=$reference \
	     OUTPUT_TYPE=ACGT \
	     MAX_TO_MERGE=999 \
	     OUTPUT=$scaffolds



## make lists of intervals
## -----------------------
GenomeAnalysisTK SplitIntervals \
		 -R $reference \
		 -L $scaffolds \
		 -O intervalLists \
		 --scatter-count 200 \
		 --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION


# ## [B]

# ## split reference by N
# ## --------------------
# picard-tools ScatterIntervalsByNs \
# 	     REFERENCE=$reference \
# 	     OUTPUT_TYPE=ACGT \
# 	     MAX_TO_MERGE=999999 \
# 	     OUTPUT=$scaffolds



# ## make lists of intervals
# ## -----------------------
# GenomeAnalysisTK SplitIntervals \
# 		 -R $reference \
# 		 -L $scaffolds \
# 		 -O intervalLists \
# 		 --scatter-count 8 \
# 		 --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION



## specify list of intervals:
## https://software.broadinstitute.org/gatk/documentation/article?id=11009
## -L intervals.list (or intervals.interval_list, or intervals.bed) when specifying a text file containing intervals
## GATK-style .list or .intervals: <chr>:<start>-<stop> (only <chr> is strictly required)

