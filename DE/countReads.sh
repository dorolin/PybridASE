#!/bin/bash

#SBATCH --mail-user=me@work
#SBATCH --mail-type=end,fail
#SBATCH --job-name="countReads"
#SBATCH --output=%x-%A-%a.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:50:00
#SBATCH --mem-per-cpu=2G
#SBATCH --partition=all


#module load vital-it/7
#module load UHTS/Analysis/subread/1.6.0

## installed a new version of the subreads package in
## $HOME/bin/subread-1.6.4-source/bin/

## run in $HOME/Style/step04_AxEx/Pax_genome/DE/

featureCounts -a $HOME/Style/Pax_genome/Peax304_final_edited.gtf -t exon -g gene_id -s 2 -o featureCounts_results.txt --byReadGroup -T 4 $HOME/Style/step02_AxEx/Pax_genome/secondPass/Aligned.sortedByCoord.out.all.bam


## consider adding -M

featureCounts -a $HOME/Style/Pax_genome/Peax304_final_edited.gtf -t exon -g gene_id -s 2 -M -o featureCounts_results_Mfull.txt --byReadGroup -T 4 $HOME/Style/step02_AxEx/Pax_genome/secondPass/Aligned.sortedByCoord.out.all.bam


featureCounts -a $HOME/Style/Pax_genome/Peax304_final_edited.gtf -t exon -g gene_id -s 2 -M --fraction -o featureCounts_results_Mfrac.txt --byReadGroup -T 4 $HOME/Style/step02_AxEx/Pax_genome/secondPass/Aligned.sortedByCoord.out.all.bam


## try what happens when setting -s 0 (default)

featureCounts -a $HOME/Style/Pax_genome/Peax304_final_edited.gtf -t exon -g gene_id -s 0 -o featureCounts_results_s0.txt --byReadGroup -T 4 $HOME/Style/step02_AxEx/Pax_genome/secondPass/Aligned.sortedByCoord.out.all.bam



## -----------------------------------------------------------------
## featureCounts

## some options:

# input-files ## aligned reads .bam/.sam

# -a <string> ## genomic features .gtf/.gff/.saf

#   .saf format:
#   GeneID Chr Start End Strand
#   497097 chr1 3204563 3207049 -
#   497097 chr1 3411783 3411982 -
#   497097 chr1 3660633 3661579 -
#   100503874 chr1 3637390 3640590 -
#   100503874 chr1 3648928 3648985 -
#   100038431 chr1 3670236 3671869 -

# -F ## annotation file format 'GTF' 'SAF'

# -f ## summarization performed at feature level (e.g., exon-level)
# -g <string> ## attribute to group features (e.g., exons) into
#             ## meta-features (e.g., genes) [default 'gene_id']
# -t <string> ## feature type for read counting [default 'exon']

# -s <int or string> ## perform strand-specific read counting:
#                  ## 0 (unstranded), 1 (stranded), 2 (reversely stranded)
# 		 ## use comma-sep list for several input files

# --byReadGroup ## count reads by read group

# -Q <int> ## minimum mapping quality score [default 0]

# --minOverlap <int> ## minimum number of overlapping bases of read
#                    ## with feature [default 1]
# --nonOverlap <int> ## maximum number of non-overlapping bases
# --fracOverlap <float> ## minimum fraction of overlapping bases [0,1]
#                       ## (better for reads with variable lengths)

# ## multi-mapping reads need NH tag and 0x000 flag (STAR does this)
# -M ## fully count every alignment of multi-mapping reads
# -M --fraction ## count each alignment fractionally
# -M --primary ## consider only primary alignments for counting

# --ignoreDup ## reads marked as duplicates will be ignored
#             ## (0x400 flag)

# -G <string> ## genome in FASTA format used to make .bam/.sam files
#             ## good to use with -J option
# -J ## count the number of reads supporting each exon-exon junction

# -T <int> ## number of threads [between 1 and 32]

# -o <string> ## name of output file


# Example: featureCounts -a annotation.gtf -t exon -g gene_id \
# -o counts.txt mapping_results1.bam mapping_results2.bam

