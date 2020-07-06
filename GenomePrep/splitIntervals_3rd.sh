#!/bin/bash

#SBATCH --mail-user=me@work
#SBATCH --mail-type=end,fail
#SBATCH --job-name="split_intervals_3rd"
#SBATCH --output=%x-%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:50:00
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=all
#SBATCH --workdir=$HOME/Style/Pax_genome/intervals/


## this script splits scaffold list from first round in smaller subsets,
## and moves the old list from second round to 'intervalLists/split_2nd/',
## and stores new lists in 'intervalLists/split_3rd/';
## a counter increments the new names accordingly


module load vital-it/7
module load UHTS/Analysis/picard-tools/2.18.11
#module load UHTS/Analysis/samtools/1.8
module load UHTS/Analysis/GenomeAnalysisTK/4.1.3.0


## run in '$HOME/Style/Pax_genome/intervals/'

## this trial: intervalLists


reference=$HOME/Style/Pax_genome/Petunia_axillaris_304.fa
## not using PexIUPAC here to avoid potential problems with not
## having only ACGT or N characters (not tested), this will work as
## long no indels were added to PexIUPAC (both .fa have the same wc)


mkdir -p intervalLists/split_3rd
cd intervalLists/split_3rd

## set start counter
CNT=582

for ((i=${CNT}; i<=582; i++)); do

    ## padded id
    PID=$(printf "%04d" ${i})

    ## the old list
    scaffolds=$HOME/Style/Pax_genome/intervals/intervalLists/${PID}-scattered.interval_list

    ## get number of sub-scaffolds in old list
    nintvl=`grep -v "^@" ${scaffolds} | grep -c "^Peax"`

    ## compute number of new splits
    splits=$(( 1 + ${nintvl} / 50 ))


    ## make lists of intervals
    ## -----------------------
    GenomeAnalysisTK SplitIntervals \
		     -R $reference \
		     -L $scaffolds \
		     -O tmp \
		     --scatter-count $splits \
		     --subdivision-mode INTERVAL_COUNT


    ## move new lists to split_3rd, and old list to split_2nd
    outcnt=`ls tmp | wc -l`

    if [[ ${outcnt} -gt 0 ]]; then
	for ((j=0; j<${outcnt}; j++)); do
	    PCNT=$(printf "%04d" ${CNT})
	    ((CNT++))
	    PID=$(printf "%04d" ${j})
	    mv tmp/${PID}-scattered.interval_list ${PCNT}-scattered.interval_list
	done
	mv ${scaffolds} ../split_2nd/
    else
	echo "ERROR: something went wrong for ${scaffolds}"
    fi

done

## ---------------------------------------------
## !!! check whether happy with the output and then:
mv -n $HOME/Style/Pax_genome/intervals/intervalLists/split_3rd/*scattered.interval_list $HOME/Style/Pax_genome/intervals/intervalLists/
## ---------------------------------------------



## specify list of intervals:
## https://software.broadinstitute.org/gatk/documentation/article?id=11009
## -L intervals.list (or intervals.interval_list, or intervals.bed) when specifying a text file containing intervals
## GATK-style .list or .intervals: <chr>:<start>-<stop> (only <chr> is strictly required)

