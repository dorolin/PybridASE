#!/bin/bash

#SBATCH --mail-user=me@work
#SBATCH --mail-type=end,fail
#SBATCH --job-name="readCorresp"
#SBATCH --output=%x-%A-%a.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5:50:00
#SBATCH --mem-per-cpu=30G
#SBATCH --array=34,35,36,37,38,39,31,32,33
#SBATCH --partition=all


## get correct array task ids for samples of interest
## grep -n "AxExMS1" $samples ## 34,35,36
## grep -n "AxExMS2" $samples ## 37,38,39
## grep -n "AxExLS" $samples ## 31,32,33

module load vital-it/7
module load UHTS/Analysis/samtools/1.8

## run in '$HOME/Style/step02_AxEx/corresp'


samples=$HOME/Style/step02_AxEx/files_AxEx.txt
filebase=$(cat $samples | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')
FILE1=$HOME/Style/step02_AxEx/Pax_genome/secondPass/${filebase}/Aligned.sortedByCoord.out.bam
FILE2=$HOME/Style/step02_AxEx/Pex_genome/secondPass/${filebase}/Aligned.sortedByCoord.out.bam

hostname; date
start_time="$(date -u +%s)"
echo "processing $filebase"

## Pax -> Pex
## ----------
$HOME/Style/scripts/readCorresp.pl \
  -f $FILE1 -s $FILE2 -o ${filebase}_aligninfo.txt

## sort by scaffold name in dictionary order, then
##  by scaffold name in numeric order, then by start
##  position in numeric order (only first entry for
##  multimappers considered for sorting)
paste -d';' <(cat ${filebase}_aligninfo.txt | cut -f2 | cut -d' ' -f2 | cut -d',' -f1 | sed 's/[0-9]\+//g' | sed 's/\./Z/g') <(cat ${filebase}_aligninfo.txt | cut -f2 | cut -d' ' -f2 | cut -d',' -f1 | sed 's/[a-zA-Z]\+//g') <(cat ${filebase}_aligninfo.txt | cut -f2 | cut -d' ' -f3 | cut -d',' -f1) <(cat ${filebase}_aligninfo.txt) >${filebase}_aligninfo_tmp.txt

cat ${filebase}_aligninfo_tmp.txt | sort -k1d,1 -k2n,2 -k3n,3 -t';' | cut -d';' -f4- >${filebase}_aligninfo_Pax.txt

rm -f ${filebase}_aligninfo_tmp.txt
rm -f ${filebase}_aligninfo.txt


## Pex -> Pax
## ----------
$HOME/Style/scripts/readCorresp.pl \
  -f $FILE2 -s $FILE1 -o ${filebase}_aligninfo.txt

## sort by scaffold name in dictionary order, then
##  by scaffold name in numeric order, then by start
##  position in numeric order (only first entry for
##  multimappers considered for sorting)
paste -d';' <(cat ${filebase}_aligninfo.txt | cut -f2 | cut -d' ' -f2 | cut -d',' -f1 | sed 's/[0-9]\+//g' | sed 's/\./Z/g') <(cat ${filebase}_aligninfo.txt | cut -f2 | cut -d' ' -f2 | cut -d',' -f1 | sed 's/[a-zA-Z]\+//g') <(cat ${filebase}_aligninfo.txt | cut -f2 | cut -d' ' -f3 | cut -d',' -f1) <(cat ${filebase}_aligninfo.txt) >${filebase}_aligninfo_tmp.txt

cat ${filebase}_aligninfo_tmp.txt | sort -k1d,1 -k2n,2 -k3n,3 -t';' | cut -d';' -f4- >${filebase}_aligninfo_Pex.txt

rm -f ${filebase}_aligninfo_tmp.txt
rm -f ${filebase}_aligninfo.txt





end_time="$(date -u +%s)"
elaps=$(dc -e "$end_time $start_time - p")
elapstime=`printf '%dh:%dm:%ds\n' $(($elaps/3600)) $(($elaps%3600/60)) $(($elaps%60))`
echo "---------------------------"
echo "readCorresp took $elapstime"
echo "---------------------------"

