#!/bin/bash

#SBATCH --mail-user=me@work
#SBATCH --mail-type=end,fail
#SBATCH --job-name="prep_oma"
#SBATCH --output=%x-%A-%a.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:50:00
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=all



module load vital-it/7

## run in [A] '$HOME/Style/OMA/Pax304-Pex304/prep'
## OR in [B] '$HOME/Style/OMA/Pax304-Pax162/prep'

dir="$HOME/Style/OMA/Pax304-Pex304/prep" ## [A]
##dir="$HOME/Style/OMA/Pax304-Pax162/prep" ## [B]
cd ${dir}


## [A]

## PEAX 304
## --------

fasta="$HOME/Style/Pax_genome/Peax304All_final.maker.proteins.fasta_renamed_ids_edited.fa"

species="PAX"

db="$HOME/Style/OMA/Pax304-Pex304/DB"

$HOME/Style/scripts/fastaAdjust.pl ${fasta} ${species} " " 1

$HOME/Style/scripts/fastaFilter.pl ${species}.fasta 10 0.2

$HOME/Style/scripts/fastaInfo.pl ${species}_good.fasta 2

cp ${species}_good.fasta ${db}/${species}.fa
sed -i "s/${species}|//g" ${db}/${species}.fa


## PEEX 304
## --------

fasta="$HOME/Style/Pex_genome/Peex304All_final.maker.proteins.fasta_renamed_ids_edited.fa"

species="PEX"

db="$HOME/Style/OMA/Pax304-Pex304/DB"

$HOME/Style/scripts/fastaAdjust.pl ${fasta} ${species} " " 1

$HOME/Style/scripts/fastaFilter.pl ${species}.fasta 10 0.2

$HOME/Style/scripts/fastaInfo.pl ${species}_good.fasta 2

cp ${species}_good.fasta ${db}/${species}.fa
sed -i "s/${species}|//g" ${db}/${species}.fa


# ## [B]

# ## PEAX 162
# ## --------

# fasta="$HOME/Style/OMA/Pax304-Pax162/prep/Peaxi162annotation_v4.PROTEIN.fasta"

# species="PAX162"

# db="$HOME/Style/OMA/Pax304-Pax162/DB"

# ## for some weird reason, most protein sequences have a "." at the end of line

# cat ${fasta} | sed 's/.$//g' > ${fasta}.edited

# $HOME/Style/scripts/fastaAdjust.pl ${fasta}.edited ${species} " " 1

# $HOME/Style/scripts/fastaFilter.pl ${species}.fasta 10 0.2

# $HOME/Style/scripts/fastaInfo.pl ${species}_good.fasta 2

# cp ${species}_good.fasta ${db}/${species}.fa
# sed -i "s/${species}|//g" ${db}/${species}.fa

# ## assuming PAX304 was already done
# scp -p $HOME/Style/OMA/Pax304-Pex304/DB/PAX.fa $HOME/Style/OMA/Pax304-Pax162/DB/PAX304.fa


