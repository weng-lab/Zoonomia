#!/bin/bash

# -- Kaili Fan
# fankaili.bio@gmail.com
# Weng Lab, UMass Chan Med
# Jun, 2022

# This script is for aggregating phyloP scores on human cCREs.

workDir=$1
scriptDir=$2
cd ${workDir}

chain_file_dir=$3
chain_file_list=$4
# chain file downloaded from
#   https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/
#   https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz100way/
cCRE_file=$5

##################
# 1. liftOver for calculating % alignment between human and other 99 species
mkdir phyloP100_alignment
sbatch ${scriptDir}Calculate_phyloP100_alignment_ratio_SLURM.sh ${workDir} ${chain_file_dir} ${chain_file_list} ${cCRE_file}


# 2. count numbers
awk '{FS=OFS="\t"}{if(NR==1){print "id","ccre","align90","align10";print $4,$5,0,0}else{print $4,$5,0,0}}' ${cCRE_file} > GRCh38-ccREs_phyloP100_alignment_count.txt
while read species0
do
    species="${species0^}"
    if [ "$species" != "Hg38" ];then
        echo ${species}
        #
        awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$6}else{if(FNR==1){print $0}else{if(a[$1]==0.9){print $1,$2,$3+1,$4}else if(a[$1]==0){print $1,$2,$3,$4+1}else{print $0}}}}' ./phyloP100_alignment/GRCh38-ccREs_liftOverTo${species}.txt GRCh38-ccREs_phyloP100_alignment_count.txt > tmp.txt
        mv tmp.txt GRCh38-ccREs_phyloP100_alignment_count.txt
    fi
done < ${chain_file_list}
