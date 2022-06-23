#!/bin/bash

# -- Kaili Fan
# fankaili.bio@gmail.com
# Weng Lab, UMass Chan Med
# Jun, 2022

# This script is for calculating %alignment for human cCREs versus 240 species.

workDir=$1
scriptDir=$2

cCRE_file=$3
PAR_region=$4
species_list=$5
human_genome=$6

cd ${workDir}

##################
# 1. remove cCREs overlapping PAR regions
intersectBed -a ${cCRE_file} -b ${PAR_region} -v  | sort -k1,1 -k2,2n > GRCh38-ccREs_rmPAR.bed


# 2. calculate alignment across all species
mkdir alignment
sbatch ${scriptDir}Calculate_alignment_on_cCRE_SLURM.sh ${species_list}
# generate matrix
paste ./alignment/*_alignment.txt > tmp.alignment.txt
echo "ccre" > tmp.txt
cut -f 4 GRCh38-ccREs_rmPAR.bed >> tmp.txt
paste tmp.txt tmp.alignment.txt > GRCh38-ccREs_alignment_matrix.txt
## count ≥90% alignment vs. ≤10% alignment
awk '{FS=OFS="\t"}{if(NR==1){print "id","align90","align10";print $1,0,0}else{print $1,0,0}}' ./alignment/Equus_asinus_alignment.tab > GRCh38-cCRE_alignment_count.txt
while read species
do
    echo ${species}
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$5}else{if(FNR==1){print $0}else{if(a[$1]>=0.9){print $1,$2+1,$3}else if(a[$1]<=0.1){print $1,$2,$3+1}else{print $0}}}}' ./alignment/${species}_alignment.tab GRCh38-cCRE_alignment_count.txt > tmp.txt
    mv tmp.txt GRCh38-cCRE_alignment_count.txt
done < ${species_list}


# 3. alignment for random regions
## generate 10k 250bp random regions
bedtools random -n 10000 -l 250 -g ${human_genome} > hg38_random_10k_250bp.bed
## calculate alignment
mkdir alignment_random
sbatch ${scriptDir}Calculate_alignment_on_randon-regions_SLURM.sh ${species_list}
paste ./alignment_random/*_alignment.txt > tmp.alignment.txt
echo "id" > tmp.txt
cut -f 4 hg38_random_10k_250bp.bed >> tmp.txt
paste tmp.txt tmp.alignment.txt > hg38_random_10k_250bp_alignment_matrix.txt
## count ≥90% alignment vs. ≤10% alignment
awk '{FS=OFS="\t"}{if(NR==1){print "id","align90","align10";print $1,0,0}else{print $1,0,0}}' ./alignment_random/Equus_asinus_alignment.tab > GRCh38-random_alignment_count.txt
while read species
do
    echo ${species}
    #
    awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$5}else{if(FNR==1){print $0}else{if(a[$1]>=0.9){print $1,$2+1,$3}else if(a[$1]<=0.1){print $1,$2,$3+1}else{print $0}}}}' ./alignment_random/${species}_alignment.tab GRCh38-random_alignment_count.txt > tmp.txt
    mv tmp.txt GRCh38-random_alignment_count.txt
done < ${species_list}



#------------------
Rscript ${scriptDir}UMAP_across_speceis.R GRCh38-ccREs_alignment_matrix.txt
Rscrip ${scriptDir}make_cCRE_heatscatter_240Species.R GRCh38-cCRE_alignment_count.txt GRCh38-random_alignment_count.txt
