#!/bin/bash

# -- Kaili Fan
# fankaili.bio@gmail.com
# Weng Lab, UMass Chan Med
# Jun, 2022

# This script is for annotating gene by related phyloP scores.

workDir=$1
scriptDir=$2
cd ${workDir}

non_ubi_PLS_file=$3
# list of non-ubi-PLSs
# File from: Fan, K., Moore, J. E., Zhang, X. O., & Weng, Z. (2021). 
# Genetic and epigenetic features of promoters with ubiquitous chromatin accessibility support ubiquitous transcription of cell-essential genes.
# Nucleic acids research, 49(10), 5705-5725."

phyloP_signal=$4
GENCODE_TSS=$5

##################
# 1. defined conserved non-ubi-PLSs
bigWigAverageOverBed ${phyloP_signal} ${non_ubi_PLS_file} non-ubi-PLS_phyloP.tab
awk '{FS=OFS="\t"}{if(NR==FNR){if($5>1){a[$4]=1}}else{if(a[$4] ){print $0}}}' non-ubi-PLS_phyloP.tab ${non_ubi_PLS_file} > conserved_non-ubi-PLS.txt


# 2. get associated genes
bedtools intersect -a ${GENCODE_TSS} -b conserved_non-ubi-PLS.txt -wa | cut -f 7 | sort -u > gene_conserved_non-ubi-PLSs.txt


# 3. get control groups
bedtools intersect -a ${GENCODE_TSS} -b ${non_ubi_PLS_file} -wa | cut -f 7 | sort -u > gene_all_non-ubi-PLSs.txt
