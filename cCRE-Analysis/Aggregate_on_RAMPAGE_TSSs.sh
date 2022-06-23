#!/bin/bash

# -- Kaili Fan
# fankaili.bio@gmail.com
# Weng Lab, UMass Chan Med
# Jun, 2022

# This script is for aggregating phyloP scores on RAMPAGE rPeaks.

RAMPAGE_rPeak=$1
# File from: Moore, J. E., Zhang, X. O., Elhajjajy, S. I., Fan, K., Pratt, H. E., Reese, F., ... & Weng, Z. (2022).
# Integration of high-resolution promoter profiling assays reveals novel, cell type–specific transcription start sites across 115 human cell and tissue types.
# Genome Research, 32(2), 389-402."
# rPeaks from supplementary files, only use TSS peaks

signal=$2 #phyloP or π
output_matrix=$3

##################
# 1. separate peaks to narrow & broad
awk '{FS=OFS="\t"}{if(($3-$2)<=9 && $6=="+"){print $0}}' ${RAMPAGE_rPeak} > narrow_rPeak_pl.bed
awk '{FS=OFS="\t"}{if(($3-$2)>9 && $6=="+"){print $0}}' ${RAMPAGE_rPeak} > broad_rPeak_pl.bed


# 2. calculate phyloP signal
n=500
# total number of bins

for peak_type in narrow broad
do
    ## cut rPeak-summit±250bp into 1bp bins
    awk -v n="$n" '{FS=OFS="\t"}{
        split($4,a,"_");
        summit=a[3];
        start=summit-250;
        for(i=1;i<=n;i++){
            print $1,int(start+(i-1)),int(start+i),peak"-"NR"_"i;1;"."
        }
    }' ${peak_type}_rPeak_pl.bed > tmp.${peak_type}_rPeak_split.bed

    ## calcualte signal for each bin
    bigWigAverageOverBed ${signal} tmp.${peak_type}_rPeak_split.bed tmp.${peak_type}_rPeak_split.tab
    cut -f 1,5 tmp.ccREs_split.tab | awk '{FS=OFS="\t"}{
        split($1,a,"_");print $1,a[1],a[2],$2
    }' | sort -k2,2 -k3,3n > tmp.${peak_type}_rPeak_signal.txt

    ## make matrix
    awk -v n="$n" '{FS=OFS="\t"}{if(NR%n==1){split($1,a,"_");printf a[1]"\t"$4"\t"}else if(NR%n!=0){printf $4"\t"}else{printf $4"\n"}}'  tmp.${peak_type}_rPeak_signal.txt > ${output_matrix}

    rm tmp.${peak_type}_rPeak_split.bed tmp.${peak_type}_rPeak_signal.txt
done
