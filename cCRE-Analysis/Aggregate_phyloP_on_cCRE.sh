#!/bin/bash

# -- Kaili Fan
# fankaili.bio@gmail.com
# Weng Lab, UMass Chan Med
# Jun, 2022

# This script is for aggregating phyloP scores on human cCREs.

cCRE_file=$1
phyloP_signal=$2
output_matrix=$3

scriptDir=$4

##################
n=500 # total number of bins

## cut cCRE-center±250bp into 1bp bins
awk -v n="$n" '{FS=OFS="\t"}{
    center=int(($2+$3)/2);
    start=center-250;
    for(i=1;i<=n;i++){
        print $1,int(start+(i-1)),int(start+i),$4"_"i;1;"."
    }
}' ${cCRE_file} > tmp.ccREs_split.bed

## calcualte signal for each bin
bigWigAverageOverBed ${phyloP_signal} tmp.ccREs_split.bed tmp.ccREs_split.tab
cut -f 1,5 tmp.ccREs_split.tab | awk '{FS=OFS="\t"}{
    split($1,a,"_");print $1,a[1],a[2],$2
}' | sort -k2,2 -k3,3n > tmp.ccREs_signal.txt

## make matrix
awk -v n="$n" '{FS=OFS="\t"}{if(NR%n==1){split($1,a,"_");printf a[1]"\t"$4"\t"}else if(NR%n!=0){printf $4"\t"}else{printf $4"\n"}}'  tmp.ccREs_signal.txt > ${output_matrix}


rm tmp.ccREs_split.bed tmp.ccREs_signal.txt
