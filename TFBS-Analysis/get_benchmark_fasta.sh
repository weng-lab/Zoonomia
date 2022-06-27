#!/bin/bash
source /home/ga42w/miniconda3/etc/profile.d/conda.sh
set -e

ACC=$1
BASE_DIR=/project/umw_zhiping_weng/andrewsg/benchmark
RESULTS_DIR=$BASE_DIR/1-Fasta
PEAK_FILE=/project/umw_zhiping_weng/andrewsg/data/cistrome/human_factor/"$ACC"_sort_peaks.narrowPeak.bed
GENOME=/project/umw_zhiping_weng/andrewsg/data/genome/hg38/hg38.fa

export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-5.3.0:$PATH
module load bedtools/2.28.0

cd $RESULTS_DIR

cat $PEAK_FILE | sort -nrk7 | grep -w '^#\|chr[1-9]\|chr[1-2][0-9]\|chr[X]' | awk 'BEGIN{FS="\t";OFS="\t"} {print $1,$2+$10-100,$2+$10+100}' | awk '$2 > 0' | head -n100000 | sort -k1,1 -k2,2n | uniq | bedtools getfasta -fi $GENOME -bed - > $ACC.pos.fasta 

fasta-shuffle-letters -kmer 2 $ACC.pos.fasta > $ACC.neg.fasta
