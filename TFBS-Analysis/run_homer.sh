#!/bin/bash
source /home/ga42w/miniconda3/etc/profile.d/conda.sh
set -e

ACC=$1
BASE_DIR=/project/umw_zhiping_weng/andrewsg/benchmark
RESULTS_DIR=$BASE_DIR/2-Motifs/HOMER/$ACC
PEAK_FILE=/project/umw_zhiping_weng/andrewsg/data/cistrome/human_factor/"$ACC"_sort_peaks.narrowPeak.bed
GENOME=/project/umw_zhiping_weng/andrewsg/data/genome/hg38/hg38.fa

rm -rf $RESULTS_DIR; mkdir -p $RESULTS_DIR; cd $RESULTS_DIR

export PATH=$PATH:$HOME/homer/bin/

START_TIME=$(date +%s)
echo "start:" $START_TIME >> $ACC.log

perl $HOME/homer/bin/findMotifs.pl $BASE_DIR/1-Fasta/$ACC.pos.fasta fasta motifResults/ -fasta $BASE_DIR/1-Fasta/$ACC.neg.fasta -len 6,12,18,24 -S 4 -p 4

STOP_TIME=$(date +%s)
echo "stop:" $STOP_TIME >> $ACC.log

touch $ACC.ok