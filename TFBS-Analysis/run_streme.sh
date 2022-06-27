#!/bin/bash
source /home/ga42w/miniconda3/etc/profile.d/conda.sh
set -e

ACC=$1
BASE_DIR=/project/umw_zhiping_weng/andrewsg/benchmark
RESULTS_DIR=$BASE_DIR/2-Motifs/STREME/$ACC
PEAK_FILE=/project/umw_zhiping_weng/andrewsg/data/cistrome/human_factor/"$ACC"_sort_peaks.narrowPeak.bed
GENOME=/project/umw_zhiping_weng/andrewsg/data/genome/hg38/hg38.fa

export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-5.3.0:$PATH

rm -rf $RESULTS_DIR; mkdir -p $RESULTS_DIR; cd $RESULTS_DIR

START_TIME=$(date +%s)
echo "start:" $START_TIME >> $ACC.log

streme -p $BASE_DIR/1-Fasta/$ACC.pos.fasta -n $BASE_DIR/1-Fasta/$ACC.neg.fasta --maxw 24 --nmotifs 16

STOP_TIME=$(date +%s)
echo "stop:" $STOP_TIME >> $ACC.log