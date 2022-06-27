#!/bin/bash
set -e

ACC=$1
BASE_DIR=/project/umw_zhiping_weng/andrewsg/benchmark
RESULTS_DIR=$BASE_DIR/2-Motifs/DREME/$ACC
PEAK_FILE=/project/umw_zhiping_weng/andrewsg/data/cistrome/human_factor/"$ACC"_sort_peaks.narrowPeak.bed
GENOME=/project/umw_zhiping_weng/andrewsg/data/genome/hg38/hg38.fa


module load perl/5.18.1
module load gcc/5.4.0
export PATH=/home/ga42w/meme/bin:/home/ga42w/meme/libexec/meme-5.3.0:$PATH

rm -rf $RESULTS_DIR; mkdir -p $RESULTS_DIR; cd $RESULTS_DIR

START_TIME=$(date +%s)
echo "start:" $START_TIME >> $ACC.log

dreme -p $BASE_DIR/1-Fasta/$ACC.pos.fasta -n $BASE_DIR/1-Fasta/$ACC.neg.fasta -m 16

STOP_TIME=$(date +%s)
echo "stop:" $STOP_TIME >> $ACC.log