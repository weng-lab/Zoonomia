#!/bin/bash
source /home/ga42w/miniconda3/etc/profile.d/conda.sh
set -e

ACC=$1
BASE_DIR=/project/umw_zhiping_weng/andrewsg/benchmark
RESULTS_DIR=$BASE_DIR/2-Motifs/ProSampler/$ACC

rm -rf $RESULTS_DIR; mkdir -p $RESULTS_DIR; cd $RESULTS_DIR

cp $BASE_DIR/1-Fasta/$ACC.pos.fasta pos.fasta
cp $BASE_DIR/1-Fasta/$ACC.neg.fasta neg.fasta

module load g++/5.4.0

START_TIME=$(date +%s)
echo "start:" $START_TIME >> $ACC.log

$HOME/ProSampler/ProSampler.unix -i pos.fasta -b neg.fasta -o $ACC

STOP_TIME=$(date +%s)
echo "stop:" $STOP_TIME >> $ACC.log

rm pos.fasta
rm neg.fasta