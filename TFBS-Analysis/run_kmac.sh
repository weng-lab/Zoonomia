#!/bin/bash
source /home/ga42w/miniconda3/etc/profile.d/conda.sh
set -e

ACC=$1
BASE_DIR=/project/umw_zhiping_weng/andrewsg/benchmark
RESULTS_DIR=$BASE_DIR/2-Motifs/KMAC/$ACC

rm -rf $RESULTS_DIR
mkdir -p $RESULTS_DIR
cd $RESULTS_DIR

START_TIME=$(date +%s)
echo "start:" $START_TIME >> $ACC.log

java -Xmx16G -jar $HOME/gem/gem.jar KMAC --pos_seq $BASE_DIR/1-Fasta/$ACC.pos.fasta --k_win 61 --k_min 5 --k_max 13 --k_top 10 --out_name $ACC --neg_seq $BASE_DIR/1-Fasta/$ACC.neg.fasta
touch $ACC.ok

STOP_TIME=$(date +%s)
echo "stop:" $STOP_TIME >> $ACC.log
