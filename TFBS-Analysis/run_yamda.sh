#!/bin/bash
source /home/ga42w/miniconda3/etc/profile.d/conda.sh
set -e

ACC=$1
BASE_DIR=/project/umw_zhiping_weng/andrewsg/benchmark
RESULTS_DIR=$BASE_DIR/2-Motifs/YAMDA/$ACC
GENOME=/project/umw_zhiping_weng/andrewsg/data/genome/hg38/hg38.fa

rm -rf $RESULTS_DIR; mkdir -p $RESULTS_DIR; cd $RESULTS_DIR

conda activate YAMDA-env
START_TIME=$(date +%s)
echo "start:" $START_TIME >> $ACC.log

python $HOME/YAMDA/run_em.py -r -e -i $BASE_DIR/1-Fasta/$ACC.pos.fasta -j $BASE_DIR/1-Fasta/$ACC.neg.fasta -oc yamda_out -n 16 

STOP_TIME=$(date +%s)
echo "stop:" $STOP_TIME >> $ACC.log

conda deactivate 

rm yamda_out/*fa