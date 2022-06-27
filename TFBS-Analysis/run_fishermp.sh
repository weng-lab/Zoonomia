#!/bin/bash
source /home/ga42w/miniconda3/etc/profile.d/conda.sh
set -e

ACC=$1
BASE_DIR=/project/umw_zhiping_weng/andrewsg/benchmark
RESULTS_DIR=$BASE_DIR/2-Motifs/FisherMP/$ACC

rm -rf $RESULTS_DIR; mkdir -p $RESULTS_DIR; cd $RESULTS_DIR
cp $BASE_DIR/1-Fasta/$ACC.pos.fasta pos.fasta
cp $BASE_DIR/1-Fasta/$ACC.neg.fasta neg.fasta

module load g++/5.4.0
export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-5.3.0:$PATH

START_TIME=$(date +%s)
echo "start:" $START_TIME >> $ACC.log

$HOME/fishermp/fishermp pos.fasta -b neg.fasta -t 4 > fishermp.out

STOP_TIME=$(date +%s)
echo "stop:" $STOP_TIME >> $ACC.log

conda activate numpy 
python $BASE_DIR/scripts/fishermp-to-meme.py fishermp.out $ACC.meme
conda deactivate 

tomtom -thresh 1 -text $ACC.meme $BASE_DIR/misc/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme > $ACC.tomtom.out

conda activate numpy
python $BASE_DIR/scripts/get_json.py $ACC.meme $ACC.tomtom.out FisherMP $ACC
conda deactivate 

rm pos.fasta
rm neg.fasta