#!/bin/bash
set -e 

export PATH=/data/common/tools/ucsc.v385/:$PATH

baseDir=/data/zusers/andrewsg/zoonomia/Roller-et-al

tf=$1
genome=$2
cat ../constrained-sites/$tf.bed | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,1000,"+"}' > $tf.tmp.bed

liftOver $tf.tmp.bed /home/andrewsg/genome/$genome/hg38To"$genome".over.chain.gz $baseDir/8-TFs/$tf.$genome.bed /dev/null
