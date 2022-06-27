#!/bin/bash
set -e

tf=$1

baseDir=/data/zusers/andrewsg/zoonomia
tmpDir=/tmp/andrewsg/$tf
# cd $tmpDir
rm -rf $tmpDir; mkdir -p $tmpDir; cd $tmpDir

cp -rf /data/zusers/andrewsg/zmotif/zmotif.py .
cp $baseDir/3-Seeds/$tf.meme .

cat $baseDir/2-Instances/$tf/$tf.final.bed | awk '$2 > 0' | \
    sort -k1,1 -k2,2n | uniq | tee $tf.pos.bed | \
    bedtools getfasta -fi /home/andrewsg/genome/hg38/hg38.fa -bed - > $tf.pos.fasta
    
singularity exec /home/andrewsg/bin/zmotif.simg python3 zmotif.py -pos_fasta $tf.pos.fasta -chrom_sizes /home/andrewsg/genome/hg38/hg38.chrom.sizes -g /home/andrewsg/genome/hg38/hg38.fa  -motif_db $tf.meme -o $tf -n 1 -w 24 -e 100 -b 32 --seed_motif $tf.meme -flank 5000

cp $tf.bed $baseDir/4-TFBSs/$tf.bed
cd $baseDir
rm -rf $tmpDir