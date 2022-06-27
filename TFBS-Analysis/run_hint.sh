#!/bin/bash
set -e

cell=$1
baseDir=/data/zusers/andrewsg/zoonomia/ATAC-seq
tmpDir=/tmp/andrewsg/$cell

rm -rf $tmpDir; mkdir -p $tmpDir; cd $tmpDir
# cd $tmpDir

echo "Copy BAM files..."
for acc in $(cat $baseDir/atac-metadata.txt | grep -w $cell | cut -f4)
do
    echo $acc
    cp $baseDir/1-BAM/$acc.bam .
done

echo "Merging BAM files..."
singularity exec $HOME/bin/atac-seq-pipeline.simg samtools merge -@ 64 $cell.merged.bam *.bam

acc=$(cat $baseDir/atac-metadata.txt | grep -w $cell | head -n1 | cut -f5)
cp $baseDir/2-Peaks/$acc.bed peaks.bed
echo "Indexing BAM file"
singularity exec $HOME/bin/atac-seq-pipeline.simg samtools index -@ 64 $cell.merged.bam

echo "Running HINT ATAC bias correction"

singularity exec /home/andrewsg/bin/rgt.simg rgt-hint tracks --bc --bigWig --organism=hg38 $cell.merged.bam peaks.bed --output-prefix=$cell
cp "$cell".bw $baseDir/3-BigWig/$cell.bc.bigWig

singularity exec /home/andrewsg/bin/rgt.simg rgt-hint tracks --raw --bigWig --organism=hg38 $cell.merged.bam peaks.bed --output-prefix=$cell
cp "$cell".bw $baseDir/3-BigWig/$cell.raw.bigWig

