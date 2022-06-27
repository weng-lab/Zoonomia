#!/bin/bash
set -e
source /home/andrewsg/miniconda3/etc/profile.d/conda.sh

baseDir=/data/zusers/andrewsg/zoonomia/Roller-et-al
dataDir=$baseDir/1-FASTQ
resultsDir=$baseDir/2-BAM
qcDir=$baseDir/3-QC

acc=$1
tmpDir=/tmp/andrewsg/$acc

rm -rf $tmpDir; mkdir -p $tmpDir
cp $dataDir/$acc.fastq.gz $tmpDir
cd $tmpDir

genome=$(cat $baseDir/metadata.txt | grep -w $acc | awk 'BEGIN{FS="\t"}{print $7; exit}')


echo $genome


singularity exec /home/andrewsg/bin/chip-seq-pipeline2.simg bwa mem -t 24 /home/andrewsg/genome/$genome/$genome.fa $acc.fastq.gz > $acc.sam

singularity exec /home/andrewsg/bin/chip-seq-pipeline2.simg samtools view -@ 24 -bS $acc.sam | samtools sort -@ 24 -o $acc.bam -
singularity exec /home/andrewsg/bin/chip-seq-pipeline2.simg samtools flagstat -@ 24 $acc.bam > $acc.flagstat

mv $acc.bam $resultsDir
mv $acc.flagstat $qcDir

cd $baseDir; rm -rf $tmpDir

