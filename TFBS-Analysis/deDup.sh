#!/bin/bash
set -e

acc=$1
baseDir=/data/zusers/andrewsg/zoonomia/Roller-et-al
dataDir=/data/zusers/andrewsg/zoonomia/Roller-et-al/2-BAM
tmpDir=/tmp/andrewsg/$acc

rm -rf $tmpDir; mkdir -p $tmpDir; cd $tmpDir
# cd $tmpDir

echo "Copying BAM file"
cp $dataDir/$acc.bam .

echo "Marking duplicates"
java -jar /home/andrewsg/bin/picard.jar MarkDuplicates\
      I=$acc.bam \
      O=$acc.dupMarked.bam \
      M=marked_dup_metrics.txt


samtools view -F 1024 -@ 4 -o $acc.dedup.bam $acc.dupMarked.bam

echo "Indexing de-duplicated BAM file"
samtools index -@ 4 $acc.dedup.bam

target=$(cat $baseDir/metadata.txt | grep -w $acc | awk 'BEGIN{FS="\t"}{print $5; exit}')
echo $target

if [ "$target" = "anti-H3K4me1" ] 
then
    N=40000000
elif [ "$target" = "input DNA" ]
then
    echo "Input"
    N=40000000
else
    N=20000000
fi

echo $N
echo "Obtaining down-sampling fraction"
frac=$( samtools idxstats $acc.dedup.bam | cut -f3 | awk -v N=$N 'BEGIN {total=0} {total += $1} END {frac=N/total; if (frac > 1) {print 1} else {print frac}}' )

echo $frac
if [ "$frac" -eq 1 ] 
then
    echo "Nothing to do"
    cp $acc.dedup.bam $acc.dedup.downSamp.bam
else
    echo "Downsampling BAM"
    samtools view -bs $frac $acc.dedup.bam > $acc.dedup.downSamp.bam
fi

cp $acc.dedup.bam $dataDir
cp $acc.dedup.downSamp.bam $dataDir

cd 
rm -rf $tmpDir