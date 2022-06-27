#!/bin/bash
set -e

baseDir=/data/zusers/andrewsg/zoonomia/Roller-et-al
resultsDir=$baseDir/7-bigWig
bamDir=$baseDir/2-BAM

indv=$1
tissue=$2
mark=$3
tmpDir=/tmp/andrewsg/"$indv"_"$tissue"_"$mark"

rm -rf $tmpDir; mkdir -p $tmpDir; cd $tmpDir
# cd $tmpDir

echo "copying ChIP BAM files"
for acc in $(cat $baseDir/metadata.txt | grep -w "$indv" | grep "$tissue" | grep "$mark" | cut -f1)
do
    echo $acc
    cp $bamDir/$acc.dedup.bam $tmpDir
done

echo "copying input BAM files"
for acc in $(cat $baseDir/metadata.txt | grep -w "input" | grep -w "$indv" | grep "$tissue" | cut -f1)
do
    echo $acc
    cp $bamDir/$acc.dedup.bam $tmpDir
done

inputBams=$(cat $baseDir/metadata.txt | grep -w "input" | grep -w "$indv" | grep "$tissue" | cut -f1  | sed 's/.*/&.dedup.bam/' | xargs)
echo $inputBams
if [ "${#inputBams[@]}" -eq 1 ]
then
    echo "Only one input BAM provided"
    mv $inputBams input.bam
else
    singularity exec /home/andrewsg/bin/chip-seq-pipeline2.simg samtools merge -@ 64 input.bam $inputBams
fi

echo "Indexing input BAM"
singularity exec /home/andrewsg/bin/chip-seq-pipeline2.simg samtools index -@ 64 input.bam

trtBams=$(cat $baseDir/metadata.txt | grep -w "$indv" | grep "$tissue" | grep "$mark" | cut -f1 | sed 's/.*/&.dedup.bam/' | xargs)
echo $trtBams
if [ ${#trtBams[@]} -eq 1 ]
then
    echo "Only one treatment BAM provided"
    mv $trtBams trt.bam
else
    singularity exec /home/andrewsg/bin/chip-seq-pipeline2.simg samtools merge -@ 64 trt.bam $trtBams
fi

echo "Indexing treatment BAM"
singularity exec /home/andrewsg/bin/chip-seq-pipeline2.simg samtools index -@ 64 trt.bam

species=$(cat $baseDir/metadata.txt | grep -w "$indv" | awk -F '\t' '{ print $2 }' | head -n1)
g=$(cat $baseDir/genome-sizes.txt | grep -w "$species" | cut -f3 | head -n1)
echo $species "Genome size="$g" bases (" $genome ")"

singularity exec /home/andrewsg/bin/deepTools.simg bamCompare -b1 trt.bam -b2 input.bam -o "$indv"_"$tissue"_"$mark".bigWig --binSize 10 --exactScaling -p 64 --effectiveGenomeSize $g -e 147 --operation ratio -v

cp "$indv"_"$tissue"_"$mark".bigWig $resultsDir