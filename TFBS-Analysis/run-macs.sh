#!/bin/bash
set -e

baseDir=/data/zusers/andrewsg/zoonomia/Roller-et-al
resultsDir=$baseDir/4-bedGraph
bamDir=$baseDir/2-BAM

indv=$1
tissue=$2
mark=$3
tmpDir=/tmp/andrewsg/"$indv"_"$tissue"_"$mark"

rm -rf $tmpDir; mkdir -p $tmpDir; cd $tmpDir

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
    singularity exec /home/andrewsg/bin/chip-seq-pipeline2.simg samtools merge -@ 24 trt.bam $trtBams
fi
echo "Indexing treatment BAM"
singularity exec /home/andrewsg/bin/chip-seq-pipeline2.simg samtools index -@ 24 trt.bam

species=$(cat $baseDir/metadata.txt | grep -w "$indv" | awk -F '\t' '{ print $2 }' | head -n1)
g=$(cat $baseDir/genome-sizes.txt | grep -w "$species" | cut -f3 | head -n1)
echo $species "Genome size="$g" bases (" $genome ")"

singularity exec /home/andrewsg/bin/chip-seq-pipeline2.simg macs2 callpeak -t trt.bam -c input.bam -B -g $g -n "$indv"_"$tissue"_"$mark" -q 0.05 --broad --broad-cutoff 0.1

cp "$indv"_"$tissue"_"$mark"_peaks.broadPeak $baseDir/5-Peaks
# cp "$indv"_"$tissue"_"$mark"_peaks.narrowPeak $baseDir/5-Peaks

singularity exec /home/andrewsg/bin/chip-seq-pipeline2.simg macs2 bdgcmp -t "$indv"_"$tissue"_"$mark"_treat_pileup.bdg -c "$indv"_"$tissue"_"$mark"_control_lambda.bdg -o "$indv"_"$tissue"_"$mark"_FE.bdg -m FE

cat "$indv"_"$tissue"_"$mark"_FE.bdg | sort -k1,1 -k2,2n --parallel=24 -S32G > $resultsDir/"$indv"_"$tissue"_"$mark".sort.FE.bg
