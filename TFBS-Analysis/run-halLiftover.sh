#!/bin/bash
set -e

baseDir=/data/zusers/andrewsg/zoonomia/Roller-et-al
tf=$1
species=$2
tmpDir=/tmp/andrewsg/"$tf"_"$species"

rm -rf $tmpDir; mkdir -p $tmpDir; cd $tmpDir
# cd $tmpDir

cat /data/zusers/andrewsg/zoonomia/7-Constrained-TFBSs-Core-Only/$tf.bed | cut -f1-6 | grep -v "chrM" > $tf.bed

split -l 5000 --additional-suffix=.bed $tf.bed tmp.

for x in $(ls tmp.*.bed)
do 
    echo "singularity exec --bind /data/zusers/andrewsg/zoonomia/alignments/:/mnt /home/andrewsg/bin/comparative-annotation-toolkit.simg halLiftover --bedType 6 /mnt/241-mammalian-2020v2.hal Homo_sapiens" $x  $species $x.$species >> cmds.sh
done

parallel -j 64 < cmds.sh
cat tmp.*.$species | sort -k1,1 -k2,2n > $baseDir/8-TFs/$tf.$species.bed
