#!/bin/bash
set -e

tf=$1

cat /data/zusers/andrewsg/zoonomia/4-TFBSs/$tf.bed | sort -k1,1 -k2,2n |\
    awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,"intergenic"}' |\
    bedtools closest -a - -b misc/exon.bed -d -t first |\
    awk 'BEGIN{FS=OFS="\t"} {if ($NF == 0) {print $1,$2,$3,$4,$5,$6,"exonic"} else {print $1,$2,$3,$4,$5,$6,$7}}' |\
    bedtools closest -a - -b misc/TSS.bed -d -t first |\
    awk 'BEGIN{FS=OFS="\t"} {if ($NF == 0) {print $1,$2,$3,$4,$5,$6,"TSS",$NF} else {print $1,$2,$3,$4,$5,$6,$7,$NF}}' |\
    bedtools closest -a - -b <(cat /data/projects/encode/Registry/V3/GRCh38/GRCh38-rDHSs.bed | sort -k1,1 -k2,2n) -d -t first |\
    awk 'BEGIN{FS=OFS="\t"} {if ($NF == 0) {print $1,$2,$3,$4,$5,$6,$7,$8,"rDHS"} else {print $1,$2,$3,$4,$5,$6,$7,$8,"no-rDHS"}}' |\
    bedtools closest -a - -b transposons.bed -d -t first |\
    awk 'BEGIN{FS=OFS="\t"} {if ($NF == 0) {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$(NF-3),$(NF-2),$(NF-1),"junk"} else {print $1,$2,$3,$4,$5,$6,$7,$8,$9,"None","None",-1.0,"junk"}}' > 5-TFBSs-Annotated/$tf.bed
    
bedtools getfasta -bed 5-TFBSs-Annotated/$tf.bed -fi /home/andrewsg/genome/hg38/hg38.fa -s | paste - - | awk '{print toupper($2)}' > 5-TFBSs-Annotated/$tf.seqs.txt

paste 5-TFBSs-Annotated/$tf.bed 5-TFBSs-Annotated/$tf.seqs.txt |\
    awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$14}' > 5-TFBSs-Annotated/$tf.tmp.bed

mv 5-TFBSs-Annotated/$tf.tmp.bed 5-TFBSs-Annotated/$tf.bed

rm 5-TFBSs-Annotated/$tf.seqs.txt