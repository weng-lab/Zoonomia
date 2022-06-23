#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=50G
#SBATCH --time=12:00:00
#SBATCH -o ./log/calculate_phyloP100_alignment_ratio_SLURM.out
#SBATCH -e ./log/calculate_phyloP100_alignment_ratio_SLURM.err
#SBATCH --partition=12hours
#SBATCH --job-name=calculate_phyloP100_alignment_ratio_SLURM
#SBATCH --array=2-100

workDir=$1
chain_file_dir=$2
chain_file_list=$3
cCRE_file=$4

cd ${workDir}

species0=`awk -v n="${SLURM_ARRAY_TASK_ID}" '{if(NR==n){print $1}}' chain_file_list`
species="${species0^}"

###############
awk '{FS=OFS="\t"}{print $0,0}' ${cCRE_file} > ./phyloP100_alignment/GRCh38-ccREs_liftOverTo${species}.txt
cp ${cCRE_file} ${species}_current_file.txt

#for i in $(seq 1 -0.01 0.01)
for i in 0.9 0.1
do
    echo ${i}
    #
    liftOver -minMatch=${i} ${species}_current_file.txt ${chainDir}hg38To${species}.over.chain tmp.${species}_liftOver.txt tmp.${species}_unMapped.txt
    #
    awk -v i="$i" '{FS=OFS="\t"}{if(NR==FNR){a[$4]=1}else{if(a[$4]){print $1,$2,$3,$4,$5,i}else{print $0}}}' tmp.${species}_liftOver.txt ./phyloP100_alignment/GRCh38-ccREs_liftOverTo${species}.txt > tmp.${species}.txt
    mv tmp.${species}.txt ./phyloP100_alignment/GRCh38-ccREs_liftOverTo${species}.txt
    #
    awk '{FS=OFS="\t"}{if(NR%2!=1){print $0}}' tmp.${species}_unMapped.txt > ${species}_current_file.txt
done


rm tmp.${species}_liftOver.txt
rm tmp.${species}_unMapped.txt
rm ${species}_current_file.txt
