#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=50G
#SBATCH --time=30:00
#SBATCH -o ./log/calculate_alignment_on_cCRE_SLURM.out
#SBATCH -e ./log/calculate_alignment_on_cCRE_SLURM.err
#SBATCH --partition=30mins
#SBATCH --job-name=calculate_alignment_on_cCRE_SLURM
#SBATCH --array=1-240%40

# -- Kaili Fan
# fankaili.bio@gmail.com
# Weng Lab, UMass Chan Med
# Jun, 2022

species_list=$1
species=`awk -v n="${SLURM_ARRAY_TASK_ID}" '{if(NR==n){print $1}}' ${species_list}`


###############
if [ ! -d /tmp/kaili/ ];then mkdir -p /tmp/kaili/; fi
cp ${workDir}signal/${species}_aligned.bw /tmp/kaili/

bigWigAverageOverBed /tmp/kaili/${species}_aligned.bw GRCh38-ccREs_rmPAR.bed /tmp/kaili/${species}_alignment.tab
#
echo ${species} > /tmp/kaili/${species}_alignment.txt
awk '{print $5}' /tmp/kaili/${species}_alignment.tab >> /tmp/kaili/${species}_alignment.txt

mv /tmp/kaili/${species}_alignment.t* ./alignment/
rm /tmp/kaili/${species}_aligned.bw
