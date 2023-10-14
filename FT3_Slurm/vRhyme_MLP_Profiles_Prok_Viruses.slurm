#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO

#SBATCH --time=00-06:00:00
#SBATCH --job-name=vRhyme_MLP_Profiles_Prok_Viruses
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=50G
#SBATCH --output=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.out
#SBATCH --error=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load cesga/2020 bowtie2 mmseqs2 samtools miniconda3/22.11.1-1

conda activate /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Conda_Envs/basic

cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/test_vRhyme

vRhyme 