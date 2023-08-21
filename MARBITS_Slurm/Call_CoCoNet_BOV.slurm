#!/bin/bash
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=3-00:00:00
#SBATCH --job-name=CoCoNet_BOV
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=100G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load coconet/1.1.0-064cb8c
source activate
export PYTHONUSERBASE=intentionally-disabled
coconet run --fasta /mnt/lustre/scratch/fcoutinho/BOV/Renamed_BOV_CheckV_Trimmed.fasta --output /mnt/lustre/scratch/fcoutinho/BOV/CoCoNet_BOV --threads 47 --bam /mnt/lustre/scratch/fcoutinho/BOV/Corrected_Abundance/*bam