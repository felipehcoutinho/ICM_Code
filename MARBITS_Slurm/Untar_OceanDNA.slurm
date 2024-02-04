#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=01-00:00:00
#SBATCH --job-name=OceanDNA_Untar
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=5G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
cd /mnt/lustre/bio/users/fcoutinho/OceanDNA_MAGs/
tar -xvf /mnt/lustre/bio/users/fcoutinho/OceanDNA_MAGs/fasta_non-representatives.tar


