#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=5-00:00:00
#SBATCH --job-name=AMG_Hunter_StG
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=150G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load diamond/2.0.7
module load python/3.8.5
module load hmmer
module load perl

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/AMG_Hunter.py --annotate True --threads 48 --cds /mnt/lustre/scratch/fcoutinho/StG/Min_Count_1k_Cluster_95ID_80Cov_GLUVAB3_VIBRANT_Genomes_Min_5kbp_rep_seq_CDS.faa