#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-06:00:00
#SBATCH --job-name=Build_MMSeqs_DB_GLUVAB3_CDS
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=400G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/GLUVAB3/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/GLUVAB3/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load mmseqs2/15c7762
#mmseqs createdb /mnt/lustre/scratch/fcoutinho/GLUVAB3/GLUVAB3_CDS.faa DB_GLUVAB3_CDS
mmseqs createindex DB_GLUVAB3_CDS tmp --split 0
