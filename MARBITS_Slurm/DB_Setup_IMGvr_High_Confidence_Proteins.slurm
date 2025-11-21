#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-12:00:00
#SBATCH --job-name=SetupDB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=75G
#SBATCH --output=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
cd /mnt/smart/scratch/vir/felipe/Databases/IMGVR/

gunzip IMGVR_all_proteins-high_confidence.faa.gz

module load mmseqs2

mmseqs createdb IMGVR_all_proteins-high_confidence.faa DBIMGVR_all_proteins-high_confidence
mmseqs createindex DBIMGVR_all_proteins-high_confidence tmp