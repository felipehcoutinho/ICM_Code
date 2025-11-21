#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=01-01:00:00
#SBATCH --job-name=MAtrix_gvSAGs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --output=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
conda activate /mnt/smart/scratch/vir/felipe/envs/basic

Rscript /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/SAG_Viruses_Interaction_Matrix_Analysis.R