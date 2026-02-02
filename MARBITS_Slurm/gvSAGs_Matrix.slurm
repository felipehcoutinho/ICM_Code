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
#SBATCH --mem=25G
#SBATCH --output=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
cd /mnt/smart/scratch/vir/felipe/gvSAGs/Read_Mapping_All_SAGs/

module load miniconda
source /home/apps/miniconda/4.5.x/etc/profile.d/conda.sh

conda activate /mnt/smart/scratch/vir/felipe/envs/basic

Rscript /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/SAG_Viruses_Interaction_Matrix_Analysis.R