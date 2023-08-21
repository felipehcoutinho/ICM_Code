#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-04:00:00
#SBATCH --job-name=Merge_VIBRANT_GLUVAB3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=END
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python
python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Index_Multi_VIBRANT_Results.py --input /mnt/lustre/scratch/fcoutinho/GLUVAB3/vibrant-results-tsv/ --out /mnt/lustre/scratch/fcoutinho/GLUVAB3/GLUVAB3_VIBRANT_Info.tsv