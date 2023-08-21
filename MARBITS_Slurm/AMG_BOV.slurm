#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-03:00:00
#SBATCH --job-name=AMG_Hunter_BOV
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=15G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load diamond/2.0.7
module load python/3.8.5
module load hmmer

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/AMG_Hunter.py --cds /mnt/lustre/scratch/fcoutinho/BOV/Renamed_BOV_CheckV_Trimmed.faa --threads 48

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/AMG_Hunter.py --annotate True --parse_only True --genome_abundance /mnt/lustre/scratch/fcoutinho/BOV/Corrected_Abundance/Clean_Names_Correct_RPKM_Abundance_Renamed_BOV_CheckV_Trimmed.tsv --cds /mnt/lustre/scratch/fcoutinho/BOV/Renamed_BOV_CheckV_Trimmed.faa