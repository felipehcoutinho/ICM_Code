#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=01-00:00:00
#SBATCH --job-name=BOV_Hmmer_KOfam
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=50G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python hmmer

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --threads 48 --call_hmmer True --hmmer_min_score 30 --hmmer_max_evalue 0.00001 --cds /mnt/lustre/scratch/fcoutinho/BOV/Filtered_BOV_Viral.faa --hmmer_db /mnt/lustre/bio/users/fcoutinho/KOfam/All_KOs.hmm
