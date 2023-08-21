#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-06:00:00
#SBATCH --job-name=RaFAH_Malaspina_Virathon
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=80G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python/3.8.5
module load hmmer/3.3
module load perl
module load R

export PATH="/mnt/lustre/bio/users/fcoutinho/Scripts/:$PATH"
python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --call_rafah True --cds /mnt/lustre/scratch/fcoutinho/Malaspina_Virus/Malaspina_CoAssembly_Viral_Scaffolds_Renamed.faa --threads 28
