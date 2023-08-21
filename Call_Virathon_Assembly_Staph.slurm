#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=01-12:00:00
#SBATCH --job-name=Assembly_Staph
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=250G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python/3.8.5
module load spades/3.13.0
module load vibrant/1.2.1

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --threads 48 --assemble True --samples_dir /mnt/lustre/scratch/fcoutinho/Staph_Phages/Clean_Reads/ --samples_extension .fastq.gz --call_vibrant True
