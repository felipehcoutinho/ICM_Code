#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-03:00:00
#SBATCH --job-name=Get_GLUVAB3_mer_Hits
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500M
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python
python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Get_Seqs.py --input_seq /mnt/lustre/scratch/fcoutinho/GLUVAB3/GLUVAB3_CDS.faa --list /mnt/lustre/scratch/fcoutinho/Malaspina_Virus/GLUVAB3_mer_Candidates.txt --output_seq /mnt/lustre/scratch/fcoutinho/Malaspina_Virus/GLUVAB3_mer_Candidates.faa