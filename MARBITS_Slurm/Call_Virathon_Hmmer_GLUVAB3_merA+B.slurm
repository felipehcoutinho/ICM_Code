#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-04:00:00
#SBATCH --job-name=Hmmer_GLUVAB3_merA+B_Virathon
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=1G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python/3.8.5
module load hmmer/3.3

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --call_hmmer True --cds /mnt/lustre/scratch/fcoutinho/GLUVAB3/GLUVAB3_CDS.faa --threads 8 --hmmer_program hmmsearch --hmmer_db /mnt/lustre/scratch/fcoutinho/Malaspina_Virus/merA+B.hmm
