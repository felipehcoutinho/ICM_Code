#!/bin/bash
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=2-00:00:00
#SBATCH --job-name=Multiphate_RefSeqVir
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=100G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
eval "$(conda shell.bash hook)"
conda activate /mnt/lustre/bio/users/fcoutinho/envs/multiphate2
python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Run_Multiphate2.py --genomes /mnt/lustre/bio/users/fcoutinho/Databases/RefSeqVir_15-07-21/Seqs/RefSeqVir_Arch+Bact_15-07-21.fasta --threads 48