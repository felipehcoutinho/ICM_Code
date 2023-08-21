#!/bin/bash
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-03:00:00
#SBATCH --job-name=CoCoNet_BOV
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=250G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
eval "$(conda shell.bash hook)"
conda activate  /home/fcoutinho/.conda/envs/coconet_dev

coconet run --fasta /mnt/lustre/scratch/fcoutinho/xlopez/RefSeq_Genomes_Phycodnaviridae.fasta --output /mnt/lustre/scratch/fcoutinho/xlopez/Test_CoCoNet/Results/ --threads 24