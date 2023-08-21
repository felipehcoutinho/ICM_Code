#!/bin/bash
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=3-00:00:00
#SBATCH --job-name=CoCoNet_TOPC
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=250G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
eval "$(conda shell.bash hook)"
conda activate /mnt/lustre/bio/users/fcoutinho/envs/coconet

coconet run --fasta /mnt/lustre/scratch/fcoutinho/TARA_Polar/Renamed_Viruses_Genomes_Filtered_All_Tara_Polar_Coassemblies_Contigs.fasta --output /mnt/lustre/scratch/fcoutinho/TARA_Polar/TOPC_CoCoNet --threads 30 --bam /mnt/lustre/scratch/fcoutinho/TARA_Polar/Abundance/BAMs/*bam