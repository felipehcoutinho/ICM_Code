#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=01-00:00:00
#SBATCH --job-name=Peatlands_VirSorter
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load virsorter/2.2.3
conda activate vs2

virsorter run --keep-original-seq --seqfile /mnt/lustre/scratch/fcoutinho/Peatlands/Peatlands_Viruses_Genomes.fasta --min-length 1000 --min-score 0.5 --jobs 24 --prep-for-dramv all
