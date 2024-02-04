#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-12:00:00
#SBATCH --job-name=Jellyfish_Malaspina_Viromes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=1G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load jellyfish/2.3.0
module load python


time python3 /mnt/lustre/bio/users/fcoutinho/Scripts/jellywrapper.py --metagenomes_extension fastq.gz --kmer 7 --threads 24 --metagenomes_dir /mnt/lustre/scratch/fcoutinho/BOV/QC_Reads/