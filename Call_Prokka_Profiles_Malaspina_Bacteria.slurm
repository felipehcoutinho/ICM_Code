#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=02-00:00:00
#SBATCH --job-name=Prokka_Profiles_Malaspina_Bacteria
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=15G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load prokka
conda activate prokka-v1.14.6

for file in *.fasta; do prokka --kingdom Bacteria --cpus 24 --noanno --outdir "${file%.fasta}" --prefix "${file%.fasta}" --locustag "${file%.fasta}" "$file" ; done ;
