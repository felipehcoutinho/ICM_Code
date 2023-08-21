#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-08:00:00
#SBATCH --job-name=Bacphlip_Malaspina_Virathon
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=3G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python/3.8.5
module load hmmer

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --threads 4 --bacphlip True --genome_files /mnt/lustre/scratch/fcoutinho/Malaspina_Virus/Malaspina_CoAssembly_Viral_Scaffolds_Renamed.fasta --index True
