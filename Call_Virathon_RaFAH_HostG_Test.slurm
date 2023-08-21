#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=2-00:00:00
#SBATCH --job-name=Virathon_RaFAH_HostG
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=250G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python/3.8.5
module load prodigal
module load perl
module load R
module load hmmer/3.3
export PATH="/mnt/lustre/bio/users/fcoutinho/Scripts/:$PATH"

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --genome_files /mnt/lustre/scratch/fcoutinho/HostG/TEST_Viral_Genomes/Test_671_PHP_Phages_Genomes.fasta --call_rafah True --threads 48

