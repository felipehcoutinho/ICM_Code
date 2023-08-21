#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-00:30:00
#SBATCH --job-name=Virathon_Hmmer_Polaribacter_pfam
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=27
#SBATCH --mem=1G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python/3.8.5
module load hmmer/3.3

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --call_hmmer True --cds /mnt/lustre/scratch/fcoutinho/Polaribacter/Polaribacter_Phage_Genomes.faa --threads 27 --hmmer_program hmmsearch --hmmer_db /mnt/lustre/repos/bio/databases/public/pfam/pfam_release_34.0/Pfam-A.hmm
