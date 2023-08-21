#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=02-00:00:00
#SBATCH --job-name=Hmmer_Rumen_pfam_Virathon
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=1G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python/3.8.5
module load hmmer/3.3

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --call_hmmer True --cds /mnt/lustre/scratch/fcoutinho/Rumen/Assembled/Rumen_Viruses_Renamed.faa --threads 48 --hmmer_program hmmsearch --hmmer_db /mnt/lustre/repos/bio/databases/public/pfam/pfam_release_34.0/Pfam-A.hmm
