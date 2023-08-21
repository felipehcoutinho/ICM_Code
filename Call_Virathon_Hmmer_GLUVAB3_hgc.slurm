#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-04:00:00
#SBATCH --job-name=hgcAB_Hmmer_GLUVAB3_Virathon
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=5G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python/3.8.5
module load hmmer/3.3

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --call_hmmer True --cds /mnt/lustre/repos/bio/databases/local/GLUVAB3/GLUVAB3_CDS.faa --threads 24 --hmmer_program hmmsearch --hmmer_db /mnt/lustre/scratch/fcoutinho/Mercury/hgcAB.hmm
