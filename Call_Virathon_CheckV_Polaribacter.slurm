#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-03:00:00
#SBATCH --job-name=Polaribacter_CheckV_Vpf_Class_Pops_Polaribacter
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=10G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=END
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python/3.8.5
module load checkv/0.7.0
module load vpf-class
module load hmmer
module load prodigal
module load blast/2.7.1

export CHECKVDB="/mnt/lustre/repos/bio/databases/public/checkv/checkv-db-v1.0/"
python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --genome_files /mnt/lustre/scratch/fcoutinho/Polaribacter/Polaribacter_Phage_Genomes_Renamed.fasta --call_checkv True --call_vpf_class True --make_pops True --threads 24