#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=01-18:00:00
#SBATCH --job-name=CheckV_Vpf_Class_RaFAH_VHMNet
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=300G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python/3.8.5
module load checkv/0.7.0
module load vpf-class
module load prodigal
module load blast/2.7.1
module load virhostmatcher-net/1.0
module load perl
module load R
module load hmmer/3.3

export PATH="/mnt/lustre/bio/users/fcoutinho/Scripts/:$PATH"
export CHECKVDB="/mnt/lustre/repos/bio/databases/public/checkv/checkv-db-v1.0/"
python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --genome_files /mnt/lustre/scratch/fcoutinho/Rumen/Assembled/Rumen_Viruses.fasta --call_checkv True --call_vpf_class True --make_pops True --call_vhm True --vhmnet_mode_short True --call_rafah True --threads 24