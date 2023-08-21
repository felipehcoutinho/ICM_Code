#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=2-00:00:00
#SBATCH --job-name=RaFAH_Shipworm_Virathon
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=100G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python/3.8.5
#module load checkv/0.7.0
#module load vpf-class
module load prodigal
#module load blast/2.7.1
#module load virhostmatcher-net/1.0
module load perl
module load R
module load hmmer/3.3

#export CHECKVDB="/mnt/lustre/repos/bio/databases/public/checkv/checkv-db-v1.0/"
export PATH="/mnt/lustre/bio/users/fcoutinho/Scripts/:$PATH"

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --genome_files /mnt/lustre/scratch/fcoutinho/Shipworm/Shipworm_Contigs+5kbp_Viruses_No_Prophage.fasta --cds /mnt/lustre/scratch/fcoutinho/Shipworm/Shipworm_Contigs+5kbp_Viruses_No_Prophage.faa --call_rafah True --threads 48

