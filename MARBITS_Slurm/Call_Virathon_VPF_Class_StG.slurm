#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=3-00:00:00
#SBATCH --job-name=VPF_Class_StG_Virathon
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=150G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python/3.8.5
module load vpf-class
module load prodigal/2.6.3
module load hmmer/3.3

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --threads 48 --call_vpf True --genome_files /mnt/lustre/scratch/fcoutinho/StG/Min_Count_1k_Cluster_95ID_80Cov_GLUVAB3_VIBRANT_Genomes_Min_5kbp_rep_seq.fasta --cds /mnt/lustre/scratch/fcoutinho/StG/Min_Count_1k_Cluster_95ID_80Cov_GLUVAB3_VIBRANT_Genomes_Min_5kbp_rep_seq_CDS.faa