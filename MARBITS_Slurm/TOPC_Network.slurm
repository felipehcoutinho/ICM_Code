#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=1-00:00:00
#SBATCH --job-name=Fastspar_TOPC_Current
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=50G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
/mnt/lustre/bio/users/fcoutinho/Fastspar/fastspar --yes --threads 48 -c /mnt/lustre/scratch/fcoutinho/StG/Networks/For_Fastspar_New_Sample_IDs_RPKM_Abundance_Renamed_Viruses_Genomes_Filtered_All_Tara_Polar_Coassemblies_Contigs.tsv -r /mnt/lustre/scratch/fcoutinho/StG/Networks/TOPC_Current_Fastspar_Correls.tsv -a /mnt/lustre/scratch/fcoutinho/StG/Networks/TOPC_Current_Fastspar_Cov.tsv