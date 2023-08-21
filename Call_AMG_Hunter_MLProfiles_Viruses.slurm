#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=02-00:00:00
#SBATCH --job-name=AMG_Hunter_MLP_Viruses
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=150G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/%J.Log.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/%J.Log.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load diamond/2.0.7
module load python/3.8.5
module load hmmer
module load perl/5.28
module load prodigal

python3 /mnt/lustre/bio/users/fcoutinho/Repos/virathon/Virathon.py --threads 23 --genome_files /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/ENA_Submission/MLP_Profiles_Euk+Prok_Virus_Genomic_Sequences.fasta --call_prodigal True

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/AMG_Hunter.py --annotate True --threads 23 --cds /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/ENA_Submission/MLP_Profiles_Euk+Prok_Virus_Genomic_Sequences.faa


