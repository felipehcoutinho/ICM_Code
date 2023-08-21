#!/bin/bash
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=24-00:00:00
#SBATCH --job-name=MaxBin_Malaspina_Profiles
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=150G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load perl
module load maxbin

for i in /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Split_Assemblies/*fasta ; do echo "Processing $i"; sid=$(echo $i | cut -d "/" -f 10 | cut -d "." -f 1) ; run_MaxBin.pl -thread 48 -contig $i -abund_list /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/MaxBin_Binning/List_Abd_Files.txt -out Bin_$sid ; done;