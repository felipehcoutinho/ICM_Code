#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=10-00:00:00
#SBATCH --job-name=QueroBins_Metabat_Malaspina_Profiles
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=50G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python
module load metabat

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/QueroBins2.py --threads 48 --min_scaffold_length 2500 --make_bins True --jgi /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/VAMB_Binning/Corrected_Profiles_Malaspina_JGI_Depth.tsv --assemblies /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Split_Assemblies/*fasta

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/QueroBins2.py --threads 48 --min_scaffold_length 2500 --make_bins True --jgi /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/VAMB_Binning/Corrected_Profiles_Malaspina_JGI_Depth.tsv --assemblies /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/CAT_Annotation/Split_Filtered_Renamed_S*fasta

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/QueroBins2.py --threads 48 --min_scaffold_length 2500 --make_bins True --jgi /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/VAMB_Binning/Corrected_Profiles_Malaspina_JGI_Depth.tsv --assemblies /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Split_Assemblies/*fasta