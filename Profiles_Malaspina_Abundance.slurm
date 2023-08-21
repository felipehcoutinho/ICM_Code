#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=05-00:00:00
#SBATCH --job-name=Profiles_Malaspina_Abundance
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=150G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python/3.8.5
module load bowtie2/2.3.4.1
module load samtools

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/QueroBins2.py --calc_abundance True --threads 48 --rename False --assemblies_files /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/VAMB_Binning/Profiles_Malaspina_Scaffolds_Filtered_2500_bp.fasta --coassembly_table /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Profiles_Malaspina_Metagenome_Info.tsv

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/QueroBins2.py --calc_abundance True --threads 48 --rename False --assemblies_files /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/VAMB_Binning/Profiles_Malaspina_Scaffolds_Filtered_2500_bp.fasta --coassembly_table /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Profiles_Malaspina_Metagenome_Info_Unfinished_Abundance.tsv