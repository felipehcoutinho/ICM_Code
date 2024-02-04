#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=01-00:00:00
#SBATCH --job-name=Index_Profiles_Malaspina_Scaffolds_MAGs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Index_Malaspina_Profiles_Scaffold_ID.py --assemblies_dir /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Filtered_Renamed_Scafolds/ --table /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/VAMB_Binning/Scaffolds_Info.tsv --bins_dir /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Original_and_Refined_Bins/ --scaffold_info Profiles_Malaspina_All_Scaffolds_Info_Redo.tsv --bin_info Profiles_Malaspina_All_Bin_Info_Redo.tsv

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Index_Malaspina_Profiles_Scaffold_ID.py --assemblies_dir /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Filtered_Renamed_Scafolds/ --table /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/VAMB_Binning/Scaffolds_Info.tsv --bins_dir /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/MAGs_Original_and_Refined/ --scaffold_info Profiles_Malaspina_MAG_Scaffolds_Info_Redo.tsv --bin_info Profiles_Malaspina_MAG_Info_Redo.tsv

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Get_Seqs.py --input_s ../Filtered_Renamed_S*fasta --list List_Non_MAG_Scaffolds_Redo.txt --output_sequences Profiles_Malaspina_Non_MAG_Scaffolds.fasta