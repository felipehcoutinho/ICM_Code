#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=01-12:00:00
#SBATCH --job-name=fastANI_Profiles_Malaspina
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load fastani

fastANI --ql /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/All_MAGs/ANI_Comps/List_Full_Path_Profiles_Malaspina_dRep_MAGs.txt --rl /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/All_MAGs/ANI_Comps/List_Full_Path_TOPC+Mlaspina_Deep_MAGs.txt --threads 24 --output fastANI_Result_Profiles_Malaspina_dRep_MAGsxTOPC_Deep_Malaspina_MAGs.tsv 