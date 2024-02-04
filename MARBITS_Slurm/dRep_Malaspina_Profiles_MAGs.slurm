#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=5-00:00:00
#SBATCH --job-name=dRep_Malaspina_Profiles_MAGs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=250G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load drep

#dRep dereplicate /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/All_MAGs/Dereplicated_Bacteria/ -p 48 -comp 50 -con 5 -l 200000 --genomeInfo /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/All_MAGs/All_MAGs_dRep_Info.csv -g /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/All_MAGs/Bacteria/*fasta 

#dRep dereplicate /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/All_MAGs/Dereplicated_Archaea/ -p 48 -comp 50 -con 5 -l 200000 --genomeInfo /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/All_MAGs/All_MAGs_dRep_Info.csv -g /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/All_MAGs/Archaea/*fasta

dRep dereplicate /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Dereplicated_MAGs_Original_and_Refined/ -p 48 -comp 50 -con 5 -l 200000 --genomeInfo /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Profiles_Malaspina_Redo_MAG_Info_for_dRep.csv -g /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/MAGs_Original_and_Refined/*fa