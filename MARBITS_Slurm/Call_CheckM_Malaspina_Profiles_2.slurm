#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=07-00:00:00
#SBATCH --job-name=CheckM_Malaspina_Profiles_2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=750G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load perl
module load checkm/1.1.3
source activate

#perl /mnt/lustre/bio/users/fcoutinho/Scripts/Call_Checkm_in_Chunks.pl --threads 12 --pplacer_threads 12 --chunk_size 1000 --directory /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Metabat_Bins_Round_1_Salvage/Refined_Bins/

#perl /mnt/lustre/bio/users/fcoutinho/Scripts/Call_Checkm_in_Chunks.pl --threads 12 --pplacer_threads 12 --chunk_size 1000 --directory /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/All_Bins/

perl /mnt/lustre/bio/users/fcoutinho/Scripts/Call_Checkm_in_Chunks.pl --threads 12 --pplacer_threads 12 --chunk_size 1000 --directory /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Refined_Bins/