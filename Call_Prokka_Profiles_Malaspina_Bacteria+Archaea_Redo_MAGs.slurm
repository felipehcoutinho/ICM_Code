#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=04-00:00:00
#SBATCH --job-name=Prokka_Profiles_Malaspina_Redo
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=15G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load prokka
conda activate prokka-v1.14.6

while read file ; do prokka --kingdom Archaea --cpus 24 --noanno --outdir Prokka_$file --prefix $file --locustag $file /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/MAGs_Original_and_Refined/$file.fa ; done < /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Prokka_All_MAGs/List_All_Archaea_MAGs.txt

while read file ; do prokka --kingdom Bacteria --cpus 24 --noanno --outdir Prokka_$file --prefix $file --locustag $file /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/MAGs_Original_and_Refined/$file.fa ; done < /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Prokka_All_MAGs/List_All_Bacteria_MAGs.txt


