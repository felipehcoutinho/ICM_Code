#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=1-00:00:00
#SBATCH --job-name=Compress_Scratch_Files
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
#cd /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/

#tar -zcvf Profiles_Malaspina_MaxBin_Results.tgz MaxBin_Binning/
#mv Profiles_Malaspina_MaxBin_Results.tgz /mnt/lustre/bio/users/fcoutinho/Profiles_Malaspina/Assemblies_Round2/MaxBin_Binning/

#tar -zcvf Profiles_Malaspina_All_Scaffolds_Abundance.tgz Scaffold_Abundances/
#mv Profiles_Malaspina_All_Scaffolds_Abundance.tgz /mnt/lustre/bio/users/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Scaffolds_Abundance/

#cd /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/
#tar -zcvf Profiles_Malaspina_MetaBat_Bins_Redo_Originals.tgz All_Bins/
#mv Profiles_Malaspina_MetaBat_Bins_Redo_Originals.tgz MetaBat_Binning_Redo/

#cd /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/
#tar -zcvf Profiles_Malaspina_MetaBat_MAGs_Redo_Originals_and_Refined.tgz MAGs_Original_and_Refined/
#mv Profiles_Malaspina_MetaBat_MAGs_Redo_Originals_and_Refined.tgz /mnt/lustre/bio/users/fcoutinho/Profiles_Malaspina/Assemblies_Round2/MetaBat_Binning_Redo/

#tar -zcvf Profiles_Malaspina_MetaBat_Bins_Redo_Refined.tgz Refined_Bins/
#mv Profiles_Malaspina_MetaBat_Bins_Redo_Refined.tgz /mnt/lustre/bio/users/fcoutinho/Profiles_Malaspina/Assemblies_Round2/MetaBat_Binning_Redo/


#cd /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Dereplicated_MAGs_Original_and_Refined/
#tar -zcvf Profiles_Malaspina_MetaBat_MAGs_Redo_Dereplicated_Genomes.tgz  dereplicated_genomes/
#mv Profiles_Malaspina_MetaBat_MAGs_Redo_Dereplicated_Genomes.tgz /mnt/lustre/bio/users/fcoutinho/Profiles_Malaspina/Assemblies_Round2/MetaBat_Binning_Redo/

cd /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/
tar -zcvf Profiles_Malaspina_MetaBat_Bins_Redo_Prokka.tgz Prokka_All_MAGs/
mv Profiles_Malaspina_MetaBat_Bins_Redo_Prokka.tgz /mnt/lustre/bio/users/fcoutinho/Profiles_Malaspina/Assemblies_Round2/MetaBat_Binning_Redo/




