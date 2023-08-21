#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=05-00:00:00
#SBATCH --job-name=DefenseFinder
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=50G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%A_%a.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.cisc.es
#SBATCH --array=1-2672%12

##############################################
module load defensefinder/1.0.9
source activate
#cut -f 43 ../Info_All_MAGs_Profiles_Malaspina_Redo.tsv | grep -v "File" | sed "s///"
CDSFILE=$(ls /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/CDS_by_MAG/*.faa | sed -n ${SLURM_ARRAY_TASK_ID}p)
defense-finder run --workers 4 --out-dir /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Defense_Finder_Results/${SLURM_ARRAY_TASK_ID}p --models-dir /mnt/lustre/repos/bio/databases/public/defensefinder/models/ $CDSFILE 

#defense-finder run --workers 1 --out-dir /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Defense_Finder_Results/TestDF_Ordered_Bin_S11.9 --models-dir /mnt/lustre/repos/bio/databases/public/defensefinder/models/ /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/CDS_by_MAG/Bin_S11.9.faa &> TestDF_Ordered_Bin_S11.9.out.txt
#defense-finder run --workers 1 --db-type unordered --out-dir /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Defense_Finder_Results/TestDF_Unordered_Bin_S11.9 --models-dir /mnt/lustre/repos/bio/databases/public/defensefinder/models/ /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/CDS_by_MAG/Bin_S11.9.faa &> TestDF_Unordered_Bin_S11.9.out.txt
#module load prokka
#conda activate prokka-v1.14.6
#prokka --kingdom Bacteria --cpus 24 --noanno --outdir Prokka_NZ_CP014270 --prefix NZ_CP014270 --locustag NZ_CP014270 /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Defense_Finder_Results/NZ_CP014270.1.fasta
#conda deactivate

#defense-finder run --workers 4 --out-dir /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Defense_Finder_Results/TestDF_Ordered_EColi_K12 --models-dir /mnt/lustre/repos/bio/databases/public/defensefinder/models/ /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Defense_Finder_Results/Prokka_NZ_CP014270/NZ_CP014270.faa &> TestDF_Ordered_NZ_CP014270.1.out.txt

#defense-finder run --workers 24 --db-type unordered --out-dir /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Defense_Finder_Results/TestDF_Unordered_EColi_K12 --models-dir /mnt/lustre/repos/bio/databases/public/defensefinder/models/ /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Defense_Finder_Results/Prokka_NZ_CP014270/NZ_CP014270.faa &> TestDF_Unordered_NZ_CP014270.1.out.txt
