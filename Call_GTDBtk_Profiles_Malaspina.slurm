#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=3-00:00:00
#SBATCH --job-name=GTDBtk_Profiles_Malaspina
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=250G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load gtdbtk/1.5
source activate

#gtdbtk classify_wf --extension fa --genome_dir /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Metabat_MAGs_Round_1/ --out_dir /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Profiles_Malaspina_Metabat_MAGs_Round_1_Output/ --pplacer_cpus 12 --cpus 48 --prefix Profiles_Malaspina_Metabat_MAGs_Round_1_GTDBtk

#gtdbtk classify_wf --extension fa --genome_dir /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Metabat_Bins_Round_1_Salvage/MAGs/ --out_dir /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Metabat_Bins_Round_1_Salvage/MAGs/Profiles_Malaspina_Metabat_Refined_MAGs_Round_1_Output/ --pplacer_cpus 12 --cpus 48 --prefix Profiles_Malaspina_Metabat_Refined_MAGs_Round_1_GTDBtk

#gtdbtk classify_wf --extension fasta --genome_dir /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/All_MAGs/Archaea/ --out_dir /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/All_MAGs/Profiles_Malaspina_All_Archaeal_MAGs_GTDBtk_Output/ --pplacer_cpus 12 --cpus 48 --prefix Profiles_Malaspina_All_Archaeal_MAGs_GTDBtk_Output

#gtdbtk classify_wf --extension fasta --genome_dir /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/All_MAGs/Bacteria/ --out_dir /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/All_MAGs/Profiles_Malaspina_All_Bacterial_MAGs_GTDBtk_Output/ --pplacer_cpus 12 --cpus 48 --prefix Profiles_Malaspina_All_Bacterial_MAGs_GTDBtk_Output

gtdbtk classify_wf --extension fa --genome_dir /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/MAGs_Original_and_Refined/ --out_dir /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/MAGs_Original_and_Refined_GTDBtk_Output/ --pplacer_cpus 12 --cpus 48 --prefix Profiles_Malaspina_MetaBat_Redo_MAGs_GTDBtk_Output



