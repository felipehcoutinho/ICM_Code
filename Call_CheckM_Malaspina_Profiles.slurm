#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=1-12:00:00
#SBATCH --job-name=CheckM_Malaspina_Profiles
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=250G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load perl
module load checkm/1.1.3
source activate

#checkm lineage_wf --tab_table --file CheckM_Bin_Info_Metabat_Bins_Round_1.tsv --threads 24 --pplacer_threads 6 --extension fa /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Metabat_Bins_Round_1/ CheckM_Results_Metabat_Bins_Round_1

#perl /mnt/lustre/bio/users/fcoutinho/Scripts/Call_Checkm_in_Chunks.pl --threads 24 --pplacer_threads 12 --extension fa --dir /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Metabat_CAT_Bins_Round_1/ &

#checkm lineage_wf -x fa -t 48 --pplacer_threads 12 --tab_table --file Chunk_6.Checkm_Results.tsv Chunk_6/ Checkm_Output_Chunk_6/

#perl /mnt/lustre/bio/users/fcoutinho/Scripts/Call_Checkm_in_Chunks.pl --threads 12 --pplacer_threads 6 --extension .fa --dir /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Metabat_CAT_Bins_Round_1/

#checkm lineage_wf -x fa -t 12 --pplacer_threads 6 --tab_table --file Refined_Bins_Checkm_Results.tsv Refined_Bins/ Checkm_Output_Refined_Bins/

perl /mnt/lustre/bio/users/fcoutinho/Scripts/Call_Checkm_in_Chunks.pl --threads 12 --pplacer_threads 6 --extension fasta --dir /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/MaxBin_Binning/MaxBin_Bins_Round_1/