#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-18:00:00
#SBATCH --job-name=Malaspina_Profiles_Counts
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%A_%a.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.csic.es
#SBATCH --array=1-5%3

##############################################
module load python

#DIR_LIST=($(</mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_MLProfiles_Large_Fraction/Malaspina_Large_Fraction_Dir_List.txt))
#DIR_LIST=($(</mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Metadata/Unique_Dir_List_Kaust_Malaspina_Profiles_MG_Large_Fraction_3-20.txt))
DIR_LIST=($(</mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Metadata/Missing_Dir_List_Kaust_Malaspina_Profiles_MG_Large_Fraction_3-20.txt))
SAMPLE_DIR=${DIR_LIST[${SLURM_ARRAY_TASK_ID}]}
SAMPLE_ID=$(echo "$SAMPLE_DIR" | cut -d "/" -f 13)
#SAMPLE_ID=$(echo "$SAMPLE_DIR" | cut -d "/" -f 12)

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Count_Reads_and_Length.py --format fastq --extension .fastq.gz --output Profiles_Malaspina_Large_Fraction_${SAMPLE_ID}_Read_Counts.tsv --directories $SAMPLE_DIR

#python3 $scripts/Count_Reads_and_Length.py --format fastq --extension fastq.gz --output Profiles_Malaspina_Large_Fraction_Read_Counts.tsv --directories /mnt/lustre/repos/bio/projects/malaspina/release1/station019/Epi/3-20/metaDNA/MP0310/rawDatas/cleanReads/ /mnt/lustre/repos/bio/projects/malaspina/release1/station019/Epi/3-20/metaDNA/MP0320/rawDatas/cleanReads/ /mnt/lustre/repos/bio/projects/malaspina/release1/station019/Epi/3-20/metaDNA/MP0322/rawDatas/cleanReads/ /mnt/lustre/repos/bio/projects/malaspina/release1/station092/Epi/3-20/metaDNA/MP1671/rawDatas/cleanReads/ /mnt/lustre/repos/bio/projects/malaspina/release1/station092/Bathy/3-20/metaDNA/MP1673-1675m/rawDatas/cleanReads/ /mnt/lustre/repos/bio/projects/malaspina/release1/station092/Meso/3-20/metaDNA/MP1677/rawDatas/cleanReads/ /mnt/lustre/repos/bio/projects/malaspina/release1/station092/Meso/3-20/metaDNA/MP1679/rawDatas/cleanReads/ /mnt/lustre/repos/bio/projects/malaspina/release1/station092/Meso/3-20/metaDNA/MP1681/rawDatas/cleanReads/ /mnt/lustre/repos/bio/projects/malaspina/release1/station092/Epi/3-20/metaDNA/MP1683/rawDatas/cleanReads/ /mnt/lustre/repos/bio/projects/malaspina/release1/station126/Bathy/3-20/metaDNA/MP2367/rawDatas/cleanReads/ /mnt/lustre/repos/bio/projects/malaspina/release1/station098/Meso/3-20/metaDNA/MP1814/rawDatas/cleanReads/

#python3 $scripts/Count_Reads_and_Length.py --format fastq --extension fastq.gz --output Profiles_Malaspina_Missing_Samples_Read_Counts.tsv --directories /mnt/lustre/repos/bio/projects/malaspina/release1/station115/Meso/3-20/metaDNA/MP2122/rawDatas/cleanReads/ /mnt/lustre/repos/bio/projects/malaspina/release1/station135/Meso/3-20/metaDNA/MP2657/rawDatas/cleanReads/ /mnt/lustre/repos/bio/projects/malaspina/release1/station046/Meso/3-20/metaDNA/MP0843/rawDatas/cleanReads/

#grep -v "Total_Reads" Profiles_Malaspina_Missing_Samples_Read_Counts.tsv >> Profiles_Malaspina_Large_Fraction_Read_Counts.tsv



