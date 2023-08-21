#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=01-00:00:00
#SBATCH --job-name=Match_Vir_MAG_Profiles_Deep_Malaspina
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python/3.8.5
module load blast/2.7.1

cd /mnt/lustre/scratch/fcoutinho/Malaspina_Virus/PHIST_Preds

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Link_Viral_Seqs_to_MAGs.py --viral_genomes /mnt/lustre/scratch/fcoutinho/Malaspina_Virus/Malaspina_CoAssembly_Viral_Scaffolds_Renamed.fasta --mags_dir /mnt/lustre/repos/bio/projects/malaspina/release1/MAGs/01_bathypelagic-MAGs-JGI/01_sequences/01_317-hq-mags/ --mags_ext fasta --threads 24 --remove_virus True

#mv No_Virus*fasta /mnt/lustre/scratch/fcoutinho/Malaspina_Virus/PHIST_Preds/No_Virus_MAGs/

python3 /mnt/lustre/bio/users/fcoutinho/PHIST/phist.py -t 24 /mnt/lustre/scratch/fcoutinho/Malaspina_Virus/Exploded_CheckV_Malaspina_CoAssembly_Viral_Scaffolds_Renamed/ /mnt/lustre/scratch/fcoutinho/Malaspina_Virus/PHIST_Preds/No_Virus_MAGs/ /mnt/lustre/scratch/fcoutinho/Malaspina_Virus/PHIST_Preds/PHIST_Output/PHIST_Kmers.csv /mnt/lustre/scratch/fcoutinho/Malaspina_Virus/PHIST_Preds/PHIST_Output/PHIST_Predictions.csv