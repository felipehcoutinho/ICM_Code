#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=21-00:00:00
#SBATCH --job-name=Diamond_Profiles_Malaspina
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=50G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load diamond/2.0.7

#while read file ; do read_r1=$(echo "$file" | cut -f 2); sample=$(echo "$file" | cut -f 1); echo "$sample $read_r1" ; done < /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Profiles_Malaspina_Metagenome_Info.tsv

while read file ; do read_r1=$(echo "$file" | cut -f 2); sample=$(echo "$file" | cut -f 1); echo "$sample $read_r1" ; diamond blastx --matrix BLOSUM45 --threads 48 --db /mnt/lustre/repos/bio/databases/public/UniProtKB/UniProt_Release_2021-04-07/uniref100.dmnd --outfmt 6 --query $read_r1 --max-target-seqs 10 --evalue 0.00001 --out $sample-X-UniRef100.blastx ; done < /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Profiles_Malaspina_Metagenome_Info.tsv




