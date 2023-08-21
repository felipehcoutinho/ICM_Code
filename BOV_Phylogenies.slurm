#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=10-00:00:00
#SBATCH --job-name=BOV_Phylogenies
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python
module load hmmer
module load muscle
module load fasttree
module load cd-hit


#ribH
#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/AMG_Hunter.py --threads 24 --kegg_phylogeny True --kegg_ko K00794 --ref /mnt/lustre/repos/bio/databases/public/UniProtKB/UniProt_Release_2021-04-07/uniref100.fasta --cds /mnt/lustre/scratch/fcoutinho/BOV/Renamed_BOV_CheckV_Trimmed.faa --derep_id 0.9

#ribBA
python3 /mnt/lustre/bio/users/fcoutinho/Scripts/AMG_Hunter.py --threads 24 --kegg_phylogeny True --kegg_ko K14652 --ref /mnt/lustre/repos/bio/databases/public/UniProtKB/UniProt_Release_2021-04-07/uniref100.fasta --cds /mnt/lustre/scratch/fcoutinho/BOV/Renamed_BOV_CheckV_Trimmed.faa --derep_id 0.9

#DHFR
python3 /mnt/lustre/bio/users/fcoutinho/Scripts/AMG_Hunter.py --threads 24 --kegg_phylogeny True --kegg_ko K00287 --ref /mnt/lustre/repos/bio/databases/public/UniProtKB/UniProt_Release_2021-04-07/uniref100.fasta --cds /mnt/lustre/scratch/fcoutinho/BOV/Renamed_BOV_CheckV_Trimmed.faa --derep_id 0.9




