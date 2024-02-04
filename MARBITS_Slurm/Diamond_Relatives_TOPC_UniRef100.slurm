#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-04:00:00
#SBATCH --job-name=TOPC_Relatives_Diamond_UniRef100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=50G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load diamond/2.0.7

diamond blastp --matrix BLOSUM45 --threads 48 --more-sensitive --db /mnt/lustre/repos/bio/databases/public/UniProtKB/UniProt_Release_2021-04-07/uniref100.dmnd --outfmt 6 --query /mnt/lustre/scratch/fcoutinho/TARA_Polar/Manual_Curation/VP_32/GLUVAB3_Relatives_TOPC_Viral_Seq_95838_CDS.faa --max-target-seqs 100 --evalue 0.001 --out /mnt/lustre/scratch/fcoutinho/TARA_Polar/Manual_Curation/VP_32/GLUVAB3_Relatives_TOPC_Viral_Seq_95838_CDSxUniRef100.blastp
