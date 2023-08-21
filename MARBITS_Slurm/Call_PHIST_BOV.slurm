#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=05-00:00:00
#SBATCH --job-name=PHIST_BOV
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=250G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load python

#python3 /mnt/lustre/bio/users/fcoutinho/PHIST/phist.py -t 47 /mnt/lustre/scratch/fcoutinho/BOV/Split_BOV_Genomes/ /mnt/lustre/bio/users/fcoutinho/OceanDNA_MAGs/Species_Rep_MAGs/ ./PHIST_Output/PHIST_Kmers.csv ./PHIST_Output/PHIST_Predictions.csv

cd /mnt/lustre/scratch/fcoutinho/BOV/PHIST_Preds_OceanDNA_No_Species_Rep_MAGs/
python3 /mnt/lustre/bio/users/fcoutinho/PHIST/phist.py -t 47 /mnt/lustre/scratch/fcoutinho/BOV/Split_BOV_Genomes/ /mnt/lustre/bio/users/fcoutinho/OceanDNA_MAGs/Non_Species_Rep_MAGs/ ./PHIST_Output/PHIST_Kmers.csv ./PHIST_Output/PHIST_Predictions.csv