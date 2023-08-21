#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-03:00:00
#SBATCH --job-name=DRAMv_MLP_Deep
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=250G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load barrnap
module load dram
source activate

DRAM-v.py annotate -i /mnt/lustre/scratch/fcoutinho/Malaspina_Virus/VirSorter2_Output/for-dramv/final-viral-combined-for-dramv.fa -v /mnt/lustre/scratch/fcoutinho/Malaspina_Virus/VirSorter2_Output/for-dramv/viral-affi-contigs-for-dramv.tab -o DRAMv_Output/ --min_contig_size 1000 --threads 24 --low_mem_mode
DRAM-v.py distill -i /mnt/lustre/scratch/fcoutinho/Malaspina_Virus/DRAMv_Output/annotations.tsv -o DRAMv_Output_Distilled


