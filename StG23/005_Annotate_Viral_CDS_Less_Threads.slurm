#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=15-00:00:00
#SBATCH --job-name=StG23_AMG_Hunter_Virus
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=150G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load diamond/2.0.7
module load python/3.8.5
module load hmmer
#module load prodigal
module load perl

cd /mnt/lustre/scratch/fcoutinho/StG_23/Viruses/Annotation/Less_Threads/

#python3 /mnt/lustre/bio/users/fcoutinho/Repos/virathon/Virathon.py --genome_files /mnt/lustre/scratch/fcoutinho/StG_23/Viruses/Genomic_Sequences/IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.fasta --call_prodigal_module True

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/AMG_Hunter.py --annotate True --threads 23 --cds /mnt/lustre/scratch/fcoutinho/StG_23/Viruses/Annotation/IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.faa

