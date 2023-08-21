#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=03-00:00:00
#SBATCH --job-name=StG_PHIST_OceanDNA
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=250G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load python
cd /mnt/lustre/scratch/fcoutinho/StG_23/Viruses/Genomic_Sequences/Split_IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives/
python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Explode_Fasta.py --assembly /mnt/lustre/scratch/fcoutinho/StG_23/Viruses/Genomic_Sequences/IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.fasta

cd /mnt/lustre/scratch/fcoutinho/StG_23/Viruses/PHIST_OceanDNA/
mkdir PHIST_Output
python3 /mnt/lustre/bio/users/fcoutinho/PHIST/phist.py -t 47 /mnt/lustre/scratch/fcoutinho/StG_23/Viruses/Genomic_Sequences/Split_IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives/ /mnt/lustre/bio/users/fcoutinho/OceanDNA_MAGs/All_MAGs/ ./PHIST_Output/PHIST_Kmers.csv ./PHIST_Output/PHIST_Predictions.csv