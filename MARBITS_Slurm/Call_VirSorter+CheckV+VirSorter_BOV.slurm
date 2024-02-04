#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=01-00:00:00
#SBATCH --job-name=BOV_VirSorter_CheckV
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=150G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
#module load virsorter
#conda activate virsorter2

#virsorter run --keep-original-seq --seqfile /mnt/lustre/scratch/fcoutinho/BOV/BOV_Filtered_Scaffolds_5kbp.fasta --include-groups dsDNAphage,ssDNA --min-length 5000 --min-score 0.5 --jobs 48 all

module load checkv/0.7.0

export CHECKVDB="/mnt/lustre/repos/bio/databases/public/checkv/checkv-db-v1.0/"
export PATH="/mnt/lustre/bio/users/fcoutinho/Scripts/:$PATH"

#checkv end_to_end /mnt/lustre/scratch/fcoutinho/BOV/final-viral-combined.fa Checkv_Filtered_BOV_Viral -t 48

#cat Checkv_Filtered_BOV_Viral/proviruses.fna Checkv_Filtered_BOV_Viral/viruses.fna > BOV_CheckV_Trimmed.fasta

#Below command for version 2.1
#virsorter run --provirus-off --prep-for-dramv --seqfile BOV_CheckV_Trimmed.fasta --include-groups dsDNAphage,ssDNA --min-length 5000 --min-score 0.5 --jobs 48 all

#Below command for version 2.2.1
#virsorter run --seqname-suffix-off --viral-gene-enrich-off --provirus-off --prep-for-dramv --seqfile BOV_CheckV_Trimmed.fasta --include-groups dsDNAphage,ssDNA --min-length 5000 --min-score 0.5 --jobs 48 all

checkv end_to_end /mnt/lustre/scratch/fcoutinho/BOV/final-viral-combined.fa Checkv_Filtered_BOV_Viral_Round2 -t 48