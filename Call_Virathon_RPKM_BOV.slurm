#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=01-12:00:00
#SBATCH --job-name=RPKM_BOV_Virathon
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=50G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load python/3.8.5
module load bowtie2/2.3.4.1
module load samtools

python3 /mnt/lustre/bio/users/fcoutinho/Repos/virathon/Virathon.py --genome_files /mnt/lustre/scratch/fcoutinho/BOV/Renamed_BOV_CheckV_Trimmed.fasta --abundance_table True --abundance_rpkm True --raw_read_table /mnt/lustre/scratch/fcoutinho/BOV/Metadata/BOV_Metadata_Edited.tsv --bowtie_k 1 --threads 47