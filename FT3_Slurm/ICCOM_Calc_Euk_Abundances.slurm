#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO

#SBATCH --time=01-12:00:00
#SBATCH --job-name=ICCOM_Calc_Euk_Abund
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=10G
#SBATCH --output=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.out
#SBATCH --error=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load cesga/2020 miniconda3/22.11.1-1 bowtie2/2.4.4 samtools/1.14

conda activate /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Conda_Envs/basic

cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/Eukaryotes/Abundance/

python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/virathon/Virathon.py --genome_files /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/ICCOMM/Eukaryotes/Sequences/Delmont_MAGs_Scaffolds.fasta --threads 47 --bowtie_k 1 --abundance_max_reads 10000000 --abundance_table True --abundance_rpkm True