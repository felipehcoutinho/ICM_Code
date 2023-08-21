#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=10-00:00:00
#SBATCH --job-name=ProkAbdStG23
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=150G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load python/3.8.5
module load bowtie2/2.3.4.1
module load samtools

python3 /mnt/lustre/bio/users/fcoutinho/Repos/virathon/Virathon.py --genome_files /mnt/lustre/bio/users/fcoutinho/OceanDNA_MAGs/Species_Rep_MAGs/*fa --merge_seqs True --out_seq_file OceanDNA_All_Species_Rep_MAGs_Scaffolds.fasta --threads 47 --bowtie_k 1 --abundance_max_reads 100000000 --abundance_table True --abundance_rpkm True --raw_read_table /mnt/lustre/scratch/fcoutinho/StG_23/Metagenomes/TARA_MGs_Paired_With_Salazar_Data_Info+Metadata.tsv