#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=05-12:00:00
#SBATCH --job-name=RPKM_ZOFE_Virathon
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load python/3.8.5
module load bowtie2/2.3.4.1
module load samtools

python3 /mnt/lustre/bio/users/fcoutinho/Repos/virathon/Virathon.py --genome_files /mnt/lustre/scratch/fcoutinho/Soil_Elena/Post_CheckV_Trimmed_True_Virus_Soil_Genomes.fasta --abundance_table True --abundance_rpkm True --raw_read_table /mnt/lustre/scratch/fcoutinho/Soil_Elena/Soil_Elena_MG_Info.tsv --bowtie_k 1 --abundance_max_reads 50000000 --threads 23