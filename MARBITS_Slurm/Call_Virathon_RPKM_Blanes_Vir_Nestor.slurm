#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=15-00:00:00
#SBATCH --job-name=rpkm_Blanes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=150G
#SBATCH --output=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es
##############################################
module load python/3.8.5
module load bowtie2/2.3.4.1
module load samtools

cd /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Abundance/

python3 /mnt/smart/users/fcoutinho/virathon/Virathon.py --genome_files /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/dsDNAphage_Blanes_virus.fasta --abundance_table True --abundance_rpkm True --raw_read_table /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Blanes_MG_Info.tsv --bowtie_k 1 --abundance_max_reads 10000000 --threads 47