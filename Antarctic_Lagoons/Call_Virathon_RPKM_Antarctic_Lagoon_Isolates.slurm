#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-12:00:00
#SBATCH --job-name=rpkm_lagoon_isolates
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=15G
#SBATCH --output=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es
##############################################
module load python/3.8.5
module load bowtie2/2.3.4.1
module load samtools

#Rename seqs with unique IDs so they can be linked to appropriate files

cd /mnt/smart/users/fcoutinho/Antarctic_Lagoons/GENOMES_Pyscrobacter_Sebas

for i in SQK*/MIFEST*/*.fna; do uid=$(echo "$i" | cut -d "/" -f 2 | sed 's/\-.*//g'); python3 /mnt/smart/users/fcoutinho/Repos/virathon/Virathon.py --genome_files $i --rename_seqs True --string_rename ${uid}_ ; mv All_Genomic.fasta ${uid}.fasta; done

#Calc abundances 

cd /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Isolate_Abundances/

python3 /mnt/smart/users/fcoutinho/Repos/virathon/Virathon.py --genome_files /mnt/smart/users/fcoutinho/Antarctic_Lagoons/GENOMES_Pyscrobacter_Sebas/*fasta --abundance_table True --abundance_rpkm True  --bowtie_k 10 --abundance_max_reads 10000000 --threads 23 --raw_read_table /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Metadata/WGS_Samples.tsv