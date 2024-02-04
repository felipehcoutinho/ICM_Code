#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=7-00:00:00
#SBATCH --job-name=RPKM_TOPC_Virathon
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=100G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python/3.8.5
module load bowtie2
module load samtools

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --threads 24 --abundance_table True --abundance_rpkm True --genome_files /mnt/lustre/scratch/fcoutinho/TARA_Polar/Renamed_Viruses_Genomes_Filtered_All_Tara_Polar_Coassemblies_Contigs.fasta --index True --metagenomes_dir /mnt/lustre/bio/shared/topc/tara.polar.mags.mapping/data/raw/ --metagenomes_extension _clean.fastq.gz

#time python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --threads 24 --abundance_table True --genome_files /mnt/lustre/scratch/fcoutinho/TARA_Polar/Renamed_Viruses_Genomes_Filtered_All_Tara_Polar_Coassemblies_Contigs.fasta --bowtiedb /mnt/lustre/scratch/fcoutinho/TARA_Polar/Abundance/Renamed_Viruses_Genomes_Filtered_All_Tara_Polar_Coassemblies_Contigs --index True --metagenomes_dir /mnt/lustre/scratch/fcoutinho/TARA_Polar/Abundance/Test_Samples/ --metagenomes_extension _clean.fastq.gz

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --threads 48 --abundance_table True --genome_files /mnt/lustre/scratch/fcoutinho/TARA_Polar/Renamed_Viruses_Genomes_Filtered_All_Tara_Polar_Coassemblies_Contigs.fasta --index True --bowtiedb /mnt/lustre/scratch/fcoutinho/TARA_Polar/Abundance/Renamed_Viruses_Genomes_Filtered_All_Tara_Polar_Coassemblies_Contigs --metagenomes_dir /mnt/lustre/scratch/fcoutinho/TARA_Polar/Abundance/Test_Samples/ --metagenomes_extension _clean.fastq.gz --abundance_max_reads 1000000


python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --threads 48 --abundance_table True --genome_files /mnt/lustre/scratch/fcoutinho/TARA_Polar/Renamed_Viruses_Genomes_Filtered_All_Tara_Polar_Coassemblies_Contigs.fasta --index True --bowtiedb /mnt/lustre/scratch/fcoutinho/TARA_Polar/Abundance/Renamed_Viruses_Genomes_Filtered_All_Tara_Polar_Coassemblies_Contigs --metagenomes_dir /mnt/lustre/bio/shared/topc/tara.polar.mags.mapping/data/raw/ --metagenomes_extension _clean.fastq.gz --abundance_max_reads 10000000