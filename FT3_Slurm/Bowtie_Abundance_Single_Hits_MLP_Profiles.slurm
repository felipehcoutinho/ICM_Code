#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO

#SBATCH --time=00-05:00:00
#SBATCH --job-name=Bowtie_Single_Hits_MLP_Profiles_Prok_Virus
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=20G
#SBATCH --output=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.out
#SBATCH --error=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load cesga/2020 miniconda3/22.11.1-1 bowtie2/2.4.4 samtools/1.14

conda activate /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Conda_Envs/basic

cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/Abundance/Single_Hits/

#python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/virathon/Virathon.py --threads 47 --genome_files /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/Sequences/Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes.fasta --bowtie_k 1 --abundance_max_reads 10000000 --abundance_table True --abundance_rpkm True --raw_read_table /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Metadata/Profiles_Malaspina_76_FL_Metagenome_Info.tsv

#python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/virathon/Virathon.py --threads 47 --genome_files /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/Sequences/Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes.fasta --bowtie_k 1 --abundance_max_reads 10000000 --abundance_table True --abundance_rpkm True --raw_read_table /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Metadata/Profiles_Malaspina_76_FL_Metagenome_Info.tsv --bowtiedb /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/Abundance/Single_Hits/All_Genomic

#python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/virathon/Virathon.py --threads 47 --genome_files /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/Sequences/Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes.fasta --bowtie_k 1 --abundance_max_reads 10000000 --abundance_table True --abundance_rpkm True --raw_read_table /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/Abundance/Single_Hits/Samples_1-5_Profiles_Malaspina_76_FL_Metagenome_Info.tsv --bowtiedb /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/Abundance/Single_Hits/All_Genomic

#python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/virathon/Virathon.py --threads 47 --genome_files /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/Sequences/Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes.fasta --bowtie_k 1 --abundance_max_reads 10000000 --abundance_table True --abundance_rpkm True --raw_read_table /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/Abundance/Single_Hits/Samples_6-76_Profiles_Malaspina_76_FL_Metagenome_Info.tsv --bowtiedb /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/Abundance/Single_Hits/All_Genomic

#Parse the files again to have the abundances in a single table:
python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/virathon/Virathon.py --genome_files /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/Sequences/Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes.fasta --threads 47 --bowtie_k 1 --abundance_max_reads 10000000 --abundance_table True --abundance_rpkm True --raw_read_table  /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Metadata/Profiles_Malaspina_76_FL_Metagenome_Info.tsv --parse_only True

#Testing
#python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/virathon/Virathon.py --threads 47 --genome_files test_genomes.fasta --bowtie_k 1 --abundance_max_reads 10000000 --abundance_table True --abundance_rpkm True --raw_read_table test_info.tsv &