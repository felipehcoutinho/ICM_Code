#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=05-00:00:00
#SBATCH --job-name=VirAbdStG23
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=100G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load python/3.8.5
module load bowtie2/2.3.4.1
module load samtools

#python3 /mnt/lustre/bio/users/fcoutinho/Repos/virathon/Virathon.py --genome /mnt/lustre/scratch/fcoutinho/StG_23/Viruses/Genomic_Sequences/IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.fasta --abundance_table True --abundance_rpkm True --raw_read_table /mnt/lustre/scratch/fcoutinho/StG_23/Metagenomes/Virome_Info.tsv

python3 /mnt/lustre/bio/users/fcoutinho/Repos/virathon/Virathon.py --genome /mnt/lustre/scratch/fcoutinho/StG_23/Viruses/Genomic_Sequences/IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.fasta --bowtiedb /mnt/lustre/scratch/fcoutinho/StG_23/Viruses/Abundance/IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives --threads 47 --bowtie_k 1 --abundance_max_reads 100000000 --abundance_table True --abundance_rpkm True --raw_read_table /mnt/lustre/scratch/fcoutinho/StG_23/Metagenomes/Virome_Info.tsv