#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=02-06:00:00
#SBATCH --job-name=rpkm_omdb_mlp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=750G
#SBATCH --output=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es
##############################################
module load python/3.8.5
module load bowtie2/2.3.4.1
module load samtools

#reads in:
# /mnt/biostore/data/fs/bio/malaspina/vertical-profiles-clean-reads/

cd /mnt/smart/scratch/vir/felipe/Databases/OMDB/

bowtie2-build --threads 47 /mnt/smart/scratch/vir/felipe/Databases/OMDB/Renamed_OMDBv2.0_SC_G_R.fa Renamed_OMDBv2.0_SC_G_R

# python3 /mnt/smart/scratch/vir/felipe/Databases/OMDB/Renamed_OMDBv2.0_SC_G_R.fa --abundance_table True --abundance_rpkm True --bowtie_k 10 --abundance_max_reads 10000000 --threads 47 --bowtiedb Renamed_OMDBv2.0_SC_G_R --raw_read_table /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Profiles_Malaspina_Metagenome_Info.tsv

#Submitted batch job 988468