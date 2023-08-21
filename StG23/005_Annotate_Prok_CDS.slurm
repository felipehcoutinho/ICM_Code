#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=15-00:00:00
#SBATCH --job-name=StG23_AMG_Hunter_Prok
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=250G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load prokka
conda activate prokka-v1.14.6

cd /mnt/lustre/bio/users/fcoutinho/OceanDNA_MAGs/CDS_All_MAGs/

while read file ; do prokka --kingdom Archaea --cpus 47 --noanno --outdir Prokka_$file --prefix $file --locustag $file /mnt/lustre/bio/users/fcoutinho/OceanDNA_MAGs/All_MAGs/$file.fa ; cat Prokka_$file/$file.faa >> /mnt/lustre/bio/users/fcoutinho/OceanDNA_MAGs/All_OCeanDNA_MAGs_PEGs.faa; done < /mnt/lustre/bio/users/fcoutinho/OceanDNA_MAGs/OceanDNA_Archaeal_MAGs_List.txt

while read file ; do prokka --kingdom Bacteria --cpus 47 --noanno --outdir Prokka_$file --prefix $file --locustag $file /mnt/lustre/bio/users/fcoutinho/OceanDNA_MAGs/All_MAGs/$file.fa ; cat Prokka_$file/$file.faa  >> /mnt/lustre/bio/users/fcoutinho/OceanDNA_MAGs/All_OCeanDNA_MAGs_PEGs.faa; done < /mnt/lustre/bio/users/fcoutinho/OceanDNA_MAGs/OceanDNA_Bacterial_MAGs_List.txt

conda deactivate

module load diamond/2.0.7
module load python/3.8.5
module load hmmer
module load prodigal

cd /mnt/lustre/scratch/fcoutinho/StG_23/Prokaryotes/Annotation/

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/AMG_Hunter.py --annotate True --threads 47 --cds  /mnt/lustre/bio/users/fcoutinho/OceanDNA_MAGs/All_OCeanDNA_MAGs_PEGs.faa

