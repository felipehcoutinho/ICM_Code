#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO

#SBATCH --time=00-05:00:00
#SBATCH --job-name=vRhyme_MLP_Profiles_Prok_Virus
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=20G
#SBATCH --output=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.out
#SBATCH --error=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load cesga/system miniconda3/22.11.1-1

conda activate /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Conda_Envs/vrhyme

cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/vRhyme_Binning/

python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Software/vRhyme/vRhyme/vRhyme -i /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/Sequences/Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes.fasta -b /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/Abundance/Single_Hits/*bam -t 47 -o /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/vRhyme_Binning/Default_Params_Single_Hits/ -g /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/Sequences/Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes.fna -p /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/Sequences/Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes.faa -l 2000