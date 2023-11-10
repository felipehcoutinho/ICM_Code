#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO

#SBATCH --time=00-06:00:00
#SBATCH --job-name=PHIST_MLP_Profiles
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=10G
#SBATCH --output=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.out
#SBATCH --error=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load cesga/2020 blast/2.13.0-Linux_x86_64 miniconda3/22.11.1-1

conda activate /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Conda_Envs/basic

export PATH="/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Software/PHIST:$PATH"

#cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/PHIST_Predictions/

#python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/virathon/Virathon.py --threads 47 --genome_files /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/Sequences/Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes.fasta --phist_host_prediction True --putative_host_genomes_directory /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Databases/OceanDNA/All_MAGs/ --remove_exact_matches True --extension_putative_host_genomes fa

#python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/virathon/Virathon.py --threads 47 --genome_files /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/Sequences/Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes.fasta --phist_host_prediction True --putative_host_genomes_directory /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Databases/OceanDNA/All_MAGs/ --remove_exact_matches True --extension_putative_host_genomes fa --parse_only True

#above command died due to exceeding disk quota of number of files

#Explode fasta of viral genomes
#cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/PHIST_Predictions/Viral_Genomes_PHIST/
#python3 $code/Explode_Fasta.py --assembly /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/Sequences/Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes.fasta
#Run phist manually
cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/PHIST_Predictions/
phist.py -t 47 Viral_Genomes_PHIST/ No_Vir_Host_Genomes/ PHIST_Output/
grep "," PHIST_Output/predictions.csv > PHIST_Output/positive_predictions.csv
#Parse otuput with Virathon
python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/virathon/Virathon.py --threads 47 --genome_files /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/Sequences/Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes.fasta --phist_host_prediction True --putative_host_genomes_directory /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Databases/OceanDNA/All_MAGs/ --remove_exact_matches True --extension_putative_host_genomes fa --parse_only True

#Testing
#cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/PHIST_Predictions/test/
#python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/virathon/Virathon.py --threads 47 --genome_files /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/Sequences/Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes.fasta --phist_host_prediction True --putative_host_genomes_directory /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/PHIST_Predictions/test/Hosts/ --remove_exact_matches True --extension_putative_host_genomes fa &

#python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/virathon/Virathon.py --threads 47 --genome_files /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/Sequences/Profiles_Malaspina_Viruses_CheckV_Trimmed_Genomes.fasta --phist_host_prediction True --putative_host_genomes_directory /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Malaspina/Profiles/Prok_Viruses/PHIST_Predictions/test/Hosts/ --remove_exact_matches True --extension_putative_host_genomes fa --parse_only True &

#rm -f DB_All_Host_Seqs.n* All*fasta *.tsv *blastn ; rm -fr No_Vir_Host_Genomes/ PHIST_Output/ Viral_Genomes_PHIST/