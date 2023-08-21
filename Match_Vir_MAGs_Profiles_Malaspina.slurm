#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=02-00:00:00
#SBATCH --job-name=Match_Vir_MAG_Profiles_Malaspina
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com
#SBATCH --constraint=supermicro

##############################################
module load python/3.8.5
module load blast/2.7.1

cd /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Link_Viral_Seqs_to_MAGs.py --viral_genomes /mnt/lustre/scratch/elopez/3_checkV_output/all_pro_viruses_checkV.fna --mags_dir /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Dereplicated_MAGs_Original_and_Refined/dereplicated_genomes/ --mags_ext fa --threads 24 --scaffold_specific True --remove_virus True

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Explode_Fasta.py --assembly /mnt/lustre/scratch/elopez/4_checkV_output/all_pro_viruses_checkV.fna

cd /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/

python3 /mnt/lustre/bio/users/fcoutinho/PHIST/phist.py -t 24 /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/Profiles_Malaspina_Viral_Genomes/  /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/No_Virus_MAGs/ /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/PHIST_Output/PHIST_Kmers.csv /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Viruses/PHIST_Output/PHIST_Predictions.csv