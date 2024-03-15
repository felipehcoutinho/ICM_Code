#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=03-00:00:00
#SBATCH --job-name=cas_MLP_MAGs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=15G
#SBATCH --output=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load python/2.7.15 hmmer/3.2

#cd /mnt/smart/users/fcoutinho/DeepCas/MLP_Profiles_cas_Search/
cd /mnt/smart/users/fcoutinho/malaspina.crispr-cas/
mkdir MLP_MAGs_CDSxcas_db.updated2019

/mnt/smart/users/fcoutinho/malaspina.crispr-cas/code/hmmscanner.sh /mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Genes_and_CDS_All_MAGs/Profiles_Malaspina_Non_MAG_Scaffolds_CDS.faa MLP_MAGs_CDSxcas_db.updated2019
