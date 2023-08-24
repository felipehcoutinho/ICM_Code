#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --time=02-23:00:00
#SBATCH --job-name=StG23_AMG_Hunter_Virus
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=100G
#SBATCH --output=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.out
#SBATCH --error=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load cesga/system miniconda3/22.11.1-1
conda activate /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Conda_Envs/basic

#python3 /mnt/lustre/bio/users/fcoutinho/Repos/virathon/Virathon.py --genome_files /mnt/lustre/scratch/fcoutinho/StG_23/Viruses/Genomic_Sequences/IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.fasta --call_prodigal_module True

cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Annotation/

# head -n 502 /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Sequences/IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.faa > short.faa
#python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/AMG_Hunter.py --kegg_db  /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Databases/KOfam/All_KOfam.hmm --kegg_json /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Databases/KOfam/ko00001.json --annotate True --threads 12 --cds /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Annotation/short.faa --skip_uniref True --skip_pfam True

python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/AMG_Hunter.py --kegg_db  /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Databases/KOfam/All_KOfam.hmm --kegg_json /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Databases/KOfam/ko00001.json --annotate True --threads 63 --cds /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Sequences/IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives_Min_Prev_50.faa --skip_uniref True --skip_pfam True

