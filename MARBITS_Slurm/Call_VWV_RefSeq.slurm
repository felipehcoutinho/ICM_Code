#!/bin/bash
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=01-00:00:00
#SBATCH --job-name=VWV_RefSeq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=100G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python
module load mmseqs2
module load muscle
export PATH="/mnt/lustre/bio/users/fcoutinho/bin:/mnt/lustre/bio/users/fcoutinho/scripts:$PATH"

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/viral_words_v0.3.py --chunks 12 --threads 24 --protein RefSeqVir_Arch+Bact_15-07-21.faa --prefix RefSeq_Round_1 --output RefSeq_Round_1 --hhsuite_query_db /mnt/lustre/bio/users/fcoutinho/Databases/PHROGs/phrogs_hhsuite_db/phrogs

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/viral_words_v0.3.py --onlyparse True --chunks 6 --threads 24 --protein RefSeqVir_Arch+Bact_15-07-21.faa --prefix RefSeq_Round_1 --output RefSeq_Round_1 --hhsuite_query_db /mnt/lustre/bio/users/fcoutinho/Databases/PHROGs/phrogs_hhsuite_db/phrogs --annotation_query_db /mnt/lustre/bio/users/fcoutinho/Databases/PHROGs/PHROGindex.tsv

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Mango_Optimize_VW_Clustering_v0.2.py --map /mnt/lustre/scratch/fcoutinho/VWV/RefSeq_Round_1/RefSeq_Round_1_OG_Map.txt --threads 6 --classification_file /mnt/lustre/scratch/fcoutinho/VWV/RefSeq_Round_1/Optimization/RefSeq_Round_1_OG_Info_No_Unk_Categ.tsv --classification_var Best_Hit_DB_Category --unknown_value NA

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Mango_Optimize_VW_Clustering_v0.2.py --map /mnt/lustre/scratch/fcoutinho/VWV/RefSeq_Round_1/RefSeq_Round_1_OG_Map.txt --threads 6 --classification_file /mnt/lustre/scratch/fcoutinho/VWV/RefSeq_Round_1/Optimization/RefSeq_Round_1_OG_Info_No_Unk_Categ.tsv --classification_var Best_Hit_DB_Annotation

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/viral_words_v0.3.py --onlyparse True --chunks 1 --threads 1 --protein RefSeqVir_Arch+Bact_15-07-21.faa --prefix RefSeq_Round_1 --output RefSeq_Round_1