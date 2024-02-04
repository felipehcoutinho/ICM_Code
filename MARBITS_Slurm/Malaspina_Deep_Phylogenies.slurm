#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=02-00:00:00
#SBATCH --job-name=Malaspina_Deep_Phylogenies
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python
module load hmmer
module load muscle
module load fasttree
module load cd-hit

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/AMG_Hunter.py --threads 48 --kegg_phylogeny True --kegg_ko K00287 --ref /mnt/lustre/repos/bio/databases/public/UniProtKB/UniProt_Release_2021-04-07/uniref100.fasta --cds /mnt/lustre/scratch/fcoutinho/Malaspina_Virus/Malaspina_CoAssembly_Viral_Scaffolds_Renamed.faa --derep_id 0.9

muscle -in Length_Filtered_Merged_NR_Matched_K00287.faa -out Aligned_Length_Filtered_Merged_NR_Matched_K00287.faa -maxhours 24 -maxiters 3

FastTreeMP -nosupport -out Tree_Length_Filtered_Merged_NR_Matched_K00287.newick Aligned_Length_Filtered_Merged_NR_Matched_K00287.faa




