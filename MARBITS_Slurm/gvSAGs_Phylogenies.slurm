#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-12:00:00
#SBATCH --job-name=Phylo_MV52
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --output=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
cd /mnt/smart/scratch/vir/felipe/gvSAGs/Phylogenies/

module load python
module load hmmer
module load muscle
module load fasttree
#module load cd-hit

#make the fasta file with all protein sequences
python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/Get_Seqs.py --fetch_all True --input /mnt/smart/scratch/vir/felipe/gvSAGs/gvSAGs.faa /mnt/smart/scratch/vir/felipe/Databases/GVDB/GVDB_genomes/GVDB_prot/*faa --matched_out gvSAG_PhyloDB.faa

#PolB
python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/AMG_Hunter.py --threads 12 --kegg_phylogeny True --kegg_ko PolB --cds /mnt/smart/scratch/vir/felipe/gvSAGs/Phylogenies/gvSAG_PhyloDB.faa --derep_id 0 --kegg_min_score 100 --ko_hmm_dir /mnt/smart/scratch/vir/felipe/Databases/GVDB/hmm/

