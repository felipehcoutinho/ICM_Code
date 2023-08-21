#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=1-00:00:00
#SBATCH --job-name=Phycodnaviridae_AMG_Hunter
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load diamond/2.0.7
module load python/3.8.5
module load hmmer
module load perl

hmmsearch -o RefSeq_Genomes_PhycodnaviridaexAll_KOs.hmmsearch --noali --cpu 24 /mnt/lustre/bio/users/fcoutinho/KOfam/All_KOs.hmm /mnt/lustre/scratch/fcoutinho/xlopez/Phylogenies/RefSeq_Genomes_Phycodnaviridae.faa

hmmsearch -o RefSeq_Genomes_PhycodnaviridaexPfam-A.hmmsearch --noali --cpu 24 /mnt/lustre/repos/bio/databases/public/pfam/pfam_release_34.0/Pfam-A.hmm /mnt/lustre/scratch/fcoutinho/xlopez/Phylogenies/RefSeq_Genomes_Phycodnaviridae.faa

#diamond blastp --matrix BLOSUM45 --threads 24 --more-sensitive --db /mnt/lustre/repos/bio/databases/public/UniProtKB/UniProt_Release_2021-04-07/uniref100.dmnd --outfmt 6 --query /mnt/lustre/scratch/fcoutinho/xlopez/Phylogenies/RefSeq_Genomes_Phycodnaviridae.faa --max-target-seqs 100 --evalue 0.001 --out RefSeq_Genomes_Phycodnaviridaexuniref100.blastp

#perl /mnt/lustre/bio/users/fcoutinho/Scripts/Enrich_m8_Table.pl --names /mnt/lustre/repos/bio/databases/public/NCBI/taxonomy/taxonomy_2021-04-06/names.dmp --nodes /mnt/lustre/repos/bio/databases/public/NCBI/taxonomy/taxonomy_2021-04-06/nodes.dmp --acc2taxid_dir /mnt/lustre/repos/bio/databases/public/UniProtKB/UniProt_Release_2021-04-07/ --description /mnt/lustre/bio/users/fcoutinho/Databases/UniRef100/UniRef100_Release_2021_04_07_Id_to_Desc.txt --input RefSeq_Genomes_Phycodnaviridaexuniref100.blastp --output Enriched_RefSeq_Genomes_Phycodnaviridaexuniref100.blastp

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/AMG_Hunter.py --annotate True --parse_only True --threads 24 --cds /mnt/lustre/scratch/fcoutinho/xlopez/Phylogenies/RefSeq_Genomes_Phycodnaviridae.faa
