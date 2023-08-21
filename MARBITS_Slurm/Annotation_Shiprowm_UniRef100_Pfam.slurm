#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=02-12:00:00
#SBATCH --job-name=Annotation_Shipworm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=150G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load diamond/2.0.7
module load perl
module load python/3.8.5
module load hmmer

#diamond blastp --matrix BLOSUM45 --threads 48 --more-sensitive --db /mnt/lustre/repos/bio/databases/public/UniProtKB/UniProt_Release_2021-04-07/uniref100.dmnd --outfmt 6 --query /mnt/lustre/scratch/fcoutinho/Shipworm/Shipworm_Contigs+5kbp_Viruses_No_Prophage.faa --max-target-seqs 100 --evalue 0.001 --out /mnt/lustre/scratch/fcoutinho/Shipworm/Shipworm_Contigs+5kbp_Viruses_No_Prophage_CDSxUniRef100.blastp

#perl /mnt/lustre/bio/users/fcoutinho/Scripts/Enrich_m8_Table.pl --names /mnt/lustre/repos/bio/databases/public/NCBI/taxonomy/taxonomy_2021-04-06/names.dmp --nodes /mnt/lustre/repos/bio/databases/public/NCBI/taxonomy/taxonomy_2021-04-06/nodes.dmp --acc2taxid_dir /mnt/lustre/repos/bio/databases/public/UniProtKB/UniProt_Release_2021-04-07/ --description /mnt/lustre/bio/users/fcoutinho/Databases/UniRef100/UniRef100_Release_2021_04_07_Id_to_Desc.txt --input /mnt/lustre/scratch/fcoutinho/Shipworm/Shipworm_Contigs+5kbp_Viruses_No_Prophage_CDSxUniRef100.blastp --output /mnt/lustre/scratch/fcoutinho/Shipworm/Enriched_Shipworm_Contigs+5kbp_Viruses_No_Prophage_CDSxUniRef100.blastp

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --genome_files /mnt/lustre/scratch/fcoutinho/Shipworm/Shipworm_Contigs+5kbp_Viruses_No_Prophage.fasta --cds /mnt/lustre/scratch/fcoutinho/Shipworm/Shipworm_Contigs+5kbp_Viruses_No_Prophage.faa --call_hmmer True --threads 48 --hmmer_program hmmsearch --hmmer_db /mnt/lustre/repos/bio/databases/public/pfam/pfam_release_34.0/Pfam-A.hmm