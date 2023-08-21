#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=01-00:00:00
#SBATCH --job-name=SAG_Vir_Analysis
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=250G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
module load python
module load perl
module load R
module load hmmer
module load prodigal
module load mmseqs2
module load diamond

#Call protein encoding genes with Prodigal (Output is in /mnt/lustre/scratch/fcoutinho/xlopez/SAG_Seqs)
python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --genome_files /mnt/lustre/scratch/fcoutinho/xlopez/SAG_Seqs/98NCLDV.fna --call_prodigal True

#Compare SAG viruses to RefSeq NCLDV viruses (Output is in /mnt/lustre/scratch/fcoutinho/xlopez/PPS_Comp)
python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --pairwise True --pps_subject_fasta /mnt/lustre/scratch/fcoutinho/xlopez/Reference_Seqs/RefSeq_NCLDV.faa --threads 24 --cds /mnt/lustre/scratch/fcoutinho/xlopez/SAG_Seqs/98NCLDV.fna.faa --genome_files /mnt/lustre/scratch/fcoutinho/xlopez/SAG_Seqs/98NCLDV.fna

#Compare SAG viruses to all GLUVAB3 viruses (Output is in /mnt/lustre/scratch/fcoutinho/xlopez/PPS_Comp)
python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --pairwise True --pps_subject_db /mnt/lustre/repos/bio/databases/local/GLUVAB3/MMSeqs2_CDS_DB/DB_GLUVAB3_CDS --threads 48 --cds /mnt/lustre/scratch/fcoutinho/xlopez/SAG_Seqs/98NCLDV.fna.faa --genome_files /mnt/lustre/scratch/fcoutinho/xlopez/SAG_Seqs/98NCLDV.fna

#Annotate SAG viral protein sequences (Output is in /mnt/lustre/scratch/fcoutinho/xlopez/Annotation)
python3 /mnt/lustre/bio/users/fcoutinho/Scripts/AMG_Hunter.py --annotate True --threads 48 --cds /mnt/lustre/scratch/fcoutinho/xlopez/SAG_Seqs/98NCLDV.fna.faa --info_cds_output CDS_Info_98NCLDV.tsv --info_genome_output Genome_Info_98NCLDV.tsv

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/AMG_Hunter.py --annotate True --threads 48 --cds /mnt/lustre/scratch/fcoutinho/xlopez/152viroph_proteins.faa --info_cds_output CDS_Info_152viroph.tsv --info_genome_output Genome_Info_152viroph.tsv

#Generate GBK files (Output is in /mnt/lustre/scratch/fcoutinho/xlopez/Genbank_Files)
perl /mnt/lustre/bio/users/fcoutinho/Scripts/Generate_GBK3.pl --skip_unk --explode --annot /mnt/lustre/scratch/fcoutinho/xlopez/Annotation/Enriched_98NCLDVproteinsxuniref100.blastp --protein /mnt/lustre/scratch/fcoutinho/xlopez/SAG_Seqs/98NCLDV.fna.faa --nucleotide /mnt/lustre/scratch/fcoutinho/xlopez/SAG_Seqs/98NCLDV.fna --gff /mnt/lustre/scratch/fcoutinho/xlopez/SAG_Seqs/98NCLDV.fna.gff

#Plot genomic map (Output is in /mnt/lustre/scratch/fcoutinho/xlopez/Annotation)
#(Notice the scape characters in the genbank file name below which ar enecessary bacause using "|" to name files is forbidden
Rscript /mnt/lustre/bio/users/fcoutinho/Scripts/Rscript Make_Genoplot.R /mnt/lustre/scratch/fcoutinho/xlopez/Genbank_Files/AH3097_RM14_A7_1_viruses_197881_cov_1994.395046\|\|full\|\|full.gbk




