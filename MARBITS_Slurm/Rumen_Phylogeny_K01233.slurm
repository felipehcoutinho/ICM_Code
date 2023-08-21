#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-12:00:00
#SBATCH --job-name=K01233_Phylogeny_Rumen
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=100G
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

cp /mnt/lustre/bio/users/fcoutinho/KOfam/profiles/K01233.hmm .
hmmpress K01233.hmm

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --threads 48 --call_hmmer True --hmmer_min_score 150 --hmmer_max_evalue 0.00001 --cds /mnt/lustre/scratch/fcoutinho/Rumen/Assembled/Rumen_Viruses_Renamed.faa --hmmer_db /mnt/lustre/scratch/fcoutinho/Rumen/Assembled/Manual_Curation/K01233.hmm

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --threads 48 --call_hmmer True --hmmer_min_score 150 --hmmer_max_evalue 0.00001 --cds /mnt/lustre/repos/bio/databases/public/UniProtKB/UniProt_Release_2021-04-07/uniref100.fasta --hmmer_db /mnt/lustre/scratch/fcoutinho/Rumen/Assembled/Manual_Curation/K01233.hmm

cut -f 4 /mnt/lustre/scratch/fcoutinho/Rumen/Assembled/Manual_Curation/OG_Pairwise_Score_Table_Rumen_Viruses_RenamedxK01233.hmmsearch.tsv | grep -v "Subject" > List_Rumen_Hits_K01233.txt

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Get_Seqs.py --input_seq /mnt/lustre/scratch/fcoutinho/Rumen/Assembled/Rumen_Viruses_Renamed.faa --output_seq Rumen_Hits_K01233.faa --list List_Rumen_Hits_K01233.txt --min 150 --max 300

cut -f 4 OG_Pairwise_Score_Table_uniref100.xK01233.hmmsearch.tsv | grep -v "Subject" > List_UniRef100_Hits_K01233.txt

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Get_Seqs.py --input_seq /mnt/lustre/repos/bio/databases/public/UniProtKB/UniProt_Release_2021-04-07/uniref100.fasta --output_seq UniRef100_Hits_K01233.faa --list List_UniRef100_Hits_K01233.txt --min 150 --max 300

cd-hit -i UniRef100_Hits_K01233.faa -o NR90_UniRef100_Hits_K01233.faa -c 0.9 -n 5 -M 100000 -T 48

cat NR90_UniRef100_Hits_K01233.faa Rumen_Hits_K01233.faa > Rumen+UniRef100_NR90_Hits_K01233.faa

muscle -maxhours 3 -in Rumen+UniRef100_NR90_Hits_K01233.faa -out Aligned_Rumen+UniRef100_NR90_Hits_K01233.faa 

FastTreeMP -nosupport -out Tree_Rumen+UniRef100_NR90_Hits_K01233.newick Aligned_Rumen+UniRef100_NR90_Hits_K01233.faa