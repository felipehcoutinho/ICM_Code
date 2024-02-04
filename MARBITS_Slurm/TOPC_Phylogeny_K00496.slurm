#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-12:00:00
#SBATCH --job-name=Phylogeny_TOPC_K00496
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=10G
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

cp /mnt/lustre/bio/users/fcoutinho/KOfam/profiles/K00496.hmm .
hmmpress K00496.hmm

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --threads 48 --call_hmmer True --hmmer_min_score 30 --hmmer_max_evalue 0.00001 --cds /mnt/lustre/scratch/fcoutinho/TARA_Polar/Renamed_Viruses_Genomes_Filtered_All_Tara_Polar_Coassemblies_Contigs.faa --hmmer_db /mnt/lustre/scratch/fcoutinho/TARA_Polar/Phylogenies/K00496/K00496.hmm

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --threads 48 --call_hmmer True --hmmer_min_score 30 --hmmer_max_evalue 0.00001 --cds /mnt/lustre/repos/bio/databases/public/UniProtKB/UniProt_Release_2021-04-07/uniref100.fasta --hmmer_db /mnt/lustre/scratch/fcoutinho/TARA_Polar/Phylogenies/K00496/K00496.hmm

cut -f 4 /mnt/lustre/scratch/fcoutinho/TARA_Polar/Phylogenies/K00496/OG_Pairwise_Score_Table_Renamed_Viruses_Genomes_Filtered_All_Tara_Polar_Coassemblies_ContigsxK00496.hmmsearch.tsv | grep -v "Subject" > List_TOPC_Hits_K00496.txt

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Get_Seqs.py --input_seq /mnt/lustre/scratch/fcoutinho/TARA_Polar/Renamed_Viruses_Genomes_Filtered_All_Tara_Polar_Coassemblies_Contigs.faa --output_seq TOPC_Hits_K00496.faa --list List_TOPC_Hits_K00496.txt --min 250 --max 450

cut -f 4 OG_Pairwise_Score_Table_uniref100.xK00496.hmmsearch.tsv | grep -v "Subject" > List_UniRef100_Hits_K00496.txt

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Get_Seqs.py --input_seq /mnt/lustre/repos/bio/databases/public/UniProtKB/UniProt_Release_2021-04-07/uniref100.fasta --output_seq UniRef100_Hits_K00496.faa --list List_UniRef100_Hits_K00496.txt --min 250 --max 450

cd-hit -i UniRef100_Hits_K00496.faa -o NR90_UniRef100_Hits_K00496.faa -c 0.9 -n 5 -M 100000 -T 48

cat NR90_UniRef100_Hits_K00496.faa TOPC_Hits_K00496.faa > TOPC+UniRef100_NR90_Hits_K00496.faa

muscle -maxhours 3 -in TOPC+UniRef100_NR90_Hits_K00496.faa -out Aligned_TOPC+UniRef100_NR90_Hits_K00496.faa 

FastTreeMP -nosupport -out Tree_TOPC+UniRef100_NR90_Hits_K00496.newick Aligned_TOPC+UniRef100_NR90_Hits_K00496.faa