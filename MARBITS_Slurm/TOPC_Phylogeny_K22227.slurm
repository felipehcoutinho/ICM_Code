#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-12:00:00
#SBATCH --job-name=Phylogeny_TOPC_K22227
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

#cp /mnt/lustre/bio/users/fcoutinho/KOfam/profiles/K22227.hmm .
#hmmpress K22227.hmm

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --threads 48 --call_hmmer True --hmmer_min_score 50 --hmmer_max_evalue 0.00001 --cds /mnt/lustre/scratch/fcoutinho/TARA_Polar/Renamed_Viruses_Genomes_Filtered_All_Tara_Polar_Coassemblies_Contigs.faa --hmmer_db /mnt/lustre/scratch/fcoutinho/TARA_Polar/Phylogenies/ahbD/K22227.hmm

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --threads 48 --call_hmmer True --hmmer_min_score 50 --hmmer_max_evalue 0.00001 --cds /mnt/lustre/repos/bio/databases/public/UniProtKB/UniProt_Release_2021-04-07/uniref100.fasta --hmmer_db /mnt/lustre/scratch/fcoutinho/TARA_Polar/Phylogenies/ahbD/K22227.hmm

cut -f 4 /mnt/lustre/scratch/fcoutinho/TARA_Polar/Phylogenies/ahbD/OG_Pairwise_Score_Table_Renamed_Viruses_Genomes_Filtered_All_Tara_Polar_Coassemblies_ContigsxK22227.hmmsearch.tsv | grep -v "Subject" > List_TOPC_Hits_K22227.txt

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Get_Seqs.py --input_seq /mnt/lustre/scratch/fcoutinho/TARA_Polar/Renamed_Viruses_Genomes_Filtered_All_Tara_Polar_Coassemblies_Contigs.faa --output_seq TOPC_Hits_K22227.faa --list List_TOPC_Hits_K22227.txt --min 300 --max 400

cut -f 4 OG_Pairwise_Score_Table_uniref100.xK22227.hmmsearch.tsv | grep -v "Subject" > List_UniRef100_Hits_K22227.txt

python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Get_Seqs.py --input_seq /mnt/lustre/repos/bio/databases/public/UniProtKB/UniProt_Release_2021-04-07/uniref100.fasta --output_seq UniRef100_Hits_K22227.faa --list List_UniRef100_Hits_K22227.txt --min 300 --max 400

cd-hit -i UniRef100_Hits_K22227.faa -o NR90_UniRef100_Hits_K22227.faa -c 0.9 -n 5 -M 100000 -T 48

cat NR90_UniRef100_Hits_K22227.faa TOPC_Hits_K22227.faa > TOPC+UniRef100_NR90_Hits_K22227.faa

muscle -maxhours 3 -in TOPC+UniRef100_NR90_Hits_K22227.faa -out Aligned_TOPC+UniRef100_NR90_Hits_K22227.faa 

FastTreeMP -nosupport -out Tree_TOPC+UniRef100_NR90_Hits_K22227.newick Aligned_TOPC+UniRef100_NR90_Hits_K22227.faa