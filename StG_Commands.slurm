#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=2-00:00:00
#SBATCH --job-name=StG
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=250G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################
#module load cd-hit/4.8.1-longSeqs
#module load mmseqs2/13-45111
#module load bowtie2/2.3.4.1
#module load samtools
#module load python
#module load perl
module load R

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Get_Seqs.py --list GLUVAB3_Passed_VIBRANT.txt --min_length 5000 --input_s /mnt/lustre/scratch/fcoutinho/StG/GLUVAB3_Genomes.fasta --output_s GLUVAB3_VIBRANT_Genomes_Min_5kbp.fasta

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Get_Seqs.py --protein True --list GLUVAB3_Passed_VIBRANT.txt --input_s /mnt/lustre/repos/bio/databases/local/GLUVAB3/GLUVAB3_CDS.faa --output_s GLUVAB3_VIBRANT_Genomes_Min_5kbp_CDS.faa

#mmseqs easy-linclust /mnt/lustre/scratch/fcoutinho/StG/GLUVAB3_VIBRANT_Genomes_Min_5kbp.fasta Cluster_95ID_80Cov_GLUVAB3_VIBRANT_Genomes_Min_5kbp mmseqs_tmp --min-seq-id 0.95 -c 0.8 --cov-mode 3 --threads 48 --max-seq-len 2473871


#bowtie2-build --threads 48 /mnt/lustre/scratch/fcoutinho/StG/Cluster_95ID_80Cov_GLUVAB3_VIBRANT_Genomes_Min_5kbp_rep_seq.fasta StG_GLUVAB3

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --abundance_table True --abundance_rpkm True --bowtiedb /mnt/lustre/scratch/fcoutinho/StG/Bowtie2_DB/StG_GLUVAB3 --metagenomes_dir /mnt/lustre/scratch/fcoutinho/StG/TARA_Viromes/ --metagenomes_extension fastq --genome_files /mnt/lustre/scratch/fcoutinho/StG/Cluster_95ID_80Cov_GLUVAB3_VIBRANT_Genomes_Min_5kbp_rep_seq.fasta --threads 48 --bowtie_mode very-fast --abundance_max_reads 10000000 --index True

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/virathon_v0.1.py --abundance_table True --abundance_rpkm True --bowtiedb /mnt/lustre/scratch/fcoutinho/StG/Bowtie2_DB/StG_GLUVAB3 --metagenomes_dir /mnt/lustre/scratch/fcoutinho/StG/TARA_Viromes/ --metagenomes_extension fastq --genome_files /mnt/lustre/scratch/fcoutinho/StG/Cluster_95ID_80Cov_GLUVAB3_VIBRANT_Genomes_Min_5kbp_rep_seq.fasta --threads 48 --bowtie_mode very-fast --abundance_max_reads 10000000 --parse_only True --abundance_min_count 1000 --index True


#module load hmmer/3.3

#python3 /mnt/lustre/bio/users/fcoutinho/Scripts/Parallel_Hmmer.py --chunks 11 --threads 4 --query /mnt/lustre/scratch/fcoutinho/StG/Min_Count_1k_Cluster_95ID_80Cov_GLUVAB3_VIBRANT_Genomes_Min_5kbp_rep_seq_CDS.faa --database /mnt/lustre/bio/users/fcoutinho/RaFAH/HP_Ranger_Model_3_Filtered_0.9_Valids.hmm --program hmmsearch

#perl /mnt/lustre/bio/users/fcoutinho/Scripts/Make_OG_Table_For_RaFAH.pl --ext faaxDB

#perl /mnt/lustre/bio/users/fcoutinho/Scripts/RaFAH_v0.2.pl --predict --threads 48 --genomexog_table_file_name /mnt/lustre/scratch/fcoutinho/StG/RaFAH/RaFAH_Genome_to_OG_Score_Min_Score_50-Max_evalue_1e-05.tsv

#cd /mnt/lustre/scratch/fcoutinho/StG/Reverse_Ecology_RF/
#Rscript /mnt/lustre/scratch/fcoutinho/StG/Reverse_Ecology_RF/RF_Reverse_Ecology.R

#cd /mnt/lustre/scratch/fcoutinho/StG/Niche_Models_ANN
Rscript /mnt/lustre/scratch/fcoutinho/StG/Niche_Models_ANN/Paralell_ANN_Niche_Models.R

#cd /mnt/lustre/scratch/fcoutinho/StG/Networks/

#Rscript /mnt/lustre/scratch/fcoutinho/StG/Networks/propr_Network.R