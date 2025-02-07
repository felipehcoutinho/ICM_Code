#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-12:00:00
#SBATCH --job-name=VP32_topc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=150G
#SBATCH --output=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
#Move to the relvant directory
# cd /mnt/smart/users/fcoutinho/TARA_Polar/VP32/

#Load the prodigla module
# module load prodigal/2.6.3

# #Run prodigal
# prodigal -q -p meta -a VP32_original_Scaffold_k117_975715.faa -d VP32_original_Scaffold_k117_975715.fna -f gff -i VP32_original_Scaffold_k117_975715.fasta -o VP32_original_Scaffold_k117_975715.gff

# module purge

###Query proteins against Pfam
# #Set up database
# cd /mnt/smart/users/fcoutinho/Databases/
# mkdir PfamA
# cd PfamA

# wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz

# gunzip Pfam-A.hmm.gz

# module load hmmer

# hmmpress Pfam-A.hmm

# #run hmmsearch
# cd /mnt/smart/users/fcoutinho/TARA_Polar/VP32/

# hmmsearch -o VP32xPfamA --tblout VP32xPfamA.tsv -T 50 --noali --cpu 23 /mnt/smart/users/fcoutinho/Databases/PfamA/Pfam-A.hmm VP32_original_Scaffold_k117_975715.faa

###Create IMGvr db to find close relatives
# module load mmseqs2/13-45111

# cd /mnt/smart/users/fcoutinho/Databases/IMGVR/

# mmseqs createdb /mnt/smart/users/fcoutinho/Databases/IMGVR/IMGVR_all_proteins-high_confidence.faa.gz IMGvrDB
# mmseqs createindex IMGvrDB tmp

# rm -fr tmp/

###Query Proteins Against IMGvr db to find close relatives
# cd /mnt/smart/users/fcoutinho/TARA_Polar/VP32/

# mmseqs easy-search VP32_original_Scaffold_k117_975715.faa /mnt/smart/users/fcoutinho/Databases/IMGVR/IMGvrDB VP32_original_Scaffold_k117_975715xIMGvr.tsv tmp --threads 47 --max-seqs 1000 --min-seq-id 0.3 --min-aln-len 30

# rm -fr tmp/

# #Calculate pairwise protein scores
# module load python
# python3 /mnt/smart/users/fcoutinho/Repos/virathon/Virathon.py --pairwise_protein_scores True --genome VP32_original_Scaffold_k117_975715.fasta --cds VP32_original_Scaffold_k117_975715.faa --pps_hits_table VP32_original_Scaffold_k117_975715xIMGvr.tsv

# #Select close relatives with min 50% AAI and 40% shared genes
# awk -F "\t" '{ if (($3 >= 50) && ($5 >= 40)) { print } }'  VP32_original_Scaffold_k117_975715xIMGvr.tsv.Pairwise_Protein_Scores.tsv | cut -f 2 > List_of_IMGvr_Close_Relatives_AAI_50_Shared_Genes_40_Scaffold_k117_975715.txt 

# ###Query proteins against uniprot using mmseqs2
# module load mmseqs2/13-45111

# mmseqs createdb /mnt/smart/users/fcoutinho/Databases/Uniprot_02-10-2024/uniprot_sprot.fasta uniprot_sprotDB
# mmseqs createindex uniprot_sprotDB tmp
# mmseqs easy-search VP32_original_Scaffold_k117_975715.faa uniprot_sprotDB VP32_original_Scaffold_k117_975715xUniprot.tsv tmp --threads 24 --max-seqs 10 --min-seq-id 0.3 --min-aln-len 30 --format-mode 0 --format-output query,theader,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits

# module purge

# #Separate target description from ID
# sed -i 's/ /'$'\t''/'  VP32_original_Scaffold_k117_975715xUniprot.tsv


###VP32 only bowtie abundance
# cd /mnt/smart/users/fcoutinho/TARA_Polar/VP32/Abundance/

# module load python/3.8.5
# module load bowtie2/2.3.4.1
# module load samtools

# python3 /mnt/smart/users/fcoutinho/Repos/virathon/Virathon.py --genome_files /mnt/smart/users/fcoutinho/TARA_Polar/VP32/VP32_original_Scaffold_k117_975715.fasta --abundance_table True --abundance_max_reads 10000000 --abundance_rpkm True --raw_read_table /mnt/smart/users/fcoutinho/TARA_Polar/VP32/Abundance/Short_TOPC_Cel_Fraction_Raw_Reads_Info.tsv --threads 23

#python3 /mnt/smart/users/fcoutinho/Repos/virathon/Virathon.py --bowtiedb All_Genomic --genome_files /mnt/smart/users/fcoutinho/TARA_Polar/VP32/VP32_original_Scaffold_k117_975715.fasta --abundance_table True --abundance_max_reads 10000000 --abundance_rpkm True --raw_read_table /mnt/smart/users/fcoutinho/TARA_Polar/VP32/Abundance/TOPC_Cel_Fraction_Raw_Reads_Info.tsv --threads 23

#python3 /mnt/smart/users/fcoutinho/Repos/virathon/Virathon.py --bowtiedb All_Genomic --genome_files /mnt/smart/users/fcoutinho/TARA_Polar/VP32/VP32_original_Scaffold_k117_975715.fasta --abundance_table True --abundance_max_reads 10000000 --abundance_rpkm True --raw_read_table /mnt/smart/users/fcoutinho/Metadata/TARA/TARA_All_Virome_File_Info+Metadata.tsv --threads 23

###Calculate abundance in additional polar MGs not from TOPC
# module load perl
# perl /mnt/smart/users/fcoutinho/Repos/ICM_Code/download_SRA.pl --list test_sra_List.txt

# while read sra_id; do echo "Processing ${i}"; //mnt/smart/shared/mlCourse/fcoutinho_data/sratoolkit.3.1.1-centos_linux64sratoolkit.2.11.3-ubuntu64/bin/prefetch $sra_id --output-directory . ; while read sra_id; do echo "Processing ${sra_id}"; echo -e "Sample\tGroup\tR1\tR2" > ${sra_id}_File_Info.tsv ; echo -e "${sra_id}\t${sra_id}\t${sra_id}_1.fastq\t${sra_id}_2.fastq" >> ${sra_id}_File_Info.tsv ; python3 /mnt/smart/users/fcoutinho/Repos/virathon/Virathon.py --genome_files /mnt/smart/users/fcoutinho/TARA_Polar/VP32/VP32_original_Scaffold_k117_975715.fasta --bowtiedb All_Genomic --abundance_table True --abundance_max_reads 10000000 --abundance_rpkm True --raw_read_table ${sra_id}_File_Info.tsv --threads 23 ; rm -f $sra_id*fastq ;done < test_sra_List.txt

# cd /mnt/smart/users/fcoutinho/TARA_Polar/VP32/Abundance/

# module load bowtie2/2.3.4.1
# module load samtools

#Not necessary tor edo the databae build step because jorge has done it before
#bowtie2-build --threads 4 /mnt/smart/users/fcoutinho/TARA_Polar/VP32/VP32_original_Scaffold_k117_975715.fasta VP_32_DB

#File SRA_ID_List.txt should contain the list of the SRA IDs compiled by Jorge and nothing else (one ID per line no extra spaces or tabulations, etc) 

# while read sra_id; do echo "Processing ${sra_id}"; /mnt/smart/shared/mlCourse/fcoutinho_data/sratoolkit.3.1.1-centos_linux64/bin/prefetch $sra_id --output-directory . ; /mnt/smart/shared/mlCourse/fcoutinho_data/sratoolkit.3.1.1-centos_linux64/bin/fasterq-dump --split-files $sra_id ; bowtie2 -x VP_32_DB -q -1 ${sra_id}_1.fastq -2 ${sra_id}_2.fastq -S ${sra_id}.sam -k 1 -u 10000000 --sensitive --no-discordant --no-mixed --no-unal --threads 24 ; samtools view -bS ${sra_id}.sam > ${sra_id}.bam ; samtools sort ${sra_id}.bam -o ${sra_id}.sorted.bam ;  samtools index ${sra_id}.sorted.bam ; samtools idxstats ${sra_id}.sorted.bam > ${sra_id}.Counts.tsv ; rm -fr $sra_id ; rm -f $sra_id*fastq ; rm -f ${sra_id}.sam ${sra_id}.bam ; done < SRA_ID_List.txt



