#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=05-00:00:00
#SBATCH --job-name=AL_Viruses
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=90G
#SBATCH --output=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
# cd /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Viruses/

# module load python

# python3 /mnt/smart/users/fcoutinho/Repos/ICM_Code/Get_Seqs.py --input_sequences /mnt/smart/users/fcoutinho/Antarctic_Lagoons/WGS_Assemblies/Antarctic_Lagoons_All_Filtered_Scaffolds_5kbp.fasta --id_split_sep _ --id_split_pos 0 --fetch_all True --group True

# module load virsorter
# conda activate /mnt/smart/apps/condaEnvs/virsorter2-2.2.4

# # virsorter config --set LOCAL_SCRATCH=$TMPDIR

# virsorter config --set HMMSEARCH_THREADS=12

# # for i in C*.fasta; do smp=$(echo "$i" | sed 's/.fasta//g') ; echo $smp; virsorter run  --db-dir --seqfile $i --exclude-lt2gene --jobs 47 --include-groups dsDNAphage --rm-tmpdir --label $smp --working-dir VS2_Output_${smp}; done

# for i in C*.fasta; do virsorter run --db-dir /mnt/smart/apps/bio/app-databases/virsorter2/db/ --seqfile $i --exclude-lt2gene --jobs 12 --include-groups dsDNAphage,NCLDV --rm-tmpdir --label "${i%.fasta}" --working-dir VS2_Output_"${i%.fasta}"; done

#gives error:
# Traceback (most recent call last):
#   File "/mnt/smart/apps/condaEnvs/virsorter2-2.2.4/lib/python3.10/site-packages/snakemake/__init__.py", line 662, in snakemake
#     success = workflow.execute(
#   File "/mnt/smart/apps/condaEnvs/virsorter2-2.2.4/lib/python3.10/site-packages/snakemake/workflow.py", line 856, in execute
#     dag.create_conda_envs(
#   File "/mnt/smart/apps/condaEnvs/virsorter2-2.2.4/lib/python3.10/site-packages/snakemake/dag.py", line 280, in create_conda_envs
#     env.create(dryrun)
#   File "/mnt/smart/apps/condaEnvs/virsorter2-2.2.4/lib/python3.10/site-packages/snakemake/deployment/conda.py", line 257, in create
#     os.makedirs(env_path, exist_ok=True)
#   File "/mnt/smart/apps/condaEnvs/virsorter2-2.2.4/lib/python3.10/os.py", line 225, in makedirs
#     mkdir(name, mode)
# PermissionError: [Errno 13] Permission denied: '/mnt/smart/apps/bio/app-databases/virsorter2/db/conda_envs/8d0693a2'



# module load python/3.8.5

# cd /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Viruses/

# python3 /mnt/smart/users/fcoutinho/Repos/ICM_Code/Get_Seqs.py --input_sequences /mnt/smart/users/fcoutinho/Antarctic_Lagoons/WGS_Assemblies/Antarctic_Lagoons_All_Filtered_Scaffolds_5kbp.fasta --id_split_sep _ --id_split_pos 0 --fetch_all True --group True

# module load hmmer
# module load prodigal
# module load vibrant

# #Command below nto working in MARBITS because the vibrant databases are stored in lustre
# for i in C*.fasta; do smp=$(echo "$i" | sed 's/.fasta//g') ; echo $smp; python3 /mnt/smart/users/fcoutinho/Repos/virathon/Virathon.py --genome_files $i --call_vibrant True --threads 4 --info_output AL_Virathon_Info_${smp}.tsv ; rm -f All_Genomic.fasta ; done


###Call genes with prodigal-gv
# cd /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Viruses/

# /mnt/smart/users/fcoutinho/Software/prodigal-gv -q -p meta -a AL_VIBRANT_Trimmed_Phage_CDS.faa -d AL_VIBRANT_Trimmed_Phage_Genes.fna -f gff -i AL_VIBRANT_Trimmed_Phage_Genomes.fasta -o AL_VIBRANT_Trimmed_Phage_Genomes.gff

###Host prediction with iPhop
cd /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Viruses/

module load iphop/iphop_env

conda activate /mnt/smart/apps/condaEnvs/iphop_env

mkdir iPhop_Results

iphop predict --num_threads 47 --min_score 75 --fa_file /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Viruses/AL_VIBRANT_Trimmed_Phage_Genomes.fasta --db_dir /mnt/smart/users/fcoutinho/Databases/iPHoP_db/Aug_2023_pub_rw/ --out_dir /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Viruses/iPhop_Results/

mkdir iPhop_Results_C1s6_Seq_182984

iphop predict --num_threads 12 --min_score 75 --fa_file /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Viruses/Genomic_C1s6_Seq_182984.fasta --db_dir /mnt/smart/users/fcoutinho/Databases/iPHoP_db/Aug_2023_pub_rw/ --out_dir /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Viruses/iPhop_Results_C1s6_Seq_182984/

###Annotation
module load python
module load hmmer
module load mmseqs2

cd /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Viruses/

python3 /mnt/smart/users/fcoutinho/Repos/ICM_Code/AMG_Hunter.py --threads 24 --cds /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Viruses/AL_VIBRANT_Trimmed_Phage_CDS.faa --info_cds_output CDS_Info_AL_Virus.tsv --info_genome_output Genome_Info_AL_Virus.tsv --pfam_annotate True --uniref_annotate True --uniref_db /mnt/smart/users/fcoutinho/Databases/Uniprot_02-10-2024/uniprot_sprot.fasta --kegg_annotate True #--parse_only True

# echo "C1s6_Seq_182984" > List_C1s6_Seq_182984.txt
# python3 /mnt/smart/users/fcoutinho/Repos/ICM_Code/Get_Seqs.py --list 

# python3 /mnt/smart/users/fcoutinho/Repos/ICM_Code/Get_Seqs.py --input AL_VIBRANT_Trimmed_Phage_Genomes.fasta --list List_C1s6_Seq_182984.txt  --matched_out Genomic_C1s6_Seq_182984.fasta

# python3 /mnt/smart/users/fcoutinho/Repos/ICM_Code/Get_Seqs.py --input AL_VIBRANT_Trimmed_Phage_CDS.faa --list List_C1s6_Seq_182984.txt  --matched_out CDS_C1s6_Seq_182984.faa --protein True

# python3 /mnt/smart/users/fcoutinho/Repos/ICM_Code/AMG_Hunter.py --threads 6 --cds /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Viruses/CDS_C1s6_Seq_182984.faa --info_cds_output CDS_Info_C1s6_Seq_182984.tsv --info_genome_output Genome_Info_C1s6_Seq_182984.tsv --pfam_annotate True --uniref_annotate True --uniref_db /mnt/smart/users/fcoutinho/Databases/Uniprot_02-10-2024/uniprot_sprot.fasta --kegg_annotate True # --parse_only True


