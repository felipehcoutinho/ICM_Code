#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=02-00:00:00
#SBATCH --job-name=AL_Viruses
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=150G
#SBATCH --output=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load python/3.8.5

cd /mnt/smart/users/fcoutinho/Antarctic_Lagoons/Viruses/

# python3 /mnt/smart/users/fcoutinho/Repos/ICM_Code/Get_Seqs.py --input_sequences /mnt/smart/users/fcoutinho/Antarctic_Lagoons/WGS_Assemblies/Antarctic_Lagoons_All_Filtered_Scaffolds_5kbp.fasta --id_split_sep _ --id_split_pos 0 --fetch_all True --group True

module load hmmer
module load prodigal
module load vibrant

#Command below nto working in MARBITS because the vibrant databases are stored in lustre
for i in C*.fasta; do smp=$(echo "$i" | sed 's/.fasta//g') ; echo $smp; python3 /mnt/smart/users/fcoutinho/Repos/virathon/Virathon.py --genome_files $i --call_vibrant True --threads 4 --info_output AL_Virathon_Info_${smp}.tsv ; rm -f All_Genomic.fasta ; done