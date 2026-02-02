#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-36:00:00
#SBATCH --job-name=cas_OceanDNAxMLP_abd
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --output=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/scratch/vir/felipe/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load mmseqs2

cd /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/cas_Abundance/

#deep
# while read i; do j=$(echo "$i" | sed 's/.*\///') ; echo "Processing $j"; mmseqs easy-search --threads 23 -s 1.0 -e 0.00001 --min-aln-len 20 --min-seq-id 0.3 --max-seqs 5 $i DB_OceanDNA_All_Species_Rep_MAGs_cas_only_CDS ${j}xDB_OceanDNA_All_Species_Rep_MAGs_cas_only_CDS.m8 mmseqs_tmp ; done < /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/cas_Abundance/Malaspina_Deep_R1_List.txt

#profiles
# while read i; do j=$(echo "$i" | sed 's/.*\///') ; echo "Processing $j"; mmseqs easy-search --threads 23 -s 1.0 -e 0.00001 --min-aln-len 20 --min-seq-id 0.3 --max-seqs 5 $i DB_OceanDNA_All_Species_Rep_MAGs_cas_only_CDS ${j}xDB_OceanDNA_All_Species_Rep_MAGs_cas_only_CDS.m8 mmseqs_tmp ; done < /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/cas_Abundance/Profiles_Malaspina_76_FL_Metagenome_R1_List.txt

while read i; do j=$(echo "$i" | sed 's/.*\///') ; echo "Processing $j"; mmseqs easy-search --threads 23 -s 1.0 -e 0.00001 --min-aln-len 20 --min-seq-id 0.3 --max-seqs 5 $i DB_OceanDNA_All_Species_Rep_MAGs_cas_only_CDS ${j}xDB_OceanDNA_All_Species_Rep_MAGs_cas_only_CDS.m8 mmseqs_tmp ; done < /mnt/smart/scratch/vir/felipe/Databases/OceanDNA/cas_Abundance/Pending_Profiles_Malaspina_76_FL_Metagenome_R1_List.txt

module load python

python3 /mnt/smart/scratch/vir/felipe/Repos/ICM_Code/Parse_M8.py --m8_files *m8 --out OceanDNA_All_CasxMLP_Partial_Profiles+Deep+Kaust_Counts.tsv --max_e 0.001 --min_bit 30 --min_ali 0 --min_id 0.3

#clean sample names on header
sed -i -E s'/(_1\.fastq\S+)|(\.pair1\S+)//g' OceanDNA_All_CasxMLP_Partial_Profiles+Deep+Kaust_Counts.tsv

#Caclulate number of sequences in each file:
# echo "Group,R1,Lines" > MLP_Profiles+Deep+Kaust_R1_Line_Count.csv 
# while read i; do j=$(echo "$i" | sed 's/.*\///') ; lcount=$(zcat $i | wc -l); echo "$i,$j,$lcount" >> MLP_Profiles+Deep+Kaust_R1_Line_Count.csv ; done < /mnt/smart/scratch/vir/felipe/DeepCas/cas1_Abundance/MLP_Profiles+Deep+Kaust_R1_List.txt

#get list of non gz files
# grep -v ".gz" /mnt/smart/scratch/vir/felipe/DeepCas/cas1_Abundance/MLP_Profiles+Deep+Kaust_R1_List.txt > /mnt/smart/scratch/vir/felipe/DeepCas/cas1_Abundance/Nong_GZ_MLP_Profiles+Deep+Kaust_R1_List.txt

# while read i; do j=$(echo "$i" | sed 's/.*\///') ; lcount=$(cat $i | wc -l); echo "$i,$j,$lcount" >> Nong_GZ_MLP_Profiles+Deep+Kaust_R1_Line_Count.csv ; done < /mnt/smart/scratch/vir/felipe/DeepCas/cas1_Abundance/Nong_GZ_MLP_Profiles+Deep+Kaust_R1_List.txt
