#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=15-00:00:00
#SBATCH --job-name=iphop_Blanes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=250G
#SBATCH --output=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
# cd /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/

# module load perl
# perl /mnt/smart/users/fcoutinho/Repos/ICM_Code/split_fasta.pl dsDNAphage_Blanes_virus+viruses_in_bins.fasta 5000

# module load iphop/iphop_env

# conda activate /mnt/smart/apps/condaEnvs/iphop_env

#Canceleld during part_17. Deleted the previously completed part fasta files so can resume with same command

# for i in Part*fasta; do smp=$(echo "$i" | sed 's/.fasta//g') ; echo $smp; mkdir iPhop_Results_${smp} ; iphop predict --min_score 75 --num_threads 23 --fa_file $i --db_dir /mnt/smart/users/fcoutinho/Databases/iPHoP_db/Aug_2023_pub_rw/ --out_dir iPhop_Results_${smp}/; rm -fr iPhop_Results_${smp}/Wdir/ ; done



###Redo for VP representatives only
cd /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/VP_Rep_iPhop/

# grep "True" /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Vpop/VPop_dsDNAphage_Blanes_virus.tsv  | cut -f 1 > List_VP_Rep.txt

# module load python perl

# python3 /mnt/smart/users/fcoutinho/Repos/ICM_Code/Get_Seqs.py --input /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/dsDNAphage_Blanes_virus+viruses_in_bins.fasta --list List_VP_Rep.txt --matched_output VP_Rep_dsDNAphage_Blanes_virus.fasta

# rm -f Unmatched_Sequences.fasta

# perl /mnt/smart/users/fcoutinho/Repos/ICM_Code/split_fasta.pl VP_Rep_dsDNAphage_Blanes_virus.fasta 11000

# module purge

module load iphop/iphop_env

conda activate /mnt/smart/apps/condaEnvs/iphop_env

for i in Part*fasta; do smp=$(echo "$i" | sed 's/.fasta//g') ; echo $smp; mkdir VP_Rep_iPhop_Results_${smp} ; iphop predict --min_score 75 --num_threads 23 --fa_file $i --db_dir /mnt/smart/users/fcoutinho/Databases/iPHoP_db/Aug_2023_pub_rw/ --out_dir VP_Rep_iPhop_Results_${smp}/; rm -fr VP_Rep_iPhop_Results_${smp}/Wdir/ ; done

head -n 1 VP_Rep_iPhop_Results_Part_0/Host_prediction_to_genus_m75.csv > VP_Rep_iPhop_Results_All.csv

cat VP_Rep_iPhop_Results_Part_*/Host_prediction_to_genus_m75.csv | grep -v "List of methods" >> VP_Rep_iPhop_Results_All.csv
