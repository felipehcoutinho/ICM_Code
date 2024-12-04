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
cd /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/

# module load perl
# perl /mnt/smart/users/fcoutinho/Repos/ICM_Code/split_fasta.pl dsDNAphage_Blanes_virus+viruses_in_bins.fasta 5000

module load iphop/iphop_env

conda activate /mnt/smart/apps/condaEnvs/iphop_env

for i in Part*fasta; do smp=$(echo "$i" | sed 's/.fasta//g') ; echo $smp; mkdir iPhop_Results_${smp} ; iphop predict --min_score 75 --num_threads 23 --fa_file $i --db_dir /mnt/smart/users/fcoutinho/Databases/iPHoP_db/Aug_2023_pub_rw/ --out_dir iPhop_Results_${smp}/; rm -fr iPhop_Results_${smp}/Wdir/ ; done