#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=00-06:00:00
#SBATCH --job-name=VirSorter_Polaribacter
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --output=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
###Virsorter
cd /mnt/smart/users/fcoutinho/Polaribacter_Phages/

# module load virsorter/2.2.3
# conda activate vs2

# cat /mnt/smart/users/fcoutinho/Polaribacter_Phages/Genomes/*fasta > All_Polaribacter_Genomes.fasta

# virsorter run --seqfile All_Polaribacter_Genomes.fasta --jobs 23 --prep-for-dramv --rm-tmpdir --min-score 0.5 --include-groups dsDNAphage all

# conda deactivate

# module purge

###CheckV
module load chekcv
source activate

export CHECKVDB="/mnt/lustre/repos/bio/databases/public/checkv/checkv-db-v1.0/"

checkv end_to_end Vibrant_Polaribacter_Virus_Genomes.fasta CheckV_Vibrant_Polaribacter_Virus_Genomes -t 24