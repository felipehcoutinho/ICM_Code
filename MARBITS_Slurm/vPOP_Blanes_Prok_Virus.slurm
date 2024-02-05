#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=03-00:00:00
#SBATCH --job-name=Vpop_Blanes_Prok_Vir
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=250G
#SBATCH --output=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/smart/users/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
#ncbi-blast/
#module load cesga/system miniconda3/22.11.1-1
module load python blast/2.7.1 prodigal
#conda activate /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Conda_Envs/basic

cd /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Vpop/

#python3 /mnt/cold02/bio/fcoutinho/virathon/Virathon.py --genome_files /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/test.fasta --call_prodigal True  --info_output /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Vpop/VPop_Test.tsv --make_pops True --threads 14 & #[1] 13148

python3 /mnt/cold02/bio/fcoutinho/virathon/Virathon.py --genome_files /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/dsDNAphage_Blanes_virus.fasta --gene_file /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/dsDNAphage_Blanes_virus.fna --info_output /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Vpop/VPop_dsDNAphage_Blanes_virus.tsv --make_pops True --threads 47