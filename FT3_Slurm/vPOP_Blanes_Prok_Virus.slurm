#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO

#SBATCH --time=00-05:00:00
#SBATCH --job-name=Vpop_Blanes_Prok_Vir
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=4G
#SBATCH --output=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.out
#SBATCH --error=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
#ncbi-blast/
module load cesga/system miniconda3/22.11.1-1
conda activate /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Conda_Envs/basic

cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Blanes_Vir_Nestor/Vpop/

#python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/virathon/Virathon.py --genome_files /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Blanes_Vir_Nestor/Vir_MAGs_Sequences/dsDNAphage_Blanes_virus.fasta --info_output /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Blanes_Vir_Nestor/Vpop/VPop_dsDNAphage_Blanes_virus.tsv --make_pops True --call_prodigal_module True --threads 47

python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/virathon/Virathon.py --genome_files /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Blanes_Vir_Nestor/Vir_MAGs_Sequences/dsDNAphage_Blanes_virus.fasta --gene_file /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Blanes_Vir_Nestor/Vir_MAGs_Sequences/dsDNAphage_Blanes_virus.fna --info_output /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Blanes_Vir_Nestor/Vpop/VPop_dsDNAphage_Blanes_virus.tsv --make_pops True --threads 47