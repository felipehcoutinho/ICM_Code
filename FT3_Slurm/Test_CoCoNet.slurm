#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO

#SBATCH --time=00-06:00:00
#SBATCH --job-name=Test_CoCoNet
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=10G
#SBATCH --output=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.out
#SBATCH --error=/mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Job_Logs/%J.err
#SBATCH --mail-type=All
#SBATCH --mail-user=fhernandes@icm.csic.es

##############################################
module load cesga/system miniconda3/22.11.1-1
conda activate /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Conda_Envs/basic

cd /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Software/CoCoNet/

#make test

coconet run --fasta tests/sim_data/assembly.fasta --bam tests/sim_data/*.bam --output output_test --threads 47 --n-train 64 --n-test 8 --batch-size 2 --min-prevalence 0 --test-ratio 0.2 --n-frags 5 --compo-neurons 8 4 --cover-neurons 8 4 --cover-kernel 2 --wsize 2 --wstep 2