#Index sequences with Virathon
python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/virathon/Virathon.py --genome /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Blanes_Vir_Nestor/Vir_MAGs_Sequences/dsDNAphage_Blanes_virus.fasta --info_output /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Blanes_Vir_Nestor/Vir_MAGs_Sequences/Seq_Info_dsDNAphage_Blanes_virus.tsv

#Run iPHOP
sbatch /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/FT3_Slurm/iPHop_Blanes_Prok_Viruses.slurm

#Split sequences per MAG assignment
conda activate /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Conda_Envs/basic
cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Blanes_Vir_Nestor/Prok_MAGs_Sequences/
python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/Get_Seqs.py --input /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Blanes_Vir_Nestor/Prok_MAGs_Sequences/bbmo_concatenated.fa --id_split_sep ";" --id_split_pos 0 --fetch_all True --group True

#Run PHIST
sbatch /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/FT3_Slurm/PHIST_Blanes_Prok_Virus.slurm