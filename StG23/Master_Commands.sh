#Subsetting viral protein file to only contain those from genomic sequences with prevalenc eof 50 or more
cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Abundance/
cut -f 1 StG23_Viruses_Min_Prev_50_Raw_Abundance.tsv | grep "IMG" > List_StG23_Viruses_Min_Prev_50_Raw_Abundance.txt

module load cesga/system miniconda3/22.11.1-1
conda activate /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/Conda_Envs/basic

cd /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Sequences/

python3 $STORE/repos/ICM_Code/Get_Seqs.py --protein True --list List_StG23_Viruses_Min_Prev_50_Raw_Abundance.txt --input_seq ../Sequences/IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.faa --matched_out IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives_Min_Prev_50.faa &

#Command below because the slurm script was queued and I had to change the name of the file
mv IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.faa Backup_IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.faa
mv IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives_Min_Prev_50.faa IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.faa