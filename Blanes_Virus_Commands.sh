#Index sequences with Virathon
python3 /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/virathon/Virathon.py --genome /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Blanes_Vir_Nestor/Vir_MAGs_Sequences/dsDNAphage_Blanes_virus.fasta --info_output /mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Blanes_Vir_Nestor/Vir_MAGs_Sequences/Seq_Info_dsDNAphage_Blanes_virus.tsv

#Run iPHOP
sbatch /mnt/netapp1/Store_CSIC/home/csic/eyg/fhc/repos/ICM_Code/FT3_Slurm/iPHop_Blanes_Prok_Viruses.slurm

#Explode protein fasta file for metabolic
cd /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Sequences/Viruses/

#mkdir Split_CDS
cd Split_CDS
module load python

python3 /mnt/smart/users/fcoutinho/ICM_Code/Get_Seqs.py --fetch_all True --group True --id_split_sep DUMMY --protein True --id_split_pos 0 --input /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Sequences/Viruses/dsDNAphage_Blanes_virus.faa

rm -fr  /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Metabolic_Output/Metabolic_Outputs_Batch_*/

cd /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Metabolic_Output
python3 /mnt/smart/users/fcoutinho/ICM_Code/Split_Files_In_Batches.py --batch_size 12000 --input /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Sequences/Viruses/Split_CDS/ --extension faa &

for i in Batch_*/*fasta ; do j=$(echo "$i" | sed 's/.fasta/.faa/'); mv $i $j ; done ;

###PHIST Clean MAGs
module load python
cd /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/PHIST_Predictions/
mkdir Split_Contig_and_Bin_Viruses
cd Split_Contig_and_Bin_Viruses
python3 /mnt/smart/users/fcoutinho/ICM_Code/Explode_Fasta.py --input /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Sequences/Prok_MAGs/final-viral-combined.fa &
python3 /mnt/smart/users/fcoutinho/ICM_Code/Explode_Fasta.py --input /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/dsDNAphage_Blanes_virus.fasta &

cat /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/Sequences/Prok_MAGs/final-viral-combined.fa /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/dsDNAphage_Blanes_virus.fasta > /mnt/smart/users/fcoutinho/Blanes_Vir_Nestor/dsDNAphage_Blanes_virus+viruses_in_bins.fasta

