module load R

Rscript /mnt/lustre/bio/users/fcoutinho/StG_23/Scripts/StG_23_IMGVR_Scaffold_Selection.R

module load python

python3 $scripts/Get_Seqs.py --list /mnt/lustre/scratch/fcoutinho/StG_23/Viruses/Genomic_Sequences/List_IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.txt --input_seq /mnt/lustre/scratch/fcoutinho/StG_23/Viruses/Genomic_Sequences/IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_Genomic_Sequences.fasta --output_seq /mnt/lustre/scratch/fcoutinho/StG_23/Viruses/Genomic_Sequences/IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.fasta  --id_split_sep "|" --id_split_pos 0 &