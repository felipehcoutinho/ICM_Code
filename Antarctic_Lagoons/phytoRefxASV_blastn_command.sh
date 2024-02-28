module load python blast/2.7.1 
makeblastdb -in PhytoRef_with_taxonomy.fasta -dbtype nucl -title PhytoRef -out PhytoRef

blastn -db PhytoRef -query /mnt/smart/users/fcoutinho/Antarctic_Lagoons/DADA_Output/All_ASV_Sequences.fasta -out Antarctic_Lagoons_ASVsxPhytoRef -outfmt 6 -evalue 0.001 -perc_identity 30 -max_target_seqs 100 -num_threads 8 &