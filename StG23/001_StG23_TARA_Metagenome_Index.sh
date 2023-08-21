#List TARA MG paired-end MetaDNA files
find /mnt/lustre/repos/bio/projects/tara/genoscope/Release_1/ -type f > TARA_Lustre_MG_Files.txt
grep "/Meta-DNA/RawDatas" TARA_Lustre_MG_Files.txt | grep "_1.fastq.gz" > TARA_MG_R1.txt
grep "/Meta-DNA/RawDatas" TARA_Lustre_MG_Files.txt | grep "_2.fastq.gz" > TARA_MG_R2.txt
#Manually edit the data and save it to /mnt/lustre/scratch/fcoutinho/StG_23/Metagenomes/TARA_Lustre_MGs_Info.tsv (non euk metagenomes)

#EUk metagenomes

