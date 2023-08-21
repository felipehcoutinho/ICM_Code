library(reshape2)

metadata<-read.table(file="/mnt/lustre/scratch/fcoutinho/StG_23/Metadata/TARA_MGs_Paired_With_Delmont_Data_Info+Metadata.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

asvdata<-read.table(file="/mnt/lustre/scratch/fcoutinho/StG_23/Eukaryotes/Abundance/metapr2_ASVs_selected_abundance_Eukaryota_2023-07-01.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

summary(asvdata)

abdata<-dcast(asvdata,asv_code~file_code,value.var="n_reads",fill=0)