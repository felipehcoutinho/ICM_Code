library(reshape2)

metadata<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/TARA_MGs_Paired_With_Delmont_Data_Info+Metadata.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

prjdata1<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Eukaryotes/Abundance/filereport_read_run_PRJEB6610.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

prjdata1<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Eukaryotes/Abundance/filereport_read_run_PRJEB4352.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

asvdata<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Eukaryotes/Abundance/metapr2_ASVs_selected_abundance_Eukaryota_2023-07-01.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

summary(asvdata)

abdata<-dcast(asvdata,file_code~asv_code,value.var="n_reads_pct",fill=0)

summary(abdata[,1:5])

abdata<-merge(abdata,metadata,by.x="file_code",by.y="ENA_Run_ID",all.x=TRUE)

summary(abdata)