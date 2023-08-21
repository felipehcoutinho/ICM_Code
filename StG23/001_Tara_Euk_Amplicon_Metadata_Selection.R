metadata1<-read.table(file="/mnt/lustre/scratch/fcoutinho/StG_23/Metadata/TARA_PANGAEA_Metadata_Pt1.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

metadata2<-read.table(file="/mnt/lustre/scratch/fcoutinho/StG_23/Metadata/TARA_PANGAEA_Metadata_Pt2.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)


enadata1<-read.table(file="/mnt/lustre/scratch/fcoutinho/StG_23/Metadata/ENA_filereport_read_run_PRJEB9737.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

summary(enadata1)

enadata2<-read.table(file="/mnt/lustre/scratch/fcoutinho/StG_23/Metadata/ENA_filereport_read_run_PRJEB6610.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

summary(enadata2)

enadata<-rbind(enadata1,enadata2)

metaprdata<-read.table(file="/mnt/lustre/scratch/fcoutinho/StG_23/Eukaryotes/Abundance/metapr2_ASVs_selected_abundance_Eukaryota_2023-07-01.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

summary(metaprdata)

fulldata<-merge(enadata,metadata1,by.x="sample_accession",by.y="Sample.ID..BioSamples.accession.number......",all.x=TRUE,suffixes=c("_ENA","_PANGAEA_1"))

fulldata<-merge(fulldata,metadata2,by.x="sample_accession",by.y="Sample.ID..BioSamples.accession.number......",all.x=TRUE,,suffixes=c("_ENA","_PANGAEA_2"))

#summary(fulldata)

submetadata<-fulldata[which(fulldata$run_accession %in% metaprdata$file_code),]

#summary(submetadata)

submetadata<-submetadata[,c(1:20,22,25:32,48,50),]

summary(submetadata)

write.table(submetadata,file="/mnt/lustre/scratch/fcoutinho/StG_23/Metadata/Preliminary_Euk_Amplicon_Sample_Metadata.tsv",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)