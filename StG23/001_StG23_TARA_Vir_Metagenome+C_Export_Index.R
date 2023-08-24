library(dplyr)
library(tidyr)

mgdata<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/Marine_Viral_Communities_Sample_Metadata.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
summary(mgdata)

guididata<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/Guidi_et_al_2016_TARA_Sample_Metadata+NPP+Carbon_Export.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

guididata$Station_Num<-gsub("^S","",guididata$code.on.figures,perl=TRUE)
guididata$Station_Num<-as.factor(gsub("[A-Z]+","",guididata$code.on.figures,perl=TRUE))
summary(guididata$Station_Num)

guididata$Sample_Type<-gsub("^S","",guididata$code.on.figures,perl=TRUE)
guididata$Sample_Type<-as.factor(gsub("[0-9]+","",guididata$Sample_Type,perl=TRUE))
guididata$Sample_Type<-as.factor(gsub("SUR","SRF",guididata$Sample_Type,perl=TRUE))
summary(guididata$Sample_Type)

guididata$Station_Alias<-as.factor(paste(guididata$Station_Num,guididata$Sample_Type,sep="_"))
summary(guididata$Station_Alias)

summary(guididata$Station_Alias %in% mgdata$Sample)

fulldata<-merge(mgdata,guididata[,which(colnames(guididata) != "Salinity")],by.x="Sample",by.y="Station_Alias",all.x=TRUE,suffixes=c("_Original","_Guidi"))

summary(fulldata)

write.table(fulldata, file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/Marine_Viral_Communities_Sample_Metadata_Updated_with_Guidi.tsv",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)