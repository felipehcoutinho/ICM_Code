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

guididata$Station_Alias_2<-as.factor(paste(guididata$station,guididata$Sample_Type,sep="_"))
summary(guididata$Station_Alias_2)

#Prokaryote MGs
mgdata<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/TARA_MGs_Paired_With_Salazar_Data_Info+Metadata.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

#mgdata$Station_Alias<-as.factor(paste(gsub("TARA_","",mgdata$Station_from_Indexing),mgdata$Layer_from_Table_W1,sep="_"))
mgdata$Station_Alias_2<-as.factor(paste(mgdata$Station_from_Indexing,mgdata$Layer_from_Table_W1,sep="_"))
summary(mgdata)

summary(guididata$Station_Alias_2 %in% mgdata$Station_Alias_2)
summary(mgdata$Station_Alias_2 %in% guididata$Station_Alias_2)

fulldata<-merge(mgdata,guididata[,which(colnames(guididata) != "Salinity")],by="Station_Alias_2",all.x=TRUE,suffixes=c("_Original","_Guidi"))

fdata_cnames<-colnames(fulldata)
fdata_cnames<- fdata_cnames[! fdata_cnames %in% c('Group')]
fdata_cnames<-c("Group",fdata_cnames)
fulldata<-fulldata[,fdata_cnames]

summary(fulldata)

write.table(fulldata, file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/TARA_MGs_Paired_With_Salazar_Data_Info+Metadata_Updated_with_Guidi.tsv",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

#Viromes
mgdata<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/Marine_Viral_Communities_Sample_Metadata.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
summary(mgdata)

summary(guididata$Station_Alias %in% mgdata$Sample)

fulldata<-merge(mgdata,guididata[,which(colnames(guididata) != "Salinity")],by.x="Sample",by.y="Station_Alias",all.x=TRUE,suffixes=c("_Original","_Guidi"))

summary(fulldata)

fulldata$Sample_ID<-fulldata$Sample
fulldata$Sample<-NULL

write.table(fulldata, file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Metadata/Marine_Viral_Communities_Sample_Metadata_Updated_with_Guidi.tsv",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

