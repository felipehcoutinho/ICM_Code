library(dplyr)
library(tidyr)

mgdata<-read.table(file="/mnt/lustre/scratch/fcoutinho/StG_23/Metagenomes/TARA_Lustre_MGs_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

mgdata$Station_Alias<-as.factor(paste(mgdata$Station,mgdata$Zone,sep=""))

deldata1<-read.table(file="/mnt/lustre/scratch/fcoutinho/StG_23/Metadata/Delmont_2022_Sample_to_ERR_ID.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
deldata1<-deldata1 %>% separate_longer_delim(All_ENA_Read_IDs, delim = ",")
colnames(deldata1)[2]<-"ENA_Run_ID"
deldata1$ENA_Run_ID<-as.factor(deldata1$ENA_Run_ID)
summary(deldata1)

deldata2<-read.table(file="/mnt/lustre/scratch/fcoutinho/StG_23/Metadata/Delmont_2022_ERR_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

summary(deldata2)

fulldata<-merge(deldata1,deldata2,by="Sample",all.x=TRUE)

fulldata$Station_Alias<-as.factor(paste(fulldata$Station,fulldata$Depth,sep=""))

summary(fulldata)

deldata3<-read.table(file="/mnt/lustre/scratch/fcoutinho/StG_23/Metadata/Delmont_2022_Station_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

summary(deldata3)

fulldata<-merge(fulldata,deldata3,by.x="Station_Alias",by.y="Station",all.x=TRUE)

summary(fulldata)

write.table(fulldata, file="/mnt/lustre/scratch/fcoutinho/StG_23/Metadata/TARA_MGs_Paired_With_Delmont_Data_Info+Metadata.tsv",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)