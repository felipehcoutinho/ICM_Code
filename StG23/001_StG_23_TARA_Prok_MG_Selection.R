library(dplyr)
library(tidyr)

mgdata<-read.table(file="/mnt/lustre/scratch/fcoutinho/StG_23/Metagenomes/TARA_Lustre_MGs_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
mgdata$Zone<-as.factor(gsub("SUR","SRF",mgdata$Zone,perl=TRUE))

#summary(mgdata)
#Read the salazar data tables that contain PANGEAE IDs linking samples to ENA run IDs
saldata1<-read.table(file="/mnt/lustre/scratch/fcoutinho/StG_23/Metadata/Salazar_2019_Table_W1.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
saldata1<-saldata1 %>% separate_longer_delim(ENA_Run_ID, delim = "|")
saldata1<-saldata1 %>% separate_longer_delim(BioSamples_ID, delim = ",")
saldata1$ENA_Run_ID<-as.factor(saldata1$ENA_Run_ID)
saldata1$BioSamples_ID<-as.factor(saldata1$BioSamples_ID)
summary(saldata1)

saldata2<-read.table(file="/mnt/lustre/scratch/fcoutinho/StG_23/Metadata/Salazar_2019_Table_W4.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
summary(saldata2)

saldata<-merge(saldata1,saldata2,by="PANGAEA.sample.id",all=TRUE,suffixes=c("_from_Table_W1","_from_Table_W4"))
summary(saldata)

fulldata<-merge(mgdata,saldata,by.x=c("Sample"),by.y=c("ENA_Run_ID"),all.x=TRUE,suffixes=c("_from_Indexing","_from_Salara_19_W1_W4_Tables"))
summary(fulldata)

#Only keep sampels listed in table W4 of Salazar et al. 2019 since those are the ones with metadata. Also, only keep samples from the 0.2-1.6 and 0.22-3 size fractions
subdata<-fulldata[which((fulldata$Size_Fraction == "0.2-1.6") | (fulldata$Size_Fraction == "0.22-3")),]
subdata<-subdata[!is.na(subdata$polar),]
#For abundance purposes samples, should be grouped by ther unique PANGAEA IDs
subdata$Group<-subdata$BioSamples_ID 
subdata<-subdata[!duplicated(subdata$PANGAEA.sample.id),]

summary(subdata)

write.table(subdata, file="/mnt/lustre/scratch/fcoutinho/StG_23/Metadata/TARA_MGs_Paired_With_Salazar_Data_Info+Metadata.tsv",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

saldata2<-merge(saldata2,saldata1[,c("PANGAEA.sample.id","BioSamples_ID")],by="PANGAEA.sample.id",all.x=TRUE,suffixes=c("_from_Table_W4","_from_Table_W1"))

#remove duplicated PANGAEA.sample.id rows in saldata2
saldata2<-saldata2[!duplicated(saldata2$PANGAEA.sample.id),]

summary(saldata2)

write.table(saldata2[,c("BioSamples_ID","PANGAEA.sample.id","Station.label","Layer","polar","lower.size.fraction","upper.size.fraction","Event.date","Latitude","Longitude","Depth.nominal","Ocean.region","Temperature","Oxygen","ChlorophyllA","Carbon.total","Salinity","Gradient.Surface.temp.SST.","Fluorescence","CO3","HCO3","Density","PO4","PAR.PC","NO3","Si","Alkalinity.total","Ammonium.5m","Depth.Mixed.Layer","Lyapunov","NO2","Depth.Min.O2","NO2NO3","Nitracline","Brunt.Väisälä","Iron.5m","Depth.Max.O2","Okubo.Weiss","Residence.time")], file="/mnt/lustre/scratch/fcoutinho/StG_23/Metadata/TARA_Salazar_19_MGs_Metadata.tsv",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)