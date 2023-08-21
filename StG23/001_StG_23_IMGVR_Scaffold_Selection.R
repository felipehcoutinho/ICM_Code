#Read in the scaffold metadata
full_data<-read.table(file="/mnt/lustre/bio/users/fcoutinho/Databases/IMG_VR_4.0/IMGVR_all_Sequence_information-high_confidence.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

#Filter full_data to only contain scaffolds with Estimated.completeness >= 75%, Estimated.completeness <= 100% and Estimated.contamination equal to 0%
filtered_data_1<-full_data[which(full_data$Estimated.completeness >= 75 & full_data$Estimated.completeness <= 100 & full_data$Estimated.contamination == 0),]
dim(filtered_data_1)

library(dplyr)
#Pick a single representative of each vOTU column based on the highest Estimated.completeness value. If there are ties, pick one at random.
filtered_data_2<-filtered_data_1 %>% group_by(vOTU) %>% slice(which.max(Estimated.completeness))
dim(filtered_data_2)
summary(filtered_data_2)
#Write the contents of the UVIG column in filtered_data_2 to the output file
write.table(filtered_data_2$UVIG, file="/mnt/lustre/scratch/fcoutinho/StG_23/Viruses/Genomic_Sequences/List_IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

#read in host info
host_info<-read.table(file="/mnt/lustre/bio/users/fcoutinho/Databases/IMG_VR_4.0/IMGVR_all_Host_information-high_confidence.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

host_info<-host_info[which(host_info$Prediction.method != ""),]

library(tidyr)
#Parse the taxonomy prediction column to extract the domain, phylum, class, order, family, genus and species of the host
host_info$Host.taxonomy.prediction<-gsub("((d__)|(p__)|(c__)|(o__)|(f__)|(g__)|(s__))","",host_info$Host.taxonomy.prediction,perl=TRUE)
host_info<-host_info %>% separate(Host.taxonomy.prediction, c("Host_Domain","Host_Phylum","Host_Class","Host_Order","Host_Family","Host_Genus","Host_Species"),sep=";")
host_info$Host_Phylum[which(host_info$Host_Phylum== "Proteobacteria")]<-host_info$Host_Class[which(host_info$Host_Phylum == "Proteobacteria")]
#summary(host_info)

#Merge filtered_data_2 with host_info based on the UVIG column
filtered_data_2<-merge(filtered_data_2,host_info,by="UVIG",all.x=TRUE)
summary(filtered_data_2)

#Print the contents of filtered_data_2 to the output file
write.table(subset(filtered_data_2,select=c("UVIG","Taxon_oid","Scaffold_oid","Ecosystem.classification","vOTU","Length","Estimated.completeness","Estimated.contamination","Taxonomic.classification","Taxonomic.classification.method","Host.taxonomy.prediction","Host.prediction.method","Host_Domain","Host_Phylum","Host_Class","Host_Order","Host_Family","Host_Genus","Host_Species")), file="/mnt/lustre/scratch/fcoutinho/StG_23/Viruses/Sequence_Info/Full_Info_IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.tsv",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)