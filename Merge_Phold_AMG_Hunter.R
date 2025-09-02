args <- commandArgs(trailingOnly = TRUE)
phold_file<-args[1]
amg_file<-args[2]
out_file<-args[3]

# phold_file<-"/mnt/smart/scratch/vir/felipe/OctoMicro/Viruses/Annotation/Prodigal-gv/OctoMicro_HQ_Candidates_phold/phold_per_cds_predictions.tsv"
phold_df<-read.table(file=phold_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

summary(phold_df)

phold_df$cds_num<-as.numeric(gsub(".+_CDS_","",phold_df$cds_id,perl=T))

start_cds_df<-aggregate(phold_df$cds_num,by=list(contig_id=phold_df$contig_id),FUN=min)

colnames(start_cds_df)[2]<-"First_CDS_Num"

phold_df<-merge(phold_df,start_cds_df,by="contig_id",all.x=TRUE)

phold_df$cds_amg_num<-(phold_df$cds_num-phold_df$First_CDS_Num)+1

phold_df$Sequence<-as.factor(paste(phold_df$contig_id,phold_df$cds_amg_num,sep="_"))

# summary(phold_df[,c(1:5,35:41)])

phold_df[,c("cds_num","First_CDS_Num","cds_amg_num")]<-NULL

# amg_file<-"/mnt/smart/scratch/vir/felipe/OctoMicro/Viruses/Annotation/Prodigal-gv/CDS_Info_OctoMicro_Virus_2.tsv"
amg_df<-read.table(file=amg_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

# summary(amg_df$Sequence %in% phold_df$Sequence)

phold_df<-merge(phold_df,amg_df,by="Sequence",all.x=TRUE)

write.table(phold_df,file=out_file,sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)