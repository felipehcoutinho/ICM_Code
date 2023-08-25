library(reshape2)

args = commandArgs(trailingOnly=TRUE)
seqinfo_file<-args[1]
importance_file<-args[2]

seqinfo_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Info_Sequences/Full_Info_IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.tsv"
importance_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/RE_RF_Models/Reverse_Ecology_Passed_RF_Importance.tsv"

seqinfo_df<-read.table(file=seqinfo_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

#seqinfo_df$UID<-as.factor(paste(seqinfo_df$UVIG,seqinfo_df$Taxon_oid,seqinfo_df$Scaffold_oid,sep="|"))
#write.table(seqinfo_df[,c("UID","UVIG","Taxon_oid","Scaffold_oid","Ecosystem.classification","vOTU","Length","Estimated.completeness","Estimated.contamination","Taxonomic.classification","Taxonomic.classification.method","Host.taxonomy.prediction","Host.prediction.method","Host_Domain_IMGVR","Host_Phylum_IMGVR","Host_Class_IMGVR","Host_Order_IMGVR","Host_Family_IMGVR","Host_Genus_IMGVR","Host_Species_IMGVR","MAG","phage_fasta_file","host_fasta_file","X.common.kmers","pvalue","adj.pvalue","Host_Domain_PHIST","Host_Phylum_PHIST","Host_Class_PHIST","Host_Order_PHIST","Host_Family_PHIST","Host_Genus_PHIST","Host_Species_PHIST")],file=seqinfo_file,sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

summary(seqinfo_df)

imp_df<-read.table(file=importance_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)
summary(imp_df)

norm_imp_df<-dcast(imp_df,Predictor~Response,value.var="Normalized_Importance",fill=NA)
summary(norm_imp_df)

norm_imp_df<-merge(norm_imp_df,seqinfo_df,by.x="Predictor",by.y="UID",all.x=TRUE)
summary(norm_imp_df)