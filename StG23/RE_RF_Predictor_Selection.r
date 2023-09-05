
# cut -f 1-8,28,34- Reverse_Ecology_Best_RF_Relative_Importance_Viruses_RPKM_All_Seqs+Seq_info.tsv | grep "Ocean\|NPP" > Sub_Reverse_Ecology_Best_RF_Relative_Importance_Viruses_RPKM_All_Seqs+Seq_info.tsv

imp_df<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/RE_RF_Models/Reverse_Ecology_Best_RF_Relative_Importance_Viruses_RPKM_All_Seqs+Seq_info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

anot_df<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Annotation/Genome_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

imp_df<-merge(imp_df,anot_df,by.x="Predictor",by.y="Sequence",all.x=TRUE)

func_df<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Annotation/CDS_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)
func_df$Scaffold<-gsub("_([0-9])+$","",func_df$Sequence,perl=TRUE)
#write.table(imp_df,file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/RE_RF_Models/Sub_Reverse_Ecology_Best_RF_Relative_Importance_Viruses_Percentage_All_Seqs+Seq_info.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)
table(imp_df$Host_Phylum_PHIST)

relevant_cols<-c(1,35,43)
quant_cutoff<-0.95
rep_var<-"NPP.8d.VGPM..mgC.m2.day."#"Mean.Flux.at.150m"
tax_level<-"Host_Phylum_PHIST"
tax_names<-c("Chloroflexota","SAR324","Nitrospinota","Nitrospinota_A","Nitrospirota","Cyanobacteria","Marinisomatota","Thermoplasmatota","Thermoproteota")#
path_name<-"Carbon fixation"#"Photosynthesis"

sub_imp_df<-imp_df[which((imp_df[[rep_var]] >=  quantile(imp_df[[rep_var]] , quant_cutoff)) & (imp_df[[tax_level]] %in% tax_names) & grepl(path_name, imp_df$KEGG_Matched_Pathways)),]

dim(sub_imp_df)

sub_imp_df[,relevant_cols]
summary(sub_imp_df[,relevant_cols])

scaff_ids<-as.vector(unique(sub_imp_df$Predictor))

sub_func_df<-func_df[which((func_df$Scaffold %in% scaff_ids) & grepl(path_name, func_df$KEGG_Best_Subject_Pathways)),]

sub_func_df

func_id<-"sdhC, frdC; succinate dehydrogenase / fumarate reductase, cytochrome b subunit"

unique(func_df[which(grepl(func_id,func_df$KEGG_Best_Subject_Function)),"KEGG_Best_Subject"])


outname<-paste("Reverse_Ecology_Best_RF_Relative_Importance_Viruses_RPKM_All+Annot_Info_Taxa_",paste(tax_names,collapse="|"),"_Response_Variable_",rep_var,"_Quantile_",quant_cutoff,"_",gsub("\\W","_",path_name),".tsv",sep="")
write.table(sub_imp_df,file=outname,sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)