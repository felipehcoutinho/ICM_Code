library(reshape2)
library(ggplot2)
library(RColorBrewer)

phylum_coloring<-c(brewer.pal(9,"Set1"),brewer.pal(4,"Pastel1"),brewer.pal(8,"Dark2")[8],brewer.pal(5,"Set2")[5],brewer.pal(8,"Accent"),brewer.pal(6,"Dark2"))
names(phylum_coloring)<-c("Actinobacteriota","Alphaproteobacteria","Cyanobacteria","Crenarchaeota","Bacteroidota","Nanoarchaeota","Verrucomicrobiota","Marinisomatota","Myxococcota","Patescibacteria","Gammaproteobacteria","Acidobacteriota","UBP7","Poribacteria","Eremiobacterota","Margulisbacteria","Gemmatimonadota","Desulfobacterota","Nitrospinota","Desulfobacterota_D","Bdellovibrionota","Hydrothermarchaeota","Chlamydiota","SAR324","Halobacteriota","Chloroflexota","Thermoplasmatota","Planctomycetota","Thermoproteota")

args = commandArgs(trailingOnly=TRUE)
seqinfo_file<-args[1]
importance_file<-args[2]
imp_info_file<-args[3]
group_var<-args[4]
min_nomr_imp<-as.numeric(args[5])

group_var<-"Host_Phylum_PHIST"
min_rel_imp<-0.01
min_rel_imp_sum<-1
seqinfo_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/Info_Sequences/Full_Info_IMGVR_Scaffolds_Comp_75-100_Max_Conta_0_vOTU_Representatives.tsv"
importance_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/RE_RF_Models/Reverse_Ecology_Best_RF_Importance.tsv"

seqinfo_df<-read.table(file=seqinfo_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)
summary(seqinfo_df)
#seqinfo_df$UID<-as.factor(paste(seqinfo_df$UVIG,seqinfo_df$Taxon_oid,seqinfo_df$Scaffold_oid,sep="|"))
#write.table(seqinfo_df[,c("UID","UVIG","Taxon_oid","Scaffold_oid","Ecosystem.classification","vOTU","Length","Estimated.completeness","Estimated.contamination","Taxonomic.classification","Taxonomic.classification.method","Host.taxonomy.prediction","Host.prediction.method","Host_Domain_IMGVR","Host_Phylum_IMGVR","Host_Class_IMGVR","Host_Order_IMGVR","Host_Family_IMGVR","Host_Genus_IMGVR","Host_Species_IMGVR","MAG","phage_fasta_file","host_fasta_file","X.common.kmers","pvalue","adj.pvalue","Host_Domain_PHIST","Host_Phylum_PHIST","Host_Class_PHIST","Host_Order_PHIST","Host_Family_PHIST","Host_Genus_PHIST","Host_Species_PHIST")],file=seqinfo_file,sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)
#seqinfo_df<-read.table(file=seqinfo_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

imp_df<-read.table(file=importance_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)
summary(imp_df)

imp_df<-merge(imp_df,seqinfo_df[,c("UID",group_var)],by.x="Predictor",by.y="UID",all.x=TRUE)

#Plot Importance Sums barplots
sum_imp_df<-aggregate(imp_df$Relative_Importance,by=list(Response=imp_df$Response,Group=imp_df[[group_var]]),FUN=sum)
colnames(sum_imp_df)[3]<-"Relative_Importance_Sum"
summary(sum_imp_df)

f_sum_imp_df<-sum_imp_df[which(sum_imp_df$Relative_Importance_Sum >= min_rel_imp_sum),]

figX<-ggplot(data=f_sum_imp_df,aes(x=Response, y=Relative_Importance_Sum, fill=Group))+
geom_bar(position="stack",stat="identity",colour="black",alpha=0.9)+ #,group=f_sum_imp_df[[group_var]]
scale_fill_manual(name="Predicted Host Taxon",values=phylum_coloring)+
theme_bw()+
theme(legend.position="top", axis.text.x = element_text(angle = 45,hjust = 1,size=5))+ #, strip.text.y.right = element_text(angle = 0, size =9),legend.title=element_text(size=9)
#facet_grid(. ~ Response,scales="free_x")+
#guides(col = guide_legend(nrow = 1))+
xlab(NULL)+
ylab("Sum of Relative Relative Importance by Predictor Group")

ggsave("/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/RE_RF_Models/Reverse_Ecology_Best_RF_Relative_Importance_Sum_by_Phylum_Barplots.pdf",plot=figX,device=cairo_pdf,width=10,height=7,pointsize=8)

#Build dataframe of relative importances by scaffold UID
norm_imp_df<-dcast(imp_df,Predictor~Response,value.var="Relative_Importance",fill=NA)
summary(norm_imp_df)

norm_imp_df<-merge(norm_imp_df,seqinfo_df,by.x="Predictor",by.y="UID",all.x=TRUE)
summary(norm_imp_df)

write.table(norm_imp_df,file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/RE_RF_Models/Reverse_Ecology_Best_RF_Relative_Importance+Seq_info.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

#Plot Filtered Importance boxplots
f_imp_df<-imp_df[which(imp_df$Relative_Importance >= min_rel_imp),]

#Replace NA is f_imp_df[[group_var]] with "Unknown"
clean_names<-as.vector(f_imp_df[[group_var]])
clean_names[is.na(clean_names)]<-"Unknown"
f_imp_df[[group_var]]<-as.factor(clean_names)
summary(f_imp_df)

figX<-ggplot()+
geom_boxplot(data=f_imp_df,aes(x=f_imp_df[[group_var]], y=log10(Relative_Importance),fill=f_imp_df[[group_var]],group=f_imp_df[[group_var]]),size=0.3,alpha=0.9)+
theme_bw()+
scale_fill_manual(name="Predicted Host Taxon",values=phylum_coloring)+
theme(legend.position="top", axis.text.x = element_text(angle = 90,hjust = 1,size=7), strip.text.y.right = element_text(angle = 0, size =9),legend.title=element_text(size=9))+
facet_grid(. ~ Response,scales="free_x")+
guides(col = guide_legend(nrow = 1))+
xlab(NULL)+
ylab("Log10(Relative Predictor Importance)")

ggsave("/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/StG23/Viruses/RE_RF_Models/Reverse_Ecology_Best_RF_Relative_Importance_by_Phylum_Boxplots.pdf",plot=figX,device=cairo_pdf,width=22,height=6,pointsize=8)