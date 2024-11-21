library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

metadata_file<-"/mnt/smart/users/fcoutinho/Antarctic_Lagoons/Metadata/Sample_Metadata.tsv"
metadata_df<-read.table(file = metadata_file, sep = "\t", header=TRUE, stringsAsFactors=T)
rownames(metadata_df)<-metadata_df$Sample_UID

asv_info_file<-"/mnt/smart/users/fcoutinho/Antarctic_Lagoons/DADA_Output/DADA2_ASVs_Info_Antarctic_Lagoons.tsv"
asv_info_df<-read.table(file = asv_info_file, sep = "\t", header=TRUE, stringsAsFactors=T)

asv_abd_file<-"/mnt/smart/users/fcoutinho/Antarctic_Lagoons/DADA_Output/DADA2_ASVs_Percentage_Abundances_Antarctic_Lagoons.tsv"
asv_perc_abd_df<-read.table(file = asv_abd_file, sep = "\t", header=TRUE, stringsAsFactors=T)
rownames(asv_perc_abd_df)<-asv_perc_abd_df[,1]
colnames(asv_perc_abd_df)[1]<-"MG_ID"

#Change to location of the file in your system
source("/mnt/smart/users/fcoutinho/Repos/ICM_Code/Microbiome_Analysis.R")

#Calculate the sum of abndances at the phylum level
phylum_abd_df<-calc_group_sums(abd_df=asv_perc_abd_df,info_df=asv_info_df,first_group_var="Phylum")

#Add metadata to the phylum abundance df
phylum_plus_meta_df<-merge(phylum_abd_df,metadata_df,by.x="Sample_UID",by.y="MG_ID",all.x=TRUE)

#Caulculate the correlation between phylum abundances and environmental variables
correl_df<-calc_pair_cor(a_vars=colnames(phylum_abd_df)[-1],b_vars=c("Depth","Ti","V","Cr","Mn","Fe","Co", "Ni", "Cu","Zn","Mo" ,"Cd", "Pb","Year"),data_df=phylum_plus_meta_df,method="spearman")

#Filter the correlations with adjusted p-values <= 0.05
f_correl_df<-correl_df[which(correl_df$Adjusted_P_Value_FDR <= 0.05),]

compl_col_grad<-rev(colorRampPalette(brewer.pal(11,"Spectral"))(n=100))

#Make the tile plots
figX<-ggplot(f_correl_df,aes(fill=Correlation_Coefficient,y=Variable_1,x=Variable_2))+
geom_tile(colour="black")+
scale_fill_gradientn(colours = compl_col_grad, limits=c(-0.75,0.75))+
theme_bw()+
labs(x="Environmental Variables",y="Phylum",fill="Spearman Correlation Coefficient")

ggsave("/mnt/smart/users/fcoutinho/Antarctic_Lagoons/DADA_Output/PhylumxEnv_Var_Heatmap_Antarctic_Lagoons.pdf",plot=figX,device=cairo_pdf,width=10,height=7,pointsize=8)