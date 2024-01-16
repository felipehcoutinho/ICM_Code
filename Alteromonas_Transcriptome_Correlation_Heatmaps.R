args = commandArgs(trailingOnly=TRUE)

#Correlation method can be changed below
cor_met<-args[1]#"pearson"#"spearman"#
correl_cutoff_values<-as.numeric(args[2])#0.7 #This are the absolute value of the correlations that will be kept for downstream analysis. The analysis is repated for each cutoff value
func_level<-as.character(args[3])#KEGG_Best_Subject_Modules
dummy<-args[4]#KEGG_Best_Subject_Modules

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(tibble)

#Only correlations with the genes below will be kept for downtream analysis
mer_gene_ids<-c("iss312_00964","iss312_00965","iss312_00966","iss312_00967","iss312_00968","iss312_00969","iss312_02735","iss312_02736","iss312_02737","iss312_02738","iss312_03666","iss312_03667","iss312_03668","iss312_03669","iss312_03670","iss312_04104","iss312_04105","iss312_04106","iss312_04107","iss312_04108","iss312_04109","iss312_04110")

#These are the subsets of samples that will be used to calculate the correlations, as it is, it does it for all control samples, all hg samples, and all samples together (hg and control)
sample_groups<-list(c("Lag_con_1","Lag_con_2","Lag_con_3","Exp_con_1","Exp_con_2","Exp_con_3","Latexp_con_1","Latexp_con_2","Latexp_con_3","Sta_con_1","Sta_con_2","Sta_con_3"),c("Lag_Hg_1","Lag_Hg_2","Lag_Hg_3","Exp_Hg_1","Exp_Hg_2","Exp_Hg_3","Latexp_Hg_1","Latexp_Hg_2","Latexp_Hg_3","Sta_Hg_1","Sta_Hg_2","Sta_Hg_3"),c("Lag_con_1","Lag_con_2","Lag_con_3","Exp_con_1","Exp_con_2","Exp_con_3","Latexp_con_1","Latexp_con_2","Latexp_con_3","Sta_con_1","Sta_con_2","Sta_con_3","Lag_Hg_1","Lag_Hg_2","Lag_Hg_3","Exp_Hg_1","Exp_Hg_2","Exp_Hg_3","Latexp_Hg_1","Latexp_Hg_2","Latexp_Hg_3","Sta_Hg_1","Sta_Hg_2","Sta_Hg_3"))

names(sample_groups)<-c("Control","Hg","All")

#For each mer gene, count the number of genes that have correlations with them (within the cutoff) according to the KEGG module to which the other gene is assigned. Account for the fact that some genes are assigned to multiple modules, and that some genes are assigned to the same module multiple times. Do the same 3 times: counting all correlations, only positive ones, and only negative ones
selections<-c("positive","negative","all")

#Color pallete of the figures
col_pal<-brewer.pal(11,"Spectral")
col_grad<-rev(colorRampPalette(col_pal)(n=100))

#Transcriptome TMM data
transc_data<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Alteromonas_Transcriptome/tmm_sub_only.csv",sep=";",dec=",",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
rownames(transc_data)<-transc_data[,1]

#AMG hunter annotation
annot_data<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Alteromonas_Transcriptome/CDS_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
colnames(annot_data)[1]<-"Gene_B"

#Add annotation t EDGER Diff Exp data
full_de_data<-rbind()

for (i in c(3:7)) {
	file_name<-paste("/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Alteromonas_Transcriptome/TS",i,".tsv",sep="")
	de_data<-read.table(file=file_name,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
	de_data$Original_File<-paste("TS",i,".tsv",sep="")
	full_de_data<-rbind(full_de_data,de_data)
}

summary(full_de_data)

full_de_data<-as.data.frame(full_de_data)

full_de_data<-merge(full_de_data,annot_data,by.x="Gene.ID",by.y="Gene_B",all.x=TRUE)

full_de_data$Original_File<-as.factor(full_de_data$Original_File)


summary(full_de_data)

write.table(full_de_data,file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Alteromonas_Transcriptome/DE_Genes+Annotation.tsv",sep="\t",append=FALSE,row.names=FALSE,quote=FALSE)

func_level<-"KEGG_Best_Subject_Modules"
#summary(full_de_data)

full_de_data$Functional_Group<-as.factor(gsub("(^\\[)|(\\]$)","",full_de_data[[func_level]],perl=TRUE))
#levels(full_de_data$Functional_Group)[match("",levels(full_de_data$Functional_Group))] <- "Unclassified"
#levels(full_de_data$Functional_Group)<-c(levels(full_de_data$Functional_Group),"Unclassified")
#full_de_data$Functional_Group<-as.factor(full_de_data$Functional_Group)
#summary(full_de_data)
#full_de_data$Functional_Group<-full_de_data$Functional_Group[is.na(full_de_data$Functional_Group)]<-"Unclassified"
#summary(full_de_data)
full_de_data$Functional_Group<-as.factor(full_de_data$Functional_Group)
#summary(full_de_data)


#full_de_data$Functional_Group[grepl("",full_de_data$Functional_Group ,perl=TRUE)]<-"Unclassified"
#full_de_data$Functional_Group[which(full_de_data$Functional_Group == "")]<-"Unclassified"
#full_de_data$Functional_Group<-as.factor(full_de_data$Functional_Group)
#summary(full_de_data)

full_de_data<-full_de_data %>% separate_longer_delim(Functional_Group, delim = "', '")
full_de_data$Functional_Group<-as.factor(gsub("\\'","",full_de_data$Functional_Group,perl=TRUE))


group_file_counts_df<-table(full_de_data[,c("Original_File","Functional_Group")])

write.table(group_file_counts_df[,c(2:ncol(group_file_counts_df))],file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Alteromonas_Transcriptome/TS_3_7_Functional_Group_Counts.tsv",sep="\t",append=FALSE,row.names=TRUE,quote=FALSE,col.names=NA)

#mer Correlations by KEGG functional categories heatmap 
#Do some edits to the functional annotation field to make it easier to parse
annot_data$Functional_Group<-annot_data[[func_level]]

#Gene x Abundance x Time Point x Module boxplots

#split the categ_var column in f_cds_info by ", " and create a new line for element of the split
summary(annot_data)

sub_annot_data<-subset(annot_data,select=c("Gene_B","Functional_Group"))
summary(sub_annot_data)

m_transc_data<-melt(transc_data,id.vars="GeneID",variable.name="Sample",value.name="Abundance")
m_transc_data<-merge(m_transc_data,sub_annot_data,by.x="GeneID",by.y="Gene_B",all.x=TRUE)
m_transc_data$Functional_Group<-as.factor(gsub("(^\\[)|(\\]$)","",m_transc_data$Functional_Group,perl=TRUE))
m_transc_data<-m_transc_data %>% separate_longer_delim(Functional_Group, delim = "', '")
m_transc_data$Functional_Group<-as.factor(gsub("\\'","",m_transc_data$Functional_Group,perl=TRUE))
m_transc_data$Functional_Group<-as.factor(m_transc_data$Functional_Group)
#m_transc_data<-m_transc_data[!is.na(m_transc_data$Functional_Group),]
m_transc_data$Functional_Group[which(m_transc_data$Functional_Group == "")]<-NA


m_transc_data$Phase<-as.factor(gsub("_(.)+$","",m_transc_data$Sample,perl=TRUE))
m_transc_data$Phase<-factor(m_transc_data$Phase,levels=c("Lag","Exp","Latexp","Sta"))

m_transc_data$Treatment<-"Control"
m_transc_data$Treatment[which(grepl("Hg",m_transc_data$Sample,perl=TRUE))]<-"Hg"
m_transc_data$Treatment<-as.factor(gsub("_(.)+$","",m_transc_data$Treatment,perl=TRUE))

summary(m_transc_data)

my_plot<-ggplot(m_transc_data,aes(x=Phase,y=Abundance))+
geom_violin()+
geom_boxplot(alpha=0)+
scale_y_log10()+
geom_signif(comparisons = list(c("Lag", "Exp"),c("Exp", "Latexp"),c("Latexp", "Sta")),map_signif_level = TRUE, test = "wilcox.test",tip_length=0, step_increase=0.1)+#c("Control","Hg") #c("Lag", "Exp"),c("Exp", "Latexp"),c("Latexp", "Sta")
theme_bw()+
theme(axis.text.x = element_text(size=12), strip.text.y.right = element_text(angle = 0, size =9))+ #angle = 45,hjust = 1,
facet_grid(Functional_Group ~  Treatment,scales="free")

outname<-paste("/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Alteromonas_Transcriptome/Alteromonas_Transcriptome_Functional_Group_",func_level,"_Boxplots.pdf",sep="")
ggsave(outname,plot=my_plot,device=cairo_pdf,width=12,height=45,pointsize=8)


####Correlations hetamaps
unique_modules_list<-unlist(strsplit(as.character(annot_data$Functional_Group),",(\\s)+'",perl=TRUE))
unique_modules_list<-gsub("'","",unique_modules_list)
unique_modules_list<-unique(unique_modules_list[!is.na(unique_modules_list)])

make_plots<-TRUE
#Iterate over each one of the sample groups
for (j in c(1:length(sample_groups))) {
	group_name<-names(sample_groups)[j]
	sample_vector<-sample_groups[[j]]

	#Calculate the correlation matrix with a subset of selected samples
	cor_matrix<-cor(t(subset(transc_data,select=sample_vector)),method=cor_met)

	#Melt the correlation matrix
	m_cor_matrix<-melt(cor_matrix)

	colnames(m_cor_matrix)<-c("Gene_A","Gene_B","Correlation_Score")

	#Remove correlations between genes that are no part of the mer operons. This also removes the redundancies.
	sub_m_cor_matrix<-m_cor_matrix[which(m_cor_matrix$Gene_A %in% mer_gene_ids),]

	for (correl_cutoff in correl_cutoff_values) {
		#Remove correlations outside the specified correlation score cutoff
		sub_m_cor_matrix<-sub_m_cor_matrix[which((sub_m_cor_matrix$Correlation_Score >= correl_cutoff) | (sub_m_cor_matrix$Correlation_Score <= -correl_cutoff)),]

		#Add the annotation data to the correlation matrix
		merged_data<-merge(sub_m_cor_matrix,annot_data,by="Gene_B",all.x=TRUE)
		
		out_table_name<-paste("Alteromonas_Transcriptome_",group_name,"_Samples_mer_Genes_Filtered_",correl_cutoff,"_",cor_met,"_Correlation_with_Annotation.tsv",sep="")
		
		#write.table(merged_data,file=out_table_name,sep="\t",append=FALSE,row.names=FALSE,quote=FALSE)

		#iterate over each selection of correlation values
		for (selection in selections) {
			#Keep only the correlactions specific to each selection
			if (selection == "positive") {
				filtered_merged_data<-merged_data[which(merged_data$Correlation_Score >= correl_cutoff),]
			} else if (selection == "negative") {
				filtered_merged_data<-merged_data[which(merged_data$Correlation_Score <= -correl_cutoff),]
			}  else if (selection == "all") {
				filtered_merged_data<-merged_data
			}

			#initiate empty data frame which will hold the counts of the correlations for each gene
			module_counts<-data.frame(matrix(nrow=length(mer_gene_ids),ncol=length(unique_modules_list)))
			colnames(module_counts)<-unique_modules_list
			rownames(module_counts)<-mer_gene_ids
			module_counts[is.na(module_counts)]<-0

			test_table_name<-paste("Test",group_name,correl_cutoff,selection,".tsv",sep="_")
			#write.table(filtered_merged_data,file=test_table_name,sep="\t",append=FALSE,row.names=FALSE,quote=FALSE)
			#These loops iterate over the correlation+annotation table and do the counting
			for (i in c(1:nrow(filtered_merged_data))) {
				Gene_B<-as.character(filtered_merged_data$Gene_A[i])
				#print("Gene_B:")
				#print(Gene_B)
				all_modules<-as.character(filtered_merged_data$Functional_Group[i])
				#print("All modules:")
				#print(all_modules)
				unique_modules<-unique(unlist(strsplit(all_modules,",(\\s)+'",perl=TRUE)))
				unique_modules<-unique(gsub("'","",unique_modules))
				unique_modules<-unique_modules[!is.na(unique_modules)]
				#print("Unique modules:")
				#print(unique_modules)
				
				for (module in unique_modules) {
					module_counts[Gene_B,module]<-module_counts[Gene_B,module] + 1
				}
			}

			#Change absolute counts to percentages (per gene)
			perc_module_counts<-(module_counts/rowSums(module_counts))*100
			perc_module_counts$Gene<-rownames(perc_module_counts)
			perc_module_counts$Total_Correlations<-rowSums(module_counts)
			perc_module_counts[is.na(perc_module_counts)]<-0
			
			m_perc_module_counts<-melt(perc_module_counts,id=c("Gene","Total_Correlations"))
			colnames(m_perc_module_counts)<-c("Gene","Total_Correlations","Module","Percentage_of_Correlations")
			
			out_name<-paste("Alteromonas_Transcriptome_",group_name,"_Samples_",selection,"_",cor_met,"_Cutoff_",correl_cutoff,"_Correlations_Module_Count.tsv",sep="")
			
			write.table(m_perc_module_counts,file=out_name,sep="\t",append=FALSE,row.names=FALSE,quote=FALSE)
			
			m_perc_module_counts$Operon<-gsub("(\\d){3}$","",gsub("iss312_","",m_perc_module_counts$Gene,perl=TRUE),perl=TRUE)
			
			#m_perc_module_counts$Gene<-factor(m_perc_module_counts$Gene,levels=mer_gene_ids)
			
			m_perc_module_counts$Module<-factor(m_perc_module_counts$Module,levels=rev(sort(as.vector(unique(m_perc_module_counts$Module)))))
	
			m_perc_module_counts$Gene<-paste(m_perc_module_counts$Gene," (",m_perc_module_counts$Total_Correlations,")",sep="")

			out_name<-paste("Alteromonas_Transcriptome_",group_name,"_Samples_",selection,"_Cutoff_",correl_cutoff,"_",cor_met,"_Correlations_",func_level,"_Count_Heatmap.pdf",sep="")

			fm_perc_module_counts<-m_perc_module_counts[which(m_perc_module_counts$Percentage_of_Correlations > 0),]
			#Make the plots
			if (make_plots == TRUE) {
				min_height<-7
				max_height<-40
				data_height<-length(unique(fm_perc_module_counts$Module))/3
				ideal_height<-max(min_height,min(data_height,max_height))
				pdf(out_name,width=10,height=ideal_height,pointsize=8)
				
				my_plot<-ggplot(fm_perc_module_counts,aes(fill=Percentage_of_Correlations,y=Module,x=Gene))+geom_tile(colour = "grey50")+theme_bw()+theme(axis.text.x = element_text(angle = 60,hjust = 1,size=8))+scale_fill_gradientn(colours =col_grad,limits = c(0,100))+theme(legend.position="top")+facet_grid(. ~ Operon, scales="free_x",space="free")

				print(my_plot)
				dev.off()
			}

		}
	}
}
