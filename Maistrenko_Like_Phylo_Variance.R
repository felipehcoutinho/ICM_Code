library("vegan")
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr, lib.loc="/mnt/lustre/bio/users/fcoutinho/Rlibs/")

###Updated MAG data
mag_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Info_All_MAGs_Profiles_Malaspina_Redo_DF_Summary+Prophage+MAG_Metadata.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

###Metagenome RPKM MAG distances
all_mag_abd_df<-read.table(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_all_MAGs/Profiles_Malaspina_All_MAGs_RPKM_Abundance.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

rownames(all_mag_abd_df)<-all_mag_abd_df[,1]
all_mag_abd_df<-all_mag_abd_df[,-1]

dist_metric<-"bray"

mag_abd_dists<-vegdist(all_mag_abd_df, method = dist_metric)

mag_abd_dists_df<-as.data.frame(as.matrix(mag_abd_dists))

#valid_mags_posit<-rownames(mag_abd_dists_df) %in% mag_data$MAG[which(mag_data$Domain == "Archaea")]
valid_mags_posit<-rownames(mag_abd_dists_df) %in% mag_data$MAG[which(mag_data$Domain == "Bacteria")]
summary(valid_mags_posit)
valid_mag_abd_dists_df<-mag_abd_dists_df[valid_mags_posit,valid_mags_posit]
dim(valid_mag_abd_dists_df)

mag_bray_pca<-princomp(valid_mag_abd_dists_df)
mag_bray_pca_loadings<-mag_bray_pca$loadings

summary(mag_bray_pca_loadings)

top10_mag_bray_pca_loadings<-as.data.frame(mag_bray_pca_loadings[,1:10])
top10_mag_bray_pca_loadings$MAG<-rownames(all_mag_abd_df)[valid_mags_posit]

colnames(top10_mag_bray_pca_loadings)<-gsub("Comp.","Abd_PC_",colnames(top10_mag_bray_pca_loadings),perl=TRUE)

summary(top10_mag_bray_pca_loadings)

###GTDBtk tree MAG phylo distance
library("ape")

#gtdb_phylo_tree<-read.tree(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/MAGs_Original_and_Refined_GTDBtk_Output/Profiles_Malaspina_MetaBat_Redo_MAGs_GTDBtk_Output.ar122.classify.tree")
gtdb_phylo_tree<-read.tree(file="/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/MAGs_Original_and_Refined_GTDBtk_Output/Profiles_Malaspina_MetaBat_Redo_MAGs_GTDBtk_Output.bac120.classify.tree")

summary(gtdb_phylo_tree)

leaves_to_remove<-gtdb_phylo_tree$tip.label[!(gtdb_phylo_tree$tip.label %in% mag_data$MAG)]

sub_gtdb_phylo_tree<-drop.tip(gtdb_phylo_tree,tip=leaves_to_remove)

summary(sub_gtdb_phylo_tree)

mag_gtdb_phylo_dists<-cophenetic.phylo(sub_gtdb_phylo_tree)

mag_gtdb_phylo_pca<-princomp(mag_gtdb_phylo_dists)
mag_gtdb_phylo_pca_loadings<-mag_gtdb_phylo_pca$loadings
top10_mag_gtdb_phylo_pca_loadings<-as.data.frame(mag_gtdb_phylo_pca_loadings[,1:10])
top10_mag_gtdb_phylo_pca_loadings$MAG<-rownames(top10_mag_gtdb_phylo_pca_loadings)

colnames(top10_mag_gtdb_phylo_pca_loadings)<-gsub("Comp.","Phylo_PC_",colnames(top10_mag_gtdb_phylo_pca_loadings),perl=TRUE)

summary(top10_mag_gtdb_phylo_pca_loadings)

###Maistrenko Phylogenetic x Environmental Variance
library("phylosignal")
library("Caper")
library("FactoMineR")
library("relaimpo", lib.loc="/mnt/lustre/bio/users/fcoutinho/Rlibs/")

mag_data_for_lm<-merge(top10_mag_gtdb_phylo_pca_loadings,subset(mag_data,select=c("MAG","Estimated_Genome_Size","Total_CDS","Total_KOs","Unique_KOs","TKO_Density","UKO_Density","CDS_Density","Domain","Phylum","Class","Order","Family","Genus","Species","Count_of_Unique_Defense_Systems","DS_Prevalence_per_Mbp","Prophage_Count","Prophage_Perc")),by="MAG",all.x=TRUE)

summary(mag_data_for_lm)

mag_data_for_lm<-merge(mag_data_for_lm,top10_mag_bray_pca_loadings,by="MAG",all.x=TRUE)

summary(mag_data_for_lm)

###Multiple models for all response variables
response_vars<-c("Estimated_Genome_Size","Total_CDS","Total_KOs","Unique_KOs","TKO_Density","UKO_Density","CDS_Density")
results_df<-rbind()

for (rvar in response_vars) {
	print(rvar)
	
	lm_fit<-lm(mag_data_for_lm[[rvar]] ~ Phylo_PC_1 + Phylo_PC_2 + Phylo_PC_3 + Phylo_PC_4 + Phylo_PC_5 + Phylo_PC_6 + Phylo_PC_7 + Phylo_PC_8 + Phylo_PC_9 + Phylo_PC_10 + Abd_PC_1 + Abd_PC_2 + Abd_PC_3 + Abd_PC_4 + Abd_PC_5 + Abd_PC_6 + Abd_PC_7 + Abd_PC_8 + Abd_PC_9 + Abd_PC_10, data=mag_data_for_lm)
	
	#print(summary(lm_fit))

	anova_lm_fit <- anova(lm_fit)
	#print(anova_lm_fit)
	anova_lm_fit_Estimate_ss<- anova_lm_fit$"Sum Sq"
	df_anova_test<-cbind(anova_lm_fit,PctExp=(anova_lm_fit_Estimate_ss/sum(anova_lm_fit_Estimate_ss)*100))

	#print(df_anova_test)
	
	df_anova_test$R_name<-rownames(df_anova_test)

	#print(" % of variance explained by predictors in ANOVA")
	#print(sum(df_anova_test$PctExp))

	car_relimpo_metrics <- calc.relimp(lm_fit, type = c("car"))
	car_relimpo_metricsDF<-data.frame(cbind(car_relimpo_metrics$car,rownames(car_relimpo_metrics)))
	car_relimpo_metricsDF$R_name<-rownames(car_relimpo_metricsDF)
	car_relimpo_metricsDF<-car_relimpo_metricsDF[car_relimpo_metricsDF$R_name!= "n_genomes_x",]
	colnames(car_relimpo_metricsDF)<-c("car","R_name")

	print(car_relimpo_metricsDF)
	print(" carscore % of variance explained by predictors using CAR metrics")
	TotVarExp<-sum(car_relimpo_metricsDF$car)*100
	VarExpByPhylo<-sum(car_relimpo_metricsDF$car[grepl("Phylo_PC_",car_relimpo_metricsDF$R_name)])*100
	VarExpByAbd<-sum(car_relimpo_metricsDF$car[grepl("Abd_PC_",car_relimpo_metricsDF$R_name)])*100
	print(TotVarExp)

	results_vec<-c(rvar,TotVarExp,VarExpByPhylo,VarExpByAbd)
	results_df<-rbind(results_df,results_vec)
}

colnames(results_df)<-c("Response_Variable","Total_Variance_Explained","Variance_Explained_by_Phylogeny","Variance_Explained_by_Abundance")
results_df<-as.data.frame(results_df)
results_df$Response_Variable<-as.factor(results_df$Response_Variable)
for (colnum in c(2:4)) {
	results_df[,colnum]<-as.numeric(results_df[,colnum])	
}
summary(results_df)

fig_6A<-ggplot(results_df,aes(x=Variance_Explained_by_Abundance,y=Variance_Explained_by_Phylogeny))+geom_smooth(method = "lm",colour="blue",se=FALSE,  formula = y ~ x)+geom_label(aes(label=Response_Variable),size=1.5)+theme_bw()+xlim(0,40)+ylim(10,80)

ggsave("Malaspina_Profiles_Figure_6B_Bacteria.pdf",plot=fig_6A,device=cairo_pdf,width=7,height=5,pointsize=8)

