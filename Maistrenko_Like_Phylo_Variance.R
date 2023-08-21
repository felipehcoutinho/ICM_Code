###Input options
mag_data_file<-"/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Info_All_MAGs_Profiles_Malaspina_Redo.tsv"
tree_file<-"/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/MAGs_Original_and_Refined_GTDBtk_Output/Profiles_Malaspina_MetaBat_Redo_MAGs_GTDBtk_Output.bac120.classify.tree"
abd_file<-"/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_dRep_MAGs/Profiles_Malaspina_dRep_MAGs_RPKM_Abundance.tsv"
level<-"Phylum"
domain<-"Bacteria"
dist_metric<-"bray"
transfo<-"none"
max_pcs<-20
min_pc_var_exp<-0.9501
response_vars<-c("Estimated_Genome_Size","Total_CDS","Total_KOs","Unique_KOs","TKO_Density","UKO_Density","CDS_Density","Count_of_Unique_Defense_Systems","DS_Prevalence_per_Mbp","Prophage_Count","Prophage_Perc")
rep_only<-TRUE

args = commandArgs(trailingOnly=TRUE)
level<-args[1]
domain<-args[2]
dist_metric<-args[3]
transfo<-args[4]
max_pcs<-as.numeric(args[5])
min_pc_var_exp<-as.numeric(args[6])
tree_file<-args[7]
abd_file<-args[8]
rep_only<-args[9]

library("vegan")
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr, lib.loc="/mnt/lustre/bio/users/fcoutinho/Rlibs/")
library("ape")
library("relaimpo", lib.loc="/mnt/lustre/bio/users/fcoutinho/Rlibs/")
###Updated MAG data
print(paste("Reading MAG metadata from: ",mag_data_file,sep=""))
mag_data<-read.table(file=mag_data_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

if (rep_only) {
	print("Removing replicate MAGs from MAG data DF")
	mag_data<-mag_data[which(mag_data$Cluster_Representative == TRUE),]
}

#function for performing PCA based on phylogenetic distances
get_phylo_pca<-function (level,taxon) {
	leaves_to_keep<-gtdb_phylo_tree$tip.label[(gtdb_phylo_tree$tip.label %in% mag_data$MAG[which(mag_data[[level]] == taxon)])]
	leaves_to_remove<-gtdb_phylo_tree$tip.label[!(gtdb_phylo_tree$tip.label %in% mag_data$MAG[which(mag_data[[level]] == taxon)])]

	#print(paste("Keeping the following",length(leaves_to_keep),"leaves from the tree:",paste(leaves_to_keep,collapse = ","),sep=" "))
	#print(paste("Removing the following",length(leaves_to_remove),"leaves from the tree:",paste(leaves_to_remove,collapse = ","),sep=" "))
	sub_gtdb_phylo_tree<-drop.tip(gtdb_phylo_tree,tip=leaves_to_remove)

	#summary(sub_gtdb_phylo_tree)
	print(paste("Calculating cophenetic distances among MAGs from ",level," ",taxon,sep=""))
	mag_gtdb_phylo_dists<-cophenetic.phylo(sub_gtdb_phylo_tree)

	print("Performing PCA on MAG cophenetic distances matrix")
	mag_gtdb_phylo_pca<-princomp(mag_gtdb_phylo_dists)
	mag_gtdb_phylo_pca_loadings<-mag_gtdb_phylo_pca$loadings
	
	print("Calculating relative abundance PCA eigenvalues")
	mag_phylo_pca_evals<-eigenvals(mag_gtdb_phylo_pca)
	evals_summary<-summary(mag_phylo_pca_evals)
	evals_CumPopExp<-evals_summary["Cumulative Proportion",]
	n_pc_for_min_var_phylo<-length(which(evals_CumPopExp <= min_pc_var_exp))
	n_pc_phylo_to_use<-min(c(n_pc_for_min_var_phylo,max_pcs))

	if (n_pc_phylo_to_use == 1) {
		top10_mag_gtdb_phylo_pca_loadings<-as.data.frame(mag_gtdb_phylo_pca_loadings[,1])
		colnames(top10_mag_gtdb_phylo_pca_loadings)<-"Phylo_PC_1"
	} else {
		top10_mag_gtdb_phylo_pca_loadings<-as.data.frame(mag_gtdb_phylo_pca_loadings[,1:n_pc_phylo_to_use])
		colnames(top10_mag_gtdb_phylo_pca_loadings)<-gsub("Comp.","Phylo_PC_",colnames(top10_mag_gtdb_phylo_pca_loadings),perl=TRUE)
	}
	
	top10_mag_gtdb_phylo_pca_loadings$MAG<-rownames(top10_mag_gtdb_phylo_pca_loadings)

	#summary(mag_phylo_pca_evals)[3,n_pc_phylo_to_use],
	print(paste("Will use ",n_pc_phylo_to_use," PCs from the phylogenetic PCA, which together explain ",round((evals_CumPopExp[n_pc_phylo_to_use]*100),digits=1)," % of the variance"),sep="")
	
	return(top10_mag_gtdb_phylo_pca_loadings)
}

#function for performing PCA based on abundance distances
get_abd_pca<-function (level,taxon) {
	valid_mags_posit<-rownames(mag_abd_dists_df) %in% mag_data$MAG[which(mag_data[[level]] == taxon)]

	#summary(valid_mags_posit)
	valid_mag_abd_dists_df<-mag_abd_dists_df[valid_mags_posit,valid_mags_posit]
	#dim(valid_mag_abd_dists_df)
	
	print(paste("Performing PCA on relative abundance distances matrix of MAGs from ",level," ",taxon,sep=""))
	mag_bray_pca<-princomp(valid_mag_abd_dists_df)
	mag_bray_pca_loadings<-mag_bray_pca$loadings

	print("Calculating relative abundance PCA eigenvalues")
	mag_abd_pca_evals<-eigenvals(mag_bray_pca)
	evals_summary<-summary(mag_abd_pca_evals)
	evals_CumPopExp<-evals_summary["Cumulative Proportion",]
	n_pc_for_min_var_abd<-length(which(evals_CumPopExp <= min_pc_var_exp))
	n_pc_abd_to_use<-min(c(n_pc_for_min_var_abd,max_pcs))
	
	if (n_pc_abd_to_use == 1) {
		top10_mag_bray_pca_loadings<-as.data.frame(mag_bray_pca_loadings[,1])
		colnames(top10_mag_bray_pca_loadings)<-"Abd_PC_1"
	} else {
		top10_mag_bray_pca_loadings<-as.data.frame(mag_bray_pca_loadings[,1:n_pc_abd_to_use])
		colnames(top10_mag_bray_pca_loadings)<-gsub("Comp.","Abd_PC_",colnames(top10_mag_bray_pca_loadings),perl=TRUE)
	}


	print(paste("Will use ",n_pc_abd_to_use," PCs from the abundance PCA, which together explain ",round((evals_CumPopExp[n_pc_abd_to_use]*100),digits=1)," % of the variance",sep=""))
	
	#top10_mag_bray_pca_loadings<-as.data.frame(mag_bray_pca_loadings[,1:n_pc_abd_to_use])
	top10_mag_bray_pca_loadings$MAG<-rownames(all_mag_abd_df)[valid_mags_posit]

	colnames(top10_mag_bray_pca_loadings)<-gsub("Comp.","Abd_PC_",colnames(top10_mag_bray_pca_loadings),perl=TRUE)
	abd_pcs<-colnames(top10_mag_bray_pca_loadings)[which(colnames(top10_mag_bray_pca_loadings) != "MAG")]
	#print(paste("Finished PCA on relative abundance distances matrix of MAGs from ",level," ",taxon,sep=""))
	return(top10_mag_bray_pca_loadings)
}

#Read phylo tree and calc distances
print(paste("Reading phylogenetic tree from: ",tree_file,sep=""))
gtdb_phylo_tree<-read.tree(file=tree_file)

###Metagenome RPKM MAG data
print(paste("Reading MAG relative abundance data: ",abd_file,sep=""))
all_mag_abd_df<-read.table(file=abd_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

rownames(all_mag_abd_df)<-all_mag_abd_df[,1]
all_mag_abd_df<-all_mag_abd_df[,-1]

#Optional abundance data transformations
if (transfo == "clr") {
	print("Performing CLR transformation of MAG relative abundance data")
	library("compositions")
	all_mag_abd_df<-clr(all_mag_abd_df)
}
print("Calculating distances among MAGs based on relative abundance data")
mag_abd_dists<-vegdist(all_mag_abd_df, method = dist_metric)
mag_abd_dists_df<-as.data.frame(as.matrix(mag_abd_dists))

results_df<-c()
taxon_mag_count<-table(mag_data[which(mag_data$Domain == domain),][[level]])
taxa_list<-names(taxon_mag_count[which(taxon_mag_count >= 10)])

####Loop over taxa
for (tx in taxa_list) {
	print(paste("Fitting models for taxon ",tx,", which is represented by ",taxon_mag_count[tx]," MAGs",sep=""))
	
	top10_mag_gtdb_phylo_pca_loadings<-get_phylo_pca(level,tx)
	phylo_pcs<-colnames(top10_mag_gtdb_phylo_pca_loadings)[which(colnames(top10_mag_gtdb_phylo_pca_loadings) != "MAG")]
	print("Phylo_PCs:")
	print(phylo_pcs)
	
	top10_mag_bray_pca_loadings<-get_abd_pca(level,tx)
	abd_pcs<-colnames(top10_mag_bray_pca_loadings)[which(colnames(top10_mag_bray_pca_loadings) != "MAG")]
	print("Abd_PCs:")
	print(abd_pcs)
	
	print("Merging data from Phylogenetic and Habitat PCA")
	mag_data_for_lm<-merge(top10_mag_gtdb_phylo_pca_loadings,subset(mag_data,select=c("MAG",response_vars,"Domain","Phylum","Class","Order","Family","Genus","Species")),by="MAG",all.x=TRUE)

	#summary(mag_data_for_lm)

	mag_data_for_lm<-merge(mag_data_for_lm,top10_mag_bray_pca_loadings,by="MAG",all.x=TRUE)

	#print(summary(mag_data_for_lm))

	###Multiple models for all response variables
	for (rvar in response_vars) {
		#print(rvar)
		print(paste("Fitting models for response variable: ",rvar,sep=""))
		mag_data_for_lm$Respone_Varibale<-mag_data_for_lm[[rvar]]
		sub_mag_data_for_lm<-subset(mag_data_for_lm,select=c("Respone_Varibale",phylo_pcs,abd_pcs))
		
		if(sd(sub_mag_data_for_lm$Respone_Varibale) == 0) {
			next
		}
		
		#print(summary(sub_mag_data_for_lm))
		lm_fit<-lm(Respone_Varibale ~ . , data=sub_mag_data_for_lm)

		#Skip vars with 0 variance in the LM coefficients
		if(sd(lm_fit$coefficients) == 0) {
			next
		}
		#print(summary(lm_fit))
		#anova_lm_fit <- anova(lm_fit)
		#print(anova_lm_fit)
		#anova_lm_fit_Estimate_ss<- anova_lm_fit$"Sum Sq"
		#df_anova_test<-cbind(anova_lm_fit,PctExp=(anova_lm_fit_Estimate_ss/sum(anova_lm_fit_Estimate_ss)*100))
		#print(df_anova_test)
		#df_anova_test$R_name<-rownames(df_anova_test)
		#print(" % of variance explained by predictors in ANOVA")
		#print(sum(df_anova_test$PctExp))

		car_relimpo_metrics <- calc.relimp(lm_fit, type = c("car"))
		car_relimpo_metricsDF<-data.frame(cbind(car_relimpo_metrics$car,rownames(car_relimpo_metrics)))
		car_relimpo_metricsDF$R_name<-rownames(car_relimpo_metricsDF)
		car_relimpo_metricsDF<-car_relimpo_metricsDF[car_relimpo_metricsDF$R_name!= "n_genomes_x",]
		colnames(car_relimpo_metricsDF)<-c("car","R_name")

		#print(car_relimpo_metricsDF)
		TotVarExp<-sum(car_relimpo_metricsDF$car)*100
		print("carscore % of total variance explained by predictors using CAR metrics")
		print(TotVarExp)
		VarExpByPhylo<-sum(car_relimpo_metricsDF$car[grepl("Phylo_PC_",car_relimpo_metricsDF$R_name)])*100
		#print("carscore % of variance explained by phylogenetic predictors using CAR metrics:")
		#print(VarExpByPhylo)
		VarExpByAbd<-sum(car_relimpo_metricsDF$car[grepl("Abd_PC_",car_relimpo_metricsDF$R_name)])*100
		#print("carscore % of variance explained by abundance predictors using CAR metrics:")
		#print(VarExpByAbd)

		results_vec<-c(level,tx,rvar,TotVarExp,VarExpByPhylo,VarExpByAbd)
		results_df<-rbind(results_df,results_vec)
	}
}

colnames(results_df)<-c("Taxonomic_Level","Taxon","Response_Variable","Total_Variance_Explained","Variance_Explained_by_Phylogeny","Variance_Explained_by_Abundance")
results_df<-as.data.frame(results_df)
results_df$Taxonomic_Level<-as.factor(results_df$Taxonomic_Level)
results_df$Taxon<-as.factor(results_df$Taxon)
results_df$Response_Variable<-as.factor(results_df$Response_Variable)
results_df$Total_Variance_Explained<-as.numeric(results_df$Total_Variance_Explained)
results_df$Variance_Explained_by_Phylogeny<-as.numeric(results_df$Variance_Explained_by_Phylogeny)
results_df$Variance_Explained_by_Abundance<-as.numeric(results_df$Variance_Explained_by_Abundance)

summary(results_df)

###Table
print("Printing results")
write.table(results_df,file=paste("Malaspina_Profiles_",domain,"_",level,"_dRep_",rep_only,"_Variance_Explained.tsv",sep=""),sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)
