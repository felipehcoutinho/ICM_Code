###Input options
abd_file<-"/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Abundance_all_MAGs/Profiles_Malaspina_All_MAGs_RPKM_Abundance.tsv"
tree_file<-"/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/MAGs_Original_and_Refined_GTDBtk_Output/Profiles_Malaspina_MetaBat_Redo_MAGs_GTDBtk_Output.bac120.classify.tree"
mag_data_file<-"/mnt/lustre/scratch/fcoutinho/Profiles_Malaspina/Assemblies_Round2/Metabat_Binning/Redo_MetaBat_Bins_Round_1/Info_All_MAGs_Profiles_Malaspina_Redo_DF_Summary+Prophage+MAG_Metadata.tsv"
level<-"Phylum"
dist_metric<-"bray"
transfo<-"none"
max_pcs<-20
min_pc_var_exp<-0.95
response_vars<-c("Estimated_Genome_Size","Total_CDS","Total_KOs","Unique_KOs","TKO_Density","UKO_Density","CDS_Density","Count_of_Unique_Defense_Systems","DS_Prevalence_per_Mbp","Prophage_Count","Prophage_Perc")

args = commandArgs(trailingOnly=TRUE)
level<-args[1]
taxon<-args[2]
dist_metric<-args[3]
transfo<-args[4]
max_pcs<-as.numeric(args[5])
min_pc_var_exp<-as.numeric(args[6])

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

#function for performing PCA based on phylogenetic distances
get_phylo_pca<-function (level,taxon) {
	leaves_to_remove<-gtdb_phylo_tree$tip.label[!(gtdb_phylo_tree$tip.label %in% mag_data$MAG[which(mag_data[[level]] == taxon)])]

	sub_gtdb_phylo_tree<-drop.tip(gtdb_phylo_tree,tip=leaves_to_remove)

	#summary(sub_gtdb_phylo_tree)
	print("Calculating cophenetic distances among MAGs")
	mag_gtdb_phylo_dists<-cophenetic.phylo(sub_gtdb_phylo_tree)

	print("Performing PCA on MAG cophenetic distances matrix")
	mag_gtdb_phylo_pca<-princomp(mag_gtdb_phylo_dists)
	mag_gtdb_phylo_pca_loadings<-mag_gtdb_phylo_pca$loadings
	
	print("Calculating relative abundance PCA eigenvalues")
	mag_phylo_pca_evals<-eigenvals(mag_gtdb_phylo_pca)
	n_pc_for_min_var_phylo<-length(which(summary(mag_phylo_pca_evals)[3,] <= min_pc_var_exp))
	n_pc_phylo_to_use<-min(c(n_pc_for_min_var_phylo,max_pcs))

	top10_mag_gtdb_phylo_pca_loadings<-as.data.frame(mag_gtdb_phylo_pca_loadings[,1:n_pc_phylo_to_use])
	top10_mag_gtdb_phylo_pca_loadings$MAG<-rownames(top10_mag_gtdb_phylo_pca_loadings)

	#summary(mag_phylo_pca_evals)[3,n_pc_phylo_to_use],
	print(paste("Will use ",n_pc_phylo_to_use," PCs from the phylogenetic PCA, which together explain "," % of the variance",sep=""))
	colnames(top10_mag_gtdb_phylo_pca_loadings)<-gsub("Comp.","Phylo_PC_",colnames(top10_mag_gtdb_phylo_pca_loadings),perl=TRUE)
	return(top10_mag_gtdb_phylo_pca_loadings)
}

#function for performing PCA based on abundance distances
get_abd_pca<-function (level,taxon) {
	valid_mags_posit<-rownames(mag_abd_dists_df) %in% mag_data$MAG[which(mag_data[[level]] == taxon)]

	#summary(valid_mags_posit)
	valid_mag_abd_dists_df<-mag_abd_dists_df[valid_mags_posit,valid_mags_posit]
	#dim(valid_mag_abd_dists_df)
	
	print("Performing PCA on MAG relative abundance distances matrix")

	mag_bray_pca<-princomp(valid_mag_abd_dists_df)
	mag_bray_pca_loadings<-mag_bray_pca$loadings

	print("Calculating relative abundance PCA eigenvalues")
	mag_abd_pca_evals<-eigenvals(mag_bray_pca)
	n_pc_for_min_var_abd<-length(which(summary(mag_abd_pca_evals)[3,] <= min_pc_var_exp))
	n_pc_abd_to_use<-min(c(n_pc_for_min_var_abd,max_pcs))
	#summary(mag_bray_pca_loadings)
	#summary(mag_abd_pca_evals)[3,n_pc_abd_to_use],
	print(paste("Will use ",n_pc_abd_to_use," PCs from the abundance PCA, which together explain "," % of the variance",sep=""))
	top10_mag_bray_pca_loadings<-as.data.frame(mag_bray_pca_loadings[,1:n_pc_abd_to_use])
	top10_mag_bray_pca_loadings$MAG<-rownames(all_mag_abd_df)[valid_mags_posit]

	colnames(top10_mag_bray_pca_loadings)<-gsub("Comp.","Abd_PC_",colnames(top10_mag_bray_pca_loadings),perl=TRUE)
	abd_pcs<-colnames(top10_mag_bray_pca_loadings)[which(colnames(top10_mag_bray_pca_loadings) != "MAG")]
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
taxon_mag_count<-table(mag_data[which(mag_data$Domain == taxon),][[level]])
taxa_list<-names(taxon_mag_count[which(taxon_mag_count >= 10)])
####Loop over taxa
for (tx in taxa_list) {
	print(paste("Fitting models for taxon: ",tx,sep=""))
	
	top10_mag_gtdb_phylo_pca_loadings<-get_phylo_pca(level,tx)
	phylo_pcs<-colnames(top10_mag_gtdb_phylo_pca_loadings)[which(colnames(top10_mag_gtdb_phylo_pca_loadings) != "MAG")]
	print(paste("Phylo_PCs: ",phylo_pcs,sep=""))
	
	top10_mag_bray_pca_loadings<-get_abd_pca(level,tx)
	abd_pcs<-colnames(top10_mag_bray_pca_loadings)[which(colnames(top10_mag_bray_pca_loadings) != "MAG")]
	print(paste("Abd_PCs: ",abd_pcs,sep=""))
	
	print("Merging data from Phylogenetic and Habitat PCA")
	mag_data_for_lm<-merge(top10_mag_gtdb_phylo_pca_loadings,subset(mag_data,select=c("MAG",response_vars,"Domain","Phylum","Class","Order","Family","Genus","Species")),by="MAG",all.x=TRUE)

	#summary(mag_data_for_lm)

	mag_data_for_lm<-merge(mag_data_for_lm,top10_mag_bray_pca_loadings,by="MAG",all.x=TRUE)

	#print(summary(mag_data_for_lm))

	###Multiple models for all response variables
	for (rvar in response_vars) {
		print(rvar)
		print(paste("Fitting models for response variable: ",rvar,sep=""))
		mag_data_for_lm$Respone_Varibale<-mag_data_for_lm[[rvar]]
		sub_mag_data_for_lm<-subset(mag_data_for_lm,select=c("Respone_Varibale",phylo_pcs,abd_pcs))
		
		if(sd(sub_mag_data_for_lm$Respone_Varibale) == 0) {
			next
		}
		
		#print(summary(sub_mag_data_for_lm))
		
		lm_fit<-lm(Respone_Varibale ~ . , data=sub_mag_data_for_lm)
		
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

		#print(car_relimpo_metricsDF)
		print("carscore % of variance explained by predictors using CAR metrics")
		TotVarExp<-sum(car_relimpo_metricsDF$car)*100
		VarExpByPhylo<-sum(car_relimpo_metricsDF$car[grepl("Phylo_PC_",car_relimpo_metricsDF$R_name)])*100
		VarExpByAbd<-sum(car_relimpo_metricsDF$car[grepl("Abd_PC_",car_relimpo_metricsDF$R_name)])*100
		print(TotVarExp)

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
write.table(results_df,file="Variance_Explained.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

library(reshape2)
library(ggplot2)
library(RColorBrewer)


###Scatterplot
fig_name<-paste("Malaspina_Profiles_PhyloxHabitat_Variance_",level,"Transf","_",transfo,"_","Distance_Metric","_",dist_metric,"_Min_Var_Exp_by_PCs_",as.character(min_pc_var_exp),"_Scatterplot.pdf",sep="")

fig1<-ggplot(results_df,aes(x=Variance_Explained_by_Abundance,y=Variance_Explained_by_Phylogeny))+geom_smooth(method = "lm",colour="darkblue",se=FALSE,  formula = y ~ x)+geom_label(aes(label=Response_Variable),size=1.5)+theme_bw()+xlab("% Variance Explained by Habitat")+ylab("% Variance Explained by Phylogeny")+xlim(0,100)+ylim(0,100)+facet_wrap(. ~ Taxon)

ggsave(fig_name,plot=fig1,device=cairo_pdf,width=20,height=20,pointsize=8)

###
results_df<-read.table(file="Variance_Explained.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
summary(results_df)


mdata<-melt(results_df,id=c("Taxonomic_Level","Taxon","Response_Variable"))
fmdata<-mdata[which(mdata$variable != "Total_Variance_Explained"),]

colnames(fmdata)<-c("Taxonomic_Level","Taxon","Response_Variable","Category","Variance_Explained")
fmdata$Category<-as.factor(gsub("Variance_Explained_by_","",fmdata$Category,perl=TRUE))

fmdata$Response_Variable<-as.factor(gsub("_"," ",fmdata$Response_Variable,perl=TRUE))

summary(fmdata)

#fig_name<-"Variance_Explained_Test_Barplot.pdf"

#Heatmap
fig_name<-paste("Malaspina_Profiles_PhyloxHabitat_Variance_",level,"Transf","_",transfo,"_","Distance_Metric","_",dist_metric,"_Min_Var_Exp_by_PCs_",as.character(min_pc_var_exp),"_Heatmap.pdf",sep="")

#fig_name<-"Variance_Explained_Test_Heatmap.pdf"

col_pal<-brewer.pal(11,"Spectral")
col_grad<-rev(colorRampPalette(col_pal)(n=100))

taxon_order<-c("Acidobacteriota","Actinobacteriota","Bacteroidota","Bdellovibrionota","Chlamydiota","Chloroflexota","Cyanobacteria","Desulfobacterota","Desulfobacterota_D","Eremiobacterota","Gemmatimonadota","Margulisbacteria","Marinisomatota","Myxococcota","Nitrospinota","Patescibacteria","Planctomycetota","Poribacteria","Alphaproteobacteria","Gammaproteobacteria","SAR324","UBP7","Verrucomicrobiota","Halobacteriota","Hydrothermarchaeota","Nanoarchaeota","Thermoplasmatota","Thermoproteota")
sub_taxon_order<-taxon_order[which(taxon_order %in% unique(fmdata$Taxon))]

summary(fmdata)

library(stringr)
fmdata$Response_Variable<-as.factor(str_wrap(fmdata$Response_Variable,width=10))

fig3<-ggplot(fmdata,aes(fill=Variance_Explained,y=Taxon,x=Category))+geom_tile(colour = "grey50")+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust = 1,size=6),axis.text.y = element_text(size=8))+scale_fill_gradientn(name="% Variance Explained",colours =col_grad,limits = c(0,100))+theme(legend.position="top")+facet_grid(. ~ Response_Variable)+ylim(rev(sub_taxon_order))

ggsave(fig_name,plot=fig3,device=cairo_pdf,width=10,height=5,pointsize=8)

###Barplot
fig_name<-paste("Malaspina_Profiles_PhyloxHabitat_Variance_",level,"Transf","_",transfo,"_","Distance_Metric","_",dist_metric,"_Min_Var_Exp_by_PCs_",as.character(min_pc_var_exp),"_Barplot.pdf",sep="")


var_coloring<-brewer.pal(9,"RdBu")[c(2,8)]
names(var_coloring)<-unique(fmdata$Category)

fig2<-ggplot(fmdata,aes(x=Taxon,y=Variance_Explained,fill=Category,group=Category))+geom_bar(position="dodge",stat="identity",colour="black",alpha=0.9,linewidth=0.1)+theme_bw()+ylab("% Variance Explained")+theme(axis.text.x = element_text(angle = 90,hjust = 1,size=9),legend.position="top")+scale_fill_manual(name="Explaining Variable",values=var_coloring)+facet_wrap(. ~ Response_Variable, nrow=3)

ggsave(fig_name,plot=fig2,device=cairo_pdf,width=9,height=10,pointsize=8)



#quit(status=0)