# require(maps)
# library(ggsignif)
# library("GGally")

# library("vegan")

library(dplyr)
library(tidyr)
library(tibble)

library(ggpubr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

###Note: for abundance matrixes, assumes that columns are taxa and that rows are samples

###Calc PCA
calc_pca<-function(pca_df=NA,scale_data=TRUE) {
    #Expects samples as rows and variables as columns. Sample names will be taken from rownames. No extra columns or non-numeric variables accepted.

    print("Removing Samples containign NA values")
    pca_df<-na.omit(pca_df)

    sdata<-pca_df
    if (scale_data == TRUE) {
        print("Scaling data")
        sdata<-as.data.frame(scale(pca_df,center=TRUE,scale=TRUE))
    } else {
        print("Data will not be scale")
    }
    
    print("Summary of PCA data:")
    print(summary(sdata))
    
    library("vegan")
    print("Calculating PCA")

    set.seed(666)
    rdadata<-rda(sdata)
    xlabel<-round((summary(rdadata)$cont$importance[2,1]*100),digits=1)
    xlabel<-paste("PC1 (",xlabel,"% explained)",sep="")

    ylabel<-round((summary(rdadata)$cont$importance[2,2]*100),digits=1)
    ylabel<-paste("PC2 (",ylabel,"% explained)",sep="")

    uscores <- data.frame(rdadata$CA$u)
    vscores <- data.frame(rdadata$CA$v)

    return(list(uscores=uscores,vscores=vscores,xlabel=xlabel,ylabel=ylabel))
}
###Split taxonomy(GTDBtk format)

split_tax<-function(data_df=NA,tax_col="classification") {

    data_df$Full_Classification<-data_df[[tax_col]]

    data_df$classification<-gsub("((d__)|(p__)|(c__)|(o__)|(f__)|(g__)|(s__))","",data_df$classification,perl=TRUE)

    data_df<-data_df %>% separate(classification, c("Domain","Phylum","Class","Order","Family","Genus","Species"),sep=";")

    #replace empty string values in columns "Domain","Phylum","Class","Order","Family","Genus","Species" with "Unclassified":
    for (level in c("Domain","Phylum","Class","Order","Family","Genus","Species")){
        data_df[[level]][(which(data_df[[level]] == ""))]<-"Unclassified"  
        data_df[[level]]<-as.factor(data_df[[level]])
    }

    return(data_df)
}

###Pairwise correlations wih FDR adjusted p-values
calc_pair_cor<-function(a_vars=NA,b_vars=NA,data_df=NA,method="spearman") {
    correls_df<-rbind()
    for (a_var in a_vars) {
        for (b_var in b_vars) {
            print(paste("Calculating correlation between",a_var,"and",b_var))
            result_cor<-cor.test(data_df[[a_var]],data_df[[b_var]],method=method,alternative="two.sided",exact=FALSE,conf.level=0.95,na.option=na.omit)
            correls_df<-rbind(c(a_var,b_var,result_cor$estimate,result_cor$p.value),correls_df)
        }
    }
    colnames(correls_df)<-c("Variable_1","Variable_2","Correlation_Coefficient","P_Value")
    
    correls_df<-as.data.frame(correls_df)
    correls_df$Correlation_Coefficient<-as.numeric(correls_df$Correlation_Coefficient)
    correls_df$P_Value<-as.numeric(correls_df$P_Value)

    correls_df$Variable_1<-as.factor(correls_df$Variable_1)
    correls_df$Variable_2<-as.factor(correls_df$Variable_2)

    correls_df$Adjusted_P_Value_FDR<-p.adjust(correls_df$P_Value,method="fdr")

    correls_df$Method<-method
    
    summary(correls_df)

    return(correls_df)
}


###pairwise mann whitney tests with FDR correction
pairwise_mw<-function(input_df=NA,meta_df=NA,group_var=NA) {
    #Expected abd_df, wide format: SAMPLES AS COLUMNS, ROWS AS variables to be tested, first column is the variable identifier
    all_results<-as.data.frame(rbind())

	rownames(input_df)<-input_df[,1] # first column names to rename rows
	input_df<-input_df[,-1] # remove first column

	rownames(meta_df)<-meta_df[,1] # first column names to rename rows
	meta_df<-meta_df[,-1] # remove first column
    
    print(paste("Metadata DF dimensions after removing first column:",nrow(meta_df),"rows,",ncol(meta_df),"columns",sep=" "))

    #make subset of the meta_df only containing sampels that are also in the abundance df
    meta_df<-meta_df[which(rownames(meta_df) %in% colnames(input_df)),]

    print(paste("Metadata DF dimensions after subsetting to match samples in abundance df:",nrow(meta_df),"rows,",ncol(meta_df),"columns",sep=" "))

	# merge with metadata
	valid_groups<-as.vector(na.omit(unique(meta_df[[group_var]])))

    print(paste("Will compare following groups:",paste(valid_groups,sep=","),sep=" "))
	
	# empty vector to save pairs of values already seen
	seen_combos<-c()
	
	for (groupA in valid_groups) {
		#make subset of input df using only group A samples
		groupA_samples<-rownames(meta_df)[which(meta_df[[group_var]] == groupA)] 
		groupA_df<-input_df[,groupA_samples]	
		for (groupB in valid_groups) {
			#Avoid testing a group of samples against itself
			if (groupA == groupB) { next }
			#Avoid testing groupB against A if A has already been tested against B
            combo<-paste(groupA,groupB,sep="_")
			if (combo %in% seen_combos) { next }
			# add to the list to no repeat the same analysis
			seen_combos<-c(seen_combos,paste(groupA,groupB,sep="_"))
            seen_combos<-c(seen_combos,paste(groupB,groupA,sep="_"))

			#make subset of input df using only group B samples
			groupB_samples<-rownames(meta_df)[which(meta_df[[group_var]] == groupB)] 
			groupB_df<-input_df[,groupB_samples]

			#Iterate over each row of the input df		
			for (rownum in 1:nrow(input_df)) {
				var<-rownames(input_df)[rownum]
				
				#get a vactor of values for the variable being tested from groupA and groupB
                # Add 1 to everything so that FC can be calculated when the mean is 0
				vals_groupA<-as.numeric(groupA_df[rownum,])+1
				vals_groupB<-as.numeric(groupB_df[rownum,])+1
				
				# perform Mann-Whitney test
				test_result<-wilcox.test(vals_groupA,vals_groupB) 
				
				# add rows with results to a dataframe (from the test, save the p-value)
                #print(c(var,groupA,groupB,mean(vals_groupA),mean(vals_groupB),mean(vals_groupA)/mean(vals_groupB),log10(mean(vals_groupA)/mean(vals_groupB)),test_result$p.value))
				all_results<-rbind(all_results,c(var,groupA,groupB,mean(vals_groupA),mean(vals_groupB),mean(vals_groupA)/mean(vals_groupB),log10(mean(vals_groupA)/mean(vals_groupB)),test_result$p.value)) 
				
			}
		}
	}

	colnames(all_results)<-c("Variable","Group_A","Group_B","Mean_A","Mean_B","Fold_Change","Log10_Fold_Change","p_value")
	all_results<-as.data.frame(all_results)
	all_results$Variable<-as.factor(all_results$Variable)
	all_results$Group_A<-as.factor(all_results$Group_A)
	all_results$Group_B<-as.factor(all_results$Group_B)
	all_results$Mean_A<-as.numeric(all_results$Mean_A)
	all_results$Mean_B<-as.numeric(all_results$Mean_B)
	all_results$Fold_Change<-as.numeric(all_results$Fold_Change)
	all_results$Log10_Fold_Change<-as.numeric(all_results$Log10_Fold_Change)
	all_results$p_value<-as.numeric(all_results$p_value)
	all_results$Adjusted_P_Value_FDR<-p.adjust(all_results$p_value,method="fdr")
	return(all_results)
}


###Diversity Metrics
calc_div<-function(abd_df=NA) {
    library("vegan")
    shannon_div<-as.data.frame(diversity(abd_df, index = "shannon"))
    shannon_div$Sample_UID<-rownames(shannon_div)
    colnames(shannon_div)[1]<-"Shannon_Diversity_Index"

    simp_div<-as.data.frame(diversity(abd_df, index = "simpson"))
    simp_div$Sample_UID<-rownames(simp_div)
    colnames(simp_div)[1]<-"Simpson_Diversity_Index"

    spec_rich<-as.data.frame(specnumber(abd_df))
    spec_rich$Sample_UID<-rownames(spec_rich)
    colnames(spec_rich)[1]<-"Richness"

    div_metrics<-merge(shannon_div,simp_div,by="Sample_UID",all.xy=TRUE)
    div_metrics<-merge(div_metrics,spec_rich,by="Sample_UID",all.xy=TRUE)

    div_metrics$Sample_UID<-as.factor(div_metrics$Sample_UID)

    return(div_metrics)
}

###NMDS
calc_nmds<-function(abd_df=NA,dist_metric="bray") {
    library("vegan")
    print(paste("Calculating",dist_metric,"distances",sep=" "))
    dists<-vegdist(abd_df, method = dist_metric)

    print("Calculating NMDS")
    set.seed(666)
    mdsresult<-metaMDS(dists,distance = dist_metric,k = 2,maxit = 999)
    data.scores<-as.data.frame(scores(mdsresult))
    data.scores$Sample_UID<-rownames(data.scores)

    nmds_data<-data.scores

    print(summary(nmds_data))

    return(nmds_data)
}

###Take community DF as input and calculate abundance stats for each variable, mean and median abundance, prevalence (non zero counts acrosss samples), and standard deviation of abundance values
calc_abund_stats<-function(abd_df=NA,exclude_cols=NA) {
    abd_df<-abd_df[,which(!(colnames(abd_df) %in% exclude_cols))]
    #calculate how many non zero values occur in each column of abd_df
    col_nz_count<-apply(abd_df,2,function(x) length(which(x > 0)))
    #calculate median values in each column of abd_df
    col_medians<-apply(abd_df,2,function(x) median(x))
    col_sd<-apply(abd_df,2,function(x) sd(x))
    col_means<-colMeans(abd_df,na.rm=TRUE)
    stats_df<-cbind(colnames(abd_df),col_nz_count,col_means,col_medians,col_sd)
    colnames(stats_df)<-c("UID","Prevalence","Mean_Abundance","Median_Abundance","SD_Abundance")
    stats_df<-as.data.frame(stats_df)
    stats_df$Prevalence<-as.numeric(stats_df$Prevalence)
    stats_df$Mean_Abundance<-as.numeric(stats_df$Mean_Abundance)
    stats_df$Median_Abundance<-as.numeric(stats_df$Median_Abundance)
    stats_df$SD_Abundance<-as.numeric(stats_df$SD_Abundance)
    stats_df$UID<-as.factor(stats_df$UID)
    return(stats_df)
}

#stats_df<-calc_abund_stats(abd_df=abd_df,exclude_cols=c(MG_ID))#,"Sample_UID","Dataset","Data_Type"
#write.table(stats_df,file=outfile,sep="\t",quote=FALSE,row.names=FALSE)


###Calc abundance sums by group
calc_group_sums<-function(abd_df=NA,info_df=NA,first_group_var=NA,debug=FALSE,transpose_abd=FALSE,info_id_var=NA,abd_id_var=NA) {
    #Expected abd_df: SAMPLES AS COLUMNS, ROWS AS TAXA, first column is the taxon identifier
    if (transpose_abd == TRUE) {
        print(paste("Transposing abundance DF"))
        abd_df<-as.data.frame(t(abd_df))
        abd_df<-as.data.frame(cbind(rownames(abd_df),abd_df))
    }

    #Unless otherwise specified in the call to the faction, assumes firt column of the info df to be sample identifiers
    if (is.na(info_id_var)) {
        info_id_var<-colnames(info_df)[1]
    } 
    print(paste("Using ",info_id_var,"as Taxon UID in info DF"))
    info_df$Taxon_UID<-info_df[[info_id_var]]

    #The grouping variable must always be defined
    if (is.na(first_group_var)) {
        print("Must define one of the columns in the info df as the grouping variable")
    } 
    print(paste("Using as group var in info df: ",first_group_var))
    info_df$Group<-info_df[[first_group_var]]

    #Unless otherwise specified in the call to the faction, assumes firt column of the abundance df to be sample identifiers
    if (is.na(abd_id_var)) {
        abd_id_var<-colnames(abd_df)[1]
    } 
    
    print(paste("Using as abundance UID var: ",abd_id_var))
    
    #Rename the previous identifier column in the abund df to ensure the Taxon_UID column is always present and to avoid adding extra columns
    if (abd_id_var != "Taxon_UID") {
        colnames(abd_df)[which(colnames(abd_df) == abd_id_var)]<-"Taxon_UID"
    }
    
    #Add the relevant group into the abundance df as an additional COLUMN
    abd_df<-merge(abd_df,info_df[,c("Taxon_UID","Group")],by="Taxon_UID",all.x=TRUE)

    #Replace NA values in the Group column by "Unclassified"
    abd_df$Group<-as.character(abd_df$Group)
    abd_df$Group[which(abd_df$Group == "NA")] <- "Unclassified"
    abd_df$Group[which(abd_df$Group == "")] <- "Unclassified"
    abd_df$Group[which(is.na(abd_df$Group))] <- "Unclassified"
    abd_df$Group<-as.factor(abd_df$Group)

    m_abd_df<-reshape2::melt(abd_df,id=c("Taxon_UID","Group"),variable.name="Sample_UID",value.name="Abundance")

    if (debug == TRUE) {
        print("Summary of m_abd_df")
        print(summary(m_abd_df))
    }

    group_abd_df<-reshape2::dcast(m_abd_df,Sample_UID ~ Group,value.var="Abundance",fun.aggregate=sum,fill=0)

    return(group_abd_df)
}

###Expand host predictions
expand_host<-function(seq_info_df,id_var,group_var,score_var,exp_vars,minimize_score) {
    library(data.table)
	# exp_vars=c("PHIST_MAG","PHIST_Phylum")
	# score_var<-"PHIST_pvalue"
	# id_var<-"Sequence" 
	# seq_info_df<-full_vir_scaff_data
	# group_var<-"Population"
	
	gdata<-seq_info_df[!is.na(seq_info_df[[group_var]]),]
	print("Dimensions of gdata after removing rows with missing columns in group var:")
	print(dim(gdata))

	gdata<-seq_info_df[!is.na(seq_info_df[[score_var]]),]
	print("Dimensions of gdata after removing rows with missing columns in score var:")
	print(dim(gdata))

	gdata<-subset(gdata,select=c(id_var,group_var,score_var,exp_vars))

	gdata$score_var<-gdata[[score_var]]
	gdata$group_var<-gdata[[group_var]]
	gdata<-as.data.table(gdata)

    if (minimize_score == TRUE) {
        print("Minimizing scores")
	    f_gdata<-gdata[gdata[, .I[which.min(score_var)], by=group_var]$V1]
    } else {
        print("Maximizing scores")
        f_gdata<-gdata[gdata[, .I[which.max(score_var)], by=group_var]$V1]
    }

	f_gdata<-as.data.frame(f_gdata)
	gdata$score_var<-NULL
	gdata$group_var<-NULL
	f_gdata$score_var<-NULL
	f_gdata$group_var<-NULL

	exp_df<-merge(seq_info_df,subset(f_gdata,select=c(group_var,exp_vars)),by=group_var,all.x=TRUE,suffixes=c("",paste("_Expanded_by_",group_var,sep="")))
	print("Summary of expanded df:")
	print(summary(exp_df))

	return(f_gdata)
}

#exp_full_vir_scaff_data<-expand_host(full_vir_scaff_data,id_var="Sequence",group_var="Population",score_var="PHIST_pvalue",exp_vars=c("PHIST_MAG","PHIST_Domain","PHIST_Phylum","PHIST_Class","PHIST_Order","PHIST_Family","PHIST_Genus","PHIST_Species"))


###Genome metrics Correlogram
cor_fn <- function(data, mapping, method="p", use="pairwise", ...){

              # grab data
              x <- eval_data_col(data, mapping$x)
              y <- eval_data_col(data, mapping$y)

              # calculate correlation
			  corr_val<-0
			  if ((length(x) > 5) & (length(y) > 5)) {
				corr_val <- cor(x, y, method=method, use="pairwise.complete.obs")
				}
              # calculate colour based on correlation value
              # Here I have set a correlation of minus one to blue, 
              # zero to white, and one to red 
              # Change this to suit: possibly extend to add as an argument of `my_fn`
              colFn <- colorRampPalette(rev(brewer.pal(9,"RdBu")), interpolate ='spline')
              fill <- colFn(100)[findInterval(corr_val, seq(-1, 1, length=100))]

              ggally_cor(data = data, size = 4, colour="black", mapping = mapping, ...) + 
                theme_void() +
                theme(panel.background = element_rect(fill=fill))
            }
			
scatter_fn <- function(data, mapping, method="p", use="pairwise", ...){
	
              # grab data
              x <- eval_data_col(data, mapping$x)
              y <- eval_data_col(data, mapping$y)

	ggally_smooth_loess(data = data, mapping=mapping, ...)+theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

scatter_fn_no_se <- function(data, mapping, method="p", use="pairwise", ...){
	
              # grab data
              x <- eval_data_col(data, mapping$x)
              y <- eval_data_col(data, mapping$y)

	ggally_smooth_loess(data = data, mapping=mapping, se=FALSE, ...)+theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

