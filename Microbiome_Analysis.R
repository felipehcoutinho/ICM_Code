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
    data.scores$Sample<-rownames(data.scores)

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
calc_group_sums<-function(abd_df=NA,info_df=NA,first_group_var=NA,debug=FALSE) {
    #Assumes firt column of the abundance df to be sample identifiers
    colnames(abd_df)[1]<-"Sample_UID"
    m_abd_df<-reshape2::melt(abd_df,id="Sample_UID",variable.name="Taxon_UID",value.name="Abundance")

    print(paste("Using as group var: ",first_group_var))
    #colnames(info_df)[colnames(info_df) == first_group_var]<-"Group"

    #Assumes firt column of the info df to be Unique taxon identifiers
    colnames(info_df)[1]<-"Taxon_UID"
    info_df$Group<-info_df[[first_group_var]]
    m_abd_df<-merge(m_abd_df,info_df[,c("Taxon_UID","Group")],by="Taxon_UID",all.x=TRUE)

    #Replace NA values in the Group column by "Unclassified"
    m_abd_df$Group<-as.character(m_abd_df$Group)
    m_abd_df$Group[which(m_abd_df$Group == "NA")] <- "Unclassified"
    m_abd_df$Group[which(m_abd_df$Group == "")] <- "Unclassified"
    m_abd_df$Group<-as.factor(m_abd_df$Group)

    if (debug == TRUE) {
        print("Summary of m_abd_df")
        print(summary(m_abd_df))
    }

    group_abd_df<-reshape2::dcast(m_abd_df,Sample_UID ~ Group,value.var="Abundance",fun.aggregate=sum,fill=0)

    return(group_abd_df)
}

###Expand host predictions
library(data.table)
expand_host<-function(seq_info_df,id_var,group_var,score_var,exp_vars) {
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

	f_gdata<-gdata[gdata[, .I[which.min(score_var)], by=group_var]$V1]

	f_gdata<-as.data.frame(f_gdata)
	gdata$score_var<-NULL
	gdata$group_var<-NULL
	f_gdata$score_var<-NULL
	f_gdata$group_var<-NULL

	exp_df<-merge(seq_info_df,subset(f_gdata,select=c(group_var,exp_vars)),by=group_var,all.x=TRUE,suffixes=c("",paste("_Expanded_by_",group_var,sep="")))

	print("Summary of expanded df:")
	print(summary(exp_df))

	return(exp_df)
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

# sub_metadata_df<-metadata_df[,num_cols]

# good_vars<-c("Temperature","Oxygen","ChlorophyllA","Carbon.total","Salinity","Gradient.Surface.temp.SST.","Fluorescence","CO3","HCO3","Density","PO4","PAR.PC","NO3","Si","Alkalinity.total","Ammonium.5m","Depth.Mixed.Layer","Lyapunov","NO2","Depth.Min.O2","NO2NO3","Nitracline","Brunt.Väisälä","Iron.5m","Depth.Max.O2","Okubo.Weiss","Residence.time","Mean.Chloro.HPLC.adjusted..mg.Chl.m3.","Mean.Flux.at.150m","NPP.8d.VGPM..mgC.m2.day.")

# plot<-ggpairs(subset(metadata_df,select=good_vars),lower = list(continuous = scatter_fn_no_se),upper = list(continuous = cor_fn),aes(alpha = 0.5))

# ggsave("Tara_Metadata_Correlogram.pdf",plot=plot,device=cairo_pdf,width=35,height=35, pointsize=8)

# mdata<-melt(subset(metadata_df,select=c("Group",good_vars)),id=c("Group"),variable.name="Variable",value.name="Value")

# plot<-ggplot(data=mdata,aes(y=Value))+geom_density()+coord_flip()+facet_wrap(Variable ~ ., ncol=6, scales="free")
# ggsave("Tara_Metadata_Desnity.pdf",plot=plot,device=cairo_pdf,width=35,height=35, pointsize=8)

