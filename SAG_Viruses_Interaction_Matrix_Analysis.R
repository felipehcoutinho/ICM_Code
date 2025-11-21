###Bipartite netwok analysis
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(bipartite)

# # trying with the data below because /mnt/smart/scratch/vir/felipe/gvSAGs/Read_Mapping_All_SAGs/df_rpkm_high_quality_namesCorrected.tsv has NAs
# abd_file<-"/mnt/smart/scratch/vir/felipe/gvSAGs/Read_Mapping_All_SAGs/Raw_Abundance_All_Genomic.tsv" # 
# # load df with SAG information 
# cell_info_file<-"/mnt/smart/scratch/vir/felipe/gvSAGs/SAG_Info.tsv"
# # load df with viral sequence info
# viral_info_file<-"/mnt/smart/scratch/vir/felipe/gvSAGs/Read_Mapping_All_SAGs/viral_fragment_data_2185ctgs.tsv"

# vinfo_df<-read.table(file=viral_info_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

# #create a column in viral info df that combine SAG_ID and viral_contig (will eb necessary to link Dfs later) 
# vinfo_df$fullname<-as.factor(paste(vinfo_df$SAG_ID,vinfo_df$viral_contig,sep="_"))

# cell_info_df<-read.table(file=cell_info_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

# read_map_df<-read.table(file=abd_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

# #melt the abd df into long format
# mdata<-melt(read_map_df,id="Sequence",variable.name="SAG_ID",value.name="Abundance")

# #Create an extra column in mdata that is the SAG ource from which the viral vontig was obtained
# mdata$Sequence_Source<-as.factor(gsub("(_K(\\d)+_NODE.*$)|(_NODE.*$)","",mdata$Sequence,perl=T))

# #Create an extra column in mdata that is the SAG source from which the reads were being mapped from
# mdata$Reads_Source<-mdata$SAG_ID

# #filter mdata to get rid of low abundances viral contig x sag read paris (more important when using the raw abundances, for RPKM this number might need to be lower)
# fmdata<-mdata#[which(mdata$Abundance > 100),]

# #count how many sequences are missing due to inconsistencies in naming of contig/SAG
# summary(unique(fmdata$Sequence_Source) %in% cell_info_df$SAG_ID)

# #list the ids of sequences that are missing due to inconsistencies in naming of contig/SAG
# unique(fmdata$Sequence_Source)[!(unique(fmdata$Sequence_Source) %in% cell_info_df$SAG_ID)]

# # add sag info for viral contigs
# fmdata<-merge(fmdata,subset(cell_info_df,select=c("SAG_ID","Species","Group","Supergroup","SAG.report")),by.x="Sequence_Source",by.y="SAG_ID",all.x=T)

# #chekc for NAs
# unique(fmdata$Sequence_Source[is.na(fmdata$SAG.report)])

# ###rename sag ids in cell_info_df to match the ids in the abundance files (specifically those that start with AH)
# rows_to_change<-grepl("^AH",cell_info_df$SAG_ID,perl=T)
# cell_info_df$SAG_ID_renamed<-as.character(cell_info_df$SAG_ID)
# cell_info_df$SAG_ID_renamed[rows_to_change]<-gsub("_RM(\\d)+.*$","",cell_info_df$SAG_ID_renamed[rows_to_change],perl=T)
# cell_info_df$SAG_ID_renamed<-as.factor(cell_info_df$SAG_ID_renamed)

# #count and list missing sequences again just in case
# summary(fmdata$Reads_Source %in% cell_info_df$SAG_ID_renamed)
# unique(fmdata$Reads_Source)[!(unique(fmdata$Reads_Source) %in% cell_info_df$SAG_ID_renamed)]

# ###Now add sag infor for the reads mapped
# fmdata<-merge(fmdata,cell_info_df[,c("SAG_ID_renamed","Group","Species","SAG.report")],by.x="Reads_Source",by.y="SAG_ID_renamed",suffixes=c("_Viral_Contig","_SAG_Reads"),all.x=T)

# #Infer the name of the genus of the SAG by removing eveyrhint after -sp or _ or a space 
# fmdata$Genus_Sequence_Source<-as.factor(gsub("(-sp.+$)|(_.+$)|(\\s).+$","",fmdata$Species_Viral_Contig,perl=T))

# fmdata$Genus_Reads_Source<-as.factor(gsub("(-sp.+$)|(_.+$)|(\\s).+$","",fmdata$Species_SAG_Reads,perl=T))

# ###List viral contigs with missing source data
# unique(fmdata[is.na(fmdata$Genus_Sequence_Source),"Sequence"])

# #ATM there are 290 scaffolds starting with GC100 with missing data: Maybe these are from the sags that were excluded

# fmdata$Matching_Taxa<-FALSE
# fmdata$Matching_Taxa[which(fmdata$Genus_Sequence_Source == fmdata$Genus_Reads_Source)]<-TRUE

# #check if the data frame has the expected contents
# # summary(fmdata[,c(1:6,18:26)])
# #calculate an interaction matrix by summing the reads mapped t the desired taxonomic level
# #have to the decide if fun.aggregate makes more sense as sum or mean/median or some other metric
# clean_df<-dcast(fmdata, Sequence ~ Genus_Reads_Source, value.var = "Abundance", fill=NaN, fun.aggregate=sum)

# rownames(clean_df)<-clean_df[,1]
# clean_df<-clean_df[,-1]


# csums<-colSums(clean_df)
# names(csums)<-colnames(clean_df)
# cs2exclude<-names(csums)[is.na(csums)]

# # clean_df<-subset(clean_df,select = -c(cs2exclude))
# # clean_df<-subset(clean_df,select = -c(AH2983,AH2996,AH3100,P2_89))

# clean_df<-clean_df[,!grepl("(NA)|(Unknown)",colnames(clean_df),perl=T)]

# table(is.na(clean_df))

#Must either remove viruses that dont map to any named taxa or call these taxa something other than NA
# clean_df<-clean_df[which(rowSums(clean_df) > 0),] 

# #Save formatted df to file
# write.table(clean_df,file="/mnt/smart/scratch/vir/felipe/gvSAGs/Read_Mapping_All_SAGs/SAG_Viral_Clean_DF.tsv",sep="\t",append=FALSE,row.names=TRUE,col.names=NA,quote=FALSE)

clean_df<-read.table(file="/mnt/smart/scratch/vir/felipe/gvSAGs/Read_Mapping_All_SAGs/SAG_Viral_Clean_DF.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,row.names=1)


###calc modules (very slow)
mod<-computeModules(clean_df,forceLPA=TRUE)

# pdf(file="/mnt/smart/scratch/vir/felipe/gvSAGs/Read_Mapping_All_SAGs/Network_Modules.pdf",width=25,height=25)
# plotModuleWeb(mod)
# dev.off()

#cant figure out how to extract module assignments from mod: https://www.rdocumentation.org/packages/bipartite/versions/2.23/topics/moduleWeb-class
lhood_val<-mod@likelihood

modules_df<-as.data.frame(mod@modules)
modorder_clean_df<-as.data.frame(mod@moduleWeb)

###Generate null models
#TODO: figure if vegan nullmodels generated with curveball cna be used for the bipartite package functions
# nullobj<-vegan::nullmodel(clean_df, method='curveball')
# vegan_null_models <-simulate(nullobj, nsim=10, seed=666)

null_models <- bipartite::nullmodel(web=clean_df, N=100, method='r2d')

#calculate null model emtrics
# H2_nulls <- sapply(null_models, function(x) H2fun(x)$H2)

#if vegan_null_models are used, gives: Error in rowSums(web) : 'x' must be an array of at least two dimensions
#H2 will be higher the more specialized the network is, i.e. in a netwrok where all virus infect all hosts equally, H2 would be 0
H2_nulls <- sapply(null_models, H2fun, H2_integer=FALSE)

#Network metrics
netlvl_metrics<-networklevel(clean_df, index=c("H2"))

#calculate the specialization level for the whole network
# h2_metric<-H2fun(clean_df,H2_integer=FALSE)

# the df containing d' values has less rows than clean_df. Some viruses are eing excluded with no warning. Perhaps this happens bc they dont have any non-zero values (i.e. are only mapped to sags for whic the taxon is NA)
d_nulls<- sapply(null_models, specieslevel, index=c("d"), level="lower")

null_metrics<-c() #as.data.frame()
# for (i in 1:ncol(d_nulls)) {
for (i in 1:length(d_nulls)) {
    sub_matrix<-as.data.frame(d_nulls[[i]])
    colnames(sub_matrix)[1]<-paste("Iteration_",i,sep="")
    null_metrics<-cbind(null_metrics,as.matrix(sub_matrix))
}

null_metrics<-as.data.frame(null_metrics)
rownames(null_metrics)<-rownames(clean_df)
null_metrics$viral_contig<-rownames(null_metrics)

mmetrics<-melt(null_metrics, id=c("viral_contig"))

summary(mmetrics)

null_d_means<-aggregate(mmetrics$value,by=list(viral_contig=mmetrics$viral_contig),FUN=mean)
null_d_sd<-aggregate(mmetrics$value,by=list(viral_contig=mmetrics$viral_contig),FUN=sd)
null_d_min<-aggregate(mmetrics$value,by=list(viral_contig=mmetrics$viral_contig),FUN=min)
null_d_max<-aggregate(mmetrics$value,by=list(viral_contig=mmetrics$viral_contig),FUN=max)

null_metrics_stats<-as.data.frame(cbind(null_d_means,null_d_sd[,2],null_d_min[,2],null_d_max[,2]))
colnames(null_metrics_stats)<-c("viral_contig","Mean","sd","Min","Max")
null_metrics_stats$viral_contig<-as.factor(null_metrics_stats$viral_contig)

#calculate the individual speclization index of each virus as defined by Bluthgen 10.1186/1472-6785-6-9
# More specialized viruses will have higher d'
#taxa that map reads to many differnt viruses will have low d' 

splvl_metrics<-specieslevel(clean_df, level="lower", index=c("d"))

splvl_metrics$viral_contig<-as.factor(rownames(splvl_metrics))

splvl_metrics<-merge(splvl_metrics,null_metrics_stats,by="viral_contig",all.x=T)

splvl_metrics$fold<-splvl_metrics$d/splvl_metrics$Mean

fold_order<-as.character(splvl_metrics$viral_contig[order(splvl_metrics$fold)])
#Below I tried to do the heatmap only for the high quality contigs , but it doe snot work because names are not matching with those in the vinfo table. Xabi's new tables will fix this
figA<-ggplot()+
#scale_y_log10()+
geom_bar(data=splvl_metrics,aes(x=viral_contig, y=d, fill=fold),position="dodge",stat="identity",colour="black",linewidth = 0)+
#geom_point(data=splvl_metrics,aes(x=viral_contig, y=Max),position="dodge",stat="identity",colour="black",alpha=0.9)+
scale_x_discrete(limits = fold_order)+
theme_bw()
#labs(title="% of total reads")+
#theme(legend.position="top",text = element_text(size = 12),axis.text.x = element_text(angle = 45,hjust = 1,size=9),legend.key.size = unit(0.5, "cm"))+facet_grid(Variable ~ . , scales="free_y")
#+scale_fill_manual(name="Taxon",values=sub_tax_coloring)
#+scale_fill_manual(name="Taxon",values=phylum_coloring)
ggsave("/mnt/smart/scratch/vir/felipe/gvSAGs/Read_Mapping_All_SAGs/d_barplots.pdf",plot=figA,device=cairo_pdf,width=45,height=10,pointsize=8)

save.image(file='/mnt/smart/scratch/vir/felipe/gvSAGs/Read_Mapping_All_SAGs/NetworkMetrics.RData')
# load('/mnt/smart/scratch/vir/felipe/gvSAGs/Read_Mapping_All_SAGs/NetworkMetrics.RData')