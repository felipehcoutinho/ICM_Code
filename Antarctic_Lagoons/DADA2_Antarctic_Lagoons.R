library(dada2)

#Declare the "path" variable with the directory containing the PostQC fastqc files
path <- "/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/Reads/Clean/"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- gsub("Clean_","",sapply(strsplit(basename(fnFs), "_L001"), `[`, 1))

# Inspect read quality profiles
pdf("/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/DADA/DADA2_Quality_Profiles_Antarctic_Lagoons_Forward.pdf")
plotQualityProfile(fnFs[1:2])
dev.off()

#Now we visualize the quality profile of the reverse reads:
pdf("/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/DADA/DADA2_Quality_Profiles_Antarctic_Lagoons_Reverse.pdf")
plotQualityProfile(fnRs[1:2])
dev.off()

# Place filtered files in filtered/ subdirectory
filtFs <- file.path("/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/DADA", "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/DADA", "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#We’ll use standard filtering parameters: maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores.
#Note: omitted the  truncLen=c(240,160) from the tutorial from the command below
#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=10)
#head(out)

#Learn the error rates
errF <- learnErrors(filtFs, multithread=47)
errR <- learnErrors(filtRs, multithread=47)

#It is always worthwhile, as a sanity check if nothing else, to visualize the estimated error rates:
pdf("/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/DADA/DADA2_Error_Rates_Antarctic_Lagoons.pdf")
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()

#We are now ready to apply the core sample inference algorithm to the filtered and trimmed sequence data.
dadaFs <- dada(filtFs, err=errF, multithread=47)
dadaRs <- dada(filtRs, err=errR, multithread=47)

#Inspecting the returned dada-class object:
dadaFs[[1]]

#Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#We can now construct an amplicon sequence variant table (ASV) table, a higher-resolution version of the OTU table produced by traditional methods.
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=47, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

taxa <- assignTaxonomy(seqtab.nochim, "/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Databases/DADA_DBs/silva_nr99_v138.1_train_set.fa.gz", multithread=47)

taxa <- addSpecies(taxa,"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Databases/DADA_DBs/silva_species_assignment_v138.1.fa.gz")

save.image(file ="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/DADA/DADA2_Antarctic_Lagoons.RData")

#load(file ="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/DADA/DADA2_Antarctic_Lagoons.RData")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

asv_info_df<-cbind(colnames(seqtab.nochim),paste("ASV_" ,seq(1:length(colnames(seqtab.nochim))),sep=""),taxa.print)
colnames(asv_info_df)[c(1,2)]<-c("Sequence","ASV_UID")
asv_info_df<-as.data.frame(asv_info_df)
summary(asv_info_df[,-1])
write.table(asv_info_df, file = "/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/DADA/DADA2_ASVs_Info_Antarctic_Lagoons.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

asv_abd_df<-seqtab.nochim
colnames(asv_abd_df)<-asv_info_df$ASV_UID
write.table(asv_abd_df, file = "/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/DADA/DADA2_ASVs_Raw_Abundances_Antarctic_Lagoons.tsv", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

#Calculate percente abundances per row
asv_perc_abd_df<-(asv_abd_df/rowSums(asv_abd_df))*100
write.table(asv_perc_abd_df, file = "/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/DADA/DADA2_ASVs_Percentage_Abundances_Antarctic_Lagoons.tsv", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

asv_perc_abd_df<-read.table(file = "/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/DADA/DADA2_ASVs_Percentage_Abundances_Antarctic_Lagoons.tsv", sep = "\t", header=TRUE)
###Phylum barplot
colnames(asv_perc_abd_df)[1]<-"Sample"
m_asv_perc_abd_df<-melt(asv_perc_abd_df,id="Sample",value.name="Abundance",variable.name="ASV_UID")
m_asv_perc_abd_df<-merge(m_asv_perc_abd_df,asv_info_df,by="ASV_UID")

summary(m_asv_perc_abd_df)

tax_abd_df<-aggregate(m_asv_perc_abd_df$Abundance, by=list(Sample=m_asv_perc_abd_df$Sample, Taxon=m_asv_perc_abd_df$Phylum), FUN=sum)
colnames(tax_abd_df)<-c("Sample","Taxon","Abundance")

summary(tax_abd_df)

#Calculate medin abundances per taxon
tax_median_abds<-aggregate(tax_abd_df$Abundance, by=list(Taxon=tax_abd_df$Taxon), FUN=median)
summary(tax_median_abds)
valid_tax<-tax_median_abds$Taxon[which(tax_median_abds$x >= 0.5)]

#Recalculate taxon abundances sums grouping small abundance phyla into "others"
tax_level<-"Phylum"
m_asv_perc_abd_df[[tax_level]]<-as.character(m_asv_perc_abd_df[[tax_level]])
m_asv_perc_abd_df[[tax_level]][!(m_asv_perc_abd_df[[tax_level]] %in% valid_tax)]<-"Others"
m_asv_perc_abd_df[[tax_level]]<-as.factor(m_asv_perc_abd_df[[tax_level]])
summary(m_asv_perc_abd_df)

tax_abd_df<-aggregate(m_asv_perc_abd_df$Abundance, by=list(Sample=m_asv_perc_abd_df$Sample, Taxon=m_asv_perc_abd_df$Phylum), FUN=sum)
colnames(tax_abd_df)<-c("Sample","Taxon","Abundance")

summary(tax_abd_df)

#Now Add core info to df
core_id<-data.frame(do.call('rbind', strsplit(as.character(tax_abd_df$Sample),'_',fixed=TRUE)))
colnames(core_id)<-c("Core","Depth_cm","Replicate")

tax_abd_df<-as.data.frame(cbind(tax_abd_df,core_id))

tax_abd_df$Depth_cm<-as.numeric(tax_abd_df$Depth_cm)*-1

summary(tax_abd_df)

#Make barplot
fig2B<-ggplot(tax_abd_df,aes(y=Abundance,x=Depth_cm,fill=Taxon))+geom_bar(position="stack",stat="identity",alpha=0.9,colour="black")+theme_bw()+coord_flip()+facet_wrap(Core ~.)

ggsave("Phylum_Barplots_Antarctic_Lagoons.pdf",plot=fig2B,device=cairo_pdf,width=15,height=10,pointsize=8)



###NMDS
dist_metric<-"bray"
dists<-vegdist(asv_perc_abd_df, method = dist_metric)

set.seed(666)
mdsresult<-metaMDS(dists,distance = dist_metric,k = 2,maxit = 999)
data.scores<-as.data.frame(scores(mdsresult))
data.scores$Sample<-rownames(data.scores)

mdata<-data.scores

core_id<-data.frame(do.call('rbind', strsplit(as.character(mdata$Sample),'_',fixed=TRUE)))
colnames(core_id)<-c("Core","Depth_cm","Replicate")
mdata<-as.data.frame(cbind(mdata,core_id))
summary(mdata)

library(RColorBrewer)
depth_col_pal<-brewer.pal(9,"RdBu")
depth_col_grad<-rev(colorRampPalette(depth_col_pal)(n=299))

mdata$Depth_cm<-as.numeric(as.character(mdata$Depth_cm))

figX<-ggplot(mdata, aes(x=NMDS1,y=NMDS2))+geom_point(size=3.5, aes(colour=Depth_cm, shape=Core))+theme_bw()+theme(text=element_text(size=16))+scale_colour_gradientn(colours = depth_col_grad, name="Depth",trans = 'reverse')
#+scale_fill_gradientn(colours = depth_col_grad, name="Depth",trans = 'reverse') ,limits = c(4000, 0),

ggsave("NMDS_Antarctic_Lagoons.pdf",plot=figX,device=cairo_pdf,width=6,height=5,pointsize=8)

