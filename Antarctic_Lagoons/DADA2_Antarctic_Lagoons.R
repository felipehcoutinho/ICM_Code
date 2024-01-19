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