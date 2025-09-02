library(RColorBrewer)
library(genoPlotR)
#library(svglite)#, lib.loc="/mnt/lustre/bio/users/fcoutinho/Rlibs/")


# gbk_file<-"/mnt/smart/scratch/vir/felipe/OctoMicro/Viruses/Annotation/Prodigal-gv/OctoMicro_HQ_Candidates_phold/Individual_Genbanks/A13_Scaffold_32211||full.genbank"
# fig_width<-15
# fig_height<-7

args = commandArgs(trailingOnly=TRUE)
gbk_file<-args[1]
fig_width<-as.numeric(args[2])
fig_height<-as.numeric(args[3])

g1_seg<-read_dna_seg_from_genbank(file=gbk_file,gene_type="arrows",extra_fields=c("function","taxon"))

# g1_seg<-add_functional_categories(g1_seg)

color_vec<-c(brewer.pal(8,"Spectral"),brewer.pal(9,"Greys")[7],brewer.pal(9,"Greys")[7])

names(color_vec)<-c("DNA, RNA and nucleotide metabolism","integration and excision","transcription regulation","moron, auxiliary metabolic gene and host takeover","head and packaging","connector","tail","lysis","other","unknown function")

mid_pos <- middle(g1_seg)
annot1 <- annotation(x1=mid_pos, text=g1_seg[["product"]])
annot1$rot<-60


g1_seg$fill<-color_vec[as.vector(g1_seg[["function"]])]
g1_seg$col<-"#000000"

annot1$text[which(annot1$text == "NA")]<-""
annot1$text[which(annot1$text == "hypothetical protein")]<-""

out_file_name<-gbk_file

out_file_name<-gsub(".*\\/","GenomePlot_",out_file_name,perl=TRUE)
out_file_name<-gsub("\\.gbk$","",out_file_name,perl=TRUE)
out_file_name<-gsub("\\W","_",out_file_name,perl=TRUE)
out_file_name<-paste(out_file_name,".svg",sep="")


#svglite(filename=out_file_name,width=fig_width,height=fig_height,pointsize=8)
#plot_gene_map(dna_segs=list(g1_seg), annotations=annot1, annotation_cex=0.8)
#dev.off()

svg(out_file_name,width=fig_width,height=fig_height,pointsize=8)
plot_gene_map(dna_segs=list(g1_seg), annotations=annot1, annotation_cex=0.8)
dev.off()


