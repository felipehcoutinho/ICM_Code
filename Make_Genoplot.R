library(RColorBrewer)
library(genoPlotR)
#library(svglite)#, lib.loc="/mnt/lustre/bio/users/fcoutinho/Rlibs/")

args = commandArgs(trailingOnly=TRUE)
gbk_file<-args[1]
fig_width<-as.numeric(args[2])
fig_height<-as.numeric(args[3])

add_functional_categories<-function(segment) {
	segment$Functional_Category<-"Others"
	
	segment$Functional_Category[grepl("virion(.*)structural",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"Structural"
	segment$Functional_Category[grepl("(tail)|(fiber)",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"Structural"
	segment$Functional_Category[grepl("capsid",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"Structural"
	segment$Functional_Category[grepl("DUF3168",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"Structural"
	segment$Functional_Category[grepl("HK97(.*)gp10",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"Structural"
	segment$Functional_Category[grepl("DNA(.*)topoisomerase",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"DNA Metabolism"
	segment$Functional_Category[grepl("DNA(.*)polymerase",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"DNA Metabolism"
	segment$Functional_Category[grepl("DNA(.*)binding",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"DNA Metabolism"
	segment$Functional_Category[grepl("integrase",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"DNA Metabolism"
	segment$Functional_Category[grepl("transcription",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"DNA Metabolism"
	segment$Functional_Category[grepl("TFIIB(.)*zinc.binding(.)*protein",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"DNA Metabolism"
	segment$Functional_Category[grepl("Deoxycytidine(.)*triphosphate(.)*deaminase",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"AMG"
	segment$Functional_Category[grepl("Ribonuclease(.)*H",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"DNA Metabolism"
	segment$Functional_Category[grepl("Guanosine(.)*methyltransferase(.)*RlmB",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"DNA Metabolism"
	segment$Functional_Category[grepl("helix(.)*turn(.)*helix",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"DNA Metabolism"
	segment$Functional_Category[grepl("helicase",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"DNA Metabolism"
	segment$Functional_Category[grepl("HNH(.)*endonuclease",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"DNA Metabolism"
	segment$Functional_Category[grepl("Ftsk_gamma",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"DNA Metabolism"
	segment$Functional_Category[grepl("Chromosome(.)*partition",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"DNA Metabolism"
	segment$Functional_Category[grepl("parB",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"DNA Metabolism"
	segment$Functional_Category[grepl("zinc(.)*finger",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"DNA Metabolism"
	segment$Functional_Category[grepl("ribityllumazine(.*)synthase",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"AMG"
	segment$Functional_Category[grepl("butanone(.*)4.phosphate(.*)synthase",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"AMG"
	segment$Functional_Category[grepl("Dihydrofolate(.)*reductase",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"AMG"
	segment$Functional_Category[grepl("thymidylate(.)*(synthase|kinase)",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"AMG"
	segment$Functional_Category[grepl("Ribonucleoside(.)*diphosphate(.)*reductase",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"AMG"
	segment$Functional_Category[grepl("dUTPase",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"AMG"
	segment$Functional_Category[grepl("Deoxyuridine(.)*triphosphate(.)*nucleotidohydrolase",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"AMG"
	segment$Functional_Category[grepl("chorismate(.)*synthase",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"AMG"
	segment$Functional_Category[grepl("(terminase)|(portal)",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"Packaging & Lysis"
	segment$Functional_Category[grepl("peptidase(.)*S24",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"Packaging & Lysis"
	segment$Functional_Category[grepl("peptidoglycan(.)*hydrolase",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"Packaging & Lysis"
	segment$Functional_Category[grepl("holin",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"Packaging & Lysis"
		segment$Functional_Category[grepl("cell(.)*wall(.)*hydrolase",segment[["function"]],ignore.case = TRUE, perl = TRUE)]<-"Packaging & Lysis"
	
	return(segment)
}

g1_seg<-read_dna_seg_from_genbank(file=gbk_file,gene_type="arrows",extra_fields=c("function","taxon"))

g1_seg<-add_functional_categories(g1_seg)


color_vec<-brewer.pal(5,"Spectral")
names(color_vec)<-c("Structural","AMG","Packaging & Lysis","DNA Metabolism","Others")

mid_pos <- middle(g1_seg)
annot1 <- annotation(x1=mid_pos, text=g1_seg[["function"]])
annot1$rot<-90


g1_seg$fill<-color_vec[as.vector(g1_seg$Functional_Category)]
g1_seg$col<-"#000000"

annot1$text[which(annot1$text == "NA")]<-""


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


