#!/usr/bin/env Rscript

#Modules avalabile in : conda activate /mnt/smart/scratch/vir/felipe/envs/basic/
####Tree decoration script
library(ggplot2)
library(ggtree)
library(tidytree)
library(ape)
library(treeio)
library(tidyr)
library(ggnewscale)
library(RColorBrewer)
library(stringr)


###Part I: read in the tree file and create the leaf info table
args <- commandArgs(trailingOnly = TRUE)
tree_file<-args[1]


load_paired_data<-function(KO) {
  #Read in the tree file from command line arguments, otherwise, manually replac ein the script to use a specific file
  #
  tree_file<-paste("/mnt/smart/scratch/vir/felipe/TARA_Polar/MV52/Aligned_Merged_NR_Matched_",KO,".faa.newick",sep="")

  tree<-ape::read.tree(tree_file)

  leaf_info_file<-paste("/mnt/smart/scratch/vir/felipe/TARA_Polar/MV52/Aligned_Merged_NR_Matched_",KO,"_Tree_Seq_Info.tsv",sep="")
  leaf_info<-read.table(file=leaf_info_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
  rownames(leaf_info)<-leaf_info$Sequence

  return(list(tree, tree_file,leaf_info))
}

out_prefix<-gsub("(^.*\\/)|(\\..+$)","",tree_file,perl=TRUE)

tree<-ape::read.tree(tree_file)
#Creat a data frame that will hold tree leaf info
leaf_info<-data.frame(tree$tip.label)
colnames(leaf_info)<-"Sequence"

#leaf_info$Alias3<-get_alias(leaf_info$Sequence)
#write.table(leaf_info,file="/mnt/smart/scratch/vir/felipe/TARA_Polar/MV52/Phylo_Uniref50/K05371_Tree_Seq_Info.tsv",sep="\t",quote=FALSE,row.names=FALSE)

#Renaming needs to be done. it seems ggtree cannot handle sequences with "|" in the name (as in IMGVR sequences)
leaf_info$Alias<-leaf_info$Sequence
#leaf_info$Alias[grepl("IMGVR",leaf_info$Sequence)]<-gsub("_\\d+$","",leaf_info$Alias[grepl("IMGVR",leaf_info$Alias)],perl=TRUE)
leaf_info$Alias[grepl("IMGVR",leaf_info$Alias)]<-gsub("\\|.+$","",leaf_info$Alias[grepl("IMGVR",leaf_info$Alias)],perl=TRUE)

leaf_info$Alias2<-leaf_info$Sequence

leaf_info$Alias2<-gsub("\\|.+_","_",leaf_info$Sequence,perl=TRUE)

leaf_info$Sequence<-as.factor(leaf_info$Sequence)
leaf_info$Alias<-as.factor(leaf_info$Alias)
leaf_info$Alias2<-as.factor(leaf_info$Alias2)

tree<-treeio::rename_taxa(tree,leaf_info,"Sequence","Alias2")

summary(leaf_info)

#This dataframe holds the imgvr and unrief lineage information
# full_info<-read.table(file="/mnt/smart/scratch/vir/felipe/Databases/MV52_PhyloDB/Minimal_Tree_Deco_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
# full_info<-read.table(file="/mnt/smart/scratch/vir/felipe/Databases/MV52_PhyloDB/Swiss_Prot_Set/Swiss_Prot_Set_Minimal_tree_Deco_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
full_info<-read.table(file="/mnt/smart/scratch/vir/felipe/Databases/MV52_PhyloDB/Swiss_Prot_Set/Swiss_Prot_Set+Extra_Seqs_Minimal_tree_Deco_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

full_info$Sequence<-as.factor(gsub("\\s.+$","",full_info$Sequence,perl=T))

#Merge with the Uniref lineage information
leaf_info<-merge(leaf_info,full_info,by.x="Alias",by.y="Sequence",all.x=TRUE)


leaf_info<-leaf_info %>% separate(Ecosystem.classification, c("Ecosystem_Level_1","Ecosystem_Level_2","Ecosystem_Level_3","Ecosystem_Level_4"),sep=";")
leaf_info$Ecosystem_Level_1<-as.factor(leaf_info$Ecosystem_Level_1)
leaf_info$Ecosystem_Level_2<-as.factor(leaf_info$Ecosystem_Level_2)
leaf_info$Ecosystem_Level_3<-as.factor(leaf_info$Ecosystem_Level_3)
leaf_info$Ecosystem_Level_4<-as.factor(leaf_info$Ecosystem_Level_4)


#Match sequences to selected taxa 
#If adding new taxa, they must follow the orger from most generic to most specific
tax_to_match<-c("Eukaryota","Archaea","Bacteria","Viruses","Alphaproteobacteria","Gammaproteobacteria","Actinomycetota","Cyanobacteriota","Bacteroidota","Marinisomatota")

leaf_info$Taxon<-NA

for (taxon in tax_to_match) {
  leaf_info$Taxon[grepl(taxon,leaf_info$Full_Lineage)]<-taxon
}

#Match sequences to selected host taxa 
leaf_info$Host_Taxon<-NA

hosts_to_match<-c("Alphaproteobacteria","Gammaproteobacteria","Actinomycetota","Cyanobacteriota","Bacteroidota","Marinisomatota")

for (taxon in hosts_to_match) {
  leaf_info$Host_Taxon[grepl(taxon,leaf_info$Full_Host)]<-taxon
}

leaf_info$Host_Taxon<-as.factor(leaf_info$Host_Taxon)
leaf_info$Taxon<-as.factor(leaf_info$Taxon)

#rownames(leaf_info)<-leaf_info$Sequence
rownames(leaf_info)<-leaf_info$Alias2


summary(leaf_info)

#write the leaf info to keep trakc of any errors
out_info_name<-paste0(out_prefix,"_Tree_Seq_Info.tsv")
write.table(leaf_info,file=out_info_name,sep="\t",quote=FALSE,row.names=FALSE)

###Part II: Start here if you already have the leaf info table and the tree file
#if you already have the leaf info table, you can read it in instead of creating it from scratch and start the script from here

# leaf_info<-read.table(file="/mnt/smart/scratch/vir/felipe/TARA_Polar/MV52/Aligned_Merged_NR_Matched_K05371_Tree_Seq_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
#leaf_info<-read.table(file="/mnt/smart/scratch/vir/felipe/TARA_Polar/MV52/Aligned_Merged_NR_Matched_K00228_Tree_Seq_Info.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
# leaf_info_file<-paste("/mnt/smart/scratch/vir/felipe/TARA_Polar/MV52/Aligned_Merged_NR_Matched_",KO,"_Tree_Seq_Info.tsv",sep="")
# leaf_info<-read.table(file=leaf_info_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
#rownames(leaf_info)<-leaf_info$Sequence

#Read tree file manually
# tree_file<-"/mnt/smart/scratch/vir/felipe/TARA_Polar/MV52/Trees_Jorge/Aligned_Merged_NR_Matched_K05371.faa.newick"
# tree<-ape::read.tree(tree_file)

###The ids in highL-ids will be highlighted in the tree with the firebrick color
highl_ids<-leaf_info$Sequence[grepl("MV_52",leaf_info$Sequence)]

tax_to_match<-c("Eukaryota","Archaea","Bacteria","Viruses","Alphaproteobacteria","Gammaproteobacteria","Actinomycetota","Cyanobacteriota","Bacteroidota","Marinisomatota")

tax_colors<-c(brewer.pal(n = 11, name = "Spectral")[c(1,3,9,11)],brewer.pal(n = 12, name = "Paired")[c(1,2,3,4,7,8)])
names(tax_colors)<-tax_to_match

eco_colors<-c(brewer.pal(n = 12, name = "Paired")[c(1,2,9,8,4,12)])
names(eco_colors)<-c("Freshwater","Marine","Non-marine Saline and Alkaline","Digestive system","Deep subsurface","Soil")

###IMPORTANT: ggtree will use the FIRST column of the leaf_info df as the identifier of the sequences in the tree. Not the rownames and not not anything else. As of now, that should be the Alias2 column

###Command below plots the trees
#Make unrooted tree preserving branch lengths
p0<-ggtree(tree,size=0.1, layout='equal_angle', aes(color=Taxon),linewidth=3,alpha=0.8) %<+% leaf_info[,-c(1,2)]+
scale_colour_manual(values=tax_colors, name="Taxon")+
geom_tippoint(aes(subset = !is.na(Host_Taxon), fill = Host_Taxon),size=7,shape = 21, na.rm = TRUE)+
scale_fill_manual(values=tax_colors, name="Host Taxon")+
#geom_point(aes(fill = Ecosystem_Level_3),size=5,shape = 21, na.rm = TRUE)+
geom_tiplab(aes(subset = isTip & label %in% highl_ids, label="MV52"), size=6, colour = 'firebrick', offset=0)+
guides(colour = guide_legend(title = "Taxon", title.theme = element_text(size = 50, face = "bold"),label.theme = element_text(size = 45), override.aes = list(size = 50, linewidth=5,starshape = NA, starstroke = 0.15), order = 1),fill = guide_legend(title = "Virus Host", title.theme = element_text(size = 50, face = "bold"),label.theme = element_text(size = 45),override.aes = list(starshape = 21, size = 10, alpha = 1, starstroke = .15), order = 2))


out2_name<-paste0(out_prefix,"_Unrooted_Decorated_Tree.pdf")
ggsave(out2_name,plot=p0,width=55,height=55,pointsize=8,limitsize = FALSE)



#### Make a rectancular tree with branch lenghts and IDs for all seqs
p1<-ggtree(tree,size=1, layout='rectangular',branch.length='none' , aes(color=Taxon)) %<+% leaf_info[,-c(1,2)]+ #layout="equal_angle"
geom_tippoint(aes(subset = !is.na(Host_Taxon), fill = Host_Taxon),size=5,shape = 21, na.rm = TRUE)+
#geom_tippoint(aes(subset = !is.na(Ecosystem_Level_3), fill = Ecosystem_Level_3),size=5,shape = 21, na.rm = TRUE)+
guides(colour = guide_legend(title = "Taxon", title.theme = element_text(size = 50, face = "bold"),label.theme = element_text(size = 45), override.aes = list(size = 50, linewidth=5,starshape = NA, starstroke = 0.15), order = 1),fill = guide_legend(title = "Virus Host", title.theme = element_text(size = 50, face = "bold"),label.theme = element_text(size = 45),override.aes = list(starshape = 21, size = 10, alpha = 1, starstroke = .15), order = 2))+

scale_colour_manual(values=tax_colors, name="Taxon")+
scale_fill_manual(values=tax_colors, name="Host Taxon")+
theme(legend.key.size=unit(5,"cm"),legend.text = element_text(size=25),legend.title = element_blank()) 


p3 <- p1 + new_scale_fill()

p4<-gheatmap(p3, leaf_info[c("Ecosystem_Level_3")], offset=3, width=0.1,colnames_angle=0, colnames_offset_y = .25)+
scale_fill_manual(values=eco_colors, name="Ecosystem")

p5 <- p4+geom_tiplab( size=3, colour = 'black', offset=0.1)

#Select rectangular tree hight based on the number of leaves
rech<-nrow(leaf_info)/15
if (rech < 30) {
  rech<-30
}

out1_name<-paste0(out_prefix,"_Decorated_Tree_Rectangular.pdf")
ggsave(out1_name,plot=p5,width=50,height=50,pointsize=8,limitsize = FALSE)


# ###Make a circular tree with no branch lenght
# p0<-ggtree(tree,size=1, layout='circular', branch.length='none',aes(color=Taxon)) %<+% leaf_info[,-1]+
# scale_colour_manual(values=tax_colors, name="Taxon")+
# geom_tippoint(aes(subset = !is.na(Host_Taxon), fill = Host_Taxon),size=5,shape = 21, na.rm = TRUE)+
# scale_fill_manual(values=tax_colors, name="Host Taxon")+
# geom_tiplab(aes(subset = isTip & label %in% highl_ids, label="MV52"), size=9, colour = 'firebrick', offset=0.5)

# # p2<-gheatmap(p0, subset(leaf_info,select=c("Taxon"))) #+
# #p2<-gheatmap(p0, leaf_info[c("Taxon")], offset=10, width=0.1)+
# # p2<-gheatmap(p0, leaf_info[c("Host_Taxon")], offset=0.5, width=0.1, colnames_angle=0, colnames_offset_y = .25)+
# # scale_fill_manual(values=tax_colors, name="Host Taxon") #+
# # #+scale_fill_viridis_d(option="D", name="Virus Host Taxon")+

# # p3 <- p2 + new_scale_fill()

# # eco_colors<-c(brewer.pal(n = 7, name = "Set2"),brewer.pal(n = 9, name = "Set1"),brewer.pal(n = 8, name = "Dark2"))[1:length(unique(leaf_info$Ecosystem_Level_3))]
# # names(eco_colors)<-unique(leaf_info$Ecosystem_Level_3)
# # eco_colors<-eco_colors[!is.na(names(eco_colors))]

# eco_colors<-c(brewer.pal(n = 12, name = "Paired")[c(1,2,9,8,4,12)])

# names(eco_colors)<-c("Freshwater","Marine","Non-marine Saline and Alkaline","Digestive system","Deep subsurface","Soil")


# p3 <- p0 + new_scale_fill()

# p4<-gheatmap(p3, leaf_info[c("Ecosystem_Level_3")], offset=1, width=0.1,colnames_angle=0, colnames_offset_y = .25)+
# scale_fill_manual(values=eco_colors, name="Ecosystem")
# #+scale_fill_viridis_d(option="C", name="Ecosystem") 
# #trying to use a second variable in the heatmpa also only works for the Taxon column
# #p2<-gheatmap(p0, leaf_info[c("Taxon","Host_Taxon","Ecosystem_Level_3")], offset=10, width=0.1) #+

# #subset(leaf_info,select=c("Ecosystem_Level_3")), aes(subset = !is.na(Ecosystem_Level_3))
# #subset(leaf_info,select=c("Taxon"))
# out3_name<-paste0(out_prefix,"_Decorated_Tree_Circular.pdf")

# ideal_dim<-nrow(leaf_info)/7

# ggsave(out3_name,plot=p4,width=55,height=55,pointsize=8,limitsize = FALSE)



