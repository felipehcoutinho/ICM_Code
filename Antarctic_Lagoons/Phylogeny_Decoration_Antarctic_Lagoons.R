####Decorate cas9 tree
library(ggplot2)
library(reshape2)
library(ggtree)
library(tidytree)
library(ape)
library(treeio)
library(RColorBrewer)

#https://yulab-smu.top/treedata-book/chapter7.html

info_df<-read.table(file="/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/DADA/DADA2_ASVs_Info_Antarctic_Lagoons.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

info_df$Sequence<-NULL

summary(info_df)

tree_file<-"/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/DADA/Aligned_All_ASV_MF.nwk"

tree<-ape::read.tree(tree_file)

# taxon_coloring<-c(brewer.pal(9,"Set1"),brewer.pal(4,"Pastel1"),brewer.pal(8,"Dark2")[8],brewer.pal(5,"Set2")[5],brewer.pal(8,"Accent"),brewer.pal(5,"Dark2"))
# names(taxon_coloring)<-c("Actinobacteriota","Alphaproteobacteria","Cyanobacteria","Crenarchaeota","Bacteroidia","Nanoarchaeota","Verrucomicrobiota","Marinisomatia","Myxococcota","Oligoflexia","Gammaproteobacteria","Acidobacteriota","UBP7","Poribacteria","Eremiobacterota","Margulisbacteria","Gemmatimonadota","Desulfobacterota","Nitrospinota","Desulfobacterota_D","Patescibacteria","Hydrothermarchaeota","Chlamydiota","SAR324","Halobacteriota","Chloroflexota","Rhodothermia","Phycisphaerae")

# valid_tax<-c("Alphaproteobacteria","Bacteroidia","Gammaproteobacteria","Marinisomatia","Oligoflexia","Phycisphaerae","Rhodothermia")

# sub_taxon_coloring<-taxon_coloring#taxon_coloring[valid_tax]

###No IDs
p1<-ggtree(tree,size=1,layout="unrooted") %<+% info_df+
geom_tippoint(aes(subset = isTip, colour = Phylum),size=5)+
#scale_colour_manual(name="Taxon",values=sub_taxon_coloring,breaks=c("Alphaproteobacteria","Gammaproteobacteria","Marinisomatia","Oligoflexia","Phycisphaerae","Bacteroidia","Rhodothermia"))+
theme(legend.key.size=unit(5,"cm"),legend.text = element_text(size=25),legend.title = element_blank()) 

ggsave("/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Antarctic_Lagoon/DADA/Tree_Aligned_All_ASV_MFs.pdf",plot=p1,width=65,height=65,pointsize=8,limitsize = FALSE)


# ###Rep IDs only
# p1<-ggtree(tree,size=1,layout="circular") %<+% info_df+
# geom_tippoint(aes(subset = isTip & label %in% info_df$Sequence[which(info_df$Source != "Reference")], colour = Class),size=5)+
# scale_colour_manual(name="Taxon",values=sub_taxon_coloring,breaks=c("Alphaproteobacteria","Gammaproteobacteria","Marinisomatia","Oligoflexia","Phycisphaerae","Bacteroidia","Rhodothermia"))+
# geom_tiplab(aes(subset = isTip & label %in% rep_ids, label=Makarova_Class), size=14, colour = 'firebrick', offset=0.2)+
# theme(legend.key.size=unit(5,"cm"),legend.text = element_text(size=25),legend.title = element_blank()) 

# ggsave("/mnt/lustre/scratch/nlsas/home/csic/eyg/fhc/Deep_Cas/MLP_cas9_IQTree_Midpoint_Rooted.pdf",plot=p1,width=55,height=55,pointsize=8,limitsize = FALSE)
