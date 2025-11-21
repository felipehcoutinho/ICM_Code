library(RColorBrewer)

#Phylum
phylum_coloring<-c(brewer.pal(9,"Set1"),brewer.pal(4,"Pastel1"),brewer.pal(8,"Dark2")[8],brewer.pal(5,"Set2")[5],brewer.pal(8,"Accent"),brewer.pal(5,"Dark2"))
names(phylum_coloring)<-c("Actinobacteriota","Alphaproteobacteria","Cyanobacteria","Crenarchaeota","Bacteroidota","Nanoarchaeota","Verrucomicrobiota","Marinisomatota","Myxococcota","Bdellovibrionota","Gammaproteobacteria","Acidobacteriota","UBP7","Poribacteria","Eremiobacterota","Margulisbacteria","Gemmatimonadota","Desulfobacterota","Nitrospinota","Desulfobacterota_D","Patescibacteria","Hydrothermarchaeota","Chlamydiota","SAR324","Halobacteriota","Chloroflexota","Thermoplasmatota","Planctomycetota")

custom_taxon_order<-c("Acidobacteriota","Actinobacteriota","Bacteroidota","Bdellovibrionota","Chlamydiota","Chloroflexota","Cyanobacteria","Desulfobacterota","Desulfobacterota_D","Eremiobacterota","Gemmatimonadota","Margulisbacteria","Marinisomatota","Myxococcota","Nitrospinota","Patescibacteria","Planctomycetota","Poribacteria","Alphaproteobacteria","Gammaproteobacteria","SAR324","UBP7","Verrucomicrobiota","Halobacteriota","Hydrothermarchaeota","Nanoarchaeota","Thermoplasmatota","Thermoproteota")

# 100 colours gradient fro the spectral palette
col_pal<-brewer.pal(11,"Spectral")
spectral_grad<-rev(colorRampPalette(col_pal)(n=100))