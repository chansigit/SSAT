#  ========================================== Cell surface protein ==========================================
path<-"C:/Users/chans/Desktop/GeneList/CellSurfaceProtein/MusMusculus_unique_peptide_count_per_protein_and_cell_type_data.csv"
table<-read.table(file=path, header = TRUE,sep="\t")
mm.cellsurfacemarker<-sort(unique(as.character(table$ENTREZ.gene.symbol)))
mm.cellsurfacemarker<-mm.cellsurfacemarker[-c(1,2)]
mm.cellsurfacemarker


path<-"C:/Users/chans/Desktop/GeneList/CellSurfaceProtein/HomoSapiens_unique_peptide_count_per_protein_and_cell_type_data.csv"
table<-read.table(file=path, header = TRUE,sep=",")
hs.cellsurfacemarker<-sort(unique(as.character(table$ENTREZ.gene.symbol)))
hs.cellsurfacemarker<-hs.cellsurfacemarker[-c(1)]
hs.cellsurfacemarker

#  ========================================== Secretory Protein ==========================================
path<-"C:/Users/chans/Desktop/GeneList/Secretome/1564645490_mus_musculus/resultsTable.csv"
table<-read.table(file=path, header = TRUE,sep=",")
mm.secretory <-sort(unique(as.character(table$Gene)))
library("clusterProfiler")
library("org.Mm.eg.db")
mm.secretory <- bitr(mm.secretory, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")$SYMBOL

path<-"C:/Users/chans/Desktop/GeneList/Secretome/1564645844_homo_sapiens/resultsTable.csv"
table<-read.table(file=path, header = TRUE,sep=",")
hs.secretory <-sort(unique(as.character(table$Gene)))
library("clusterProfiler")
library("org.Hs.eg.db")
hs.secretory <- bitr(hs.secretory, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")$SYMBOL

#  ========================================== TF ==========================================
path<-"C:/Users/chans/Desktop/GeneList/TF/Mus_musculus_TF.txt"
table<-read.table(file=path, header = TRUE,sep="\t")
mm.tf <-sort(unique(as.character(table$Symbol)))
mm.tf

path<-"C:/Users/chans/Desktop/GeneList/TF/Homo_sapiens_TF.txt"
table<-read.table(file=path, header = TRUE,sep="\t")
hs.tf <-sort(unique(as.character(table$Symbol)))
hs.tf



write.table(mm.tf,                sep="\t" ,file="mm.tf.csv")
write.table(hs.tf,                sep="\t" ,file="hs.tf.csv")
write.table(mm.secretory,         sep="\t" ,file="mm.secretory.csv")
write.table(hs.secretory,         sep="\t" ,file="hs.secretory.csv")
write.table(mm.cellsurfacemarker, sep="\t" ,file="mm.cellsurfacemarker.csv")
write.table(hs.cellsurfacemarker, sep="\t" ,file="hs.cellsurfacemarker.csv")

save(mm.tf,                 file="mm.tf.rda",        compress = TRUE)
save(hs.tf,                 file="hs.tf.rda",        compress = TRUE)
save(mm.secretory,          file="mm.secretory.rda", compress = TRUE)
save(hs.secretory,          file="hs.secretory.rda", compress = TRUE)
save(mm.cellsurfacemarker,  file="mm.cellsurfacemarker.rda", compress = TRUE)
save(hs.cellsurfacemarker,  file="hs.cellsurfacemarker.rda", compress = TRUE)
