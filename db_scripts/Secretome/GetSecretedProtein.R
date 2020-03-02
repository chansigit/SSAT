#  ================================= Auxilary symbol conversion function =====================================
ConvertHumanGeneListToMM <- function(x){
  
  # Load human ensembl attributes
  human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  # Load mouse ensembl attributes
  mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  # Link both datasets and retrieve mouse genes from the human genes
  genes.list = biomaRt::getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows = T)
  
  # Get unique names of genes (in case gene names are duplicated)
  mouse.gene.list <- unique(genes.list[, 2])
  
  # # Print the first 6 genes found to the screen
  # print(head(mouse.gene.list))
  return(mouse.gene.list)
}


#  ========================================== Build human gene list ==========================================
## VerSeDa - human
path<-"VerSeDa/1564645844_homo_sapiens/resultsTable.csv"
table<-read.table(file=path, header = TRUE,sep=",")
hs.secretory1 <-sort(unique(as.character(table$Gene)))
library("clusterProfiler")
library("org.Hs.eg.db")
hs.secretory1 <- bitr(hs.secretory1, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")$SYMBOL

## HGNC - Chemokine and Interleukin - human
path<-"HGNC_Chemokine_Interleukin/group-483.csv"
table<-read.table(file=path, header = TRUE,sep=",")
hs.secretory2<-sort(unique(as.character(table$`Approved.symbol`)))

path<-"HGNC_Chemokine_Interleukin/group-601.csv"
table<-read.table(file=path, header = TRUE,sep=",")
hs.secretory3<-sort(unique(as.character(table$`Approved.symbol`)))

hs.secretory <- sort(unique(c(hs.secretory1, hs.secretory2, hs.secretory3))) 



#  ========================================= Build mouse gene list ==========================================
path<-"VerSeDa/1564645490_mus_musculus/resultsTable.csv"
table<-read.table(file=path, header = TRUE,sep=",")
mm.secretory1 <-sort(unique(as.character(table$Gene)))
library("clusterProfiler")
library("org.Mm.eg.db")
mm.secretory1 <- bitr(mm.secretory1, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")$SYMBOL

## HUMAN SYMBOLS
mm.secretory <- sort( unique(c(mm.secretory1, ConvertHumanGeneListToMM(hs.secretory)) ) )



#  ========================================== Output ==========================================
print(length(hs.secretory))
print(length(mm.secretory))
write.csv(hs.secretory, file="hs.secretory.csv")
write.csv(mm.secretory, file="mm.secretory.csv")
save(hs.secretory, file="hs.secretory.rda", compress=T)
save(mm.secretory, file="mm.secretory.rda", compress=T)