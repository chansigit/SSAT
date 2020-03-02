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
## Surfaceome - human
path<-"Surfaceome_human/table_S3_surfaceome.csv"
table<-read.table(file=path, header = TRUE,sep=",")
hs.cellsurfacemarker1<-sort(unique(as.character(table$Uniprot.gene)))
#hs.cellsurfacemarker1

## HGNC CD molecules- human
path<-"HGNC_CDMolecules/group-471.csv"
table<-read.table(file=path, header = TRUE,sep=",")
hs.cellsurfacemarker2<-sort(unique(as.character(table$`Approved.symbol`)))
#hs.cellsurfacemarker2

## HUMAN SYMBOLS
hs.cellsurfacemarker <- sort( unique(c(hs.cellsurfacemarker1, hs.cellsurfacemarker2) ) )



#  ========================================== Build mosue gene list ==========================================
## CellSurfaceProteinAtlas - mouse
path<-"CellSurfaceProteinAtlas/MusMusculus_unique_peptide_count_per_protein_and_cell_type_data.csv"
table<-read.table(file=path, header = TRUE,sep="\t")
glist<-sort(unique(as.character(table$ENTREZ.gene.symbol)))
mm.cellsurfacemarker1<-glist[-c(1,2)]
#mm.cellsurfacemarker1

## MOUSE SYMBOLS
mm.cellsurfacemarker <- sort(unique(c(mm.cellsurfacemarker1, ConvertHumanGeneListToMM(hs.cellsurfacemarker))))



#  ========================================== Output ==========================================
print(length(hs.cellsurfacemarker))
print(length(mm.cellsurfacemarker))
write.csv(hs.cellsurfacemarker, file="hs.cellsurfacemarker.csv")
write.csv(mm.cellsurfacemarker, file="mm.cellsurfacemarker.csv")
save(hs.cellsurfacemarker, file="hs.cellsurfacemarker.rda", compress=T)
save(mm.cellsurfacemarker, file="mm.cellsurfacemarker.rda", compress=T)