#load(url("https://github.com/chansigit/SSAT/raw/master/mm.cellsurfacemarker.rda"))
#load(url("https://github.com/chansigit/SSAT/raw/master/mm.secretory.rda"))
#load(url("https://github.com/chansigit/SSAT/raw/master/mm.tf.rda"))

annotate.genelist<- function(markers, tf, surface, secretory, ...){
    # annotate secretory protein
    is.secretory <-as.character(markers$gene %in% secretory)
    is.secretory <-replace(is.secretory, is.secretory=="TRUE",  "Secretory")
    is.secretory <-replace(is.secretory, is.secretory=="FALSE", "")
    
    # annotate surface protein
    is.surface <-as.character(markers$gene %in% surface)
    is.surface <-replace(is.surface, is.surface=="TRUE",  "Surface")
    is.surface <-replace(is.surface, is.surface=="FALSE", "")
    
    # annotate tf
    is.tf <-as.character(markers$gene %in% tf)
    is.tf <-replace(is.tf, is.tf=="TRUE",  "TF")
    is.tf <-replace(is.tf, is.tf=="FALSE", "")
    
    markers$is.secretory <- is.secretory
    markers$is.surface   <- is.surface
    markers$is.tf        <- is.tf
    
    markers
}

# use this function manually
## annotate.genelist(seurat.deg.markers, tf=mm.tf, surface=mm.cellsurfacemarker, secretory=mm.secretory)

# load this function from github
# source("https://raw.github.com/chansigit/SSAT/master/annotate.genelist.R")
