# SSAT: Secretory protein, Surface protein, And TF 

Gene annotation tool with a curated list of Secretory protein, Surface protein, And Transcription factors

---

## Use

#### Annotate Seurat DEG table



```R
# load gene list from github
load(url("https://github.com/chansigit/SSAT/raw/master/mm.cellsurfacemarker.rda"))
load(url("https://github.com/chansigit/SSAT/raw/master/mm.secretory.rda"))
load(url("https://github.com/chansigit/SSAT/raw/master/mm.tf.rda"))

# load annotation function from github
source("https://raw.github.com/chansigit/SSAT/master/annotate.genelist.R")

# use this function manually
annotated.seurat.deg.markers<- annotate.genelist(seurat.deg.markers, tf=mm.tf, surface=mm.cellsurfacemarker, secretory=mm.secretory)
```

## News

2020-03-03: Surface Marker List updated, 3413 mouse genes, 2846 human genes

2019-08-01: Database established