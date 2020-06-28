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
seurat.deg.markers$gene  <-rownames(seurat.deg.markers) # if markers are generated with FindMarkers
annotated.seurat.deg.markers<- annotate.genelist(seurat.deg.markers, tf=mm.tf, surface=mm.cellsurfacemarker, secretory=mm.secretory)
```

#### Annotate Scanpy DEG table
```python
import urllib
#response = urllib.request.urlopen("https://raw.githubusercontent.com/chansigit/SSAT/master/mm.ribo.genes.tsv")
#mm_ribo_genes = [eval(str.strip(str(line, encoding = "utf-8"))) for line in response]
mm_ribo_genes = [x for x in filter((lambda g: g.startswith("Rpl") or g.startswith("Rps")),  adata.var_names)]

response = urllib.request.urlopen("https://raw.githubusercontent.com/chansigit/SSAT/master/mm.tf.csv")
mm_tf_genes = [eval(str.strip(str(line, encoding = "utf-8"))) for line in response]
del mm_tf_genes[0]

response = urllib.request.urlopen("https://raw.githubusercontent.com/chansigit/SSAT/master/mm.secretory.csv")
mm_secretory_genes = [eval(str.strip(str(line, encoding = "utf-8"))) for line in response]
del mm_secretory_genes[0]
mm_secretory_genes= [x[1] for x in mm_secretory_genes]

response = urllib.request.urlopen("https://raw.githubusercontent.com/chansigit/SSAT/master/mm.cellsurfacemarker.csv")
mm_cellsurfacemarker_genes = [eval(str.strip(str(line, encoding = "utf-8"))) for line in response]
del mm_cellsurfacemarker_genes[0]
mm_cellsurfacemarker_genes= [x[1] for x in mm_cellsurfacemarker_genes]

####################################################################################
de_results = adata.uns['rank_genes_groups']
cluster_names = de_results['names'].dtype.names

avg_logFC_cutoff = 0.75
padj_cutoff = 0.005
filter_ribo_genes = True

de = pd.DataFrame()

for cluster in cluster_names:
    gene_column=de_results['names'][cluster]
    padj_column=de_results['pvals_adj'][cluster]
    lnfc_column=de_results['logfoldchanges'][cluster]
    de_df = pd.DataFrame({"names":de_results['names'][cluster], 
                          "cluster":str(cluster),
                          "type":"",
                          "pvals":de_results['pvals'][cluster],
                          "pvals_adj":de_results['pvals_adj'][cluster],
                          "logfoldchanges":de_results['logfoldchanges'][cluster],
                          "scores":de_results['scores'][cluster] })
    #de_df['cluster']=str(cluster)
    
    de_df = de_df[ de_df["logfoldchanges"]>avg_logFC_cutoff ]
    de_df = de_df[ de_df["pvals_adj"]<padj_cutoff ]
    de_df = de_df.sort_values('logfoldchanges',ascending=False)
    if filter_ribo_genes:
        de_df = de_df[~de_df["names"].isin( mm_ribo_genes) ]
    else:
        de_df.loc[de_df["names"].isin( mm_ribo_genes), "type"]="ribo"
    
    de_df.loc[de_df["names"].isin( mm_tf_genes), "type"]="TF"
    de_df.loc[de_df["names"].isin( mm_secretory_genes), "type"]="Secreted"
    de_df.loc[de_df["names"].isin( mm_cellsurfacemarker_genes), "type"]="Surface"
    
    #print(de_df)
    de= de.append(de_df, ignore_index=True)
    
de
```

## News

2020-03-03: Surface Marker List updated (3413 mouse genes, 2846 human genes); Secreted protein gene updated (713 human genes, 918 mouse genes). 

2019-08-01: Database established
