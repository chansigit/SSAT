import urllib
response = urllib.request.urlopen("https://raw.githubusercontent.com/chansigit/SSAT/master/mm.ribo.genes.tsv")
mm_ribo_genes = [eval(str.strip(str(line, encoding = "utf-8"))) for line in response]
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
