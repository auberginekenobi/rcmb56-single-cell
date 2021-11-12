---
#title: "CopyKAT analysis with scRNA rcmb56 ht data"
#output: html_notebook
---
  
#Step 1l installation

library(devtools)
install_github("navinlabcode/copykat")

#Step 2 preparing readcount input file
## Use the seurat .Rda object that we generated from the seurat analysis?
library(Seurat)
file <- "/Users/AirSunita/Desktop/scRNA_project/seurat_project/rcmb56_ht/filtered_feature_bc_matrix/"
raw <- Read10X(data.dir = file)
rna_counts <- raw$`Gene Expression`
raw <- CreateSeuratObject(counts = rna_counts, project = "copycat.test", min.cells = 0, min.features = 0)
exp.rawdata <- as.matrix(raw@assays$RNA@counts)
write.table(exp.rawdata, file="exp.rawdata.txt", sep="\t", quote = FALSE, row.names = TRUE)

#Step 3 running copykat

library(copykat)
copykat.test <- copykat(rawmat=exp.rawdata, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, sam.name="test", distance="euclidean", norm.cell.names="", n.cores=4,output.seg="FALSE")



