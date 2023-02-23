library (readr)
library (argparser)
library(ArchR)
library(pheatmap)
library (Seurat)
library (tidyr)
library (dplyr)
library (chromVARmotifs)
library (gtools)

projHeme2 <- readRDS('/media/ps/ /outdir/0.rds/projHeme2-reDim.rds')
setwd('/media/ps/ /outdir2/')
projHeme2 <- addHarmony(
     ArchRProj = projHeme2,
     reducedDims = "IterativeLSI",
     name = "Harmony",
     groupBy = "Sample"
 )
## Dimensionality Reduction After Harmony
 projHeme2 <- addUMAP(
     ArchRProj = projHeme2, 
     reducedDims = "Harmony", 
     name = "UMAPHarmony", 
     nNeighbors = 30, 
     minDist = 0.5, 
     metric = "cosine"
 )
 p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
 plotPDF(p1, name = file.path(outdir_dim,paste0(prefix,".UMPA_Harmony_samples.pdf")), ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
 p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
 plotPDF(p2, name = file.path(outdir_dim,paste0(prefix,".UMAP_Harmony_Cluster.pdf")), ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

