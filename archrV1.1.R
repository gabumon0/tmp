library (readr)
library (argparser)
library(ArchR)
library(pheatmap)
library (Seurat)
library (tidyr)
library (dplyr)
library (chromVARmotifs)
library (gtools)


## 物种还可以是果蝇
argv <- arg_parser('')
argv <- add_argument(argv,"--fragment", help="the cell gene express file list, split by ,")
argv <- add_argument (argv,"--spname",help="the samples name list ,split by ,")
argv <- add_argument (argv,"--threads",help="threads used",default = 1)
argv <- add_argument (argv,"--species",help="hg19 or hg38 or mm10 or mm9",default = "hg38")
argv <- add_argument (argv,"--outdir",help="path of outdir")
argv <- add_argument (argv,"--prefix",help="output file prefix")
argv <- add_argument (argv,"--method",help="ArchR cluster method: Seurat or scran,default = Seurat",default="Seurat")
argv <- add_argument (argv,"--resolution", help="the cluster resolution,default=2.5",default=2.5)
argv <- add_argument(argv,"--dim", help="use cluster number 1:30",default=30)
argv <- add_argument (argv,"--scRNA", help ="path of scRNA rds, or no",default = "no")
argv <- add_argument (argv,"--testMethod",help = "wilcoxon, ttest or binomial",default = "wilcoxon")
argv <- add_argument (argv,"--peakcells", help="call peaks for all cluster with more than 50 cells",default=50)
argv <- add_argument (argv,"--filterdoublet",help = "filter doublet or not",default = "F")
argv <- add_argument(argv,"--amuletFile",help="The full path of MultipletBarcodes_01.txt for each sample, file order should match the order of inSample, split by ,", default = "F")
#argv <- add_argument (argv,"--object",help="1:basic analysis,2:Dimensionality Reduction,3:scRNA-seq based annotation,4:Identify peaks and plot marker peak,5:Motif analysis,6:Co-accessibility and Peak2GeneLinkage,7:Trajectory Analysis")



argv <- parse_args(argv)
if (argv$amuletFile == "T") {
    amulet_File<-unlist(strsplit(argv$amuletFile,split=","))
}

outdir <- argv$outdir
dir.create (outdir,showWarnings=T)
method <- argv$method
resol <- as.numeric(argv$resolution)
cnum <- as.numeric(argv$dim)
scRNA <- argv$scRNA
testmethod <- argv$testMethod
prefix <- argv$prefix
peakfilter <- as.numeric(argv$peakcells)
outrds <- file.path (outdir,"0.rds")
dir.create (outrds,showWarnings = T)
filterdoublet <- argv$filterdoublet
set.seed (1111)
#obj <- as.numeric(argv$object)
# wd <- getwd()
## Setting a Genome and GeneAnnotation
if (argv$species == "hg38") {
    addArchRGenome("hg38")
} else if (argv$species == "hg19") {
    addArchRGenome("hg19")
} else if (argv$species == "mm10") {
    addArchRGenome("mm10")
} else if (argv$species == "mm9") {
    addArchRGenome("mm9")
}




# ### ignore the chromosome prefixes
# #addArchRChrPrefix(chrPrefix = FALSE)

# ## Setting default number of Parallel threads.
addArchRThreads(threads = as.numeric(argv$threads))   
print ("Creating Arrow Files")
#### Creating Arrow Files
### Per-cell Quality Control (TSS enrichment score greater than 4 and more than 1000 unique nuclear fragments)
inputFiles <- c()
if (grepl(',',argv$fragment)) {
    inputFiles <- unlist(strsplit(argv$fragment,split=","))
}else {
    inputFiles <- c(argv$fragment)
}
samplename <- c()
if (grepl(",",argv$spname)) {
    samplename <- unlist(strsplit(argv$spname,split=','))
}else {
    samplename <- c(argv$spname)
}

names (inputFiles) <- samplename
print (inputFiles)

ArrowFiles <- createArrowFiles(
inputFiles = inputFiles,
sampleNames = names(inputFiles),
filterTSS = 4, #Dont set this too high because you can always increase later
filterFrags = 1000, 
addTileMat = TRUE,
addGeneScoreMat = TRUE
)
ArrowFiles

## Inferring scATAC-seq Doublets (数据需不需要去doublet需要根据R2 确认)
### This is likely a good place to start but we encourage all users to inspect the pre- and post-doublet removal data to understand how doublet removal is affecting the cells
### ArchR reports the R2 value for the UMAP projection for each Arrow file.
### If these R2 values are much lower (i.e. less than 0.9), skipping doublet prediction.
doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)

### Creating an ArchRProject
projHeme1 <- ArchRProject(
ArrowFiles = ArrowFiles, 
outputDirectory = outdir,
copyArrows = TRUE  #This is recommened so that if you modify the Arrow files you have an original copy for later usage.##看一下修改为F结果是怎么样的
)   
projHeme1
#### which data matrices are available within the ArchRProject
getAvailableMatrices(projHeme1)
###########画图
### Plotting QC metrics
#### log10(Unique Fragments) vs TSS enrichment score
df <- getCellColData(projHeme1, select = c("log10(nFrags)", "TSSEnrichment"))
p <- ggPoint(
    x = df[,1], 
    y = df[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
#plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = projHeme1, addDOC = FALSE)
dir.create (file.path(outdir,"Plots"),showWarnings=T) 
dir.create (file.path(outdir,"Plots","0.QCmetrics"),showWarnings=T) 
dip <- file.path (outdir,"Plots","0.QCmetrics")
plotPDF(p, name = file.path("0.QCmetrics","TSS-vs-Frags.pdf"), ArchRProj = projHeme1, addDOC = FALSE)
png(file.path(dip,"TSS-vs-Frags.png"))
p
dev.off()
# #### Make a ridge plot for each sample for the TSS enrichment scores.
# p1 <- plotGroups(
#     ArchRProj = projHeme1, 
#     groupBy = "Sample", 
#     colorBy = "cellColData", 
#     name = "TSSEnrichment",
#     plotAs = "ridges"
# )
# plotPDF(p1, name = file.path("0.QCmetrics","sample-TSS-ridgePlot.pd"), ArchRProj = projHeme1, addDOC = FALSE,width=5,height=5)

#### Make a violin plot for each sample for the TSS enrichment scores
p2 <- plotGroups(
    ArchRProj = projHeme1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
)
plotPDF(p2, name = file.path("0.QCmetrics","sample-TSS-VlnPlot.pdf"), ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)
png(file.path(dip,"sample-TSS-VlnPlot.png"))
p2
dev.off()
#### Make a ridge plot for each sample for the log10(unique nuclear fragments)
# p3 <- plotGroups(
#     ArchRProj = projHeme1, 
#     groupBy = "Sample", 
#     colorBy = "cellColData", 
#     name = "log10(nFrags)",
#     plotAs = "ridges"
# )
# plotPDF(p3, name = file.path("0.QCmetrics","sample-Frags-ridgePlot.pdf"), ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)

#### Make a violin plot for each sample for the log10(unique nuclear fragments).
p4 <- plotGroups(
    ArchRProj = projHeme1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
)
plotPDF(p4, name = file.path("0.QCmetrics","sample-Frags-VlnPlot.pdf"), ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)
png(file.path(dip,"sample-Frags-VlnPlot.png"))
p4
dev.off()
#### Plotting Sample Fragment Size Distribution and TSS Enrichment Profiles.
p5 <- plotFragmentSizes(ArchRProj = projHeme1)
plotPDF(p5, name = file.path("0.QCmetrics","Sample-FragSizes.pdf"), ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)
png(file.path(dip,"Sample-FragSizes.png"))
p5
dev.off()
p6 <- plotTSSEnrichment(ArchRProj = projHeme1)
plotPDF(p6, name = file.path("0.QCmetrics","Sample-TSSProfile.pdf"), ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)
png(file.path(dip,"Sample-TSSProfile.png"))
p6
dev.off()
#### Saving and Loading an ArchRProject  
saveArchRProject(ArchRProj = projHeme1, load = FALSE)
file.copy (file.path(outdir,"Save-ArchR-Project.rds"),file.path(outrds,"projHeme1.rds"),recursive = TRUE)
# ### Filtering Doublets from an ArchRProject  (根据前面决定要不要filter Doublet)
# projHeme2 <- filterDoublets(projHeme1)
if (filterdoublet == "T"){
    #AMULET Doublet
    sample_Amulet_doublet<-c()
    for(i in 1:length(samplename)){
    file<-read.table(amulet_File[i]) 
    colnames(file)<-"barcode"   
    file$barcode_new<-paste0(samplename[i],"#",file$barcode)
    sample_Amulet_doublet<-c(sample_Amulet_doublet,file$barcode_new)   #AMULET doublet 
    }
    write.table(sample_Amulet_doublet,paste0(outdir,"/",prefix,".allSample_Amulet_doublet.txt"),row.names=F,col.names = FALSE,quote=F) 
    print("Amulet_doublet file done")
    #ArchR Doublet
    projHeme_0 <- projHeme1 
    projHeme_0 <- filterDoublets(projHeme_0)  #ArchR_Doublet barcode
    ArchR_doublet<-setdiff(rownames(projHeme1@cellColData),rownames(projHeme_0@cellColData)) #doublet
    write.table(ArchR_doublet,paste0(outdir,"/",prefix,".allSample_ArchR_doublet.txt"),row.names=F,col.names = FALSE,quote=F)
    print("ArchR_doublet file done")
    #doublet
    Co_doublet<-intersect(ArchR_doublet,sample_Amulet_doublet)
    write.table(Co_doublet,paste0(outdir,"/",prefix,".allSample_CoDoublet.txt"),row.names=F,col.names = FALSE,quote=F)
    print("Co_doublet file done")
    #projHeme_0 CoDoublet
    projHeme1@cellColData <- projHeme1@cellColData[!rownames(projHeme1@cellColData) %in% Co_doublet, , drop = FALSE]
    saveArchRProject(ArchRProj = projHeme1, load = FALSE)
    file.copy (file.path(outdir,"Save-ArchR-Project.rds"),file.path(outrds,"projHeme1.rds"),recursive = TRUE)
    
}else{
    saveArchRProject(ArchRProj = projHeme1, load = FALSE)
    file.copy (file.path(outdir,"Save-ArchR-Project.rds"),file.path(outrds,"projHeme1.rds"),recursive = TRUE)
}

# if (filterdoublet == "T"){
#     projHeme_0 <- projHeme1
#     projHeme1 <- filterDoublets(projHeme1)
#     saveArchRProject(ArchRProj = projHeme1, load = FALSE)
#     file.copy (file.path(outdir,"Save-ArchR-Project.rds"),file.path(outrds,"projHeme1.rds"))
#     #得到 doublet的barcode
#     ArchR_doublet<-setdiff(rownames(projHeme_0@cellColData),rownames(projHeme1@cellColData))
#     write.table(file.path(outdir,"Plots/0.QCmetrics"),paste0(prefix,".ArchR_doublet.txt"),row.names=F,col.names=F,quote=F)
#     print("ArchR_doublet filtered")
# }else{
#     saveArchRProject(ArchRProj = projHeme1, load = FALSE)
#     file.copy (file.path(outdir,"Save-ArchR-Project.rds"),file.path(outrds,"projHeme1.rds"))
# }

#####################
##### Dimensionality Reduction with ArchR
outdir_dim <- "1.DimReduction"
dir.create (file.path(outdir,"Plots",outdir_dim),showWarnings=T)
projHeme2 <- addIterativeLSI(
    ArchRProj = projHeme1,
    useMatrix = "TileMatrix",
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list(     #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:cnum
)

#### Cluster with ArchR 
if (method == "Seurat") {
    projHeme2 <- addClusters(
        input = projHeme2,
        reducedDims = "IterativeLSI",
        method = "Seurat",
        name = "Clusters",
        resolution = resol     ### 分群resolution 可改
    )

    ## UMAP plot
    projHeme2 <- addUMAP(
        ArchRProj = projHeme2, 
        reducedDims = "IterativeLSI", 
        name = "UMAP", 
        nNeighbors = 30, 
        minDist = 0.5, 
        metric = "cosine"
    )
    p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")   ## colorBy 指定使用的matrix的名字（cellColData 类型与Seurat中的metadata），name 指定具体的列
    plotPDF(p1, name = file.path(outdir_dim,paste0(prefix,".UMAP_IterativeLSI_samples.pdf")), ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
    png (file.path(outdir,"Plots",outdir_dim,paste0(prefix,".UMAP_IterativeLSI_samples.png")))
    print (p1)
    dev.off()
    p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
    plotPDF(p2, name =file.path(outdir_dim,paste0(prefix,".UMAP_IterativeLSI_Cluster.pdf")), ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
    png (file.path(outdir,"Plots",outdir_dim,paste0(prefix,".UMAP_IterativeLSI_Cluster.png")))
    print (p2)
    dev.off()
    ## tSNE plot
    projHeme2 <- addTSNE(
        ArchRProj = projHeme2, 
        reducedDims = "IterativeLSI", 
        name = "TSNE", 
        perplexity = 30
    )
    p3 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
    #plotPDF(p3, name = paste0(prefix,".tSNE_IterativeLSI_samples.pdf"), ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
    plotPDF(p3, name = file.path(outdir_dim,paste0(prefix,".tSNE_IterativeLSI_samples.pdf")), ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
    png (file.path(outdir,"Plots",outdir_dim,paste0(prefix,".tSNE_IterativeLSI_samples.png")))
    print (p3)
    dev.off()

    p4 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "TSNE")
    #plotPDF(p4, name = paste0(prefix,".tSNE_IterativeLSI_Cluster.pdf"), ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
    plotPDF(p4, name = file.path(outdir_dim,paste0(prefix,".tSNE_IterativeLSI_Cluster.pdf")), ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
    png (file.path(outdir,"Plots",outdir_dim,paste0(prefix,".tSNE_IterativeLSI_Cluster.png")))
    print (p4)
    dev.off()

    ## cluster prop plot
    pal <- paletteDiscrete(values = projHeme2$Clusters)
    cluster <- table (projHeme2$Clusters,projHeme2$Sample)
    write.table (cluster,file.path(outdir,"Plots",outdir_dim,paste0(prefix,".CellsPerCluster.xls")),col.names = NA,sep='\t',quote = F)
    freq_table <- prop.table (x=table (projHeme2$Clusters,projHeme2$Sample),margin = 2)
    write.table (freq_table,file.path(outdir,"Plots",outdir_dim,paste0(prefix,".CellsPropPerCluster.xls")),col.names = NA,sep='\t',quote =F)
    pdf (file.path(outdir,"Plots",outdir_dim,paste0(prefix,".CellsPropPer.pdf")))
    par(mar=c(9,4,4,2)+0.1)
    barplot(height=freq_table,width = 5,xlim=c(1,75),col =pal,legend = rownames(freq_table),args.legend = list(x = "right"),las=2,xlab="")
    dev.off()
    png (file.path(outdir,"Plots",outdir_dim,paste0(prefix,".CellsPropPer.png")))
    par(mar=c(9,4,4,2)+0.1)
    barplot(height=freq_table,width = 5,xlim=c(1,75),col =pal,legend = rownames(freq_table),args.legend = list(x = "right"),las=2,xlab="")
    dev.off()       

}else {
    projHeme2 <- addClusters(
        input = projHeme2,
        reducedDims = "IterativeLSI",
        method = "scran",
        name = "ScranClusters",
        k = 15       ### 可改
    )
    ## umap plot
    projHeme2 <- addUMAP(
        ArchRProj = projHeme2, 
        reducedDims = "IterativeLSI", 
        name = "UMAP", 
        nNeighbors = 30, 
        minDist = 0.5, 
        metric = "cosine"
    )
    p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
    plotPDF(p1, name = file.path(outdir_dim,paste0(prefix,".UMAP_IterativeLSI_samples.pdf")), ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
    png (file.path(outdir,"Plots",outdir_dim,paste0(prefix,".UMAP_IterativeLSI_samples.png")))
    print (p1)
    dev.off()
    p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
    plotPDF(p2, name = file.path(outdir_dim,paste0(prefix,".UMAP_IterativeLSI_Cluster.pdf")), ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
    png (file.path(outdir,"Plots",outdir_dim,paste0(prefix,".UMAP_IterativeLSI_Cluster.png")))
    print (p2)
    dev.off()
    ## tSNE plot
    projHeme2 <- addTSNE(
        ArchRProj = projHeme2, 
        reducedDims = "IterativeLSI", 
        name = "TSNE", 
        perplexity = 30
    )
    p3 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
    plotPDF(p3, name = file.path(outdir_dim,paste0(prefix,".tSNE_IterativeLSI_samples.pdf")), ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
    png (file.path(outdir,"Plots",outdir_dim,paste0(prefix,".tSNE_IterativeLSI_samples.png")))
    print (p3)
    dev.off()
    p4 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "TSNE")
    plotPDF(p4, name = file.path(outdir_dim,paste0(prefix,".tSNE_IterativeLSI_Cluster.pdf")), ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
    png (file.path(outdir,"Plots",outdir_dim,paste0(prefix,".tSNE_IterativeLSI_Cluster.png")))
    print (p4)
    dev.off()    
    # cM <- confusionMatrix(paste0(projHeme2$Clusters), paste0(projHeme2$Sample))
    # cluster_table <- as.data.frame (cM) 
    # write.table (cluster_table,file.path(outdir_dim,paste0 (prefix,".CellsPerCluster.xls")),col.names=NA,sep="\t",quote=F) 
    ## cluster prop plot
    pal <- paletteDiscrete(values = projHeme2$Clusters)
    cluster <- table (projHeme2$Clusters,projHeme2$Sample)
    write.table (cluster,file.path(outdir,"Plots",outdir_dim,paste0(prefix,".CellsPerCluster.xls")),col.names = NA,sep='\t',quote = F)
    freq_table <- prop.table (x=table (projHeme2$Clusters,projHeme2$Sample),margin = 2)
    write.table (freq_table,file.path(outdir,"Plots",outdir_dim,paste0(prefix,".CellsPropPerCluster.xls")),col.names = NA,sep='\t',quote =F)
    pdf (file.path(outdir,"Plots",outdir_dim,paste0(prefix,".CellsPropPer.pdf")))
    par(mar=c(9,4,4,2)+0.1)
    barplot(height=freq_table,width = 5,xlim=c(1,75),col =pal,legend = rownames(freq_table),args.legend = list(x = "right"),las=2,xlab="")
    dev.off()
    png (file.path(outdir,"Plots",outdir_dim,paste0(prefix,".CellsPropPer.png")))
    par(mar=c(9,4,4,2)+0.1)
    barplot(height=freq_table,width = 5,xlim=c(1,75),col =pal,legend = rownames(freq_table),args.legend = list(x = "right"),las=2,xlab="")
    dev.off()       
}
saveArchRProject(ArchRProj = projHeme2, load = FALSE)
file.copy (file.path(outdir,"Save-ArchR-Project.rds"),file.path(outrds,"projHeme2-reDim.rds"),recursive = TRUE)
###############去批次的这部分可以单独拆出来写一个脚本
# #### Batch Effect Correction wtih Harmony (需要确认要不要去批次) (输出去批次的结果，根据结果判断需不需要去批次)
# projHeme2 <- addHarmony(
#     ArchRProj = projHeme2,
#     reducedDims = "IterativeLSI",
#     name = "Harmony",
#     groupBy = "Sample"
# )
## Dimensionality Reduction After Harmony
# projHeme2 <- addUMAP(
#     ArchRProj = projHeme2, 
#     reducedDims = "Harmony", 
#     name = "UMAPHarmony", 
#     nNeighbors = 30, 
#     minDist = 0.5, 
#     metric = "cosine"
# )
# p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
# plotPDF(p1, name = file.path(outdir_dim,paste0(prefix,".UMPA_Harmony_samples.pdf")), ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
# p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
# plotPDF(p2, name = file.path(outdir_dim,paste0(prefix,".UMAP_Harmony_Cluster.pdf")), ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
# projHeme2 <- addTSNE(
#     ArchRProj = projHeme2, 
#     reducedDims = "Harmony", 
#     name = "TSNEHarmony", 
#     perplexity = 30
# )
# p3 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "TSNEHarmony")
# plotPDF(p3, name = file.path(outdir_dim,paste0(prefix,".tSNE_Harmony_samples.pdf")), ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
# p4 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "TSNEHarmony")
# plotPDF(p4, name = file.path(outdir_dim,paste0(prefix,".tSNE_Harmony_Cluster.pdf")), ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
# saveArchRProject(ArchRProj = projHeme2, outputDirectory = "Save-ProjHeme2", load = FALSE)
############################



### Defining Cluster Identity with scRNA-seq (Unconstrained Integration)
outdir_cellanno <- "3.Cellannotation"
dir.create (file.path(outdir,"Plots",outdir_cellanno),showWarnings=T)
dianno <- file.path(outdir,"Plots",outdir_cellanno)
seRNA <- readRDS (scRNA)
seRNA$barcode <- rownames (seRNA@meta.data)
seRNA$BioClassification <- as.character (seRNA@active.ident)
projHeme3 <- addGeneIntegrationMatrix(
    ArchRProj = projHeme2, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",     ## 如果需要去批次的话，这里要改成Harmony
    seRNA = seRNA,
    addToArrow = TRUE,
    groupRNA = "BioClassification",
    nameCell = "predictedCell_Un", 
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un",
    force = TRUE
)
## pal <- paletteDiscrete(values = seRNA$BioClassification) ## 这里可以修改注释颜色
color.use<-c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","DarkGreen","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080","#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49","#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")
source('/Public/Script/shouhou/SCRIPT/Seurat_Monocle_modify/color_protocol.R')
color.use <- c(color_protocol, color.use)
pal <- color.use[1:length(levels(seRNA))]
names (pal) <- levels(seRNA)


p1 <- plotEmbedding(
    projHeme3, 
    colorBy = "cellColData", 
    name = "predictedGroup_Un", 
    pal = pal
)
#### 得到降维每个cluster对应的细胞类型 Labeling scATAC-seq clusters with scRNA-seq information
#### cluster 方式进行注释
# cM <- confusionMatrix(projHeme3$Clusters, projHeme3$predictedGroup_Un)
# labelOld <- rownames(cM)
# labelOld
# labelNew <- colnames(cM)[apply(cM, 1, which.max)]
# labelNew
# remapClust <- unique (seRNA$BioClassification)
# names (remapClust) <- unique (seRNA$BioClassification)
# remapClust <- remapClust[names(remapClust) %in% labelNew]
# labelNew2 <- mapLabels(labelNew, oldLabels = names(remapClust), newLabels = remapClust)
# labelNew2
# projHeme3$Clusters2 <- mapLabels(projHeme3$Clusters, newLabels = labelNew2, oldLabels = labelOld)
### cellbase 方式进行注释

## 筛选细胞
print (table(projHeme3$predictedGroup_Un))
print ("### filter celltype numer <50###")
celltype <- names(table(projHeme3$predictedGroup_Un))[table(projHeme3$predictedGroup_Un)>50]
print (celltype)
idxTypes <- BiocGenerics::which(projHeme3$predictedGroup_Un %in% celltype)
cellsTypes <- projHeme3$cellNames[idxTypes]
projHeme3 <- projHeme3[cellsTypes, ]

projHeme3$Clusters2 <- projHeme3$predictedGroup_Un
#### saveRDS 
saveArchRProject(ArchRProj = projHeme3, load = FALSE)
file.copy (file.path(outdir,"Save-ArchR-Project.rds"),file.path(outrds,"projHeme3.anno.rds"),recursive = TRUE)
## get bigwig file
bg <-  getGroupBW(
    ArchRProj = projHeme3,
    groupBy = "Clusters2",
    normMethod = "ReadsInTSS",
    tileSize = 100,
    maxCells = 1000,
    ceiling = 4,
    verbose = TRUE,
    threads = 1,
    logFile = createLogFile("getGroupBW")
)

## celltype color
clustcol <- pal[which (names(pal) %in% unique(projHeme3$Clusters2))]
clustcol.barplot <- rev (clustcol)
anno.meta <- projHeme3@cellColData
anno.meta$Clusters2 <- factor (anno.meta$Clusters2,levels = names (clustcol.barplot))
anno.meta$Sample <- factor (anno.meta$Sample,levels = mixedsort (unique (anno.meta$Sample)))
## sample color
color.sample <- color.use[1:length(unique (projHeme3$Sample))]
names (color.sample) <-levels (anno.meta$Sample)

#### UMAP plot 
p1 <- plotEmbedding(projHeme3, colorBy = "cellColData", name = "Clusters2",baseSize=8,size = 0.2,quantileCut = c(0.025,0.975),embedding = 'UMAP')
plotPDF(p1, name = file.path(outdir_cellanno,paste0(prefix,"-UMAP-RNA-Integration.pdf")), ArchRProj = projHeme3, addDOC = FALSE)
p2 <- plotEmbedding(projHeme3, colorBy = "cellColData", name = "Sample",baseSize=8,size = 0.2,quantileCut = c(0.025,0.975),embedding = 'UMAP')
plotPDF(p2, name = file.path(outdir_cellanno,paste0(prefix,"-sample-UMAP-RNA-Integration.pdf")), ArchRProj = projHeme3, addDOC = FALSE)
## uMAP replot
umap.file <- projHeme3@embeddings$UMAP$df
colnames (umap.file) <- c("UMAP_1","UMAP_2")
write.table (umap.file,file.path(dianno,paste0(prefix,".UMAP_Dimension.file.xls")),col.names = NA,sep='\t',quote = F)
data.umap <- data.frame(barcode = rownames (projHeme3),celltype = anno.meta$Clusters2,umap.file)


p3 <- ggplot (data.umap,mapping = aes (x=UMAP_1,y=UMAP_2,color=celltype)) + geom_point(size = 0.2)  + scale_color_manual(values = clustcol) + 
        theme_classic() + theme(legend.title=element_blank(),legend.text = element_text (size = 10)) +
        guides(color = guide_legend(override.aes = list(size = 2)))    
pdf (file.path(dianno,paste0(prefix,".labumap.integration.pdf")))
p3
dev.off()
png (file.path(dianno,paste0(prefix,".labumap.integration.png")))
p3
dev.off()
## umap sample replot

sample.umap <- data.frame (barcode = rownames (projHeme3),sample = anno.meta$Sample,umap.file)
p4 <- ggplot (sample.umap,mapping = aes (x=UMAP_1,y=UMAP_2,color=sample)) + geom_point(size = 0.2)  + scale_color_manual(values = color.sample) + 
        theme_classic() + theme(legend.title=element_blank(),legend.text = element_text (size = 10)) +
        guides(color = guide_legend(override.aes = list(size = 2)))

pdf (file.path(dianno,paste0(prefix,".labumap.sample.pdf")))
p4
dev.off()
png (file.path(dianno,paste0(prefix,".labumap.sample.png")))
p4
dev.off()

## TSNE plot 
tsne.file <- projHeme3@embeddings$TSNE$df
colnames (tsne.file) <- c("TSNE_1","TSNE_2")
write.table (tsne.file,file.path(dianno,paste0(prefix,".TSNE_Dimension.file.xls")),col.names = NA,sep='\t',quote = F)
# data <- data.frame (barcode = rownames (projHeme3),celltype = projHeme3@cellColData$Clusters2)
# data.tsne <- bind_cols (data,tsne.file)
# clustcol <- pal[unique (data.tsne$celltype)]
data.tsne <- data.frame(barcode = rownames (projHeme3),celltype = anno.meta$Clusters2,tsne.file)
p3 <- ggplot (data.tsne,mapping = aes (x=TSNE_1,y=TSNE_2,color=celltype)) + geom_point(size = 0.2)  + scale_color_manual(values = clustcol) + 
        theme_classic() + theme(legend.title=element_blank(),legend.text = element_text (size = 10)) +
        guides(color = guide_legend(override.aes = list(size = 2)))
pdf (file.path(dianno,paste0(prefix,".labtsne.integration.pdf")))
p3
dev.off()
png (file.path(dianno,paste0(prefix,".labtsne.integration.png")))
p3
dev.off()
## tsne sample replot
# data.sample <- data.frame (barcode = rownames (projHeme3),sample = projHeme3@cellColData$Sample)
# sample.tsne <- bind_cols (data.sample,tsne.file)
sample.tsne <- data.frame (barcode = rownames (projHeme3),sample = anno.meta$Sample,tsne.file)
p4 <- ggplot (sample.tsne,mapping = aes (x=TSNE_1,y=TSNE_2,color=sample)) + geom_point(size = 0.2)  + scale_color_manual(values = color.sample) + 
        theme_classic() + theme(legend.title=element_blank(),legend.text = element_text (size = 10)) +
        guides(color = guide_legend(override.aes = list(size = 2)))

pdf (file.path(dianno,paste0(prefix,".labtsne.sample.pdf")))
p4
dev.off()
png (file.path(dianno,paste0(prefix,".labtsne.sample.png")))
p4
dev.off()
## cluster prop plot
# data <- data.frame (celltype = names (pal),pal) 
# data <- data[unique (projHeme3$Clusters2),]
# data <- data[order(data$celltype),]

cluster <- table (anno.meta$Clusters2,anno.meta$Sample)
write.table (cluster,file.path(dianno,paste0(prefix,".CellsPerCellType.xls")),col.names = NA,sep='\t',quote = F)
# freq_table <- prop.table (x=table (projHeme3$Clusters2,projHeme3$Sample),margin = 2)
freq_table <- prop.table (x=table (anno.meta$Clusters2,anno.meta$Sample),margin = 2)
write.table (freq_table,file.path(dianno,paste0(prefix,".CellsPropPerCellTypes.xls")),col.names = NA,sep='\t',quote =F)
table <- melt (freq_table)
colnames (table) <- c("celltype","sample","value")
p2 <- ggplot (table,aes (x=sample,y=value,fill = celltype)) + geom_bar (stat = 'identity',size = 0.3,color = "black",width = 0.7) +
theme_classic() +
scale_fill_manual (values = clustcol.barplot) +
xlab ("") + ylab("") + theme (legend.title = element_blank(),legend.text=element_text(size=12)) +
theme (axis.text= element_text(color = "black",size = 12),
axis.text.x = element_text(angle = 90,hjust = 1),
axis.line.x = element_line (color ="transparent"),
axis.ticks.x = element_blank(),aspect.ratio = 1) +
scale_y_continuous (expand = c(0,0))
png (file.path(dianno,paste0(prefix,".CellsPropPerCellTypes.png")))
p2
dev.off()
pdf (file.path(dianno,paste0(prefix,".CellsPropPerCellTypes.pdf")))
p2
dev.off()

# pdf (file.path(dianno,paste0(prefix,".CellsPropPerCellTypes.pdf")))
# par(mar=c(9,4,4,2)+0.1)
# barplot(height=freq_table,width = 3,xlim=c(1,75),col =data$pal,legend = rownames(freq_table),args.legend = list(x = "right"),las=2,xlab="")
# dev.off()
# png (file.path(dianno,paste0(prefix,".CellsPropPerCellTypes.png")))
# par(mar=c(9,4,4,2)+0.1)
# barplot(height=freq_table,width = 3,xlim=c(1,75),col =data$pal,legend = rownames(freq_table),args.legend = list(x = "right"),las=2,xlab="")
# dev.off()    


# cM <- confusionMatrix(paste0(projHeme3$Clusters2), paste0(projHeme3$Sample))
# freq_table <- as.data.frame (cM) 
# write.table (freq_table,file.path(outdir,"Plots",outdir_cellanno,paste0 (prefix,".CellsPerCellType.xls")),col.names=NA,sep="\t",quote=F)

## annotation result cp 
annors.out <- file.path(dianno,prefix)
dir.create (annors.out)
file.copy (file.path(dianno,paste0(prefix,".labumap.integration.pdf")),annors.out,recursive = TRUE)
file.copy (file.path(dianno,paste0(prefix,".labumap.integration.png")),annors.out,recursive = TRUE)
file.copy (file.path(dianno,paste0(prefix,".CellsPropPerCellTypes.pdf")),annors.out,recursive = TRUE)
file.copy (file.path(dianno,paste0(prefix,".CellsPropPerCellTypes.png")),annors.out,recursive = TRUE)
file.copy (file.path(dianno,paste0(prefix,".CellsPerCellType.xls")),annors.out,recursive = TRUE)
file.copy (file.path(dianno,paste0(prefix,".labumap.sample.pdf")),annors.out,recursive = TRUE)
file.copy (file.path(dianno,paste0(prefix,".labumap.sample.png")),annors.out,recursive = TRUE)
file.copy ("/SGRNJ03/PiplineTest01/Software_test/zhengzhanye/ArchR_V1.0/readme/annotation.pdf",annors.out,recursive = TRUE)
## 得到降维坐标
## projHeme3@embeddings$UMAP@listData$df

### Identify peaks
projHeme4 <- addGroupCoverages(ArchRProj = projHeme3, groupBy = "Clusters2")
### Calling Peaks with ArchR
### 细胞数目过少时，不能找marker peak,默认最少50个细胞
print ("filter cell number > 50,can not find marker peak under this threshold")

pathToMacs2 <- findMacs2()
projHeme4 <- addReproduciblePeakSet(
    ArchRProj = projHeme4, 
    groupBy = "Clusters2", 
    pathToMacs2 = pathToMacs2,
    minCells = 50
)
getPeakSet(projHeme4)
## 得到各个细胞类型鉴定的peak文件 (确认一下这里得到的是不是每一种细胞类型的peak)
df1 <- DataFrame (iranges = getPeakSet(projHeme4)) 
peak <- data.frame (celltype = rownames (df1),data.frame (getPeakSet(projHeme4)))
dio_macscallpeak <- file.path(outdir,"Plots","4.1.MacsCallPeak")
dir.create (dio_macscallpeak)
for (ct in unique (peak$celltype)) {
    rs <- filter (peak,celltype == ct)
    rs <- select (rs,-N)
    rs <- select (rs,-c(score,replicateScoreQuantile,groupScoreQuantile,Reproducibility,GroupReplicate,distToTSS,nearestTSS,GC,idx))
    write.table (rs,file.path (dio_macscallpeak,paste0(ct,".macsCallpeak.xls")),sep='\t',quote =F,row.names = F,col.names = T)
}

### Add Peak Matrix
projHeme4 <- addPeakMatrix(projHeme4) 
saveArchRProject(ArchRProj = projHeme4, load = FALSE)
file.copy (file.path(outdir,"Save-ArchR-Project.rds"),file.path(outrds,"projHeme4.peak.rds"),recursive = TRUE)

#### peak type pieplot
meta <- data.frame (projHeme4@peakSet)
data <- meta %>% group_by (peakType) %>% summarise (n=n()) %>% as.data.frame()
colnames (data) <- c("peaktype","num")
data$prop <- data$num/sum(data$num)
data <- data[order(data$prop),]
data$peaktype <- factor (data$peaktype,levels = data$peaktype)
write.table (data,file.path(dio_macscallpeak,"scATAC.peaktype.xls"),sep='\t',quote = F,row.names = F,col.names = T)
p <- ggplot (data,mapping = aes(x = 'Content',y = prop, fill = peaktype))+
      geom_bar(stat = 'identity', position = 'stack',color = "white")+
      coord_polar(theta = 'y')  +  theme_minimal() +
      theme (axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())+
      theme (panel.grid = element_blank())+
      theme (legend.title = element_blank()) +
      scale_fill_brewer(palette="PuBu")
pdf (file.path(dio_macscallpeak,"scATAC.peakType.pieplot.pdf"),height = 4,width = 4)
print (p)
dev.off()
png (file.path(dio_macscallpeak,"scATAC.peakType.pieplot.png"),units = 'cm',res = 300,height = 9,width = 9)
print (p)
dev.off()
## add readme 
file.copy ("/SGRNJ03/PiplineTest01/Software_test/zhengzhanye/ArchR_V1.0/readme/MacsCallPeak.readme.pdf",dio_macscallpeak,recursive = TRUE)
#### Identifying Marker Peaks with ArchR
#Our scRNA labels 
## table(projHeme4$Clusters2,projHeme4$Sample)每种细胞类型在每个样本中的细胞数
table(projHeme4$Clusters2)
# ### 细胞数目过少时，不能找marker peak,默认最少50个细胞
# print ("filter cell number > 50,can not find marker peak under this threshold")
# ### subset the project to keep all cells corresponding to a specific celltype: 
# celltype <- names(table(projHeme4$Clusters2))[table(projHeme4$Clusters2)>peakfilter]
# idxTypes <- BiocGenerics::which(projHeme4$Clusters2 %in% celltype)
# cellsTypes <- projHeme4$cellNames[idxTypes]
# projHeme4 <- projHeme4[cellsTypes, ]


# markersPeaks <- getMarkerFeatures(
#     ArchRProj = projHeme4, 
#     useMatrix = "PeakMatrix", 
#     groupBy = "Clusters2",
#     bias = c("TSSEnrichment", "log10(nFrags)"),
#     testMethod = "wilcoxon",
#     binarize = TRUE
# )
#testMethod = "wilcoxon" ## "wilcoxon" "ttest", and "binomial".
if (testmethod != "binomial") {
    markersPeaks <- getMarkerFeatures(
        ArchRProj = projHeme4, 
        useMatrix = "PeakMatrix", 
        groupBy = "Clusters2",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        testMethod = testmethod,
        binarize = TRUE
    )   
}else{
    markersPeaks <- getMarkerFeatures(
        ArchRProj = projHeme4, 
        useMatrix = "PeakMatrix", 
        groupBy = "Clusters2",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        testMethod = "binomial",
        binarize = TRUE
    )
}

## markergene
markersGS <- getMarkerFeatures(
    ArchRProj = projHeme4, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters2",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
mm <- as.data.frame (markerList) 
markerGenes <- mm$name

### markerpeak 
dio_markerpeak <- file.path(outdir,"Plots","4.2.Markerpeak")
dir.create (dio_markerpeak)
markerPeaksList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1") ###FDR <= 0.1 & Log2FC >= 1
for (ctp in names (markerPeaksList)) {
    print (ctp)
    markerpeak <- data.frame (markerPeaksList[ctp])
    markerpeak <- select (markerpeak,-c("group","group_name"))
    #markerpeak <- select (markerpeak,c("idx","seqnames","start","end","Log2FC","FDR","MeanDiff"))
    markerpeak <- select (markerpeak,c("seqnames","start","end","Log2FC","FDR","MeanDiff"))
    if (nrow (markerpeak) > 1) {
        write.table (markerpeak,file.path(dio_markerpeak,paste0(ctp,".markerpeak.xls")),sep='\t',quote=F,row.names = F,col.names = T)
        ### ### Marker Peak MA and Volcano Plots
        pma <- markerPlot(seMarker = markersPeaks, name = ctp, cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "MA")
        pdf (file.path(dio_markerpeak,paste0(ctp,".MarkerPeak.MA.pdf")))
        print (pma)
        dev.off()
        png (file.path(dio_markerpeak,paste0(ctp,".MarkerPeak.MA.png")))
        print (pma)
        dev.off()

        pv <- markerPlot(seMarker = markersPeaks, name = ctp, cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
        pdf (file.path(dio_markerpeak,paste0(ctp,".MarkerPeak.Volcano.pdf")))
        print (pv)
        dev.off()
        png (file.path(dio_markerpeak,paste0(ctp,".MarkerPeak.Volcano.png")))
        print (pv)
        dev.off()
        ### Marker Peak  in Browser Tracks
       
        p <- plotBrowserTrack(
            ArchRProj = projHeme4, 
            groupBy = "Clusters2", 
            geneSymbol = markerGenes[4],
            features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)[ctp], 
            upstream = 50000,
            downstream = 50000
        )
        pdf (file.path(dio_markerpeak,paste0(ctp,".MarkerPeak.Tracks.pdf")))
        print (grid::grid.draw(p[[markerGenes[4]]]))
        dev.off()
        png (file.path(dio_markerpeak,paste0(ctp,".MarkerPeak.Tracks.png")))
        print (grid::grid.draw(p[[markerGenes[4]]]))
        dev.off()

    }  
}
# markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE) ##return a GRangesList object


### Marker Peak Heatmaps
if (ncol(markersPeaks) > 2) {
    heatmapPeaks <- markerHeatmap(
    seMarker = markersPeaks, 
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
    transpose = FALSE
    )
    pdf (file.path(dio_markerpeak,"MarkerPeak.Heatmap.pdf"),width = 7, height = 20)
    draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "right")
    dev.off()
    png (file.path(dio_markerpeak,"MarkerPeak.Heatmap.png"),height = 950)
    heatmapPeaks
    dev.off()
}

## add readme
file.copy ("/SGRNJ03/PiplineTest01/Software_test/zhengzhanye/ArchR_V1.0/readme/Markerpeak.readme.pdf",dio_markerpeak,recursive = TRUE)

############Motif Enrichment in Marker Peaks
dir.create (file.path(outdir,"Plots","5.MarkerpeakMotifEnrichment"),showWarnings=T)
outdir_motif <- file.path(outdir,"Plots","5.MarkerpeakMotifEnrichment")

projHeme5 <- addMotifAnnotations(ArchRProj = projHeme4, motifSet = "JASPAR2020", name = "Motif")
projHeme5 <- addImputeWeights(projHeme5)

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = projHeme5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
heatmapEM.table <- plotEnrichHeatmap(enrichMotifs, n = 1000, transpose = FALSE,returnMatrix = TRUE)
heatmapEM.table <- data.frame (motiftpm = rownames (heatmapEM.table),heatmapEM.table)
#ht <- separate (heatmapEM.table,motiftpm,into=c("motif"),sep= " ")
ht <- separate (heatmapEM.table,motiftpm,into=c("motif"),sep= "_")
write.table (ht,file.path(outdir_motif,"MarkerPeakEnrichMotif.xls"),sep='\t',quote=F,row.names =F,col.names = T)
## 重新画图
heatmapEM.plot <- plotEnrichHeatmap(enrichMotifs, n = 10, transpose = FALSE,returnMatrix = TRUE)
heatmapEM.plot <- data.frame (motiftpm = rownames (heatmapEM.plot),heatmapEM.plot)
# ht.plot <- separate (heatmapEM.plot,motiftpm,into=c("motif"),sep= " ")
ht.plot <- separate (heatmapEM.plot,motiftpm,into=c("motif"),sep= "_")
rownames (ht.plot) <- ht.plot$motif
write.table (ht.plot,file.path(outdir_motif,"MarkerPeakEnrichMotif.Plot.xls"),sep='\t',quote=F,row.names =F,col.names = T)
ht.plot <- ht.plot[,-1]
ht.plot <- as.matrix (ht.plot)
pdf (file.path(outdir_motif,"MarkerPeakEnrichMotif.pdf"),height = 13)
print (pheatmap (ht.plot,color =  paletteContinuous(set = "comet", n = 100),cluster_row = F))
dev.off()
png (file.path(outdir_motif,"MarkerPeakEnrichMotif.png"),height = 800)
print (pheatmap (ht.plot,color =  paletteContinuous(set = "comet", n = 100),cluster_row = F))
dev.off()
# ## 自带画heatmap 方法
# heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = FALSE)
# pdf (file.path(outdir_motif,"MarkerPeakEnrichMotif.pdf"))
# ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
# dev.off()
# png (file.path(outdir_motif,"MarkerPeakEnrichMotif.png"))
# ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
# dev.off()
# ##################################
# ###  Pairwise Testing Between Groups 比较不同细胞类型下的差异peak （这部分可以拿出来单独做，老师指定相应细胞类型在做这个比较展示）
# #outdir_diffpeak <- "5.ComparePeakPlot"
# dir.create (file.path(outdir,"Plots","5.ComparePeakPlot"))
# outdir_diffpeak <- file.path(outdir,"Plots","5.ComparePeakPlot")
# ids <- combn (cellnames,2) #random  select two celltype
# cell1 <- ids[1,1]
# cell2 <- ids[2,1]
# print (cell1)
# print (cell2)
# pdf (file.path(outdir_diffpeak,paste0(cell1,"-vs-",cell2,"")))
# markerTest <- getMarkerFeatures(
#     ArchRProj = projHeme4, 
#     useMatrix = "PeakMatrix",
#     groupBy = "Clusters2",
#     testMethod = "wilcoxon",
#     #testMethod = "binomial",
#     binarize = FALSE, 
#     bias = c("TSSEnrichment", "log10(nFrags)"),
#     useGroups = cell1,
#     bgdGroups = cell2
# )
# ### MA和小提琴图展示
# # pma <- markerPlot(seMarker = markerTest, name = cell1, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
# # plotPDF(pma, name = file.path(outdir_diffpeak,paste0(cell1,"-vs-",cell2,"-Markers-MA.pdf")), width = 5, height = 5, ArchRProj = projHeme4, addDOC = FALSE)
# # pv <- markerPlot(seMarker = markerTest, name = cell2, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
# # plotPDF(pv, name = file.path(outdir_diffpeak,paste0(cell1,"-vs-",cell2,"-Volcano-MA.pdf")), width = 5, height = 5, ArchRProj = projHeme4, addDOC = FALSE)
##########################################################

# ###  Motif Enrichment in Differential Peaks (基于以上两种细胞类型差异peak 做motif比较)
# outdir_diffmotif <- "5.DiffpeakMotifEnrichment"
# dir.create (file.path(outdir,"Plots",outdir_diffmotif),showWarnings=T)
# ### compare up enrichment
# motifsUp <- peakAnnoEnrichment(
#     seMarker = markerTest,
#     ArchRProj = projHeme4,
#     peakAnnotation = "Motif",
#     cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
# )
# df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
# df <- df[order(df$mlog10Padj, decreasing = TRUE),]
# df$rank <- seq_len(nrow(df))
# write.table (df,file.path (outdir,"Plots",outdir_diffmotif,paste0(cell1,"-vs-",cell2,"-Markers-Motifs-Up-Enriched.xls")),row.names=F,col.names=T,sep='\t',quote=F)
# ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
# geom_point(size = 1) +
# ggrepel::geom_label_repel(
#         data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
#         size = 1.5,
#         nudge_x = 2,
#         color = "black"
# ) + theme_ArchR() + 
# ylab("-log10(P-adj) Motif Enrichment") + 
# xlab("Rank Sorted TFs Enriched") +
# scale_color_gradientn(colors = paletteContinuous(set = "comet")
# )
# plotPDF(ggUp, name = file.path(outdir_diffmotif,paste0(cell1,"-vs-",cell2,"-Markers-Motifs-UP-Enriched.pdf")), width = 5, height = 5, ArchRProj = projHeme4, addDOC = FALSE)
# ## compare down enrichment
# motifsDo <- peakAnnoEnrichment(
#     seMarker = markerTest,
#     ArchRProj = projHeme4,
#     peakAnnotation = "Motif",
#     cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
# )
# df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
# df <- df[order(df$mlog10Padj, decreasing = TRUE),]
# df$rank <- seq_len(nrow(df))
# write.table (df,file.path (outdir,"Plots",outdir_diffmotif,paste0(cell1,"-vs-",cell2,"-Markers-Motifs-Down-Enriched.xls")),row.names=F,col.names=T,sep='\t',quote=F)
# ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
# geom_point(size = 1) +
# ggrepel::geom_label_repel(
#         data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
#         size = 1.5,
#         nudge_x = 2,
#         color = "black"
# ) + theme_ArchR() + 
# ylab("-log10(FDR) Motif Enrichment") +
# xlab("Rank Sorted TFs Enriched") +
# scale_color_gradientn(colors = paletteContinuous(set = "comet")
# )
# plotPDF(ggDo, name = file.path(outdir_diffmotif,paste0(cell1,"-vs-",cell2,"-Markers-Motifs-Down-Enriched.pdf")), width = 5, height = 5, ArchRProj = projHeme4, addDOC = FALSE)
# saveArchRProject(ArchRProj = projHeme4, outputDirectory = "Save-ProjHeme4", load = FALSE)
##################################################

## add readme
file.copy ("/SGRNJ03/PiplineTest01/Software_test/zhengzhanye/ArchR_V1.0/readme/MarkerpeakMotifEnrichment.readme.pdf",outdir_motif,recursive = TRUE)

### Chapter 13 ChromVAR Deviatons Enrichment with ArchR
dir.create (file.path(outdir,"Plots","6.ChromVAR.Variable.motif"),showWarnings=T)
outdir_chromvar <- file.path(outdir,"Plots","6.ChromVAR.Variable.motif")

if("Motif" %ni% names(projHeme5@peakAnnotation)){
    projHeme5 <- addMotifAnnotations(ArchRProj = projHeme5, motifSet = "JASPAR2020", name = "Motif")
}

projHeme5 <- addBgdPeaks(projHeme5)
projHeme5 <- addDeviationsMatrix(
    ArchRProj = projHeme5, 
    peakAnnotation = "Motif",
    force = TRUE,
    matrixName = "MotifMatrix"
)
plotVarDev <- getVarDeviations(projHeme5, name = "MotifMatrix", plot = TRUE)
pdf (file.path(outdir_chromvar,"Variable-Motif-Deviation-Scores.pdf"))
plotVarDev
dev.off()
png (file.path(outdir_chromvar,"Variable-Motif-Deviation-Scores.png"))
plotVarDev
dev.off()
## extract a subset of motifs for downstream analysis #(select top5 motif)
varmotif <- getVarDeviations(projHeme5, name = "MotifMatrix", plot = FALSE) %>% as.data.frame ()
varmotif <- select (varmotif,-c("seqnames","idx"))
write.table (varmotif,file.path(outdir_chromvar,"ALL.Var-motif-Deviations.xls"),row.names=F,col.names=T,sep='\t',quote=F)
varmotif <- separate(varmotif,name,into=c("motif_name"),sep="_",remove=F)
motifs <- head (varmotif$motif_name,2)  
#motifs <- c("GATA1", "CEBPA", "EBF1", "IRF4", "TBX21", "PAX5")
markerMotifs <- getFeatures(projHeme5, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

### MotifMatrix contains seqnames for both z-scores and deviations, shown above by “z:” and “deviations:”.
###  remove z:SREBF1_22
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
# markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
# markerMotifs
## plot the distribution of chromVAR deviation scores for each cluster. 
for (motif in markerMotifs) {
    marker <- unlist(strsplit(motif,':'))[2]
    p <- plotGroups(ArchRProj = projHeme5, 
        groupBy = "Clusters2", 
        colorBy = "MotifMatrix", 
        name = motif,
        imputeWeights = getImputeWeights(projHeme5)
    )
    pdf (file.path(outdir_chromvar,paste0(marker,".Plot-Groups-Deviations-w-Imputation.pdf")),width = 5,height = 5)
    print (p)
    dev.off()
    png (file.path(outdir_chromvar,paste0(marker,".Plot-Groups-Deviations-w-Imputation.png")))
    print (p)
    dev.off()
}


### 1.overlay the z-scores on our UMAP embedding as we’ve done previously for gene scores.
for (motif in sort(markerMotifs)) {
    p <- plotEmbedding(
        ArchRProj = projHeme5,
        colorBy = "MotifMatrix", 
        name = motif, 
        embedding = "UMAP",
        imputeWeights = getImputeWeights(projHeme5)
    )
    pdf (file.path(outdir_chromvar,paste0(motif,".motif.umap.pdf")),width = 5,height = 5)
    print (p)
    dev.off()
    png (file.path(outdir_chromvar,paste0(motif,".motif.umap.png")))
    print (p)
    dev.off()
}

## 2. overlay the gene scores for each of these TFs on the UMAP embedding
markerRNA <- getFeatures(projHeme5, select = paste(motifs, collapse="|"), useMatrix = "GeneScoreMatrix")
markerRNA <- markerRNA[markerRNA %in% motifs]
markerRNA
for (gene in sort(markerRNA)) {
    p <- plotEmbedding(
        ArchRProj = projHeme5, 
        colorBy = "GeneScoreMatrix", 
        name = gene, 
        embedding = "UMAP",
        imputeWeights = getImputeWeights(projHeme5)
    )
    pdf (file.path(outdir_chromvar,paste0(gene,".motif.link.genescore.umap.pdf")),width = 5,height = 5)
    print (p)
    dev.off()
    png (file.path(outdir_chromvar,paste0(gene,".motif.link.genescore.umap.png")))
    print (p)
    dev.off()  
}


# ## plot the linked gene expression for each of these TFs on the UMAP embedding 
markerRNA <- getFeatures(projHeme5, select = paste(motifs, collapse="|"), useMatrix = "GeneIntegrationMatrix")
markerRNA <- markerRNA[markerRNA %in% motifs]
markerRNA
for (gene in sort(markerRNA)) {
    p <- plotEmbedding(
        ArchRProj = projHeme5, 
        colorBy = "GeneIntegrationMatrix", 
        name = gene, 
        embedding = "UMAP",
        continuousSet = "blueYellow",
        imputeWeights = getImputeWeights(projHeme5)
    )
    pdf (file.path(outdir_chromvar,paste0(gene,".motif.link.geneIntegrationScore.pdf")),width = 5,height = 5)
    print (p)
    dev.off()
    png (file.path(outdir_chromvar,paste0(gene,".motif.link.geneIntegrationScore.png")))
    print (p)
    dev.off()  
}
saveArchRProject(ArchRProj = projHeme5, load = FALSE)
file.copy (file.path(outdir,"Save-ArchR-Project.rds"),file.path(outrds,"projHeme5.motif.rds"),recursive = TRUE)

## add readme
file.copy ("/SGRNJ03/PiplineTest01/Software_test/zhengzhanye/ArchR_V1.0/readme/ChromVAR.Variable.motif.readme.pdf",outdir_chromvar,recursive = TRUE)

######################################################################
## 这部分内容待定
# ##  Identification of Positive TF-Regulators
# outdir_positivetf <- "8.Positive_TF_Regulators"
# dir.create (file.path(outdir,"Plots",outdir_positivetf),showWarnings=T)
# ## Step 1. Identify Deviant TF Motifs
# seGroupMotif <- getGroupSE(ArchRProj = projHeme5, useMatrix = "MotifMatrix", groupBy = "Clusters2")
# ## subset this SummarizedExperiment to just the deviation z-scores.
# seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

# rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
# rowMaxs(assay(seZ) - assay(seZ)[,x])
# }) %>% Reduce("cbind", .) %>% rowMaxs

# ## Step 2. Identify Correlated TF Motifs and TF Gene Score/Expression
# corGSM_MM <- correlateMatrices(
#     ArchRProj = projHeme5,
#     useMatrix1 = "GeneScoreMatrix",
#     useMatrix2 = "MotifMatrix",
#     reducedDims = "IterativeLSI"
# )

# corGIM_MM <- correlateMatrices(
#     ArchRProj = projHeme5,
#     useMatrix1 = "GeneIntegrationMatrix",
#     useMatrix2 = "MotifMatrix",
#     reducedDims = "IterativeLSI"
# )

# ## Step 3. Add Maximum Delta Deviation to the Correlation Data Frame
# corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
# corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]

# ### Step 4. Identify Positive TF Regulators
# corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
# corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
# corGSM_MM$TFRegulator <- "NO"
# corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
# sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])

# Positive_TF_Regulators <- data.frame (positive_tf_regulators=sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1]))
# write.table (Positive_TF_Regulators,file.path(outdir,"Plots",outdir_positivetf,"Positive_TF_Regulators.xls"),quote=F,sep='\t',row.names=F,col.names=F)

# pdf (file.path(outdir,"Plots",outdir_positivetf,"Positive_TF_Regulators.pdf"),onefile = FALSE)
# p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
# geom_point() + 
# theme_ArchR() +
# geom_vline(xintercept = 0, lty = "dashed") + 
# scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
# xlab("Correlation To Gene Score") +
# ylab("Max TF Motif Delta") +
# scale_y_continuous(
#     expand = c(0,0), 
#     limits = c(0, max(corGSM_MM$maxDelta)*1.05)
# )
# data <- data.frame(corGSM_MM)
# data$label=ifelse(data$TFRegulator=="YES" & data$maxDelta > 7 ,data$MotifMatrix_matchName,"")
# p1 <- p+geom_text_repel(data = data, aes(x = cor, 
#                                 y = maxDelta, 
#                                 label = label),
#                         size = 3,box.padding = unit(0.5, "lines"),
#                         point.padding = unit(0.8, "lines"), 
#                         segment.color = "black", 
#                         show.legend = FALSE)

# print (p1)
# dev.off()
# saveArchRProject(ArchRProj = projHeme5, outputDirectory = "Save-ProjHeme5", load = FALSE)

######################################3

# ### Chapter 14 Footprinting with ArchR
# outdir_footprint <- "7.Footprint"
# dir.create (file.path(outdir,"Plots",outdir_footprint),showWarnings=T)
# motifPositions <- getPositions(projHeme5)
# ### subset this GRangesList to a few TF motifs that we are interested in
# markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
# # markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
# markerMotifs
# projHeme5 <- addGroupCoverages(ArchRProj = projHeme5, groupBy = "Clusters2",force = TRUE)
# ### compute footprints for the subset of marker motifs use getFootprints
# seFoot <- getFootprints(
#   ArchRProj = projHeme5, 
#   positions = motifPositions[markerMotifs], 
#   groupBy = "Clusters2"
# )

# ### 14.2 Normalization of Footprints for Tn5 Bias (两种normMethod 选择其中一种进行)
# ## 14.2.1 Subtracting the Tn5 Bias
# ## By default, these plots will be saved in the outputDirectory of the ArchRProject.
# plotFootprints(
#   seFoot = seFoot,
#   ArchRProj = projHeme5, 
#   normMethod = "Subtract",
#   plotName = "Footprints-Subtract-Bias.pdf",#file.path(outdir_footprint,"Footprints-Subtract-Bias.pdf"),
#   addDOC = FALSE,
#   smoothWindow = 5
# )
# ## 14.2.2 Dividing by the Tn5 Bias
# plotFootprints(
#   seFoot = seFoot,
#   ArchRProj = projHeme5, 
#   normMethod = "Divide",
#   plotName = "Footprints-Divide-Bias.pdf",##  file.path(outdir_footprint,"Footprints-Divide-Bias.pdf")
#   addDOC = FALSE,
#   smoothWindow = 5
# )


### Co-accessibility with ArchR
dir.create (file.path (outdir,"Plots","7.Co-accessibility"))
outdir_coaccess <- file.path (outdir,"Plots","7.Co-accessibility")
projHeme5 <- addCoAccessibility(
    ArchRProj = projHeme5,
    reducedDims = "IterativeLSI"
)
### retrieve this co-accessibility information from the ArchRProject  ( correlation column gives the numeric correlation of the accessibility between those two peaks)
cA <- getCoAccessibility(
    ArchRProj = projHeme5,
    corCutOff = 0.5,
    resolution = 1000,
    returnLoops = FALSE
)
allpeak.tpm <- DataFrame (iranges = metadata(cA)[[1]])
allpeak <- data.frame (celltype = rownames (allpeak.tpm),data.frame (allpeak.tpm$iranges)) 
allpeak <- data.frame (allpeak,peak = paste0(allpeak$seqnames,":",allpeak$start,"-",allpeak$end))%>% select (celltype,peak)

coaccess <- data.frame (cA) %>% select (queryHits,subjectHits,seqnames,correlation)
# copeak <- data.frame (query = allpeak[coaccess$queryHits,],subject = allpeak[coaccess$subjectHits,],coaccess[,c("seqnames","correlation")])
copeak <- data.frame (query = allpeak[coaccess$queryHits,],subject = allpeak[coaccess$subjectHits,],coaccess)
## 只保留queryHits subjectHits seqnames correlation这四列，然后根据Hits 和seqnames 两列从peak文件中提取peak位置信息，然后再保存这个表
write.table (copeak,file.path(outdir_coaccess,"co-accessibility.xls"),col.names=T,row.names=F,sep='\t',quote=F)

### Plotting browser tracks of Co-accessibility
markersGS <- getMarkerFeatures(
    ArchRProj = projHeme5, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters2",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
mm <- as.data.frame (markerList) 
write.table (mm,file.path(outdir_coaccess,"celltype.markerGeneList.xls"),col.names=T,row.names = F,sep='\t',quote=F)

percelltype_top1 <- group_by (mm,group) %>% sample_n (,size=1)
markerGenes <- percelltype_top1$name[1:5]
for (genes in markerGenes) {
    p <- plotBrowserTrack(
        ArchRProj = projHeme5, 
        groupBy = "Clusters2",      #geneSymbol = markerGenes, 
        geneSymbol = genes,
        upstream = 50000,
        downstream = 50000,
        loops = getCoAccessibility(projHeme5)
    )
    plotPDF(plotList = p, 
        name = file.path ("7.Co-accessibility",paste0 ("Plot-Tracks-",genes,"-with-CoAccessibility.pdf")),#name = paste0 ("Plot-Tracks-",genes,"-with-CoAccessibility.pdf"), 
        ArchRProj = projHeme5, 
        addDOC = FALSE, width = 7, height = 5
    )
}

## add readme
file.copy ("/SGRNJ03/PiplineTest01/Software_test/zhengzhanye/ArchR_V1.0/readme/Co-accessibility.readme.pdf",outdir_coaccess,recursive = TRUE)
###############################
### Peak2GeneLinkage with ArchR
###############################
dir.create (file.path(outdir,"Plots","8.Peak2GeneLinkage"),showWarnings=T)
outdir_peak2gene <- file.path(outdir,"Plots","8.Peak2GeneLinkage")
projHeme5 <- addPeak2GeneLinks(
    ArchRProj = projHeme5,
    reducedDims = "IterativeLSI"
)

### Plotting a heatmap of peak-to-gene links 
pdf (file.path (outdir_peak2gene,"Peak2GeneLinks-heatmap.pdf"))
plotPeak2GeneHeatmap(ArchRProj = projHeme5, groupBy = "Clusters2",k = 5,returnMatrices = FALSE)
dev.off()
# get heatmap 对应聚类结果
p <- plotPeak2GeneHeatmap(ArchRProj = projHeme5, groupBy = "Clusters2",k = 5,returnMatrices = TRUE)
hp <- data.frame (p$Peak2GeneLinks,kmeansId = p$ATAC$kmeansId)
hp <- hp[order(hp$kmeansId),]
write.table (hp,file.path(outdir_peak2gene,"peak2GeneLinks.xls"),col.names=T,row.names=F,sep='\t',quote=F)

## 这部分得到结果fp2g 与上述hp相同 (hp 中有kmeanID 的信息)
# p2g <- getPeak2GeneLinks(
#     ArchRProj = projHeme5,
#     corCutOff = 0.45,
#     resolution = 1,
#     returnLoops = FALSE
# )
# ## 提取idxATAC 和idxRNA  对应peak和gene
# mATAC <- readRDS(metadata(p2g)$seATAC)[p2g$idxATAC, ]
# mRNA <- readRDS(metadata(p2g)$seRNA)[p2g$idxRNA, ]
# p2g$peak <- paste0(rowRanges(mATAC))
# p2g$gene <- rowData(mRNA)$name
# fp2g <- data.frame (p2g)
# write.table (fp2g,file.path(outdir_peak2gene,"peak2GeneLinks.xls"),col.names=T,row.names=F,sep='\t',quote=F)
## Plotting browser tracks with peak-to-gene links outdir_peak2gene (## 按照Correlation进行排序,筛选前几个gene展示)
p2g.genes <- hp[order(hp$Correlation),] %>%filter (gene != "NA") %>% select (gene) %>% unique()
p2g.genes.use <- p2g.genes[1:5,]   
for (gene in p2g.genes.use) {
    p <- plotBrowserTrack(
        ArchRProj = projHeme5, 
        groupBy = "Clusters2", 
        geneSymbol = gene, 
        upstream = 50000,
        downstream = 50000,
        loops = getPeak2GeneLinks(projHeme5)
    )
    plotPDF(plotList = p, 
        name = file.path ("8.Peak2GeneLinkage",paste0("Tracks-",gene,"-Peak2GeneLinks.pdf")), #name = paste0("Plot-Tracks-",gene,"-with-Peak2GeneLinks.pdf"),
        ArchRProj = projHeme5, 
        addDOC = FALSE, width = 7, height = 5

    )
}

# saveArchRProject(ArchRProj = projHeme5)
# file.copy (file.path(outdir,"Save-ArchR-Project.rds"),file.path(outrds,"projHeme6.coaccess.peak2gene.rds"),recursive = TRUE)
## add readme
file.copy ("/SGRNJ03/PiplineTest01/Software_test/zhengzhanye/ArchR_V1.0/readme/Peak2GeneLinkage.readme.pdf",outdir_peak2gene,recursive = TRUE)

#####
###################################################
## Trajectory Analysis with ArchR
###################################################
dir.create (file.path(outdir,"Plots","9.Trajectory"),showWarnings=T)
outdir_trajectory <- file.path(outdir,"Plots","9.Trajectory")
## p2 <- plotEmbedding(ArchRProj = projHeme5, colorBy = "cellColData", name = "Clusters2", embedding = "UMAP")##只是为了展示降维结果
## Pseudo-time UMAPs and individual feature plots
#trajectory <- c("Progenitor", "GMP", "Mono")
projHeme5 <- addImputeWeights(projHeme5)
trajectory <- unique(projHeme5$Clusters2)[c(1,2,3)]
trajectory
### call this trajectory “MyeloidU”. Cells that are not part of the trajectory are labeled with NA.
## “MyeloidU” that stores the pseudo-time value for each cell in the trajectory.
projHeme5 <- addTrajectory(
    ArchRProj = projHeme5, 
    name = "MyeloidU", 
    groupBy = "Clusters2",
    trajectory = trajectory, 
    embedding = "UMAP", 
    force = TRUE
)
### exclude cells with NA values because these are not part of the trajectory.
#head(projHeme5$MyeloidU[!is.na(projHeme5$MyeloidU)])
p <- plotTrajectory(projHeme5, trajectory = "MyeloidU", colorBy = "cellColData", name = "MyeloidU")
pdf (file.path(outdir_trajectory,"Plot-Selected-Traj-UMAP.pdf"),width = 5,height = 5)
p[[1]]
dev.off()
png (file.path(outdir_trajectory,"Plot-Selected-Traj-UMAP.png"),width = 300,height = 300)
p[[1]]
dev.off()


###############################
## Pseudo-time heatmaps
###############################
## create these pseudo-time heatmaps for motifs by passing the corresponding matrix to the useMatrix parameter
trajMM  <- getTrajectory(ArchRProj = projHeme5, name = "MyeloidU", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"),labelTop = 10,returnMatrix = FALSE) ##labelTop 图中展示几个label
pdf (file.path(outdir_trajectory,"motifs-Traj-Heatmaps.pdf"))
print (p1)
dev.off()
png (file.path(outdir_trajectory,"motifs-Traj-Heatmaps.png"))
print (p1)
dev.off()
p1data <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"),labelTop = 10,returnMatrix = TRUE)
write.table (p1data,file.path(outdir_trajectory,"motifs-Traj-Heatmaps.xls"),row.names = T,col.names = T,sep='\t',quote=F)
## create these pseudo-time heatmaps gene scores
trajGSM <- getTrajectory(ArchRProj = projHeme5, name = "MyeloidU", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- trajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"),labelTop = 10,returnMatrix = FALSE)
pdf (file.path(outdir_trajectory,"GeneScore-Traj-Heatmaps.pdf"))
print (p2)
dev.off()
png (file.path(outdir_trajectory,"GeneScore-Traj-Heatmaps.png"))
print (p2)
dev.off()
p2data <- plotTrajectoryHeatmap(trajGSM, pal = paletteContinuous(set = "solarExtra"),labelTop = 10,returnMatrix = TRUE)
write.table (p2data,file.path(outdir_trajectory,"GeneScore-Traj-Heatmaps.xls"),row.names = T,col.names = T,sep='\t',quote=F)

## create these pseudo-time heatmaps for gene expression
trajGIM <- getTrajectory(ArchRProj = projHeme5, name = "MyeloidU", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"),labelTop = 10,returnMatrix = FALSE)
pdf (file.path(outdir_trajectory,"GeneIntegration-Traj-Heatmaps.pdf"))
print (p3)
dev.off()
png (file.path(outdir_trajectory,"GeneIntegration-Traj-Heatmaps.png"))
print (p3)
dev.off()
p3data <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"),labelTop = 10,returnMatrix = TRUE)
write.table (p3data,file.path(outdir_trajectory,"GeneIntegration-Traj-Heatmaps.xls"),row.names = T,col.names = T,sep='\t',quote=F)
## ## create these pseudo-time heatmaps for peak accessibility
trajPM  <- getTrajectory(ArchRProj = projHeme5, name = "MyeloidU", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"),labelTop = 10,returnMatrix = FALSE)
pdf (file.path(outdir_trajectory,"PeakMatrix-Traj-Heatmaps.pdf"))
print (p4)
dev.off()
png (file.path(outdir_trajectory,"PeakMatrix-Traj-Heatmaps.png"))
print (p4)
dev.off()
p4data <- plotTrajectoryHeatmap(trajPM,  pal = paletteContinuous(set = "blueYellow"),labelTop = 10,returnMatrix = TRUE)
write.table (p4data,file.path(outdir_trajectory,"PeakMatrix-Traj-Heatmaps.xls"),row.names = T,col.names = T,sep='\t',quote=F)
saveArchRProject(ArchRProj = projHeme5)
file.copy (file.path(outdir,"Save-ArchR-Project.rds"),file.path(outrds,"projHeme7.traj.rds"),recursive = TRUE)
#saveArchRProject(ArchRProj = projHeme5, outputDirectory = "Save-Co-accessibility", load = FALSE)



########################################################################
## gene expression within the cells that are relevant to our trajectory
########################################################################
gene <- strsplit (rownames (p2data)[1],":")[[1]][2]
## GeneScoreMatrix 
p1 <- plotTrajectory(projHeme5, trajectory = "MyeloidU", colorBy = "GeneScoreMatrix", name =gene, continuousSet = "horizonExtra")
pdf (file.path (outdir_trajectory,paste0(gene,".pseudo-time.UMAP.GeneScoreMatrix.pdf")),width = 5,height = 5)
p1[[1]]
dev.off()
png (file.path (outdir_trajectory,paste0(gene,".pseudo-time.UMAP.GeneScoreMatrix.png")),width = 300,height = 300)
p1[[1]]
dev.off()
pdf (file.path (outdir_trajectory,paste0(gene,".pseudo-time.GeneScoreMatrix.pdf")),width = 5,height = 5)
p1[[2]]
dev.off()
png (file.path (outdir_trajectory,paste0(gene,".pseudo-time.GeneScoreMatrix.png")),width = 300,height = 300)
p1[[2]]
dev.off()
## GeneIntegrationMatrix
gene <- strsplit (rownames (p3data)[1],":")[[1]][2]
p2 <- plotTrajectory(projHeme5, trajectory = "MyeloidU", colorBy = "GeneIntegrationMatrix", name =gene, continuousSet = "blueYellow")
pdf (file.path (outdir_trajectory,paste0(gene,".pseudo-time.UMAP.GeneIntegrationMatrix.pdf")),width = 5,height = 5)
p2[[1]]
dev.off()
png (file.path (outdir_trajectory,paste0(gene,".pseudo-time.UMAP.GeneIntegrationMatrix.png")),width = 300,height = 300)
p2[[1]]
dev.off()
pdf (file.path (outdir_trajectory,paste0(gene,".pseudo-time.GeneIntegrationMatrix.pdf")),width = 5,height = 5)
p2[[2]]
dev.off()
png (file.path (outdir_trajectory,paste0(gene,".pseudo-time.GeneIntegrationMatrix.png")),width = 300,height = 300)
p2[[2]]
dev.off()

## add readme
file.copy ("/SGRNJ03/PiplineTest01/Software_test/zhengzhanye/ArchR_V1.0/readme/Monocle3.readme.pdf",outdir_trajectory,recursive = TRUE)
