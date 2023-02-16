
setwd(outdir)
outdir_cellanno <- "3.Cellannotation"
dianno <- file.path(outdir,"Plots",outdir_cellanno)
projHeme3 <- readRDS('0.rds/projHeme4.peak.rds')

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

annors.out <- file.path(dianno,prefix)
dir.create (annors.out)
file.copy (file.path(dianno,paste0(prefix,".labumap.integration.pdf")),annors.out,recursive = TRUE)
file.copy (file.path(dianno,paste0(prefix,".labumap.integration.png")),annors.out,recursive = TRUE)
file.copy (file.path(dianno,paste0(prefix,".CellsPropPerCellTypes.pdf")),annors.out,recursive = TRUE)
file.copy (file.path(dianno,paste0(prefix,".CellsPropPerCellTypes.png")),annors.out,recursive = TRUE)
file.copy (file.path(dianno,paste0(prefix,".CellsPerCellType.xls")),annors.out,recursive = TRUE)
file.copy (file.path(dianno,paste0(prefix,".labumap.sample.pdf")),annors.out,recursive = TRUE)
file.copy (file.path(dianno,paste0(prefix,".labumap.sample.png")),annors.out,recursive = TRUE)

