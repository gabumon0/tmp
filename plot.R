
setwd(outdir)
outdir_cellanno <- "3.Cellannotation"
dianno <- file.path(outdir,"Plots",outdir_cellanno)
anno.meta <- projHeme3@cellColData
anno.meta$Clusters2 <- as.factor (anno.meta$Clusters2)
anno.meta$Sample <- factor (anno.meta$Sample,levels = mixedsort (unique (anno.meta$Sample)))

projHeme3 <- readRDS('0.rds/projHeme4.peak.rds')

data.umap <- data.frame(barcode = projHeme3$cellNames,celltype = anno.meta$Clusters2,umap.file)


p3 <- ggplot (data.umap,mapping = aes (x=UMAP_1,y=UMAP_2,color=celltype)) + geom_point(size = 0.2)  + scale_color_manual(values = clustcol) + 
        theme_classic() + theme(legend.title=element_blank(),legend.text = element_text (size = 10)) +
        guides(color = guide_legend(override.aes = list(size = 2)))    
pdf (file.path(dianno,paste0(prefix,".labumap.integration.pdf")))
p3
dev.off()
png (file.path(dianno,paste0(prefix,".labumap.integration.png")))
p3
dev.off()

tsne.file <- projHeme3@embeddings$TSNE$df
colnames (tsne.file) <- c("TSNE_1","TSNE_2")
write.table (tsne.file,file.path(dianno,paste0(prefix,".TSNE_Dimension.file.xls")),col.names = NA,sep='\t',quote = F)
data.tsne <- data.frame(barcode = projHeme3$cellNames,celltype = anno.meta$Clusters2,tsne.file)
p3 <- ggplot (data.tsne,mapping = aes (x=TSNE_1,y=TSNE_2,color=celltype)) + geom_point(size = 0.2)  + scale_color_manual(values = clustcol) + 
        theme_classic() + theme(legend.title=element_blank(),legend.text = element_text (size = 10)) +
        guides(color = guide_legend(override.aes = list(size = 2)))
pdf (file.path(dianno,paste0(prefix,".labtsne.integration.pdf")))
p3
dev.off()
png (file.path(dianno,paste0(prefix,".labtsne.integration.png")))
p3
dev.off()


cluster <- table (anno.meta$Clusters2,anno.meta$Sample)
write.table (cluster,file.path(dianno,paste0(prefix,".CellsPerCellType.xls")),col.names = NA,sep='\t',quote = F)
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

