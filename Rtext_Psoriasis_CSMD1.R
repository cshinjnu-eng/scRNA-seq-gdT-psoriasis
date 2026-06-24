#######1.Initialize#################
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
setwd("F:\\03_Rtest\\5.psoriasis/Monocle/")
library(Seurat)
##options(Seurat.object.assay.version = 'v5')
library(usethis)
library(devtools)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(batchelor)
#remotes::install_github('satijalab/seurat-wrappers@community-vignette')
library(SeuratWrappers)
library(magrittr)
library(tidyverse)
library(clusterProfiler)
library(GO.db)
library(org.Hs.eg.db)
library(DOSE)
library(DoubletFinder)
library(SingleR)
library(celldex)
library(Rcpp)
library(harmony)
library(pheatmap)
#install.packages("scCustomize")
library(scCustomize)
library(RColorBrewer)
dir.create("SeuratV4")
install.packages('Seurat', repos = c('https://satijalab.r-universe.dev'), lib = "SeuratV4")
.libPaths(c("SeuratV4", .libPaths()))
library(Seurat)
packageVersion("Seurat")

set.seed(42)
#Read data
#Method 1
obj_dir <- c("NPBMC_gdt1/","NPBMC_gdt2/","NPBMC_gdt3/","PBMC_gdt1/","PBMC_gdt2/","PBMC_gdt3/")
names(obj_dir) = c("NPBMC1","NPBMC2","NPBMC3","PBMC1","PBMC2","PBMC3")

counts <- Read10X(data.dir = obj_dir)
scRNA1 = CreateSeuratObject(counts, min.cells=1)

dim(scRNA1) 
table(scRNA1@meta.data$orig.ident) 
scRNAlist <- SplitObject(scRNA1, split.by = "orig.ident")
#scRNAlist


##Method 2
obj = SplitObject(scRNA1, split.by = "orig.ident")
obj_rm=list()
doublets_plot = list()
pc.num = 1:20#
dir.create("SingleCell_QC")
##########################2.Remove double cells#############################
RemoveDoublets <-function(
    object,
    doublet.rate,
    pN=0.25,
    pc.num=1:30
){
  sweep.res.list <- paramSweep(object, PCs = pc.num, sct = F)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  homotypic.prop <- modelHomotypic(object$seurat_clusters)
  nExp_poi <- round(doublet.rate*ncol(object)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  seu.scored <- doubletFinder(object, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, 
                              nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
  cname <-colnames(seu.scored[[]])
  DF<-cname[grep('^DF',cname)]
  seu.scored[["doublet"]] <- as.numeric(seu.scored[[DF]]=="Doublet")
  seu.removed <- subset(seu.scored, subset = doublet != 1)
  p1 <- DimPlot(seu.scored, group.by = DF)
  res.list <- list("plot"=p1, "obj"=seu.removed)
  return(res.list)
}

for( i in names(obj)){
  obj[[i]] <- NormalizeData(obj[[i]])
  obj[[i]] <- FindVariableFeatures(obj[[i]], selection.method = "vst", nfeatures = 2000)
  obj[[i]] <- ScaleData(obj[[i]])
  obj[[i]] <- RunPCA(obj[[i]])
  obj[[i]] <- RunUMAP(obj[[i]], dims = 1:20)
  obj[[i]] <- FindNeighbors(obj[[i]], dims = pc.num) %>% FindClusters(resolution = 0.3)#

  cell_count <- ncol(obj[[i]])
  if ( cell_count < 750 ){
    doublet_rate = 0.004
  }else if ( cell_count %in% 750:1499 ) {
    doublet_rate = 0.008
  }else if ( cell_count %in% 1500:2499 ) {
    doublet_rate = 0.016
  }else if ( cell_count %in% 2500:3499 ) {
    doublet_rate = 0.023
  }else if ( cell_count %in% 3500:4499 ) {
    doublet_rate = 0.031
  }else if ( cell_count %in% 4500:5499 ) {
    doublet_rate = 0.039
  }else if ( cell_count %in% 5500:6499 ) {
    doublet_rate = 0.046
  }else if ( cell_count %in% 6500:7499 ) {
    doublet_rate = 0.054
  }else if ( cell_count %in% 7500:8499 ) {
    doublet_rate = 0.061
  }else if ( cell_count %in% 8500:9499 ) {
    doublet_rate = 0.069
  }else if ( cell_count %in% 9500:10499 ) {
    doublet_rate = 0.076
  }else if ( cell_count >= 10500 ) {
    doublet_rate = 0.1
  }
  tmp <- RemoveDoublets(obj[[i]], doublet.rate=doublet_rate, pc.num=pc.num)
  obj_rm[[i]] <- tmp$obj
  doublets_plot[[i]] <- tmp$plot
}

p1=doublets_plot[["NPBMC1"]]
p2=doublets_plot[["NPBMC2"]]
p3=doublets_plot[["NPBMC3"]]
p4=doublets_plot[["PBMC1"]]
p5=doublets_plot[["PBMC2"]]
p6=doublets_plot[["PBMC3"]]

dir.create("Doublets")
ggsave("Doublets/NPBMC1_before_MNN.pdf", plot = p1, width = 6, height = 5,dpi = 300) 
ggsave("Doublets/NPBMC1_before_MNN.png", plot = p1, width = 6, height = 5,dpi = 600)
ggsave("Doublets/NPBMC2_before_MNN.pdf", plot = p2, width = 6, height = 5,dpi = 300) 
ggsave("Doublets/NPBMC2_before_MNN.png", plot = p2, width = 6, height = 5,dpi = 600)
ggsave("Doublets/NPBMC3_before_MNN.pdf", plot = p3, width = 6, height = 5,dpi = 300) 
ggsave("Doublets/NPBMC3_before_MNN.png", plot = p3, width = 6, height = 5,dpi = 600)
ggsave("Doublets/PBMC1_before_MNN.pdf", plot = p4, width = 6, height = 5,dpi = 300) 
ggsave("Doublets/PBMC1_before_MNN.png", plot = p4, width = 6, height = 5,dpi = 600)
ggsave("Doublets/PBMC2_before_MNN.pdf", plot = p5, width = 6, height = 5,dpi = 300) 
ggsave("Doublets/PBMC2_before_MNN.png", plot = p5, width = 6, height = 5,dpi = 600)
ggsave("Doublets/PBMC3_before_MNN.pdf", plot = p6, width = 6, height = 5,dpi = 300) 
ggsave("Doublets/PBMC3_before_MNN.png", plot = p6, width = 6, height = 5,dpi = 600)


#######3.Standard Quality Control#####
scRNA <- merge(obj_rm[[1]], y=c(obj_rm[[2]],obj_rm[[3]],obj_rm[[4]],obj_rm[[5]],obj_rm[[6]]))
Idents(scRNA) <- 'orig.ident'
scRNA[["percent.MT"]] = PercentageFeatureSet(scRNA, 
                                             pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) 
beforeQC_vlnplot = VlnPlot(scRNA, 
                           features = c("nFeature_RNA", 
                                        "nCount_RNA", 
                                        "percent.MT",
                                        "percent.HB"), 
                           ncol = 2, 
                           pt.size = 0.01)

dir.create("SingleCell_QC")
ggsave("SingleCell_QC/BeforeQC_nFeature_nCount_percent.mt_percent.HB_vlnplot.pdf", 
       plot = beforeQC_vlnplot,width = 15, height = 15,dpi = 300)
ggsave("SingleCell_QC/BeforeQC_nFeature_nCount_percent.mt_percent.HB_vlnplot.png",
       plot = beforeQC_vlnplot,width = 15, height = 15,dpi = 600)

summary(scRNA@meta.data[,c("nFeature_RNA","nCount_RNA","percent.MT","percent.HB")])

plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.MT")
plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.HB")

pearplot <- CombinePlots(plots = list(plot1, plot2, plot3), nrow=1, legend="none") 
ggsave("SingleCell_QC/pearplot_before_qc.pdf", plot = pearplot, width = 15, height = 5,dpi = 300) 
ggsave("SingleCell_QC/pearplot_before_qc.png", plot = pearplot, width = 15, height = 5,dpi = 600)


minGene=200
maxGene=8000 
pctMT=10
pctHB=3

scRNA = subset(scRNA, 
               subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.MT < 10 & percent.HB < 3)

afterQC_vlnplot = VlnPlot(scRNA, 
                          features = c("nFeature_RNA", 
                                       # "nCount_RNA",
                                       "percent.MT",
                                       "percent.HB"), 
                          ncol = 2, 
                          pt.size = 0.01)

ggsave("SingleCell_QC/afterQC_nFeature_nCount_percent.mt_percent.HB_vlnplot.pdf", plot = afterQC_vlnplot
       ,width = 15, height = 15,dpi = 300)
ggsave("SingleCell_QC/afterQC_nFeature_nCount_percent.mt_percent.HB_vlnplot.png", plot = afterQC_vlnplot
       ,width = 15, height = 15,dpi = 600)


scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst")
scRNA <- ScaleData(scRNA, features = VariableFeatures(scRNA))



##########4.Dimensionality Reduction Clustering###########

scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA)) 
plot1 <- DimPlot(scRNA, reduction = "pca", group.by="orig.ident")
plot2 <- ElbowPlot(scRNA, ndims=30, reduction="pca") 

plot3 <- plot1+plot2
dir.create("Clustering")
ggsave("Clustering/pca.png", plot = plot3, width = 8, height = 4)
plot3

scRNA <- FindNeighbors(scRNA, dims = 1:15) #KNN + SNN
scRNA <- FindClusters(scRNA, resolution = 0.5)#Louvain
table(scRNA@meta.data$seurat_clusters)
metadata <- scRNA@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'Clustering/cell_cluster.csv',row.names = F, quote = F)
#UMAP
scRNA <- RunUMAP(scRNA, dims = 1:15)
#group_by_cluster
plot3 = DimPlot(scRNA, reduction = "umap", label=T)
ggsave("Clustering/preUMAP.png",plot3,width = 7, height = 7, dpi = 600)
ggsave("Clustering/preUMAP.pdf",plot3,width = 7, height = 7, dpi = 300)
plot3
plot4 = DimPlot(scRNA, reduction = "umap", group.by='orig.ident')
plot5 = DimPlot(scRNA, reduction = "umap", split.by='orig.ident')

plotc <- plot4+plot5 
ggsave("Clustering/preUMAP_cluster_sample.png", plot = plotc, width = 10, height = 5)
plotc
saveRDS(scRNA,"Step3_Clustering.rds")

#################5.Remove batch effects##############


scRNAlist <- SplitObject(scRNA, split.by = "orig.ident")
scRNAlist <- lapply(scRNAlist, FUN = function(x) NormalizeData(x))
scRNAlist <- lapply(scRNAlist, FUN = function(x) FindVariableFeatures(x))
scRNA_mnn <- RunFastMNN(object.list = scRNAlist)
scRNA_mnn <- FindVariableFeatures(scRNA_mnn)
scRNA_mnn <- RunUMAP(scRNA_mnn, reduction = "mnn", dims = 1:15)
scRNA_mnn <- FindNeighbors(scRNA_mnn, reduction = "mnn", dims = 1:15)#SNN + KNN
scRNA_mnn <- FindClusters(scRNA_mnn)#Louvain
p1 <- DimPlot(scRNA_mnn, group.by = "orig.ident", pt.size=0.1) + 
  ggtitle("Integrated by fastMNN")
p2 <- DimPlot(scRNA, group.by="orig.ident", pt.size=0.1) + 
  ggtitle("No integrated")
p = p1 + p2 + plot_layout(guides='collect')
saveRDS(scRNA_mnn,"Step4_MNN.rds")#
table(scRNA_mnn@meta.data$orig.ident)
group <- data.table::fread("group.csv", header = TRUE)
head(group)
metadata <- FetchData(scRNA_mnn, "orig.ident")
metadata$cell_id <- rownames(metadata)
metadata <- left_join(x = metadata, y = group, by = "orig.ident")
rownames(metadata) <- metadata$cell_id
scRNA_mnn <- AddMetaData(scRNA_mnn, metadata = metadata)
table(scRNA_mnn@meta.data$group)

####over


#########6.Basic Analysis and Cell Type Identification###############

mycolor1<-c(
  "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
  "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
  "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93")

mycolor2<-c(
  "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F","#FF6A6A","#7FFFD4", "#AB82FF","#90EE90",
  "#00CD00","#008B8B","#6495ED","#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00","#000000","#708090")

mycolor3<-c(
  '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', 
  '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', 
  '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', 
  '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')

mycolor4=c(
  "#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
  "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
  "#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
  "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

p1 <- DimPlot(scRNA_mnn, group.by = "orig.ident", pt.size=0.1,shuffle = T)+scale_color_manual(values = mycolor1)
p1 
p3 <- DimPlot(scRNA_mnn, pt.size=0.1,label=T)+scale_color_manual(values = mycolor1)
p3
p4 <- DimPlot(scRNA_mnn, split.by="orig.ident", pt.size=0.1, label = T,shuffle = T)+scale_color_manual(values = mycolor1)
p4
p = p3 + p4 + plot_layout(guides='collect')
p
ggsave('Clustering/origident.png', p1, width=8, height=6, dpi=600)
ggsave('Clustering/origident.pdf', p1, width=8, height=6, dpi=300)
ggsave('Clustering/cluster.png', p3, width=8, height=6, dpi=600)
ggsave('Clustering/cluster.pdf', p3, width=8, height=6, dpi=300)
ggsave('Clustering/sample.png', p4, width=24, height=3, dpi=600)
ggsave('Clustering/sample.pdf', p4, width=24, height=3, dpi=300)
ggsave('Clustering/fastMNN.png', p, width=15, height=3, dpi=600)

p5 <- DimPlot(scRNA_mnn, group.by = "group", pt.size=0.1,shuffle = T)+scale_color_manual(values = mycolor1)
p5 
p6 <- DimPlot(scRNA_mnn, split.by="group", pt.size=0.1, label = T,shuffle = T)+scale_color_manual(values = mycolor1)
p6
p7 = p3 + p6 + plot_layout(guides='collect')
p7
ggsave('Clustering/group.png', p5, width=8, height=6, dpi=600)
ggsave('Clustering/group.pdf', p5, width=8, height=6, dpi=300)
ggsave('Clustering/splited_group.png', p6, width=12, height=6, dpi=600)
ggsave('Clustering/splited_group.pdf', p6, width=12, height=6, dpi=300)
ggsave('Clustering/combined_group.png', p7, width=15, height=4, dpi=600)
ggsave('Clustering/combined_group.pdf', p7, width=15, height=4, dpi=300)

sample_table <- as.data.frame(table(scRNA_mnn@meta.data$orig.ident,scRNA_mnn@meta.data$seurat_clusters))
names(sample_table) <- c("orig.ident","seurat_clusters","CellNumber")
colour = mycolor1
plot_sample_1 <-ggplot(sample_table,aes(x=orig.ident,weight=CellNumber,fill=seurat_clusters))+
  geom_bar(position="fill")+
  scale_fill_manual(values=colour) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")+RotatedAxis()
plot_sample_1

ggsave('Clustering/cluster_rate_bysamples.pdf', plot_sample_1, width=5, height=6, dpi=300)
ggsave('Clustering/cluster_rate_bysamples.png', plot_sample_1, width=5, height=6, dpi=600)
write.table(sample_table, 
            "Clustering/cluster_rate_bysamples.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")

group_table <- as.data.frame(table(scRNA_mnn@meta.data$group,scRNA_mnn@meta.data$seurat_clusters))
names(group_table) <- c("group","seurat_clusters","CellNumber")
colour = mycolor1
plot_sample_2<-ggplot(group_table,aes(x=group,weight=CellNumber,fill=seurat_clusters))+
  geom_bar(position="fill")+
  scale_fill_manual(values=colour) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")+RotatedAxis()
plot_sample_2

ggsave('Clustering/cluster_rate_bygroups.pdf', plot_sample_2, width=4, height=7, dpi=300)
ggsave('Clustering/cluster_rate_bygroups.png', plot_sample_2, width=4, height=7, dpi=600)
write.table(group_table, 
            "Clustering/cluster_rate_bygroups.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")

dir.create("Marker")
all.markers = FindAllMarkers(scRNA_mnn, 
                             min.pct = 0.25, 
                             logfc.threshold = 0.25, 
                             only.pos = TRUE)
head(all.markers)

top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top5 = all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.table(all.markers, 
            "Marker/all_Markers_of_each_clusters.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")
write.table(top10, 
            "Marker/top10_Markers_of_each_clusters.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")

scRNA_mnn <- ScaleData(scRNA_mnn, features = row.names(scRNA_mnn))
heatmap_plot1 = DoHeatmap(object = scRNA_mnn, label = F , 
                          features = as.character(unique(top10$gene)),   
                          group.by = "seurat_clusters",  
                          assay = "RNA",  
                          group.colors = mycolor1)+ 
  scale_fill_gradientn(colors = c("navy","white","firebrick3"))+ 
  theme(axis.text.y = element_text(size = 4))

heatmap_plot2 = DoHeatmap(object = scRNA_mnn, label = F , 
                          features = as.character(unique(top10$gene)),   
                          group.by = "orig.ident",  
                          assay = "RNA",  
                          group.colors =mycolor2)+ 
  scale_fill_gradientn(colors = c("navy","white","firebrick3"))+ 
  theme(axis.text.y = element_text(size = 4))

ggsave("Marker/top10_marker_of_each_cluster_heatmap.pdf", width = 6, height = 8,dpi = 300,
       plot = heatmap_plot1)
ggsave("Marker/top10_marker_of_each_cluster_heatmap.png", width = 6, height = 8, dpi = 600,
       plot = heatmap_plot1)
ggsave("Marker/top10_marker_of_each_sample_heatmap.pdf", width = 6, height = 8,dpi = 300,
       plot = heatmap_plot2)
ggsave("Marker/top10_marker_of_each_sample_heatmap.png", width = 6, height = 8, dpi = 600,
       plot = heatmap_plot2)

VlnPlot1 = VlnPlot(scRNA_mnn, 
                   cols = mycolor1, 
                   features = "LAG3", 
                   pt.size = 0.01, 
                   group.by = "seurat_clusters")
VlnPlot2 = VlnPlot(scRNA_mnn, 
                   cols = mycolor2, 
                   features = "AC105402.3", 
                   pt.size = 0.01, 
                   group.by = "group")
VlnPlot3 = VlnPlot(scRNA_mnn, 
                   cols = mycolor1, 
                   features = c("AC105402.3"), 
                   pt.size = 0.01, 
                   split.plot = T,
                   split.by = "group")
VlnPlot4 = VlnPlot(scRNA_mnn, 
                   cols = mycolor2, 
                   features = "STEAP1B", 
                   pt.size = 0.01, 
                   group.by = "orig.ident")
VlnPlot5 = VlnPlot(scRNA_mnn, 
                   cols = mycolor1, 
                   features = c("STEAP1B"), 
                   pt.size = 0.01, 
                   split.plot = T,
                   split.by = "orig.ident")

dir.create("./vlnplot")
extra_gene = read.table("genelist.txt", sep = "\t", header = T)
formated_gene = CaseMatch(search = as.vector(extra_gene$gene),
                          match = rownames(x = scRNA_mnn))
for (gene in as.vector(formated_gene)) {
  vlnplot = VlnPlot(scRNA_mnn, 
                    features = gene, 
                    pt.size = 0, 
                    split.plot = T,
                    split.by = "orig.ident")
  ggsave(paste0("vlnplot/",gene,"_vlnplot.pdf"), width = 10.11, height = 3.36, plot = vlnplot)
  ggsave(paste0("vlnplot/",gene,"_vlnplot.png"), width = 10.11, height = 3.36, plot = vlnplot)
}

FeaturePlot(scRNA_mnn, 
            features = c("DUSP1"), 
            label = F,
            pt.size = 0.4,
            split.by = "group",
            reduction = "umap",
) + theme(legend.position = "right") & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) 

FeaturePlot_scCustom(scRNA_mnn, 
                     features = "DUSP1", 
                     label = F,
                     pt.size = 0.4,
                     split.by = "orig.ident",
                     reduction = "umap")
FeaturePlot_scCustom(scRNA_mnn, 
                     features = "DUSP1", 
                     label = F,
                     pt.size = 0.4,
                     split.by = "group",
                     reduction = "umap")
FeaturePlot_scCustom(scRNA_mnn, features = c("CD4","FOXP3","IL2RA"), pt.size = 1.5,split.by = "orig.ident" )
FeaturePlot_scCustom(scRNA_mnn, features = "Ngp", pt.size = 2,split.by = "orig.ident" )
FeaturePlot(scRNA_mnn, 
            features = c("CD4","FOXP3"), pt.size = 0.4, blend = TRUE)

dir.create("./featureplot_kelimei")
extra_gene = read.table("genelist_kelimei.txt", sep = "\t", header = T)
formated_gene = CaseMatch(search = as.vector(extra_gene$gene),
                          match = rownames(scRNA_mnn))
for (gene in as.vector(formated_gene)) {
  featureplot = FeaturePlot_scCustom(scRNA_mnn, 
                                     features = gene, 
                                     label = F,
                                     split.by = "group",
                                     pt.size = 0.1,
                                     reduction = "umap")
  ggsave(paste0("featureplot_kelimei/",gene,"featureplot.pdf"), width = 10, height = 4,plot = featureplot,dpi = 300)
  ggsave(paste0("featureplot_kelimei/",gene,"featureplot.png"), width = 10, height = 4,plot = featureplot,dpi = 600)
}

dir.create("Marker/featureplot/")
for (gene in as.vector(top1$gene)) {
  featureplot = FeaturePlot(scRNA_mnn, 
                            features = gene, 
                            pt.size = 0.5, 
                            cols = c("lightgrey", "#DE1F1F"))
  ggsave(paste0("Marker/featureplot/",gene,"_featureplot.pdf"), plot = featureplot)
  ggsave(paste0("Marker/featureplot/",gene,"_featureplot.png"), plot = featureplot)
}

Idents(scRNA_mnn) <- 'group'
Diff_exp = FindAllMarkers(scRNA_mnn, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25, 
                          only.pos = TRUE)
head(Diff_exp)
top10 = Diff_exp %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top20 = Diff_exp %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

dir.create("Diffexp")
write.table(Diff_exp, 
            "Diffexp/group_psoriasis-vs-control-diff.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")
write.table(top10, 
            "Diffexp/group_psoriasis-vs-control-top10-diff.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")
write.table(top20, 
            "Diffexp/group_psoriasis-vs-control-top20-diff.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")
heatmap_plot_top10 = DoHeatmap(object = scRNA_mnn, label = F , 
                               features = as.character(unique(top10$gene)),   
                               group.by = "group",  
                               assay = "RNA",  
                               group.colors = mycolor2)+ 
  scale_fill_gradientn(colors = c("navy","white","firebrick3"))+ 
  theme(axis.text.y = element_text(size = 6))
ggsave("Diffexp/top10_marker_of_each_sample_heatmap.pdf", width = 6, height = 8, dpi = 300,
       plot = heatmap_plot_top10)
ggsave("Diffexp/top10_marker_of_each_sample_heatmap.png", width = 6, height = 8, dpi = 600,
       plot = heatmap_plot_top10)
heatmap_plot_top20 = DoHeatmap(object = scRNA_mnn, label = F , 
                               features = as.character(unique(top20$gene)),   
                               group.by = "group",  
                               assay = "RNA",  
                               group.colors = mycolor2)+ 
  scale_fill_gradientn(colors = c("navy","white","firebrick3"))+ 
  theme(axis.text.y = element_text(size = 6)) 
ggsave("Diffexp/top20_marker_of_each_sample_heatmap.pdf", width = 6, height = 8, dpi = 300,
       plot = heatmap_plot_top20)
ggsave("Diffexp/top20_marker_of_each_sample_heatmap.png", width = 6, height = 8, dpi = 600,
       plot = heatmap_plot_top20)

genes_symbol <- as.character(Diff_exp$gene)
eg = bitr(genes_symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
id = as.character(eg[,2])
dir.create("enrichment")
ego <- enrichGO(gene = id,
                OrgDb = org.Hs.eg.db,
                ont = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
GO_dot <- dotplot(ego)
GO_bar <- barplot(ego)
res_plot <- CombinePlots(list(GO_dot,GO_bar), nrow=1)
ggsave("enrichment/GO_results.pdf", plot=res_plot, width = 15,height = 5,dpi = 300)
ggsave("enrichment/GO_results.png", plot=res_plot, width = 15,height = 5,dpi = 600)
ego_all <- enrichGO(gene = id,
                    OrgDb = org.Hs.eg.db,
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE)

dotplot(ego_all,split = "ONTOLOGY") + facet_grid(ONTOLOGY~.,scales = "free")
barplot(ego_all,split = "ONTOLOGY")+ facet_grid(ONTOLOGY~.,scales = "free")

ggsave("enrichment/GO_ALL_results.pdf", plot=dotplot, width = 6,height = 15,dpi = 300)
ggsave("enrichment/GO_ALL_results.png", plot=dotplot, width = 6,height = 15,dpi = 600)

ggsave("enrichment/GO_ALL_results.pdf", plot=barplot, width = 6,height = 15,dpi = 300)
ggsave("enrichment/GO_ALL_results.png", plot=barplot, width = 6,height = 15,dpi = 600)

Idents(scRNA_mnn)="celltype"
levels(Idents(scRNA_mnn))
cell.markers <- FindAllMarkers(object = scRNA_mnn, 
                               only.pos = FALSE, 
                               test.use = "wilcox", 
                               slot = "data", 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25
)
table(scRNA_mnn@meta.data$celltype,scRNA_mnn@meta.data$seurat_clusters)

Treg_degs = FindMarkers( scRNA_mnn, 
                         logfc.threshold = 0.25,
                         min.pct = 0.1, 
                         only.pos = FALSE,
                         ident.1 = "Treg", ident.2 = "Tconv") %>% 
  mutate( gene = rownames(.) )

Treg_degs_fil = Treg_degs %>% 
  filter( pct.1 > 0.1 & p_val_adj < 0.05 ) %>% 
  filter( abs( avg_log2FC ) > 0.5 )

library(ggrepel) 
colnames(Treg_degs_fil)
table(abs(Treg_degs_fil$avg_log2FC) > 0.5)
plotdt = Treg_degs_fil %>% 
  mutate(gene = ifelse(abs(avg_log2FC) >= 0.5, gene, NA))
volcanoplot= 
  ggplot(plotdt, aes(x = avg_log2FC, y = -log10(p_val_adj), 
                     size = pct.1, 
                     color = avg_log2FC)) +
  geom_point() + 
  ggtitle(label = "aHE vs bCD", subtitle = F) + 
  geom_text_repel(aes(label = gene), size = 3, color = "black") + 
  theme_bw() + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5), 
    plot.background = element_rect(fill = "transparent",colour = NA) 
  ) +
  scale_color_gradient2(low = "olivedrab", high = "salmon2", 
                        mid = "grey", midpoint = 0) + 
  scale_size(range = c(1,3))
dir.create("subTreg_degs")
ggsave("subTreg_degs/volcanoplot.pdf", width = 7, height = 6,dpi = 300,
       plot = volcanoplot)
ggsave("subTreg_degs/volcanoplot.png", width = 7, height = 6, dpi = 600,
       plot = volcanoplot)

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
Treg_degs_fil$gene <- rownames(Treg_degs_fil)
ids=bitr(Treg_degs_fil$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
Treg_degs_fil=merge(Treg_degs_fil,ids,by.x='gene',by.y='SYMBOL')
head(Treg_degs_fil)
Treg_degs_fil <- Treg_degs_fil[order(Treg_degs_fil$avg_log2FC,decreasing = T),]
Treg_degs_list <- as.numeric(Treg_degs_fil$avg_log2FC)
names(Treg_degs_list) <- Treg_degs_fil$ENTREZID
head(Treg_degs_list)

Treg_de <- names(Treg_degs_list)[abs(Treg_degs_list) > 0.5]
head(Treg_de)
Treg_ego <- enrichGO(Treg_de, OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)
head(Treg_ego)

GO = dotplot(Treg_ego, showCategory=10,title="Treg vs Tconv-GO")

Treg_ekg <- enrichKEGG(gene= Treg_de, organism = "hsa",pvalueCutoff = 0.05)
head(Treg_ekg)
dotplot(Treg_ekg, showCategory=10,title="aHE vs bCD KEGG")
KEGG = barplot(Treg_ekg, showCategory=10,title="Treg vs Tconv-KEGG")

ggsave("subenrichment/GO.pdf", plot=GO, width = 7,height = 5,dpi = 300)
ggsave("subenrichment/GO.png", plot=GO, width = 7,height = 5,dpi = 600)

ggsave("subenrichment/KEGG.pdf", plot=KEGG, width = 7,height = 5,dpi = 300)
ggsave("subenrichment/KEGG.png", plot=KEGG, width = 7,height = 5,dpi = 600)

Treg_degs$gene <- rownames(Treg_degs)
ids=bitr(Treg_degs$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
Treg_degs=merge(Treg_degs,ids,by.x='gene',by.y='SYMBOL')
head(Treg_degs)

Treg_degs <- Treg_degs[order(Treg_degs$avg_log2FC,decreasing = T),]
Treg.markers_list <- as.numeric(Treg_degs$avg_log2FC)
names(Treg.markers_list) <- Treg_degs$ENTREZID
head(Treg.markers_list)

Treg_gsekg <- gseKEGG(Treg.markers_list,organism = "hsa",pvalueCutoff = 0.05)
head(Treg_gsekg)
Treg_gsekg_arrange <- arrange(Treg_gsekg,desc(abs(NES)))
head(Treg_gsekg_arrange)
library(clusterProfiler)

library(enrichplot)
library(ggupset)

color <- c("#f7ca64", "#43a5bf", "#86c697", "#a670d6", "#ef998a")
gsekp1 <- gseaplot2(Treg_gsekg_arrange, 1:4, color = color, pvalue_table=F, base_size=14)
gsekp2 <- upsetplot(Treg_gsekg_arrange, n=5)
cowplot::plot_grid(gsekp1, gsekp2, rel_widths=c(1, .6), labels=c("A", "B"))

dir.create("SingleR")
Idents(scRNA_mnn) <- 'seurat_clusters'
load("HumanPrimaryCellAtlas_hpca.se_human.RData")
norm_count = GetAssayData(scRNA_mnn, 
                          slot="data")
pred<- SingleR(test = norm_count, 
               ref = hpca.se, 
               labels = hpca.se$label.main)
scRNA_mnn$singleR_celltype = pred$labels

celltype_plot = DimPlot(scRNA_mnn, 
                        reduction = "umap", 
                        group.by = "singleR_celltype")
ggsave("SingleR/celltype_plot.pdf", celltype_plot)
ggsave("SingleR/celltype_plot.png", celltype_plot)
celltype_plot

score_heatmap <-  plotScoreHeatmap(pred, clusters = scRNA_mnn@meta.data$seurat_clusters, order.by = "clusters")
ggsave("SingleR/score_heatmap.pdf", score_heatmap, width = 8, height = 10)
ggsave("SingleR/score_heatmap.png", score_heatmap, width = 8, height = 10)
score_heatmap

tab <- table(cluster=scRNA_mnn@meta.data$seurat_clusters, label=pred$labels) 
cluster_type <- pheatmap::pheatmap(log10(tab+10))
ggsave("SingleR/cluster_type.pdf", cluster_type)
ggsave("SingleR/cluster_type.png", cluster_type)
cluster_type

scRNA_mnn <- RenameIdents(scRNA_mnn, '0'='T_cells','1'='B_cell','2'='T_cells','3'='B_cell','4'='T_cells','5'='T_cells','6'='T_cells','7'='Monocyte','8'='B_cell','9'='Monocyte','10'='B_cell','11'='Pre-B_cell_CD34-','12'='CMP')
scRNA_mnn[["celltype"]] <- Idents(scRNA_mnn)
head(scRNA_mnn@meta.data)

saveRDS(scRNA_mnn, "Step5_CellTyping.rds")

p1 <- DimPlot(scRNA_mnn, group.by = "seurat_clusters")
p2 <- DimPlot(scRNA_mnn, group.by = "celltype")
p3 = p1+p2
ggsave("SingleR/celltype_contrast.pdf",p3,width = 12, height = 6)
ggsave("SingleR/celltype_contrast.png",p3, width = 12, height = 6)
p3

DimPlot(scRNA_mnn, group.by = "celltype", split.by = "orig.ident")

sessionInfo()

set.seed(42)
scRNA_mnn=readRDS("Step4_MNN.rds") 
mycolor1<-c(
  "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
  "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
  "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93")

mycolor2<-c(
  "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F","#FF6A6A","#7FFFD4", "#AB82FF","#90EE90",
  "#00CD00","#008B8B","#6495ED","#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00","#000000","#708090")

mycolor3<-c(
  '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', 
  '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', 
  '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', 
  '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')

mycolor4=c(
  "#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
  "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
  "#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
  "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

table(scRNA_mnn@meta.data$orig.ident)
group <- data.table::fread("group.csv", header = TRUE)
head(group)
metadata <- FetchData(scRNA_mnn, "orig.ident")
metadata$cell_id <- rownames(metadata)
metadata <- left_join(x = metadata, y = group, by = "orig.ident")
rownames(metadata) <- metadata$cell_id
scRNA_mnn <- AddMetaData(scRNA_mnn, metadata = metadata)
table(scRNA_mnn@meta.data$group)
saveRDS(scRNA_mnn,"Step5_MNN_group.rds")

scRNA_mnn=readRDS("Step5_MNN_group.rds")

sc = scRNA_mnn

gene = c("CSMD1",
         "TRDC","TRGC1","TRGC2","CD3D","CD3E","CD3G",
         "TRDV1","TRDV2","TRGV9","TRGV2",
         "IL17A","IL17F","IL23R","CCR6","RORC","CCR2","CXCR6",
         "GZMK","GZMB","PRF1","GNLY","NKG7",
         "FOXP1","STAT1","NFATC3",
         "TRAV1-2","TRAV10","TRAJ18","KLRC1",
         "TRAC","TRBC1","TRBC2",
         "IFNG-AS1","SLC4A10","FOS",
         "FGFBP2","MYOM2","TRGC1",
         "CSMD1","A2M","MAP2K5",
         "KLRF1","SGCD","LYST",
         "TNFAIP3","CD74","PIK3R1",
         "LEF1","NELL2","BACH2",
         "GNLY","ITGB1","LGALS1",
         "MT-ND3","RUNX1","CHST11",
         "LEF1", "TCF7", "SELL", "CCR7",
         "FOS","IFNG-AS1",
         "GZMB", "FGFBP2","CX3CR1",
         "CSMD1","RORA-AS1",
         "KLRF1", "NCR1","SH2D1B",
         "TNFAIP3","CD74",
         "ITGB1","GNLY",
         "ZEB2"
)
missing_genes <- c()

for (gene in gene) {
  if (!(gene %in% rownames(sc))) {
    message(paste("Gene", gene, "not found, skipping..."))
    missing_genes <- c(missing_genes, gene)
    next
  }
  
  featureplot = FeaturePlot_scCustom(seurat_object = sc, features = gene)
  
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/fanxiu/FeaturePlot/", gene, "_featureplot.pdf"), 
         plot = featureplot, width = 8, height = 7)
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/fanxiu/FeaturePlot/", gene, "_featureplot.png"), 
         plot = featureplot, width = 8, height = 7)
}

if (length(missing_genes) > 0) {
  cat("\nThe following genes were NOT found and skipped:\n")
  print(missing_genes)
  cat("\n")
} else {
  cat("All genes are present in the Seurat object.\n")
}

gene = c("CSMD1",
         "TRDC","TRGC1","TRGC2","CD3D","CD3E","CD3G",
         "TRDV1","TRDV2","TRGV9","TRGV2",
         "IL17A","IL17F","IL23R","CCR6","RORC","CCR2","CXCR6",
         "GZMK","GZMB","PRF1","GNLY","NKG7",
         "FOXP1","STAT1","NFATC3",
         "TRAV1-2","TRAV10","TRAJ18","KLRC1",
         "TRAC","TRBC1","TRBC2",
         "IFNG-AS1","SLC4A10","FOS",
         "FGFBP2","MYOM2","TRGC1",
         "CSMD1","A2M","MAP2K5",
         "KLRF1","SGCD","LYST",
         "TNFAIP3","CD74","PIK3R1",
         "LEF1","NELL2","BACH2",
         "GNLY","ITGB1","LGALS1",
         "MT-ND3","RUNX1","CHST11",
         "LEF1", "TCF7", "SELL", "CCR7",
         "FOS","IFNG-AS1",
         "GZMB", "FGFBP2","CX3CR1",
         "CSMD1","RORA-AS1",
         "KLRF1", "NCR1","SH2D1B",
         "TNFAIP3","CD74",
         "ITGB1","GNLY",
         "ZEB2"
)
missing_genes <- c()

for (gene in gene) {
  if (!(gene %in% rownames(sc))) {
    message(paste("Gene", gene, "not found, skipping..."))
    missing_genes <- c(missing_genes, gene)
    next
  }
  
  featureplot = FeaturePlot_scCustom(seurat_object = sc, features = gene, split.by = "group")
  
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/fanxiu/FeaturePlot/", gene, "_featureplot_group.pdf"), 
         plot = featureplot, width = 16, height = 7)
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/fanxiu/FeaturePlot/", gene, "_featureplot_group.png"), 
         plot = featureplot, width = 16, height = 7)
}

if (length(missing_genes) > 0) {
  cat("\nThe following genes were NOT found and skipped:\n")
  print(missing_genes)
  cat("\n")
} else {
  cat("All genes are present in the Seurat object.\n")
}

gene = c("CSMD1",
         "TRDC","TRGC1","TRGC2","CD3D","CD3E","CD3G",
         "TRDV1","TRDV2","TRGV9","TRGV2",
         "IL17A","IL17F","IL23R","CCR6","RORC","CCR2","CXCR6",
         "GZMK","GZMB","PRF1","GNLY","NKG7",
         "FOXP1","STAT1","NFATC3",
         "TRAV1-2","TRAV10","TRAJ18","KLRC1",
         "TRAC","TRBC1","TRBC2",
         "IFNG-AS1","SLC4A10","FOS",
         "FGFBP2","MYOM2","TRGC1",
         "CSMD1","A2M","MAP2K5",
         "KLRF1","SGCD","LYST",
         "TNFAIP3","CD74","PIK3R1",
         "LEF1","NELL2","BACH2",
         "GNLY","ITGB1","LGALS1",
         "MT-ND3","RUNX1","CHST11",
         "LEF1", "TCF7", "SELL", "CCR7",
         "FOS","IFNG-AS1",
         "GZMB", "FGFBP2","CX3CR1",
         "CSMD1","RORA-AS1",
         "KLRF1", "NCR1","SH2D1B",
         "TNFAIP3","CD74",
         "ITGB1","GNLY",
         "ZEB2"
)
missing_genes <- c()

for (gene in gene) {
  if (!(gene %in% rownames(sc))) {
    message(paste("Gene", gene, "not found, skipping..."))
    missing_genes <- c(missing_genes, gene)
    next
  }
  
  featureplot = FeaturePlot_scCustom(seurat_object = sc, features = gene,split.by = "orig.ident")
  
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/fanxiu/FeaturePlot/", gene, "_featureplot_sample.pdf"), 
         plot = featureplot, width = 42, height = 7)
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/fanxiu/FeaturePlot/", gene, "_featureplot_sample.png"), 
         plot = featureplot, width = 42, height = 7)
}

if (length(missing_genes) > 0) {
  cat("\nThe following genes were NOT found and skipped:\n")
  print(missing_genes)
  cat("\n")
} else {
  cat("All genes are present in the Seurat object.\n")
}

gene = c("CSMD1",
         "TRDC","TRGC1","TRGC2","CD3D","CD3E","CD3G",
         "TRDV1","TRDV2","TRGV9","TRGV2",
         "IL17A","IL17F","IL23R","CCR6","RORC","CCR2","CXCR6",
         "GZMK","GZMB","PRF1","GNLY","NKG7",
         "FOXP1","STAT1","NFATC3",
         "TRAV1-2","TRAV10","TRAJ18","KLRC1",
         "TRAC","TRBC1","TRBC2",
         "IFNG-AS1","SLC4A10","FOS",
         "FGFBP2","MYOM2","TRGC1",
         "CSMD1","A2M","MAP2K5",
         "KLRF1","SGCD","LYST",
         "TNFAIP3","CD74","PIK3R1",
         "LEF1","NELL2","BACH2",
         "GNLY","ITGB1","LGALS1",
         "MT-ND3","RUNX1","CHST11",
         "LEF1", "TCF7", "SELL", "CCR7",
         "FOS","IFNG-AS1",
         "GZMB", "FGFBP2","CX3CR1",
         "CSMD1","RORA-AS1",
         "KLRF1", "NCR1","SH2D1B",
         "TNFAIP3","CD74",
         "ITGB1","GNLY",
         "ZEB2"
)
missing_genes <- c()

for (gene in gene) {
  if (!(gene %in% rownames(sc))) {
    message(paste("Gene", gene, "not found, skipping..."))
    missing_genes <- c(missing_genes, gene)
    next
  }
  
  featureplot = FeaturePlot_scCustom(seurat_object = sc, features = gene,split.by = "celltype", num_columns = 3)
  
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/fanxiu/FeaturePlot/", gene, "_featureplot_celltype.pdf"), 
         plot = featureplot, width = 24, height = 21)
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/fanxiu/FeaturePlot/", gene, "_featureplot_celltype.png"), 
         plot = featureplot, width = 24, height = 21)
}

if (length(missing_genes) > 0) {
  cat("\nThe following genes were NOT found and skipped:\n")
  print(missing_genes)
  cat("\n")
} else {
  cat("All genes are present in the Seurat object.\n")
}

DimPlot(scRNA_mnn, pt.size=0.1,label=T,group.by = "seurat_clusters")+scale_color_manual(values = mycolor1)
scRNA_mnn@meta.data$celltype = "NA"
unique(scRNA_mnn@meta.data$celltype)
scRNA_mnn@meta.data[which(scRNA_mnn@meta.data$seurat_clusters %in% c(11,12)), "celltype"] = "Terminal-branch"
scRNA_mnn@meta.data[which(scRNA_mnn@meta.data$seurat_clusters %in% c(0,1)), "celltype"] = "Cytotoxic effector"
scRNA_mnn@meta.data[which(scRNA_mnn@meta.data$seurat_clusters %in% c(5)), "celltype"] = "NK-like"
scRNA_mnn@meta.data[which(scRNA_mnn@meta.data$seurat_clusters %in% c(10)), "celltype"] = "Cytotoxic memory-like"
scRNA_mnn@meta.data[which(scRNA_mnn@meta.data$seurat_clusters %in% c(9)), "celltype"] = "Naive-like"
scRNA_mnn@meta.data[which(scRNA_mnn@meta.data$seurat_clusters %in% c(7)), "celltype"] = "Activation-regulated"
scRNA_mnn@meta.data[which(scRNA_mnn@meta.data$seurat_clusters %in% c(2,3,4)), "celltype"] = "Transitional activated"
scRNA_mnn@meta.data[which(scRNA_mnn@meta.data$seurat_clusters %in% c(6,8)), "celltype"] = "CSMD1+"
p5=DimPlot(scRNA_mnn, pt.size=0.1,label=T,label.size = 3.5,group.by = "celltype")+scale_color_manual(values = mycolor4)
p5
p6=DimPlot(scRNA_mnn, pt.size=0.1,label=F,label.size = 3.5,group.by = "celltype",split.by="orig.ident", ncol=3)+scale_color_manual(values = mycolor4)
p6
p7=DimPlot(scRNA_mnn, pt.size=0.1,label=F,label.size = 3.5,group.by = "celltype",split.by="group", ncol=3)+scale_color_manual(values = mycolor4)
p7
ggsave('Clustering/celltype.png', p5, width=9, height=6, dpi=600)
ggsave('Clustering/celltype.pdf', p5, width=9, height=6, dpi=300)
ggsave('Clustering/celltype_sample.png', p6, width=11, height=6, dpi=600)
ggsave('Clustering/celltype_sample.pdf', p6, width=11, height=6, dpi=300)
ggsave('Clustering/celltype_group.png', p7, width=10, height=4, dpi=600)
ggsave('Clustering/celltype_group.pdf', p7, width=10, height=4, dpi=300)

saveRDS(scRNA_mnn,"/thinker/3.tangjiale/shinanxi/rds/rds/Overall_shinanxi_New.rds")

sample_table <- as.data.frame(table(scRNA_mnn@meta.data$orig.ident,scRNA_mnn@meta.data$celltype))
names(sample_table) <- c("orig.ident","celltype","CellNumber")
colour = mycolor4
plot_sample_1 <-ggplot(sample_table,aes(x=orig.ident,weight=CellNumber,fill=celltype))+
  geom_bar(position="fill")+
  scale_fill_manual(values=colour) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")+RotatedAxis()
plot_sample_1

ggsave('Clustering/celltype_rate_bysamples.pdf', plot_sample_1, width=5, height=6, dpi=300)
ggsave('Clustering/celltype_rate_bysamples.png', plot_sample_1, width=5, height=6, dpi=600)
write.table(sample_table, 
            "Clustering/celltype_rate_bysamples.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")

group_table <- as.data.frame(table(scRNA_mnn@meta.data$group,scRNA_mnn@meta.data$celltype))
names(group_table) <- c("group","celltype","CellNumber")
colour = mycolor4
plot_sample_2<-ggplot(group_table,aes(x=group,weight=CellNumber,fill=celltype))+
  geom_bar(position="fill")+
  scale_fill_manual(values=colour) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")+RotatedAxis()
plot_sample_2

ggsave('Clustering/celltype_rate_bygroups.pdf', plot_sample_2, width=4, height=7, dpi=300)
ggsave('Clustering/celltype_rate_bygroups.png', plot_sample_2, width=4, height=7, dpi=600)
write.table(group_table, 
            "Clustering/celltype_rate_bygroups.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")

scRNA = readRDS("/thinker/3.tangjiale/shinanxi/rds/rds/Overall_shinanxi_New.rds")
table(scRNA$seurat_clusters)
scRNA@meta.data$clusters = as.integer(scRNA@meta.data$seurat_clusters)
scRNA@meta.data$clusters = factor(scRNA@meta.data$clusters,levels = as.character(unique(sort(scRNA@meta.data$clusters))))
table(scRNA$clusters)

Idents(scRNA) <- "celltype"
levels(Idents(scRNA))
table(scRNA$celltype)
My_levels <- c("Transitional activated","Cytotoxic effector","CSMD1+","NK-like","Activation-regulated","Naive-like","Cytotoxic memory-like","Terminal-branch")
Idents(scRNA) <- factor(Idents(scRNA), levels= My_levels)
scRNA@meta.data$celltype <- factor(Idents(scRNA), levels= My_levels)
table(Idents(scRNA))
saveRDS(scRNA,"/thinker/3.tangjiale/shinanxi/rds/rds/Overall_shinanxi_New.rds")




#单细胞基础转录组基础流程
library(Seurat)#########SeruatV4和v5共存的安装方案
##稍微进阶一些，我们也不需要卸载V5，也可以做到使用V4，只需要在安装的时候，指定一个新的安装路径
#dir.create("SeuratV4")
# 然后安装的时候，指定安装目录
#install.packages('Seurat', repos = c('https://satijalab.r-universe.dev'), lib = "SeuratV4")
##在使用的时候，.libPaths配置好路径，使得library加载的时候，能够先加载到V4。
.libPaths(c("SeuratV4", .libPaths()))
library(Seurat)
packageVersion("Seurat")##################注意，这里SeuratV4/V5版本会影响FinderMaker等函数的结果
#install.packages("remotes")
#remotes::install_github("satijalab/seurat-data")
library(SeuratData)
library(patchwork)
library(ggplot2)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("batchelor")
library(batchelor)
#remotes::install_github('satijalab/seurat-wrappers@community-vignette')
library(SeuratWrappers)
library(magrittr)
library(tidyverse)
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(GO.db)
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
library(DOSE)
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
#BiocManager::install("SingleR")
library(SingleR)
#BiocManager::install("celldex")
library(celldex)
library(harmony)
library(pheatmap)


###########7.Cell cycle analysis###########
```r
dir.create("/thinker/3.tangjiale/shinanxi/cellcycle")

library(Seurat)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(MatrixGenerics)
library(matrixStats)
library(SeuratObject)
library(SingleCellExperiment)
library(scran)
library(org.Hs.eg.db)
library(AnnotationDbi)

dir.create("/thinker/3.tangjiale/shinanxi/cellcycle/plot")
seurat_obj = scRNA
sce <- as.SingleCellExperiment(seurat_obj)
human_cycle_markers <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
gene_symbols <- rownames(sce)
ensembl_ids <- mapIds(org.Hs.eg.db, keys = gene_symbols, column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
head(ensembl_ids)
rownames(sce) <- ensembl_ids
valid_genes <- !is.na(ensembl_ids)
sce <- sce[valid_genes, ]
cycle_pred <- cyclone(sce, pairs = human_cycle_markers)
saveRDS(cycle_pred,file = "/thinker/3.tangjiale/shinanxi/cellcycle/Overall_cycle_pred.rds")
head(cycle_pred$phases)
table(cycle_pred$phases)
seurat_obj$G1_score <- cycle_pred$scores$G1
seurat_obj$S_score <- cycle_pred$scores$S
seurat_obj$G2M_score <- cycle_pred$scores$G2M
seurat_obj$predicted_phase <- cycle_pred$phases
saveRDS(seurat_obj,file = "/thinker/3.tangjiale/shinanxi/cellcycle/Overall_cycle.rds")

mycolor2<-c('#32CD32', '#0000CD', '#FF3030')

library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(viridis)

plot_data <- data.frame(
  G1 = cycle_pred$scores$G1,
  G2M = cycle_pred$scores$G2M,
  S = cycle_pred$scores$S,
  phase = cycle_pred$phases
)

scatter_plot <- ggplot(plot_data, 
                       aes(x = G1, y = G2M, color = phase)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(
    name = "Cell Cycle Phase",
    values = c("G1" = '#32CD32', "S" =  "#0000CD","G2M" = "#FF3030")) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
  labs(title = "Cell Cycle Phase Prediction Scores",
       x = "G1 Score", y = "G2/M Score") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90"),
    aspect.ratio = 1
  )
dir.create("/thinker/3.tangjiale/shinanxi/cellcycle/plot")
print(scatter_plot)
ggsave('/thinker/3.tangjiale/shinanxi/cellcycle/plot/scatter_plot.png', scatter_plot, width=7.5, height=7.5, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/cellcycle/plot/scatter_plot.pdf', scatter_plot, width=7.5, height=7.5, dpi=300)

phase_percent <- prop.table(table(cycle_pred$phases)) * 100
phase_labels <- paste0(names(phase_percent), " (", round(phase_percent, 1), "%)")

scatter_plot <- ggplot(plot_data, 
                       aes(x = G1, y = G2M, color = phase)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(
    name = "Cell Cycle Phase",
    values = c("G1" = '#32CD32', "S" =  "#0000CD","G2M" = "#FF3030"),
    labels = phase_labels
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
  labs(title = "Cell Cycle Phase Prediction Scores",
       subtitle = paste("Total cells:", sum(!is.na(cycle_pred$phases))),
       x = "G1 Score", y = "G2/M Score") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90"),
    aspect.ratio = 1
  ) +
  annotate("text", x = max(plot_data$G1)*0.9, y = max(plot_data$G2M)*0.95,
           label = paste("G1:", round(phase_percent["G1"], 1), "%\n",
                         "S:", round(phase_percent["S"], 1), "%\n",
                         "G2M:", round(phase_percent["G2M"], 1), "%"),
           color = "black", size = 4, hjust = 1)

print(scatter_plot)
ggsave('/thinker/3.tangjiale/shinanxi/cellcycle/plot/scatter2_plot.png', scatter_plot, width=7.5, height=7.5, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/cellcycle/plot/scatter2_plot.pdf', scatter_plot, width=7.5, height=7.5, dpi=300)

p1 = DimPlot(seurat_obj, reduction = "umap", group.by = "predicted_phase")+scale_color_manual(values = c("G1" = '#32CD32', "S" =  "#0000CD","G2M" = "#FF3030"))
p2 = DimPlot(seurat_obj, reduction = "umap", group.by = "predicted_phase",split.by = "celltype",ncol = 2)+scale_color_manual(values = c("G1" = '#32CD32', "S" =  "#0000CD","G2M" = "#FF3030"))
p3 = DimPlot(seurat_obj, reduction = "umap", group.by = "predicted_phase",split.by = "group",ncol = 2)+scale_color_manual(values = c("G1" = '#32CD32', "S" =  "#0000CD","G2M" = "#FF3030"))
ggsave('/thinker/3.tangjiale/shinanxi/cellcycle/plot/dimplot.png', p1, width=10, height=8, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/cellcycle/plot/dimplot.pdf', p1, width=10, height=8, dpi=300)
ggsave('/thinker/3.tangjiale/shinanxi/cellcycle/plot/dimplot_bycelltype.png', p2, width=10, height=13, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/cellcycle/plot/dimplot_bycelltype.pdf', p2, width=10, height=13, dpi=300)
ggsave('/thinker/3.tangjiale/shinanxi/cellcycle/plot/dimplot_bygroup.png', p3, width=10, height=5, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/cellcycle/plot/dimplot_bygroup.pdf', p3, width=10, height=5, dpi=300)

scRNA = seurat_obj
scRNA$predicted_phase <- factor(scRNA$predicted_phase, levels = c("G1", "S", "G2M"))
p1 = CellStatPlot(scRNA, stat.by = "predicted_phase", group.by = "celltype", label = T,plot_type = "trend",palcolor = c("G1" = '#32CD32', "S" =  "#0000CD","G2M" = "#FF3030"),alpha = 1)
p2 = CellStatPlot(scRNA, stat.by = "predicted_phase", group.by = "celltype", label = F,plot_type = "trend",palcolor = c("G1" = '#32CD32', "S" =  "#0000CD","G2M" = "#FF3030"),alpha = 1)
p3 = CellStatPlot(scRNA, stat.by = "predicted_phase", group.by = "group", label = T,plot_type = "trend",palcolor = c("G1" = '#32CD32', "S" =  "#0000CD","G2M" = "#FF3030"),alpha = 1)
p4 = CellStatPlot(scRNA, stat.by = "predicted_phase", group.by = "group", label = F,plot_type = "trend",palcolor = c("G1" = '#32CD32', "S" =  "#0000CD","G2M" = "#FF3030"),alpha = 1)
p5 = CellStatPlot(scRNA, stat.by = "predicted_phase", group.by = "orig.ident", label = T,plot_type = "trend",palcolor = c("G1" = '#32CD32', "S" =  "#0000CD","G2M" = "#FF3030"),alpha = 1)
p6 = CellStatPlot(scRNA, stat.by = "predicted_phase", group.by = "orig.ident", label = F,plot_type = "trend",palcolor = c("G1" = '#32CD32', "S" =  "#0000CD","G2M" = "#FF3030"),alpha = 1)
ggsave('/thinker/3.tangjiale/shinanxi/cellcycle/plot/dotplot_bycelltype.png', p1, width=13, height=12, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/cellcycle/plot/dotplot_bycelltype.pdf', p1, width=13, height=12, dpi=300)
ggsave('/thinker/3.tangjiale/shinanxi/cellcycle/plot/dotplot_bycelltype_nl.png', p2, width=13, height=12, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/cellcycle/plot/dotplot_bycelltype_nl.pdf', p2, width=13, height=12, dpi=300)

ggsave('/thinker/3.tangjiale/shinanxi/cellcycle/plot/dotplot_bygroup.png', p3, width=6, height=12, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/cellcycle/plot/dotplot_bygroup.pdf', p3, width=6, height=12, dpi=300)
ggsave('/thinker/3.tangjiale/shinanxi/cellcycle/plot/dotplot_bygroup_nl.png', p4, width=6, height=12, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/cellcycle/plot/dotplot_bygroup_nl.pdf', p4, width=6, height=12, dpi=300)

ggsave('/thinker/3.tangjiale/shinanxi/cellcycle/plot/dotplot_bysam.png', p5, width=12, height=12, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/cellcycle/plot/dotplot_bysam.pdf', p5, width=12, height=12, dpi=300)
ggsave('/thinker/3.tangjiale/shinanxi/cellcycle/plot/dotplot_bysam_nl.png', p6, width=12, height=12, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/cellcycle/plot/dotplot_bysam_nl.pdf', p6, width=12, height=12, dpi=300)
```


###########8.Basic Map Plotting#######

```r
MH = scRNA
MH$sub_type <- paste0("C",MH$clusters)

df_data <- MH@meta.data %>%
  distinct(celltype, sub_type,clusters) %>%
  dplyr::arrange(clusters) 

df_data$cellnumber <- prop.table(table(MH$clusters)) * 100

df_data$celltype <- factor(df_data$celltype, levels = c("Transitional activated","Cytotoxic effector","CSMD1+","NK-like","Activation-regulated","Naive-like","Cytotoxic memory-like","Terminal-branch"))

df_data <- df_data[order(df_data$celltype), ]

df_data$num <- seq(1:nrow(df_data))

df_data$cell_name <- paste0(df_data$celltype,"_", df_data$sub_type)

df_data$cell_name <- paste0(df_data$num," ", df_data$cell_name )
table(df_data$celltype)

group_colors <- c("Transitional activated" = "#FF34B3",
                  "Cytotoxic effector"  = "#00F5FF",
                  "CSMD1+" = "#BC8F8F",
                  "NK-like" = "#ADFF2F",
                  "Activation-regulated" = "#FFFF02",
                  "Naive-like" = "#00CD00",
                  "Cytotoxic memory-like" = "#FF6A6A",
                  "Terminal-branch" = "#7FFFD4") 

df_data$color <- group_colors[df_data$celltype]

markers = c("IL7R",
            "FGFBP2",
            "CSMD1",
            "KLRF1",
            "GZMK",
            "LEF1",
            "ITGB1",
            "ZEB2")

table(df_data$celltype)

seg_idx <-  cumsum(table(df_data$celltype))

seq_lwd <- rep(0.35, length(seg_idx))

circos.clear()
group_size <- table(df_data$celltype)

circos.par(start.degree = 270, cell.padding = c(0, 0, 0, 0), 
           gap.after = c(rep(2, length(group_size)-1),15),
           circle.margin = c(0.2, 0.2, 0.2, 0.2))

circos.initialize(factors = df_data$celltype,
                  xlim = cbind(0, group_size))

circos.track(
  ylim = c(0, 1), 
  bg.border = NA, 
  track.height = 0.01,
  
  panel.fun = function(x, y) {
    
    sector_index = get.cell.meta.data("sector.index") 
    group_size = group_size[sector_index] 
    
    for (i in 1:group_size) {
      circos.text(
        x = i - 0.5, 
        y = 0.5, 
        labels = df_data$cell_name[df_data$celltype == sector_index][i], 
        col= df_data$color[df_data$celltype == sector_index][i],
        
        font = 2,
        facing = "reverse.clockwise", 
        niceFacing = TRUE,
        adj = c(1, 0.5), 
        cex = 0.5)
    }
  }
)

circos.track(ylim = c(0, 1),
             bg.border = NA, 
             track.height = 0.06,
             bg.col=group_colors,
             
             
             panel.fun=function(x, y) {
               
               xlim = get.cell.meta.data("xlim") 
               ylim = get.cell.meta.data("ylim")
               
               sector.index = get.cell.meta.data("sector.index")
               
               
             })

MH$sub_type <- factor(MH$sub_type, levels = df_data$sub_type)
MH$orig.ident <- factor(MH$orig.ident, levels = c("NPBMC1","NPBMC2","NPBMC3","PBMC1","PBMC2","PBMC3"))

df_data$NPBMC1 <- prop.table(table(MH$orig.ident,MH$sub_type),2)[1,] 
df_data$NPBMC2 <- prop.table(table(MH$orig.ident,MH$sub_type),2)[2,]
df_data$NPBMC3 <- prop.table(table(MH$orig.ident,MH$sub_type),2)[3,]
df_data$PBMC1 <- prop.table(table(MH$orig.ident,MH$sub_type),2)[4,]
df_data$PBMC2 <- prop.table(table(MH$orig.ident,MH$sub_type),2)[5,]
df_data$PBMC3 <- prop.table(table(MH$orig.ident,MH$sub_type),2)[6,]

sample <- c("NPBMC1","NPBMC2","NPBMC3","PBMC1","PBMC2","PBMC3")
sample_cols <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#98FB98","#1E90FF") %>% setNames(., sample)

circos.track(
  
  ylim = c(0,1),
  bg.border = NA, 
  track.height = 0.06,
  
  panel.fun = function(x, y) {
    
    sector_index = get.cell.meta.data("sector.index")
    group_data = df_data[df_data$celltype == sector_index, ]
    
    for (i in 1:nrow(group_data)) {
      circos.rect(
        xleft = i - 0.9, 
        xright = i - 0.1,
        ybottom = c(0, cumsum(as.vector(group_data[i, sample[1:(length(sample)-1)]]))),
        ytop = cumsum(as.vector(group_data[i, sample])),
        col = sample_cols, 
        border = NA
      )
    }
  }
)

MH$sub_type <- factor(MH$sub_type, levels = df_data$sub_type)
MH$group <- factor(MH$group , levels = c("Control","Psoriasis"))

df_data$Control <- prop.table(table(MH$group,MH$sub_type),2)[1,] 
df_data$Psoriasis <- prop.table(table(MH$group,MH$sub_type),2)[2,]

cyl <- c("Control","Psoriasis")
cyl_cols <- c("#FFE4B5","#800080") %>% setNames(., cyl)

circos.track(
  
  ylim = c(0,1),
  bg.border = NA, 
  track.height = 0.06,
  
  panel.fun = function(x, y) {
    
    sector_index = get.cell.meta.data("sector.index")
    group_data = df_data[df_data$celltype == sector_index, ]
    
    for (i in 1:nrow(group_data)) {
      circos.rect(
        xleft = i - 0.9, 
        xright = i - 0.1,
        ybottom = c(0, cumsum(as.vector(group_data[i, cyl[1:(length(cyl)-1)]]))),
        ytop = cumsum(as.vector(group_data[i, cyl])),
        col = cyl_cols, 
        border = NA
      )
    }
  }
)

DefaultAssay(MH) <- "RNA"

Idents(MH) <- "sub_type"

p  = DotPlot(MH, features = markers, cols = c("lightgrey", "#FF0000"))

dot_data <- p$data
colnames(dot_data)[4] <- 'sub_type'
dot_data$ypos <- rep(seq(1:length(markers)),13)

dot_data <- merge(dot_data, df_data, by = 'sub_type')
dot_data$sub_type <- factor(dot_data$sub_type, levels = unique(p$data$id))
dot_data <- dot_data[order(dot_data$sub_type), ]

p_color <-  ggplot_build(p)$data[[1]]

dot_data$dot_color <- p_color$colour

seg_data <- data.frame("celltype" = c("Transitional activated" ,
                                      "Cytotoxic effector",
                                      "CSMD1+",
                                      "NK-like",
                                      "Activation-regulated",
                                      "Naive-like",
                                      "Cytotoxic memory-like",
                                      "Terminal-branch" ) 
                       ,
                       "marker_genes" = markers,
                       "seg_y" = 1:length(markers))

circos.track(
  ylim = c(0, length(markers)+1),
  bg.border = "black", 
  track.height = 0.25, 
  
  panel.fun = function(x, y) {
    
    sector_index = get.cell.meta.data("sector.index")
    
    seg_sub <- seg_data[seg_data$celltype == sector_index,]
    
    cell_subtypes = unique(dot_data$sub_type[dot_data$celltype == sector_index])
    
    circos.segments(x0 = 0, 
                    y0 = seg_sub$seg_y, 
                    x1 = length(cell_subtypes), 
                    y1 = seg_sub$seg_y, 
                    lty=1, lwd=0.15)
    
    for (i in 1:length(cell_subtypes)) {
      
      subtype_data <- dot_data[dot_data$sub_type == cell_subtypes[i], ]
      
      circos.points(
        x = i - 0.5, 
        y = subtype_data$ypos,
        pch = 16,
        col = subtype_data$dot_color,
        cex = subtype_data$pct.exp / 60
      )}
    
    if(length(seg_sub$seg_y) >1){
      
      y1 <- c(seg_sub$seg_y[min(1:length(seg_sub$seg_y))]-1,seg_sub$seg_y[max(1:length(seg_sub$seg_y))]+1)
      
    }else{
      
      y1 <- seg_sub$seg_y
    }
    
    circos.text(x = length(cell_subtypes)/2,
                y = y1, 
                font = 3, 
                cex = 0.3, 
                labels = seg_sub$marker_genes, 
                col="black",
                facing="bending.inside", 
                niceFacing = TRUE) 
    
  }
)

circos.track(
  
  ylim = c(0,max(df_data$cellnumber)),
  bg.border = NA, 
  track.height = 0.1,
  
  panel.fun = function(x, y) {
    
    sector_index = get.cell.meta.data("sector.index")
    cellnumber_data = df_data[df_data$celltype == sector_index, ]
    
    circos.segments(x0 = nrow(cellnumber_data), 
                    y0 = 0, 
                    x1 = nrow(cellnumber_data), 
                    y1 = 100, 
                    lty=2, lwd=0.35)
    
    for (i in 1:nrow(cellnumber_data)) {
      circos.rect(
        xleft = i - 0.8, 
        xright = i - 0.2,
        ybottom =0,
        ytop = cellnumber_data$cellnumber[i],
        col = "#A4804C", 
        border = NA
      )
    }
  }
)

dir.create("/thinker/3.tangjiale/shinanxi/base")
setwd("/thinker/3.tangjiale/shinanxi/base")

scRNA = MH

p1 = DimPlot(scRNA, group.by = "celltype",pt.size=0.5,label=T,raster = F)+scale_color_manual(values = mycolor2_celltype)
p4 = DimPlot(scRNA, group.by = "celltype",pt.size=0.5,label=F,raster = F)+scale_color_manual(values = mycolor2_celltype)
p2 = DimPlot(scRNA, group.by = "celltype",split.by="orig.ident", pt.size=0.5, label = F,shuffle = T, ncol = 3,raster = F)+scale_color_manual(values = mycolor2_celltype)
p3 = DimPlot(scRNA, group.by = "celltype",split.by="orig.ident", pt.size=0.5, label = F,shuffle = T, ncol = 3,raster = F)+scale_color_manual(values = mycolor2_celltype)+ NoLegend()
p5 = DimPlot(scRNA, group.by = "celltype",split.by="group", pt.size=0.5, label = F,shuffle = T, ncol = 2,raster = F)+scale_color_manual(values = mycolor2_celltype)
p6 = DimPlot(scRNA, group.by = "celltype",split.by="group", pt.size=0.5, label = F,shuffle = T, ncol = 2,raster = F)+scale_color_manual(values = mycolor2_celltype)+ NoLegend()

ggsave('./UMAPcelltype.png', p1, width=10, height=8, dpi=600)
ggsave('./UMAPcelltype.pdf', p1, width=10, height=8, dpi=300)
ggsave('./UMAPcelltype_nl.png', p4, width=10, height=8, dpi=600)
ggsave('./UMAPcelltype_nl.pdf', p4, width=10, height=8, dpi=300)
ggsave('./UMAPcelltype_bysample.png', p2, width=18, height=8, dpi=600)
ggsave('./UMAPcelltype_bysample.pdf', p2, width=18, height=8, dpi=300)
ggsave('./UMAPcelltype_nlbysample.png', p3, width=18, height=8, dpi=600)
ggsave('./UMAPcelltype_nlbysample.pdf', p3, width=18, height=8, dpi=300)
ggsave('./UMAPcelltype_bygroup.png', p5, width=18, height=8, dpi=600)
ggsave('./UMAPcelltype_bygroup.pdf', p5, width=18, height=8, dpi=300)
ggsave('./UMAPcelltype_nlbygroup.png', p6, width=18, height=8, dpi=600)
ggsave('./UMAPcelltype_nlbygroup.pdf', p6, width=18, height=8, dpi=300)

p1 = DimPlot(scRNA, group.by = "clusters",pt.size=0.5,label=F,raster = F)+scale_color_manual(values = mycolor1_clusters)
ggsave('./clusters.png', p1, width=10, height=8, dpi=600)
ggsave('./clusters.pdf', p1, width=10, height=8, dpi=300)
p2 = DimPlot(scRNA, group.by = "group",pt.size=0.5,label=F,raster = F)+scale_color_manual(values = mycolor3_group)
ggsave('./group.png', p2, width=10, height=8, dpi=600)
ggsave('./group.pdf', p2, width=10, height=8, dpi=300)

library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggnetwork)
library(dplyr)
library(ggforce)

library(S4Vectors)

cellType<- levels(MH@meta.data$celltype)
Sample <- unique(MH@meta.data$orig.ident)

{
  mycolor2_celltype_new <- list(
    discrete = mycolor2_celltype
  )
  mycolor3_samples_new <- list(
    discrete = mycolor3_samples
  )
  color_assignments <- setNames(
    c(mycolor2_celltype_new$discrete[1:length(cellType)], mycolor3_samples_new$discrete[1:length(Sample)]),
    c(cellType,Sample)
  )
}

{
  data <- MH@meta.data %>%
    group_by(celltype,orig.ident) %>%
    tally() %>%
    ungroup() %>%
    gather_set_data(1:2) %>%
    dplyr::mutate(
      x = factor(x, levels = unique(x)),
      y = factor(y, levels = unique(y))
    )
  
  DataFrame(data)
  
  data_labels <- tibble(
    group = c(
      rep('celltype', length(cellType)),
      rep('orig.ident', length(Sample))
    )
  ) %>%
    mutate(
      hjust = ifelse(group == 'celltype', 1, 0),
      nudge_x = ifelse(group == 'celltype', -0.1, 0.1)
    )
  
  DataFrame(data_labels)
  
}

p1 <- ggplot(data, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = celltype), alpha = 0.75, axis.width = 0.15) +
  geom_parallel_sets_axes(aes(fill = y), color = 'black', axis.width = 0.1) +
  geom_text(
    aes(y = n, split = y), stat = 'parallel_sets_axes', fontface = 'bold',
    hjust = data_labels$hjust, nudge_x = data_labels$nudge_x
  ) +
  scale_x_discrete(labels = c('celltype','sample')) +
  scale_fill_manual(values = color_assignments) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.title = element_blank(),
    axis.text.x = element_text(face = 'bold', colour = 'black', size = 15),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )

p1
ggsave('./sangji_celltypesample.png', p1, width=12, height=12, dpi=600)
ggsave('./sangji_celltypesample.pdf', p1, width=12, height=12, dpi=300)

cellType<- levels(MH@meta.data$clusters)
Sample <- levels(MH@meta.data$celltype)

{
  mycolor2_celltype_new <- list(
    discrete = mycolor1_clusters
  )
  mycolor3_samples_new <- list(
    discrete = mycolor2_celltype
  )
  color_assignments <- setNames(
    c(mycolor2_celltype_new$discrete[1:length(cellType)], mycolor3_samples_new$discrete[1:length(Sample)]),
    c(cellType,Sample)
  )
}

{
  data <- MH@meta.data %>%
    group_by(clusters,celltype) %>%
    tally() %>%
    ungroup() %>%
    gather_set_data(1:2) %>%
    dplyr::mutate(
      x = factor(x, levels = unique(x)),
      y = factor(y, levels = levels(y))
    )
  
  DataFrame(data)
  
  data_labels <- tibble(
    group = c(
      rep('clusters', length(cellType)),
      rep('celltype', length(Sample))
    )
  ) %>%
    mutate(
      hjust = ifelse(group == 'clusters', 1, 0),
      nudge_x = ifelse(group == 'clusters', -0.1, 0.1)
    )
  
  DataFrame(data_labels)
  
}

p1 <- ggplot(data, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = clusters), alpha = 0.75, axis.width = 0.15) +
  geom_parallel_sets_axes(aes(fill = y), color = 'black', axis.width = 0.1) +
  geom_text(
    aes(y = n, split = y), stat = 'parallel_sets_axes', fontface = 'bold',
    hjust = data_labels$hjust, nudge_x = data_labels$nudge_x
  ) +
  scale_x_discrete(labels = c('clusters','celltype')) +
  scale_fill_manual(values = color_assignments) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.title = element_blank(),
    axis.text.x = element_text(face = 'bold', colour = 'black', size = 15),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )

p1
ggsave('./sangji_clusterscelltype.png', p1, width=12, height=12, dpi=600)
ggsave('./sangji_clusterscelltype.pdf', p1, width=12, height=12, dpi=300)

cellType<- unique(MH@meta.data$orig.ident)
Sample <- unique(MH@meta.data$group)

{
  mycolor2_celltype_new <- list(
    discrete = mycolor3_samples
  )
  mycolor3_samples_new <- list(
    discrete = mycolor3_group
  )
  color_assignments <- setNames(
    c(mycolor2_celltype_new$discrete[1:length(cellType)], mycolor3_samples_new$discrete[1:length(Sample)]),
    c(cellType,Sample)
  )
}

{
  data <- MH@meta.data %>%
    group_by(orig.ident,group) %>%
    tally() %>%
    ungroup() %>%
    gather_set_data(1:2) %>%
    dplyr::mutate(
      x = factor(x, levels = unique(x)),
      y = factor(y, levels = unique(y))
    )
  
  DataFrame(data)
  
  data_labels <- tibble(
    group = c(
      rep('orig.ident', length(cellType)),
      rep('group', length(Sample))
    )
  ) %>%
    mutate(
      hjust = ifelse(group == 'orig.ident', 1, 0),
      nudge_x = ifelse(group == 'orig.ident', -0.1, 0.1)
    )
  
  DataFrame(data_labels)
  
}

p1 <- ggplot(data, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = orig.ident), alpha = 0.75, axis.width = 0.15) +
  geom_parallel_sets_axes(aes(fill = y), color = 'black', axis.width = 0.1) +
  geom_text(
    aes(y = n, split = y), stat = 'parallel_sets_axes', fontface = 'bold',
    hjust = data_labels$hjust, nudge_x = data_labels$nudge_x
  ) +
  scale_x_discrete(labels = c('orig.ident','group')) +
  scale_fill_manual(values = color_assignments) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.title = element_blank(),
    axis.text.x = element_text(face = 'bold', colour = 'black', size = 15),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )

p1
ggsave('./sangji_samplegroup.png', p1, width=12, height=12, dpi=600)
ggsave('./sangji_samplegroup.pdf', p1, width=12, height=12, dpi=300)

library(ggplot2)
library(tidyverse)
library(reshape2)
library(Seurat)
library(tidyr)
library(dplyr)

library(SCP)
library(BiocParallel)
library(ggplot2)
library(gtable)
register(MulticoreParam(workers = 8, progressbar = TRUE))
library(Seurat)
library(SeuratData)
scRNA = MH
p1 = CellStatPlot(scRNA, stat.by = "celltype", group.by = "orig.ident", label = T,plot_type = "trend",palcolor = mycolor2_celltype,alpha = 1)
p2 = CellStatPlot(scRNA, stat.by = "celltype", group.by = "orig.ident", label = F,plot_type = "trend",palcolor = mycolor2_celltype,alpha = 1)

ggsave('./celltype_bysample.png', p1, width=10, height=8, dpi=600)
ggsave('./celltype_bysample.pdf', p1, width=10, height=8, dpi=300)
ggsave('./pp_celltype_bysample.png', p2, width=10, height=8, dpi=600)
ggsave('./pp_celltype_bysample.pdf', p2, width=10, height=8, dpi=300)

color = c("Control" = "#FFE4B5",
          "Psoriasis" = "#800080"
)

p1= CellStatPlot(scRNA, stat.by = "group", group.by = "celltype", palcolor = mycolor3_group, stat_type = "count", position = "dodge")
ggsave('./celltypecompare.png', p1, width=10, height=8, dpi=600)
ggsave('./celltypecompare.pdf', p1, width=10, height=8, dpi=300)  

p1 = CellStatPlot(scRNA, stat.by = "celltype", group.by = "orig.ident", label = T,plot_type = "area",palcolor = mycolor2_celltype,alpha = 1)
p2 = CellStatPlot(scRNA, stat.by = "celltype", group.by = "orig.ident", label = F,plot_type = "area",palcolor = mycolor2_celltype,alpha = 1)
ggsave('./celltype1_bysample.png', p1, width=10, height=8, dpi=600)
ggsave('./celltype1_bysample.pdf', p1, width=10, height=8, dpi=300)
ggsave('./pp_celltype1_bysample.png', p2, width=10, height=8, dpi=600)
ggsave('./pp_celltype1_bysample.pdf', p2, width=10, height=8, dpi=300)

p1 = CellStatPlot(scRNA, stat.by = "celltype", group.by = "group", label = T,plot_type = "trend",palcolor = mycolor2_celltype,alpha = 1)
p2 = CellStatPlot(scRNA, stat.by = "celltype", group.by = "group", label = F,plot_type = "trend",palcolor = mycolor2_celltype,alpha = 1)

ggsave('./celltype_bygroup.png', p1, width=4, height=8, dpi=600)
ggsave('./celltype_bygroup.pdf', p1, width=4, height=8, dpi=300)
ggsave('./pp_celltype_bygroup.png', p2, width=4, height=8, dpi=600)
ggsave('./pp_celltype_bygroup.pdf', p2, width=4, height=8, dpi=300)

p1 = CellStatPlot(scRNA, stat.by = "celltype", group.by = "group", label = T,plot_type = "area",palcolor = mycolor2_celltype,alpha = 1)
p2 = CellStatPlot(scRNA, stat.by = "celltype", group.by = "group", label = F,plot_type = "area",palcolor = mycolor2_celltype,alpha = 1)
ggsave('./celltype1_bygroup.png', p1, width=4, height=8, dpi=600)
ggsave('./celltype1_bygroup.pdf', p1, width=4, height=8, dpi=300)
ggsave('./pp_celltype1_bygroup.png', p2, width=4, height=8, dpi=600)
ggsave('./pp_celltype1_bygroup.pdf', p2, width=4, height=8, dpi=300)

p1 = CellStatPlot(scRNA, stat.by = "clusters", group.by = "orig.ident", label = T,plot_type = "trend",palcolor = mycolor1_clusters,alpha = 1)
p2 = CellStatPlot(scRNA, stat.by = "clusters", group.by = "orig.ident", label = F,plot_type = "trend",palcolor = mycolor1_clusters,alpha = 1)

ggsave('./clusters_bysample.png', p1, width=10, height=8, dpi=600)
ggsave('./clusters_bysample.pdf', p1, width=10, height=8, dpi=300)
ggsave('./pp_clusters_bysample.png', p2, width=10, height=8, dpi=600)
ggsave('./pp_clusters_bysample.pdf', p2, width=10, height=8, dpi=300)

p1 = CellStatPlot(scRNA, stat.by = "clusters", group.by = "orig.ident", label = T,plot_type = "area",palcolor = mycolor1_clusters,alpha = 1)
p2 = CellStatPlot(scRNA, stat.by = "clusters", group.by = "orig.ident", label = F,plot_type = "area",palcolor = mycolor1_clusters,alpha = 1)
ggsave('./clusters1_bysample.png', p1, width=10, height=8, dpi=600)
ggsave('./clusters1_bysample.pdf', p1, width=10, height=8, dpi=300)
ggsave('./pp_clusters1_bysample.png', p2, width=10, height=8, dpi=600)
ggsave('./pp_clusters1_bysample.pdf', p2, width=10, height=8, dpi=300)

p1 = CellStatPlot(scRNA, stat.by = "clusters", group.by = "group", label = T,plot_type = "trend",palcolor = mycolor1_clusters,alpha = 1)
p2 = CellStatPlot(scRNA, stat.by = "clusters", group.by = "group", label = F,plot_type = "trend",palcolor = mycolor1_clusters,alpha = 1)

ggsave('./clusters_bygroup.png', p1, width=4, height=8, dpi=600)
ggsave('./clusters_bygroup.pdf', p1, width=4, height=8, dpi=300)
ggsave('./pp_clusters_bygroup.png', p2, width=4, height=8, dpi=600)
ggsave('./pp_clusters_bygroup.pdf', p2, width=4, height=8, dpi=300)

p1 = CellStatPlot(scRNA, stat.by = "clusters", group.by = "group", label = T,plot_type = "area",palcolor = mycolor1_clusters,alpha = 1)
p2 = CellStatPlot(scRNA, stat.by = "clusters", group.by = "group", label = F,plot_type = "area",palcolor = mycolor1_clusters,alpha = 1)
ggsave('./clusters1_bygroup.png', p1, width=4, height=8, dpi=600)
ggsave('./clusters1_bygroup.pdf', p1, width=4, height=8, dpi=300)
ggsave('./pp_clusters1_bygroup.png', p2, width=4, height=8, dpi=600)
ggsave('./pp_clusters1_bygroup.pdf', p2, width=4, height=8, dpi=300)

mycolor1_clusters<-c(
  "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
  "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
  "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93")
mycolor2_celltype <-c(
  "#FF34B3","#00F5FF","#BC8F8F","#ADFF2F","#FFFF02",
  "#00CD00","#FF6A6A","#7FFFD4", "#AB82FF")

Idents(data) <- 'celltype'
levels(Idents(data))

table(seurat_obj@meta.data$celltype)

scRNA1 <- subset(seurat_obj, idents = 'Terminal-branch')
DimPlot(scRNA1)

mycolor3_samples=c(
  "#DC143C","#0000FF","#20B2AA","#FFA500","#98FB98","#1E90FF")
mycolor3_group=c(
  "#FFE4B5","#800080")
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)
.libPaths(c("/thinker/3.tangjiale/SeuratV4/", .libPaths()))
.libPaths(c("/thinker/3.tangjiale/02_package/", .libPaths()))

library(Seurat)

library(SeuratData)
library(patchwork)
library(ggplot2)

library(batchelor)

library(SeuratWrappers)
library(magrittr)
library(tidyverse)

library(clusterProfiler)
library(GO.db)

library(org.Hs.eg.db)
library(DOSE)

library(DoubletFinder)

library(SingleR)

library(celldex)
library(harmony)
library(pheatmap)

scRNA=readRDS("/thinker/3.tangjiale/shinanxi/rds/Overall_shinanxi_New.rds") 
Idents(scRNA) <- 'celltype'
levels(Idents(scRNA))

dir.create("/thinker/3.tangjiale/shinanxi/Marker/")
dir.create("/thinker/3.tangjiale/shinanxi/Marker/plot")

all.markers = FindAllMarkers(scRNA, 
                             min.pct = 0.25, 
                             logfc.threshold = 0.25, 
                             only.pos = TRUE)
head(all.markers)

top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 
top5 = all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top15 = all.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)

write.table(all.markers, 
            "/thinker/3.tangjiale/shinanxi/Marker/all_Markers_of_each_celltype.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")
write.table(top10, 
            "/thinker/3.tangjiale/shinanxi/Marker/top10_Markers_of_each_celltype.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")
write.table(top5, 
            "/thinker/3.tangjiale/shinanxi/Marker/top5_Markers_of_each_celltype.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")
write.table(top15, 
            "/thinker/3.tangjiale/shinanxi/Marker/top15_Markers_of_each_celltype.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")

library(SCP)
mycolor2_celltype <-c(
  "#FF34B3","#00F5FF","#BC8F8F","#ADFF2F","#FFFF02",
  "#00CD00","#FF6A6A","#7FFFD4", "#AB82FF")
p = FeatureStatPlot(scRNA, stat.by = top5$gene, group.by = "celltype", bg.by = "celltype", add_box = TRUE, stack = TRUE,
                    palcolor = mycolor2_celltype,bg_palcolor = "white",flip = T)

gene_cell_exp <- AverageExpression(scRNA,
                                   features = top10$gene,
                                   group.by = 'celltype',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)

library(ComplexHeatmap)

df <- data.frame(colnames(gene_cell_exp))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c(  "Transitional activated"="#FF34B3",
                                                    "Cytotoxic effector"="#00F5FF",
                                                    "CSMD1+"="#BC8F8F",
                                                    "NK-like"="#ADFF2F",
                                                    "Activation-regulated"="#FFFF02",
                                                    "Naive-like"="#00CD00",
                                                    "Cytotoxic memory-like"= "#FF6A6A",
                                                    "Terminal-branch" = "#7FFFD4")))

marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
p = Heatmap(marker_exp,
            cluster_rows = F,
            cluster_columns = F,
            show_column_names = F,
            show_row_names = T,
            column_title = NULL,
            heatmap_legend_param = list(
              title=' '),
            col = colorRampPalette(c("#0000EF","black","#FDFE00"))(100),
            border = 'black',
            rect_gp = gpar(col = "black", lwd = 1),
            row_names_gp = gpar(fontsize = 10),
            column_names_gp = gpar(fontsize = 10),
            top_annotation = top_anno)
pdf("/thinker/3.tangjiale/shinanxi/Marker/plot/Heatmap_top10_260427.pdf", width = 5, height = 16)
print(p)
dev.off()
png("/thinker/3.tangjiale/shinanxi/Marker/plot/Heatmap_top10_260427.png", width = 5, height = 16, res = 600, units = "in")
print(p)
dev.off()

marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
p = Heatmap(marker_exp,
            cluster_rows = F,
            cluster_columns = F,
            show_column_names = F,
            show_row_names = T,
            column_title = NULL,
            heatmap_legend_param = list(
              title=' '),
            col = colorRampPalette(c("#426daa","white","#ae443e"))(100),
            border = 'black',
            rect_gp = gpar(col = "black", lwd = 1),
            row_names_gp = gpar(fontsize = 10),
            column_names_gp = gpar(fontsize = 10),
            top_annotation = top_anno)
pdf("/thinker/3.tangjiale/shinanxi/Marker/plot/Heatmap1.pdf", width = 5, height = 8)
print(p)
dev.off()
png("/thinker/3.tangjiale/shinanxi/Marker/plot/Heatmap1.png", width = 5, height = 8, res = 600, units = "in")
print(p)
dev.off()

KS_plot_density <- function(obj,
                            marker,
                            dim=c("TSNE","UMAP"),
                            size,
                            ncol=NULL
){
  require(ggplot2)
  require(ggrastr)
  require(Seurat)
  
  cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494'))
  warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c'))
  mypalette <- c(rev(cold(11)), warm(10))
  
  if(dim=="TSNE"){
    
    xtitle = "tSNE1"
    ytitle = "tSNE2"
    
  }
  
  if(dim=="UMAP"){
    
    xtitle = "UMAP1"
    ytitle = "UMAP2"
  }
  
  
  if(length(marker)==1){
    
    plot <- FeaturePlot(obj, features = marker)
    data <- plot$data
    
    
    if(dim=="TSNE"){
      
      colnames(data)<- c("x","y","ident","gene")
      
    }
    
    if(dim=="UMAP"){
      
      colnames(data)<- c("x","y","ident","gene")
    }
    
    
    p <- ggplot(data, aes(x, y)) +
      geom_point_rast(shape = 21, stroke=0.25,
                      aes(colour=gene, 
                          fill=gene), size = size) +
      geom_density_2d(data=data[data$gene>0,], 
                      aes(x=x, y=y), 
                      bins = 5, colour="black") +
      scale_fill_gradientn(colours = mypalette)+
      scale_colour_gradientn(colours = mypalette)+
      theme_bw()+ggtitle(marker)+
      labs(x=xtitle, y=ytitle)+
      theme(
        plot.title = element_text(size=12, face="bold.italic", hjust = 0),
        axis.text=element_text(size=8, colour = "black"),
        axis.title=element_text(size=12),
        legend.text = element_text(size =10),
        legend.title=element_blank(),
        aspect.ratio=1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    return(p)
    
  }else{
    
    gene_list <- list()
    
    
    
    for (i in 1:length(marker)) {
      plot <- FeaturePlot(obj, features = marker[i])
      data <- plot$data
      
      
      if(dim=="TSNE"){
        
        colnames(data)<- c("x","y","ident","gene")
      }
      
      if(dim=="UMAP"){
        
        colnames(data)<- c("x","y","ident","gene")
      }
      
      gene_list[[i]] <- data
      names(gene_list) <- marker[i]
    }
    
    plot_list <- list()
    
    
    for (i in 1:length(marker)) {
      
      p <- ggplot(gene_list[[i]], aes(x, y)) +
        geom_point_rast(shape = 21, stroke=0.25,
                        aes(colour=gene, 
                            fill=gene), size = size) +
        geom_density_2d(data=gene_list[[i]][gene_list[[i]]$gene>0,], 
                        aes(x=x, y=y), 
                        bins = 5, colour="black") +
        scale_fill_gradientn(colours = mypalette)+
        scale_colour_gradientn(colours = mypalette)+
        theme_bw()+ggtitle(marker[i])+
        labs(x=xtitle, y=ytitle)+
        theme(
          plot.title = element_text(size=12, face="bold.italic", hjust = 0),
          axis.text=element_text(size=8, colour = "black"),
          axis.title=element_text(size=12),
          legend.text = element_text(size =10),
          legend.title=element_blank(),
          aspect.ratio=1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      
      plot_list[[i]] <- p
    }
    
    
    Seurat::CombinePlots(plot_list, ncol = ncol)
    
    
  }
  
  
}

top5
dir.create("/thinker/3.tangjiale/shinanxi/Marker/plot/featureplot/")
for (gene in top5$gene ) {
  featureplot = KS_plot_density(obj=scRNA, 
                                marker= gene,
                                dim = "UMAP", size =1)
  
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/Marker/plot/featureplot/",gene,"_featureplot.pdf"), plot = featureplot,width = 8, height = 7)
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/Marker/plot/featureplot/",gene,"_featureplot.png"), plot = featureplot,width = 8, height = 7)
}

dir.create("/thinker/3.tangjiale/shinanxi/Marker/plot/vlnplot")
for (gene in top5$gene) {
  vlnplot = VlnPlot(scRNA, 
                    features = gene, 
                    pt.size = 0, 
                    group.by = "orig.ident",
                    cols = mycolor3_samples)
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/Marker/plot/vlnplot/",gene,"_samplevlnplot.pdf"), plot = vlnplot,width = 6, height = 5)
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/Marker/plot/vlnplot/",gene,"_samplevlnplot.png"), plot = vlnplot,width = 6, height = 5)
}
for (gene in top5$gene) {
  vlnplot = VlnPlot(scRNA, 
                    features = gene, 
                    pt.size = 0, 
                    group.by = "celltype",
                    cols = mycolor2_celltype)
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/Marker/plot/vlnplot/",gene,"_celltypevlnplot.pdf"), plot = vlnplot,width = 8, height = 5)
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/Marker/plot/vlnplot/",gene,"_celltypevlnplot.png"), plot = vlnplot,width = 8, height = 5)
}

for (gene in top5$gene) {
  vlnplot = VlnPlot(scRNA, 
                    features = gene, 
                    pt.size = 0, 
                    group.by = "group",
                    cols = mycolor2_celltype)
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/Marker/plot/vlnplot/",gene,"_groupvlnplot.pdf"), plot = vlnplot,width = 4, height = 5)
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/Marker/plot/vlnplot/",gene,"_groupvlnplot.png"), plot = vlnplot,width = 4, height = 5)
}

library(dplyr)
library(dendextend)
top5

gene_cell_exp <- AverageExpression(scRNA,
                                   features = top5$gene,
                                   group.by = 'celltype',
                                   slot = 'data') 
gene_cell_exp
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
rownames(gene_cell_exp)
top5$key

key_genes <- sapply(strsplit(top5$key, "_"), function(x) tail(x, 1))
exp_genes <- rownames(gene_cell_exp)
common_genes <- intersect(key_genes, exp_genes)
top5_filtered <- top5[key_genes %in% common_genes, ]
gene_cell_exp_filtered <- gene_cell_exp[exp_genes %in% common_genes, ]

key_genes_filtered <- sapply(strsplit(top5_filtered$key, "_"), function(x) tail(x, 1))
mapping_vec <- top5_filtered$key
names(mapping_vec) <- key_genes_filtered
mapping_vec <- mapping_vec[!duplicated(names(mapping_vec))]  

rownames(gene_cell_exp_filtered) <- mapping_vec[rownames(gene_cell_exp_filtered)]

a = sapply(strsplit(rownames(gene_cell_exp_filtered), "_"), function(x) tail(x, 1))
b = rownames(gene_cell_exp)

common_genes1 <- setdiff(a, b)
common_genes2 <- setdiff(b, a)
common_genes1
common_genes2

rownames(gene_cell_exp_filtered) <- gsub("_", "-", rownames(gene_cell_exp_filtered), fixed = TRUE)
rownames(gene_cell_exp_filtered)

hc <- as.dist((1- cor(t(gene_cell_exp_filtered))) / 2) %>% 
  hclust(method = "ward.D") %>%
  as.dendrogram() %>% 
  set("labels_cex", 0.5)

circlize_dendrogram(hc,
                    labels_track_height = 0.6,
                    dend_track_height = 0.2)

prefixes <- gsub("-.*", "", rownames(gene_cell_exp_filtered))

cell_color <- c(
  "Transitional activated"="#FF34B3",
  "Cytotoxic effector"="#00F5FF",
  "CSMD1+"="#BC8F8F",
  "NK-like"="#ADFF2F",
  "Activation-regulated"="#FFFF02",
  "Naive-like"="#00CD00",
  "Cytotoxic memory-like"= "#FF6A6A",
  "Terminal-branch" = "#7FFFD4")
cell_colors <- cell_color[prefixes]
cell_colors[is.na(cell_colors)] <-"gray"

cell_colors <-as.data.frame(cell_colors)
cell_colors$celltypes <- rownames(gene_cell_exp_filtered)
rownames(cell_colors) <- rownames(gene_cell_exp_filtered)
cell_colors <- cell_colors[rownames(gene_cell_exp_filtered)[order.dendrogram(hc)],]
labels_colors(hc) <- cell_colors$cell_colors

pdf("/thinker/3.tangjiale/shinanxi/Marker/plot/celltypetop5_tree_circle.pdf", width = 6, height = 6)
circlize_dendrogram(hc,
                    labels_track_height = 0.6,
                    dend_track_height = 0.2)
dev.off()
png("/thinker/3.tangjiale/shinanxi/Marker/plot/celltypetop5_tree_circle.png", width = 6, height = 6, res = 600, units = "in")
circlize_dendrogram(hc,
                    labels_track_height = 0.6,
                    dend_track_height = 0.2)
dev.off()

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(scCustomize)
library(ggrepel)

Idents(scRNA) <- 'group'
levels(Idents(scRNA))

Diff_exp = FindAllMarkers(scRNA, 
                          logfc.threshold = 0, 
                          only.pos = F)

write.csv(Diff_exp, "/thinker/3.tangjiale/shinanxi/Diff_exp_V4/Diff_exp.csv")
write.table(Diff_exp, "/thinker/3.tangjiale/shinanxi/Diff_exp_V4/Diff_exp.xlsx")
class(Diff_exp)
head(Diff_exp)

Diff_exp <- Diff_exp %>%
  rename(group = cluster)
head(Diff_exp)

cell_DEs <- mutate(Diff_exp, sig = ifelse(p_val_adj < 0.05 & avg_log2FC > 0.2 | p_val_adj < 0.05 & avg_log2FC < -0.2, "Sig", "NS"))

cols <- data.frame(group=unique(scRNA@meta.data$group),
                   color = mycolor3_group[1:length(unique(scRNA@meta.data$group))])

cell_DEs <- full_join(cell_DEs, cols, by = "group")

cell_DEs <-group_by(cell_DEs, group)

sigs<-dplyr::filter(cell_DEs, sig == "Sig")
top_5<-sigs %>% group_by(group) %>% slice_max(order_by = avg_log2FC, n = 5)
bot_5<-sigs %>% group_by(group) %>% slice_min(order_by = avg_log2FC, n = 5)

write.csv(top_5,"/thinker/3.tangjiale/shinanxi/Diff_exp_V4/top_5.csv")
write.table(top_5,"/thinker/3.tangjiale/shinanxi/Diff_exp_V4/top_5.xlsx")
write.csv(bot_5,"/thinker/3.tangjiale/shinanxi/Diff_exp_V4/bot_5.csv")
write.table(bot_5,"/thinker/3.tangjiale/shinanxi/Diff_exp_V4/bot_5.xlsx")
five<-dplyr::full_join(top_5, bot_5)

cell_DEs$group<-factor(cell_DEs$group, levels = c("Control" ,"Psoriasis"))
cell_DEs_sig<-dplyr::filter(cell_DEs, sig!= "NS")
col_for_plot<-cell_DEs_sig$color

p1 = ggplot(cell_DEs, aes(x = group, y = avg_log2FC)) +
  geom_jitter_rast(data=cell_DEs[cell_DEs$sig == "NS",],
                   color="grey",
                   width = 0.2, 
                   height = 0.0, 
                   alpha = .25, 
                   shape = 16) +
  geom_jitter_rast(data=cell_DEs[cell_DEs$sig != "NS",],
                   color=col_for_plot,
                   width = 0.2, 
                   height = 0.0, 
                   shape = 16) + 
  theme_bw() +
  theme(axis.text.x = element_text(color = 'black',size = 10), 
        legend.position = "none", 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_text(color = 'black', size = 8))+ 
  geom_text_repel(data = five,label = five$gene, 
                  fontface = "italic", size = 2,
                  max.overlaps = Inf) 
dir.create('/thinker/3.tangjiale/shinanxi/Diff_exp_V4/plot')
ggsave('/thinker/3.tangjiale/shinanxi/Diff_exp_V4/plot/top5bot5.png', p1, width=6, height=8, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/Diff_exp_V4/plot/top5bot5.pdf', p1, width=6, height=8, dpi=300)

dir.create('/thinker/3.tangjiale/shinanxi/Diff_exp_V4/plot/top5_featureplot/')
dir.create('/thinker/3.tangjiale/shinanxi/Diff_exp_V4/plot/top5_vlnplot/')
dir.create('/thinker/3.tangjiale/shinanxi/Diff_exp_V4/plot/bot5_featureplot/')
dir.create('/thinker/3.tangjiale/shinanxi/Diff_exp_V4/plot/bot5_vlnplot/')

for (gene in as.vector(top_5$gene)) {
  p = FeaturePlot_scCustom(scRNA, 
                           features = gene, 
                           label = F,
                           pt.size = 0.3,
                           reduction = "umap",num_columns = 2,split.by = "group")
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/Diff_exp_V4/plot/top5_featureplot/",gene,"_groupfeatureplot.pdf"), plot = p,width=15, height=7, dpi=300)
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/Diff_exp_V4/plot/top5_featureplot/",gene,"_groupfeatureplot.png"), plot = p,width=15, height=7, dpi=600)
}

for (gene in as.vector(c(top_5$gene))) {
  vlnplot = VlnPlot(scRNA, 
                    features = gene, 
                    pt.size = 0, 
                    group.by = "celltype",
                    cols = mycolor2_celltype)
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/Diff_exp_V4/plot/top5_vlnplot/",gene,"_celltypevlnplot.pdf"), plot = vlnplot,width=8, height=5, dpi=300)
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/Diff_exp_V4/plot/top5_vlnplot/",gene,"_celltypevlnplot.png"), plot = vlnplot,width=8, height=5, dpi=600)
}

for (gene in as.vector(c(top_5$gene))) {
  vlnplot = VlnPlot(scRNA, 
                    features = gene, 
                    pt.size = 0, 
                    group.by = "group",
                    cols = mycolor3_group)
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/Diff_exp_V4/plot/top5_vlnplot/",gene,"_groupvlnplot.pdf"), plot = vlnplot,width=8, height=5, dpi=300)
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/Diff_exp_V4/plot/top5_vlnplot/",gene,"_groupvlnplot.png"), plot = vlnplot,width=8, height=5, dpi=600)
}

for (gene in as.vector(bot_5$gene)) {
  p = FeaturePlot_scCustom(scRNA, 
                           features = gene, 
                           label = F,
                           pt.size = 0.3,
                           reduction = "umap",num_columns = 2,split.by = "group")
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/Diff_exp_V4/plot/bot5_featureplot/",gene,"_groupfeatureplot.pdf"), plot = p,width=15, height=7, dpi=300)
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/Diff_exp_V4/plot/bot5_featureplot/",gene,"_groupfeatureplot.png"), plot = p,width=15, height=7, dpi=600)
}

for (gene in as.vector(c(bot_5$gene))) {
  vlnplot = VlnPlot(scRNA, 
                    features = gene, 
                    pt.size = 0, 
                    group.by = "celltype",
                    cols = mycolor2_celltype)
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/Diff_exp_V4/plot/bot5_vlnplot/",gene,"_celltypevlnplot.pdf"), plot = vlnplot,width=8, height=5, dpi=300)
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/Diff_exp_V4/plot/bot5_vlnplot/",gene,"_celltypevlnplot.png"), plot = vlnplot,width=8, height=5, dpi=600)
}

for (gene in as.vector(c(bot_5$gene))) {
  vlnplot = VlnPlot(scRNA, 
                    features = gene, 
                    pt.size = 0, 
                    group.by = "group",
                    cols = mycolor3_group)
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/Diff_exp_V4/plot/bot5_vlnplot/",gene,"_groupvlnplot.pdf"), plot = vlnplot,width=8, height=5, dpi=300)
  ggsave(paste0("/thinker/3.tangjiale/shinanxi/Diff_exp_V4/plot/bot5_vlnplot/",gene,"_groupvlnplot.png"), plot = vlnplot,width=8, height=5, dpi=600)
}
```















###############9.Pathway Enrichment Analysis####################
```r
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
sce.all.int = scRNA
Idents(sce.all.int) <- 'celltype'
levels(Idents(sce.all.int))
sce.markers<-FindAllMarkers(sce.all.int,
                            only.pos=TRUE,
                            min.pct=0.2,
                            return.thresh=0.01
)
head(sce.markers)
dir.create("/thinker/3.tangjiale/shinanxi/rich_V4")
write.csv(sce.markers, "/thinker/3.tangjiale/shinanxi/rich_V4/markers.csv")
write.table(sce.markers, "/thinker/3.tangjiale/shinanxi/rich_V4/markers.xlsx")
library(dplyr)
top10<-sce.markers%>%
  group_by(cluster)%>%
  dplyr::filter(avg_log2FC>1)%>%
  slice_head(n=10)%>%
  ungroup()
write.csv(top10, "/thinker/3.tangjiale/shinanxi/rich_V4/top10.csv")
write.table(top10, "/thinker/3.tangjiale/shinanxi/rich_V4/top10.xlsx")
Symbol<-mapIds(get("org.Hs.eg.db"),keys=sce.markers$gene,keytype="SYMBOL",column="ENTREZID")
head(Symbol)
ids<-bitr(sce.markers$gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
head(ids)

data<-merge(sce.markers,ids,by.x="gene",by.y="SYMBOL")
head(data)
gcSample<-split(data$ENTREZID,data$cluster)
gcSample

xx_kegg<-compareCluster(gcSample,fun="enrichKEGG",
                        organism="hsa",
                        pvalueCutoff=1,
                        qvalueCutoff=1)

res<-xx_kegg@compareClusterResult
head(res)

enrich<-res%>%
  group_by(Cluster)%>%
  top_n(n=5,wt=-pvalue)

dt<-enrich
dt<-dt[order(dt$pvalue,dt$Cluster,decreasing=F),]

dt$Description<-factor(dt$Description,levels=unique(dt$Description))
colnames(dt)

library(ggforce)
head(dt)
table(dt$Cluster)


colors1 <- c(
  "Transitional activated"="#FF34B3",
  "Cytotoxic effector"="#00F5FF",
  "CSMD1+"="#BC8F8F",
  "NK-like"="#ADFF2F",
  "Activation-regulated"="#FFFF02",
  "Naive-like"="#00CD00",
  "Cytotoxic memory-like"= "#FF6A6A",
  "Terminal-branch" = "#7FFFD4")

p <- ggplot(dt) +
  geom_link(aes(x = 0, y = Description,
                xend = -log10(pvalue), yend = Description,
                alpha = after_stat(index),
                color = Cluster,
                size = after_stat(index)), 
            n = 500, show.legend = T) 
p1 <- p +
  geom_point(aes(x = -log10(pvalue),y = Description), color = "black", fill = "white",size = 6,shape = 21) +
  geom_text(aes(x = -log10(pvalue), y = Description), label=dt$Count, size=3, nudge_x=0.05) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face ="bold.italic"),
        axis.line = element_line(color = "black", size = 0.6),
        axis.text = element_text(face = "bold"),
        axis.title = element_text(size = 13)
  ) +
  xlab("-Log10 Pvalue") + ylab("") + 
  scale_color_manual(values = colors1)

p2 <- p1 + 
  facet_wrap(~Cluster,scales = "free",ncol = 2)

dir.create("/thinker/3.tangjiale/shinanxi/rich_V4/plot")
ggsave('/thinker/3.tangjiale/shinanxi/rich_V4/plot/KEGG_celltype.png', p1, width=15, height=10, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/rich_V4/plot/KEGG_celltype.pdf', p1, width=15, height=10, dpi=300)
ggsave('/thinker/3.tangjiale/shinanxi/rich_V4/plot/KEGG_facet_celltype.png', p2, width=15, height=12, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/rich_V4/plot/KEGG_facet_celltype.pdf', p2, width=15, height=12, dpi=300)

xx_go<-compareCluster(gcSample,fun="enrichGO",OrgDb="org.Hs.eg.db",ont="BP",pvalueCutoff=1,qvalueCutoff=1)

res<-xx_go@compareClusterResult

enrich<-res%>%
  group_by(Cluster)%>%
  top_n(n=5,wt=-pvalue)

dt<-enrich
dt<-dt[order(dt$pvalue,dt$Cluster,decreasing=F),]

dt$Description<-factor(dt$Description,levels=unique(dt$Description))
colnames(dt)

library(ggforce)
head(dt)
table(dt$Cluster)
colors1 <- c(
  "Transitional activated"="#FF34B3",
  "Cytotoxic effector"="#00F5FF",
  "CSMD1+"="#BC8F8F",
  "NK-like"="#ADFF2F",
  "Activation-regulated"="#FFFF02",
  "Naive-like"="#00CD00",
  "Cytotoxic memory-like"= "#FF6A6A",
  "Terminal-branch" = "#7FFFD4")
p <- ggplot(dt) +
  geom_link(aes(x = 0, y = Description,
                xend = -log10(pvalue), yend = Description,
                alpha = after_stat(index),
                color = Cluster,
                size = after_stat(index)), 
            n = 500, show.legend = T) 
p1 <- p +
  geom_point(aes(x = -log10(pvalue),y = Description), color = "black", fill = "white",size = 6,shape = 21) +
  geom_text(aes(x = -log10(pvalue), y = Description), label=dt$Count, size=3, nudge_x=0.05) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face ="bold.italic"),
        axis.line = element_line(color = "black", size = 0.6),
        axis.text = element_text(face = "bold"),
        axis.title = element_text(size = 13)
  ) +
  xlab("-Log10 Pvalue") + ylab("") + 
  scale_color_manual(values = colors1)

p2 <- p1 + 
  facet_wrap(~Cluster,scales = "free",ncol = 2)
ggsave('/thinker/3.tangjiale/shinanxi/rich_V4/plot/GO_celltype.png', p1, width=15, height=10, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/rich_V4/plot/GO_celltype.pdf', p1, width=15, height=10, dpi=300)
ggsave('/thinker/3.tangjiale/shinanxi/rich_V4/plot/GO_facet_celltype.png', p2, width=20, height=12, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/rich_V4/plot/GO_facet_celltype.pdf', p2, width=20, height=12, dpi=300)

sce.all.int = scRNA
Idents(sce.all.int) <- 'group'
levels(Idents(sce.all.int))
sce.markers<-FindAllMarkers(sce.all.int,
                            only.pos=TRUE,
                            min.pct=0.2,
                            return.thresh=0.01
)
write.csv(sce.markers, "/thinker/3.tangjiale/shinanxi/rich_V4/groupDiff_exp.csv")
write.table(sce.markers, "/thinker/3.tangjiale/shinanxi/rich_V4/groupDiff_exp.xlsx")
head(sce.markers)
library(dplyr)
top10<-sce.markers%>%
  group_by(cluster)%>%
  dplyr::filter(avg_log2FC>1)%>%
  slice_head(n=10)%>%
  ungroup()
write.csv(top10, "/thinker/3.tangjiale/shinanxi/rich_V4/top10Diff_exp.csv")
write.table(top10, "/thinker/3.tangjiale/shinanxi/rich_V4/top10Diff_exp.xlsx")
Symbol<-mapIds(get("org.Hs.eg.db"),keys=sce.markers$gene,keytype="SYMBOL",column="ENTREZID")
head(Symbol)
ids<-bitr(sce.markers$gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
head(ids)

data<-merge(sce.markers,ids,by.x="gene",by.y="SYMBOL")
head(data)
gcSample<-split(data$ENTREZID,data$cluster)
gcSample

xx_kegg<-compareCluster(gcSample,fun="enrichKEGG",
                        organism="hsa",
                        pvalueCutoff=1,
                        qvalueCutoff=1)


res<-xx_kegg@compareClusterResult
head(res)

enrich<-res%>%
  group_by(Cluster)%>%
  top_n(n=10,wt=-pvalue)

dt<-enrich
dt<-dt[order(dt$pvalue,dt$Cluster,decreasing=F),]

dt$Description<-factor(dt$Description,levels=unique(dt$Description))
colnames(dt)

library(ggforce)
head(dt)
table(dt$Cluster)

colors1 <-c(
  'Control' =  "#FFE4B5",  
  'Psoriasis' = "#800080")
p <- ggplot(dt) +
  geom_link(aes(x = 0, y = Description,
                xend = -log10(pvalue), yend = Description,
                alpha = after_stat(index),
                color = Cluster,
                size = after_stat(index)), 
            n = 500, show.legend = T) 
p1 <- p +
  geom_point(aes(x = -log10(pvalue),y = Description), color = "black", fill = "white",size = 6,shape = 21) +
  geom_text(aes(x = -log10(pvalue), y = Description), label=dt$Count, size=3, nudge_x=0.05) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face ="bold.italic"),
        axis.line = element_line(color = "black", size = 0.6),
        axis.text = element_text(face = "bold"),
        axis.title = element_text(size = 13)
  ) +
  xlab("-Log10 Pvalue") + ylab("") + 
  scale_color_manual(values = colors1)

p2 <- p1 + 
  facet_wrap(~Cluster,scales = "free",ncol = 2)
ggsave('/thinker/3.tangjiale/shinanxi/rich_V4/plot/KEGG_ori.ident.png', p1, width=14, height=8, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/rich_V4/plot/KEGG_ori.ident.pdf', p1, width=14, height=8, dpi=300)
ggsave('/thinker/3.tangjiale/shinanxi/rich_V4/plot/KEGG_facet_ori.ident.png', p2, width=17, height=9, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/rich_V4/plot/KEGG_facet_ori.ident.pdf', p2, width=17, height=9, dpi=300)

xx_go<-compareCluster(gcSample,fun="enrichGO",OrgDb="org.Hs.eg.db",ont="BP",pvalueCutoff=1,qvalueCutoff=1)
res<-xx_go@compareClusterResult

enrich<-res%>%
  group_by(Cluster)%>%
  top_n(n=10,wt=-pvalue)

dt<-enrich
dt<-dt[order(dt$pvalue,dt$Cluster,decreasing=F),]

dt$Description<-factor(dt$Description,levels=unique(dt$Description))
colnames(dt)

library(ggforce)
colors1 <-c(
  'Control' =  "#FFE4B5",  
  'Psoriasis' = "#800080")
p <- ggplot(dt) +
  geom_link(aes(x = 0, y = Description,
                xend = -log10(pvalue), yend = Description,
                alpha = after_stat(index),
                color = Cluster,
                size = after_stat(index)), 
            n = 500, show.legend = T) 
p1 <- p +
  geom_point(aes(x = -log10(pvalue),y = Description), color = "black", fill = "white",size = 6,shape = 21) +
  geom_text(aes(x = -log10(pvalue), y = Description), label=dt$Count, size=3, nudge_x=0.05) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face ="bold.italic"),
        axis.line = element_line(color = "black", size = 0.6),
        axis.text = element_text(face = "bold"),
        axis.title = element_text(size = 13)
  ) +
  xlab("-Log10 Pvalue") + ylab("") + 
  scale_color_manual(values = colors1)

p2 <- p1 + 
  facet_wrap(~Cluster,scales = "free",ncol = 2)
ggsave('/thinker/3.tangjiale/shinanxi/rich_V4/plot/GO_ori.ident.png', p1, width=16, height=8, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/rich_V4/plot/GO_ori.ident.pdf', p1, width=16, height=8, dpi=300)
ggsave('/thinker/3.tangjiale/shinanxi/rich_V4/plot/GO_facet_ori.ident.png', p2, width=15, height=4, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/rich_V4/plot/GO_facet_ori.ident.pdf', p2, width=15, height=4, dpi=300)

mycolor1_clusters<-c(
  "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
  "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
  "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93")
mycolor2_celltype <-c(
  "#FF34B3","#00F5FF","#BC8F8F","#ADFF2F","#FFFF02",
  "#00CD00","#FF6A6A","#7FFFD4", "#AB82FF")
mycolor3_samples=c(
  "#DC143C","#0000FF","#20B2AA","#FFA500","#98FB98","#1E90FF")
mycolor3_group=c(
  "#FFE4B5","#800080")


scRNA = readRDS("/thinker/3.tangjiale/shinanxi/rds/Overall_shinanxi_New.rds")


setwd("/thinker/3.tangjiale/shinanxi/base")
p1 = DimPlot(scRNA, group.by = "clusters",pt.size=0.5,label=F,raster = F)+scale_color_manual(values = mycolor1_clusters)
ggsave('./clusters.png', p1, width=10, height=8, dpi=600)
ggsave('./clusters.pdf', p1, width=10, height=8, dpi=300)
p2 = DimPlot(scRNA, group.by = "group",pt.size=0.5,label=F,raster = F)+scale_color_manual(values = mycolor3_group)
ggsave('./group.png', p2, width=10, height=8, dpi=600)
ggsave('./group.pdf', p2, width=10, height=8, dpi=300)
.libPaths()
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggnetwork)
library(dplyr)
library(ggforce)
library(S4Vectors)
library(ggplot2)
cellType<- levels(MH@meta.data$celltype)
Sample <- unique(MH@meta.data$orig.ident)
{
  mycolor2_celltype_new <- list(
    discrete = mycolor2_celltype
  )
  mycolor3_samples_new <- list(
    discrete = mycolor3_samples
  )
  color_assignments <- setNames(
    c(mycolor2_celltype_new$discrete[1:length(cellType)], mycolor3_samples_new$discrete[1:length(Sample)]),
    c(cellType,Sample)
  )
}

{
  data <- MH@meta.data %>%
    group_by(celltype,orig.ident) %>%
    tally() %>%
    ungroup() %>%
    gather_set_data(1:2) %>%
    dplyr::mutate(
      x = factor(x, levels = unique(x)),
      y = factor(y, levels = unique(y))
    )
  
  DataFrame(data)
  
  data_labels <- tibble(
    group = c(
      rep('celltype', length(cellType)),
      rep('orig.ident', length(Sample))
    )
  ) %>%
    mutate(
      hjust = ifelse(group == 'celltype', 1, 0),
      nudge_x = ifelse(group == 'celltype', -0.1, 0.1)
    )
  
  DataFrame(data_labels)
  
}


p1 <- ggplot(data, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = celltype), alpha = 0.75, axis.width = 0.15) +
  geom_parallel_sets_axes(aes(fill = y), color = 'black', axis.width = 0.1) +
  geom_text(
    aes(y = n, split = y), stat = 'parallel_sets_axes', fontface = 'bold',
    hjust = data_labels$hjust, nudge_x = data_labels$nudge_x
  ) +
  scale_x_discrete(labels = c('celltype','sample')) +
  scale_fill_manual(values = color_assignments) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.title = element_blank(),
    axis.text.x = element_text(face = 'bold', colour = 'black', size = 15),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )

p1
ggsave('./sangji_celltypesample.png', p1, width=12, height=12, dpi=600)
ggsave('./sangji_celltypesample.pdf', p1, width=12, height=12, dpi=300)

cellType<- levels(MH@meta.data$clusters)
Sample <- levels(MH@meta.data$celltype)
{
  mycolor2_celltype_new <- list(
    discrete = mycolor1_clusters
  )
  mycolor3_samples_new <- list(
    discrete = mycolor2_celltype
  )
  color_assignments <- setNames(
    c(mycolor2_celltype_new$discrete[1:length(cellType)], mycolor3_samples_new$discrete[1:length(Sample)]),
    c(cellType,Sample)
  )
}

{
  data <- MH@meta.data %>%
    group_by(clusters,celltype) %>%
    tally() %>%
    ungroup() %>%
    gather_set_data(1:2) %>%
    dplyr::mutate(
      x = factor(x, levels = unique(x)),
      y = factor(y, levels = levels(y))
    )
  
  DataFrame(data)
  
  data_labels <- tibble(
    group = c(
      rep('clusters', length(cellType)),
      rep('celltype', length(Sample))
    )
  ) %>%
    mutate(
      hjust = ifelse(group == 'clusters', 1, 0),
      nudge_x = ifelse(group == 'clusters', -0.1, 0.1)
    )
  
  DataFrame(data_labels)
  
}


p1 <- ggplot(data, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = clusters), alpha = 0.75, axis.width = 0.15) +
  geom_parallel_sets_axes(aes(fill = y), color = 'black', axis.width = 0.1) +
  geom_text(
    aes(y = n, split = y), stat = 'parallel_sets_axes', fontface = 'bold',
    hjust = data_labels$hjust, nudge_x = data_labels$nudge_x
  ) +
  scale_x_discrete(labels = c('clusters','celltype')) +
  scale_fill_manual(values = color_assignments) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.title = element_blank(),
    axis.text.x = element_text(face = 'bold', colour = 'black', size = 15),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )

p1
ggsave('./sangji_clusterscelltype.png', p1, width=12, height=12, dpi=600)
ggsave('./sangji_clusterscelltype.pdf', p1, width=12, height=12, dpi=300)

cellType<- unique(MH@meta.data$orig.ident)
Sample <- unique(MH@meta.data$group)
{
  mycolor2_celltype_new <- list(
    discrete = mycolor3_samples
  )
  mycolor3_samples_new <- list(
    discrete = mycolor3_group
  )
  color_assignments <- setNames(
    c(mycolor2_celltype_new$discrete[1:length(cellType)], mycolor3_samples_new$discrete[1:length(Sample)]),
    c(cellType,Sample)
  )
}

{
  data <- MH@meta.data %>%
    group_by(orig.ident,group) %>%
    tally() %>%
    ungroup() %>%
    gather_set_data(1:2) %>%
    dplyr::mutate(
      x = factor(x, levels = unique(x)),
      y = factor(y, levels = unique(y))
    )
  
  DataFrame(data)
  
  data_labels <- tibble(
    group = c(
      rep('orig.ident', length(cellType)),
      rep('group', length(Sample))
    )
  ) %>%
    mutate(
      hjust = ifelse(group == 'orig.ident', 1, 0),
      nudge_x = ifelse(group == 'orig.ident', -0.1, 0.1)
    )
  
  DataFrame(data_labels)
  
}


p1 <- ggplot(data, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = orig.ident), alpha = 0.75, axis.width = 0.15) +
  geom_parallel_sets_axes(aes(fill = y), color = 'black', axis.width = 0.1) +
  geom_text(
    aes(y = n, split = y), stat = 'parallel_sets_axes', fontface = 'bold',
    hjust = data_labels$hjust, nudge_x = data_labels$nudge_x
  ) +
  scale_x_discrete(labels = c('orig.ident','group')) +
  scale_fill_manual(values = color_assignments) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.title = element_blank(),
    axis.text.x = element_text(face = 'bold', colour = 'black', size = 15),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )

p1
ggsave('./sangji_samplegroup.png', p1, width=12, height=12, dpi=600)
ggsave('./sangji_samplegroup.pdf', p1, width=12, height=12, dpi=300)

ta= CellStatPlot(scRNA, stat.by = "orig.ident", group.by = "celltype", palcolor = mycolor3_samples, stat_type = "count", position = "dodge", label = TRUE)

library(Seurat)
library(tidyverse)
library(ggpubr)

metadata <- seurat_obj@meta.data %>%
  select(celltype, group)

plot_data <- metadata %>%
  count(celltype, group) %>%
  group_by(celltype) %>%
  mutate(
    mean_count = mean(n),
    sd_count = sd(n)
  ) %>% 
  ungroup()

p <- ggplot(plot_data, aes(x = celltype, y = n, fill = group)) +
  geom_col(
    position = position_dodge(0.8),
    width = 0.7,
    alpha = 0.7
  ) +
  geom_errorbar(
    aes(ymin = mean_count - sd_count, ymax = mean_count + sd_count),
    position = position_dodge(0.8),
    width = 0.2,
    color = "black",
    linewidth = 0.5
  ) +
  scale_fill_manual(values = mycolor3_group, name = "Group") +
  scale_color_manual(values = mycolor3_group, name = "Group")+
  labs(
    x = "cellype",
    y = "count",
    fill = "Group"
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plot_data_sample <- seurat_obj@meta.data %>%
  count(celltype, orig.ident, group)

plot_data_group <- plot_data_sample %>%
  group_by(celltype, group) %>%
  summarise(
    mean_count = mean(n),
    sd_count = sd(n)
  ) %>%
  ungroup()

mycolor3_group=c(
  "#FFE4B5","#800080")

p <- ggplot() +
  geom_col(
    data = plot_data_group,
    aes(x = celltype, y = mean_count, fill = group),
    position = position_dodge(0.8), width = 0.7, alpha = 0.7
  ) +
  geom_errorbar(
    data = plot_data_group,
    aes(x = celltype, 
        ymin = mean_count - sd_count, 
        ymax = mean_count + sd_count,
        group = group),
    position = position_dodge(0.8), width = 0.2, 
    color = "black", linewidth = 0.5
  ) +
  geom_point(
    data = plot_data_sample,
    aes(x = celltype, y = n, color = group, group = group),
    position = position_jitterdodge(
      jitter.width = 0.1,
      dodge.width = 0.8
    ),
    size = 2.5, alpha = 0.9
  ) +
  scale_fill_manual(values = mycolor3_group, name = "Group") +
  scale_color_manual(values = mycolor3_group, name = "Group") +
  labs(x = "Cell Type", y = "Cell Count") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(face = "bold", size = 12),
    legend.position = "top"
  )
```

#######10.Cellchat###############
```r
devtools::install_github("immunogenomics/presto")
devtools::install_github("jinworks/CellChat")
library(CellChat)
library(Seurat)

sceV5 <- readRDS("/thinker/3.tangjiale/shinanxi/rds/Overall_shinanxi_New.rds")
Overall = sceV5
unique(sceV5$group)

Control <- subset(sceV5,group=='Control')
Psoriasis <- subset(sceV5,group=='Psoriasis')


Overall_input <- GetAssayData(Overall, layer = 'data')
Overall_meta <- Overall@meta.data[,c("orig.ident","celltype")]
colnames(Overall_meta) <-  c("group","labels")
identical(colnames(Overall_input),rownames(Overall_meta)) 

Overall.cellchat <- createCellChat(object = Overall_input, meta = Overall_meta, group.by = "labels")

levels(Overall.cellchat@idents) 
groupSize <- as.numeric(table(Overall.cellchat@idents)) 

CellChatDB <- CellChatDB.human 

Overall.cellchat@DB <- CellChatDB


Overall.cellchat <- subsetData(Overall.cellchat) 
future::plan("multisession", workers = 30) 
Overall.cellchat <- identifyOverExpressedGenes(Overall.cellchat)
options(future.globals.maxSize = 2000 * 1024^2) 
Overall.cellchat <- identifyOverExpressedInteractions(Overall.cellchat)

groupSize
current_levels <- levels(Overall.cellchat@idents)
print(current_levels)

if (!is.factor(Overall.cellchat@idents)) {
  Overall.cellchat@idents <- factor(Overall.cellchat@idents)
}

Overall.cellchat@idents <- droplevels(Overall.cellchat@idents)

new_levels <- levels(Overall.cellchat@idents)
print(new_levels)

options(future.globals.maxSize = 4000 * 1024^2)  

future::plan("multisession", workers = 32)

Overall.cellchat <- computeCommunProb(Overall.cellchat, type = "triMean")
dir.create("/thinker/3.tangjiale/shinanxi/CellchatKS")
saveRDS(Overall.cellchat,"/thinker/3.tangjiale/shinanxi/CellchatKS/Overall.cellchat.rds")
save.image("/thinker/3.tangjiale/shinanxi/CellchatKS/Overall.cellchat.Rda")



Overall.cellchat <- filterCommunication(Overall.cellchat, min.cells = 10)

Overall.cellchat <- computeCommunProbPathway(Overall.cellchat)

Overall.cellchat <- aggregateNet(Overall.cellchat)

df.net <- subsetCommunication(Overall.cellchat) 



mycolor1_clusters<-c(
  "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
  "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
  "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93")
mycolor2_celltype <-c(
  "#FF34B3","#00F5FF","#BC8F8F","#ADFF2F","#FFFF02",
  "#00CD00","#FF6A6A","#7FFFD4", "#AB82FF")
mycolor3_samples=c(
  "#DC143C","#0000FF","#20B2AA","#FFA500",,"#98FB98","#1E90FF")







table(sceV5$celltype)
colors1 <- c(
  "Transitional activated"="#FF34B3",
  "Cytotoxic effector"="#00F5FF",
  "CSMD1+"="#BC8F8F",
  "NK-like"="#ADFF2F",
  "Activation-regulated"="#FFFF02",
  "Naive-like"="#00CD00",
  "Cytotoxic memory-like"= "#FF6A6A",
  "Terminal-branch" = "#7FFFD4")

groupSize <- as.numeric(table(Overall.cellchat@idents)) 
p1 = netVisual_circle(Overall.cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Overall_Number_of_interactions",color.use = colors1)
p2 = netVisual_circle(Overall.cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Overall_Interaction_weights_strength",color.use = colors1)
dir.create("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/")
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Number_of_interactions.pdf", width = 12, height = 12)
print(p1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Number_of_interactions.png", width = 12, height = 12, res = 600, units = "in")
print(p1)
dev.off()

pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Interaction_weights_strength.pdf", width = 12, height = 12)
print(p2)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Interaction_weights_strength.png", width = 12, height = 12, res = 600, units = "in")
print(p2)
dev.off()
p3 = pheatmap::pheatmap(Overall.cellchat@net$count, border_color = "black", 
                        cluster_cols = F, fontsize = 10, cluster_rows = F,
                        display_numbers = T,number_color="black",number_format = "%.0f")
p4 = pheatmap::pheatmap(Overall.cellchat@net$count, border_color = "black", 
                        cluster_cols = T, fontsize = 10, cluster_rows = T,
                        display_numbers = T,number_color="black",number_format = "%.0f")

p5 = pheatmap::pheatmap(Overall.cellchat@net$weight, border_color = "black", 
                        cluster_cols = F, fontsize = 10, cluster_rows = F,
                        display_numbers = F,number_color="black",number_format = "%.0f")
p6 = pheatmap::pheatmap(Overall.cellchat@net$weight, border_color = "black", 
                        cluster_cols = T, fontsize = 10, cluster_rows = T,
                        display_numbers = F,number_color="black",number_format = "%.0f")
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Number_of_heatmap.png", width = 10, height = 10,dpi = 600,
       plot = p3)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Number_of_heatmap.pdf", width = 10, height = 10, dpi = 300,
       plot = p3)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Number_of_clusterheatmap.png", width = 10, height = 10,dpi = 600,
       plot = p4)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Number_of_clusterheatmap.pdf", width = 10, height = 10, dpi = 300,
       plot = p4)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_of_heatmap.png", width = 10, height = 10,dpi = 600,
       plot = p5)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_of_heatmap.pdf", width = 10, height = 10, dpi = 300,
       plot = p5)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_of_clusterheatmap.png", width = 10, height = 10,dpi = 600,
       plot = p6)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_of_clusterheatmap.pdf", width = 10, height = 10, dpi = 300,
       plot = p6)
table(sceV5$celltype)
celltype_order <- c("Transitional activated","Cytotoxic effector","CSMD1+","NK-like","Activation-regulated","Naive-like","Cytotoxic memory-like","Terminal-branch")
color.use <- colors1

mat <- as.data.frame(Overall.cellchat@net$count)
mat <- mat[celltype_order,]
mat <- mat[,celltype_order] %>% as.matrix()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_number_single_interactions.pdf", width = 12, height = 12)
par(mfrow = c(3,3), xpd=TRUE,mar = c(1.5, 1.5, 1.5, 1.5))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, 
                   weight.scale = T, arrow.size=0.05,
                   arrow.width=1, edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   color.use = color.use)
}
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_number_single_interactions.png", width = 12, height = 12, res = 600, units = "in")
par(mfrow = c(3,3), xpd=TRUE,mar = c(1.5, 1.5, 1.5, 1.5))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, 
                   weight.scale = T, arrow.size=0.05,
                   arrow.width=1, edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   color.use = color.use)
}
dev.off()


mat <- as.data.frame(Overall.cellchat@net$weight)
mat <- mat[celltype_order,]
mat <- mat[,celltype_order] %>% as.matrix()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_single_interactions.pdf", width = 12, height = 12)
par(mfrow = c(3,3), xpd=TRUE,mar = c(1.5, 1.5, 1.5, 1.5))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, 
                   weight.scale = T, arrow.size=0.05,
                   arrow.width=1, edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   color.use = color.use)
}
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_single_interactions.png", width = 12, height = 12, res = 600, units = "in")
par(mfrow = c(3,3), xpd=TRUE,mar = c(1.5, 1.5, 1.5, 1.5))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, 
                   weight.scale = T, arrow.size=0.05,
                   arrow.width=1, edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   color.use = color.use)
}
dev.off()

Overall.cellchat@netP$pathways


p1_1 = netVisual_aggregate(Overall.cellchat, signaling = "CLEC",color.use = colors1)
p1_2 = netVisual_aggregate(Overall.cellchat, signaling = "MHC-I",color.use = colors1)
p1_3 = netVisual_aggregate(Overall.cellchat, signaling = "MIF",color.use = colors1)
p1_4 = netVisual_aggregate(Overall.cellchat, signaling = "CD99",color.use = colors1)
p1_5 = netVisual_aggregate(Overall.cellchat, signaling = "ADGRE",color.use = colors1)
p1_6 = netVisual_aggregate(Overall.cellchat, signaling = "TGFb",color.use = colors1)

pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_CLEC_interactions.pdf", width = 12, height = 12)
print(p1_1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_CLEC_interactions.png", width = 12, height = 12, res = 600, units = "in")
print(p1_1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_MHC-I_interactions.pdf", width = 12, height = 12)
print(p1_2)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_MHC-I_interactions.png", width = 12, height = 12, res = 600, units = "in")
print(p1_2)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_MIF_interactions.pdf", width = 12, height = 12)
print(p1_3)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_MIF_interactions.png", width = 12, height = 12, res = 600, units = "in")
print(p1_3)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_CD99_interactions.pdf", width = 12, height = 12)
print(p1_4)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_CD99_interactions.png", width = 12, height = 12, res = 600, units = "in")
print(p1_4)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_ADGRE_interactions.pdf", width = 12, height = 12)
print(p1_5)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_ADGRE_interactions.png", width = 12, height = 12, res = 600, units = "in")
print(p1_5)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_TGFb_interactions.pdf", width = 12, height = 12)
print(p1_6)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_TGFb_interactions.png", width = 12, height = 12, res = 600, units = "in")
print(p1_6)
dev.off()
p2_1 = netVisual_aggregate(Overall.cellchat, signaling = "CLEC", layout = "chord",color.use = colors1)
p2_2 = netVisual_aggregate(Overall.cellchat, signaling = "MHC-I", layout = "chord",color.use = colors1)
p2_3 = netVisual_aggregate(Overall.cellchat, signaling = "MIF", layout = "chord",color.use = colors1)
p2_4 = netVisual_aggregate(Overall.cellchat, signaling = "CD99", layout = "chord",color.use = colors1)
p2_5 = netVisual_aggregate(Overall.cellchat, signaling = "ADGRE", layout = "chord",color.use = colors1)
p2_6 = netVisual_aggregate(Overall.cellchat, signaling = "TGFb", layout = "chord",color.use = colors1)

pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_CLEC_cirplot.pdf", width = 12, height = 12)
print(p2_1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_CLEC_cirplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_MHC-I_cirplot.pdf", width = 12, height = 12)
print(p2_2)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_MHC-I_cirplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_2)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_MIF_cirplot.pdf", width = 12, height = 12)
print(p2_3)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_MIF_cirplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_3)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_CD99_cirplot.pdf", width = 12, height = 12)
print(p2_4)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_CD99_cirplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_4)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_ADGRE_cirplot.pdf", width = 12, height = 12)
print(p2_5)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_ADGRE_cirplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_5)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_TGFb_cirplot.pdf", width = 12, height = 12)
print(p2_6)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_TGFb_cirplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_6)
dev.off()
p3_1 = netVisual_heatmap(Overall.cellchat, signaling = "CLEC",color.heatmap = "Reds",color.use = colors1)
p3_2 = netVisual_heatmap(Overall.cellchat, signaling = "MHC-I",color.heatmap = "Reds",color.use = colors1)
p3_3 = netVisual_heatmap(Overall.cellchat, signaling = "MIF",color.heatmap = "Reds",color.use = colors1)
p3_4 = netVisual_heatmap(Overall.cellchat, signaling = "CD99",color.heatmap = "Reds",color.use = colors1)
p3_5 = netVisual_heatmap(Overall.cellchat, signaling = "ADGRE",color.heatmap = "Reds",color.use = colors1)
p3_6 = netVisual_heatmap(Overall.cellchat, signaling = "TGFb",color.heatmap = "Reds",color.use = colors1)
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_CLEC_heatmap.pdf", width = 12, height = 12)
print(p3_1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_CLEC_heatmap.png", width = 12, height = 12, res = 600, units = "in")
print(p3_1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_MHC-I_heatmap.pdf", width = 12, height = 12)
print(p3_2)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_MHC-I_heatmap.png", width = 12, height = 12, res = 600, units = "in")
print(p3_2)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_MIF_heatmap.pdf", width = 12, height = 12)
print(p3_3)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_MIF_heatmap.png", width = 12, height = 12, res = 600, units = "in")
print(p3_3)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_CD99_heatmap.pdf", width = 12, height = 12)
print(p3_4)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_CD99_heatmap.png", width = 12, height = 12, res = 600, units = "in")
print(p3_4)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_ADGRE_heatmap.pdf", width = 12, height = 12)
print(p3_5)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_ADGRE_heatmap.png", width = 12, height = 12, res = 600, units = "in")
print(p3_5)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_TGFb_heatmap.pdf", width = 12, height = 12)
print(p3_6)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_TGFb_heatmap.png", width = 12, height = 12, res = 600, units = "in")
print(p3_6)
dev.off()
p4 = netVisual_heatmap(Overall.cellchat,color.heatmap = "Reds",color.use = colors1)
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_all_heatmap.pdf", width = 12, height = 12)
print(p4)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_weight_all_heatmap.png", width = 12, height = 12, res = 600, units = "in")
print(p4)
dev.off()


p2_1 = netAnalysis_contribution(Overall.cellchat, signaling = "CLEC")
p2_2 = netAnalysis_contribution(Overall.cellchat, signaling = "MHC-I")
p2_3 = netAnalysis_contribution(Overall.cellchat, signaling = "MIF")
p2_4 = netAnalysis_contribution(Overall.cellchat, signaling = "CD99")
p2_5 = netAnalysis_contribution(Overall.cellchat, signaling = "ADGRE")
p2_6 = netAnalysis_contribution(Overall.cellchat, signaling = "TGFb")

pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_CLEC_dotplot.pdf", width = 12, height = 12)
print(p2_1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_CLEC_dotplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_MHC-I_dotplot.pdf", width = 12, height = 12)
print(p2_2)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_MHC-I_dotplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_2)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_MIF_dotplot.pdf", width = 12, height = 12)
print(p2_3)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_MIF_dotplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_3)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_CD99_dotplot.pdf", width = 12, height = 12)
print(p2_4)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_CD99_dotplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_4)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_ADGRE_dotplot.pdf", width = 12, height = 12)
print(p2_5)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_ADGRE_dotplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_5)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_TGFb_dotplot.pdf", width = 12, height = 12)
print(p2_6)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_TGFb_dotplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_6)
dev.off()
pairLR1 <- extractEnrichedLR(Overall.cellchat, signaling = "CLEC", geneLR.return = FALSE)
pairLR2 <- extractEnrichedLR(Overall.cellchat, signaling = "MHC-I", geneLR.return = FALSE)
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_CLEC2D_KLRB1_interactions.pdf", width = 12, height = 12)
netVisual_individual(Overall.cellchat, signaling = "CLEC", pairLR.use = "CLEC2D_KLRB1", layout = "circle",color.use = colors1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_CLEC2D_KLRB1_interactions.png", width = 12, height = 12, res = 600, units = "in")
netVisual_individual(Overall.cellchat, signaling = "CLEC", pairLR.use = "CLEC2D_KLRB1", layout = "circle",color.use = colors1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_CLEC2C_KLRB1_interactions.pdf", width = 12, height = 12)
netVisual_individual(Overall.cellchat, signaling = "CLEC", pairLR.use = "CLEC2C_KLRB1", layout = "circle",color.use = colors1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_CLEC2C_KLRB1_interactions.png", width = 12, height = 12, res = 600, units = "in")
netVisual_individual(Overall.cellchat, signaling = "CLEC", pairLR.use = "CLEC2C_KLRB1", layout = "circle",color.use = colors1)
dev.off()




pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_HLA-E_KLRK1_interactions.pdf", width = 12, height = 12)
netVisual_individual(Overall.cellchat, signaling = "MHC-I", pairLR.use = "HLA-E_KLRK1", layout = "circle",color.use = colors1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_HLA-E_KLRK1_interactions.png", width = 12, height = 12, res = 600, units = "in")
netVisual_individual(Overall.cellchat, signaling = "MHC-I", pairLR.use = "HLA-E_KLRK1", layout = "circle",color.use = colors1)
dev.off()


p = netVisual_bubble(Overall.cellchat, remove.isolate = FALSE)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_all_sender_receiver.png", width = 15, height = 10,dpi = 600,
       plot = p)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_all_sender_receiver.pdf", width = 15, height = 10, dpi = 300,
       plot = p)



p2 = netVisual_bubble(Overall.cellchat, remove.isolate = FALSE,sources.use = "Cytotoxic memory-like")
p3 = netVisual_bubble(Overall.cellchat, remove.isolate = FALSE,targets.use = "Cytotoxic memory-like")
p4 = netVisual_bubble(Overall.cellchat, remove.isolate = FALSE,sources.use = "Activation-regulated")
p5 = netVisual_bubble(Overall.cellchat, remove.isolate = FALSE,targets.use = "Activation-regulated")
p6 = netVisual_bubble(Overall.cellchat, remove.isolate = FALSE,sources.use = "CSMD1+")
p7 = netVisual_bubble(Overall.cellchat, remove.isolate = FALSE,targets.use = "CSMD1+")
p8 = netVisual_bubble(Overall.cellchat, remove.isolate = FALSE,sources.use = "Transitional activated")
p9 = netVisual_bubble(Overall.cellchat, remove.isolate = FALSE,targets.use = "Transitional activated")
p10 = netVisual_bubble(Overall.cellchat, remove.isolate = FALSE,sources.use = "Naive-like")
p11 = netVisual_bubble(Overall.cellchat, remove.isolate = FALSE,targets.use = "Naive-like")
p12 = netVisual_bubble(Overall.cellchat, remove.isolate = FALSE,sources.use = "Cytotoxic effector")
p13 = netVisual_bubble(Overall.cellchat, remove.isolate = FALSE,targets.use = "Cytotoxic effector")
p14 = netVisual_bubble(Overall.cellchat, remove.isolate = FALSE,sources.use = "NK-like")
p15 = netVisual_bubble(Overall.cellchat, remove.isolate = FALSE,targets.use = "NK-like")
p16 = netVisual_bubble(Overall.cellchat, remove.isolate = FALSE,sources.use = "Terminal-branch")
p17 = netVisual_bubble(Overall.cellchat, remove.isolate = FALSE,targets.use = "Terminal-branch")
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Cytotoxic memory-like_sender.png", width = 10, height = 10,dpi = 600,
       plot = p2)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Cytotoxic memory-like_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p2)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Cytotoxic memory-like_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p3)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Cytotoxic memory-like_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p3)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Activation-regulated_sender.png", width = 10, height = 10,dpi = 600,
       plot = p4)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Activation-regulated_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p4)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Activation-regulated_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p5)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Activation-regulated_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p5)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_CSMD1+_sender.png", width = 10, height = 10,dpi = 600,
       plot = p6)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_CSMD1+_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p6)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_CSMD1+_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p7)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_CSMD1+_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p7)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Transitional activated_sender.png", width = 10, height = 10,dpi = 600,
       plot = p8)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Transitional activated_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p8)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Transitional activated_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p9)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Transitional activated_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p9)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Naive-like_sender.png", width = 10, height = 10,dpi = 600,
       plot = p10)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Naive-like_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p10)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Naive-like_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p11)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Naive-like_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p11)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Cytotoxic effector_sender.png", width = 10, height = 10,dpi = 600,
       plot = p12)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Cytotoxic effector_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p12)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Cytotoxic effector_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p13)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Cytotoxic effector_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p13)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_NK-like_sender.png", width = 10, height = 10,dpi = 600,
       plot = p14)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_NK-like_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p14)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_NK-like_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p15)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_NK-like_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p15)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Terminal-branch_sender.png", width = 10, height = 10,dpi = 600,
       plot = p16)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Terminal-branchs_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p16)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Terminal-branch_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p17)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_Terminal-branch_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p17)


Overall.cellchat <- netAnalysis_computeCentrality(Overall.cellchat, slot.name = "netP")
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_CLEC_charactor.png", width = 10, height = 8, res = 600, units = "in")
netAnalysis_signalingRole_network(Overall.cellchat, signaling = "CLEC", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_CLEC_charactor.pdf", width = 10, height = 8)
netAnalysis_signalingRole_network(Overall.cellchat, signaling = "CLEC", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_MHC-I_charactor.png", width = 10, height = 8, res = 600, units = "in")
netAnalysis_signalingRole_network(Overall.cellchat, signaling = "MHC-I", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_MHC-I_charactor.pdf", width = 10, height = 8)
netAnalysis_signalingRole_network(Overall.cellchat, signaling = "MHC-I", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_MIF_charactor.png", width = 10, height = 8, res = 600, units = "in")
netAnalysis_signalingRole_network(Overall.cellchat, signaling = "MIF", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_MIF_charactor.pdf", width = 10, height = 8)
netAnalysis_signalingRole_network(Overall.cellchat, signaling = "MIF", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_CD99_charactor.png", width = 10, height = 8, res = 600, units = "in")
netAnalysis_signalingRole_network(Overall.cellchat, signaling = "CD99", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_CD99_charactor.pdf", width = 10, height = 8)
netAnalysis_signalingRole_network(Overall.cellchat, signaling = "CD99", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_ADGRE_charactor.png", width = 10, height = 8, res = 600, units = "in")
netAnalysis_signalingRole_network(Overall.cellchat, signaling = "ADGRE", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_ADGRE_charactor.pdf", width = 10, height = 8)
netAnalysis_signalingRole_network(Overall.cellchat, signaling = "ADGRE", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_TGFb_charactor.png", width = 10, height = 8, res = 600, units = "in")
netAnalysis_signalingRole_network(Overall.cellchat, signaling = "TGFb", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_TGFb_charactor.pdf", width = 10, height = 8)
netAnalysis_signalingRole_network(Overall.cellchat, signaling = "TGFb", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
gg1 <- netAnalysis_signalingRole_scatter(Overall.cellchat,color.use = colors1)
gg2_1 <- netAnalysis_signalingRole_scatter(Overall.cellchat, signaling = "CLEC",color.use = colors1)
gg2_2 <- netAnalysis_signalingRole_scatter(Overall.cellchat, signaling = "MHC-I",color.use = colors1)
gg2_3 <- netAnalysis_signalingRole_scatter(Overall.cellchat, signaling = "MIF",color.use = colors1)
gg2_4 <- netAnalysis_signalingRole_scatter(Overall.cellchat, signaling = "CD99",color.use = colors1)
gg2_5 <- netAnalysis_signalingRole_scatter(Overall.cellchat, signaling = "ADGRE",color.use = colors1)
gg2_6 <- netAnalysis_signalingRole_scatter(Overall.cellchat, signaling = "TGFb",color.use = colors1)
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_all_signalingRole.png", width = 8, height = 8, res = 600, units = "in")
print(gg1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_all_signalingRole.pdf", width = 8, height = 8)
print(gg1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_CLEC_signalingRole.png", width = 8, height = 8, res = 600, units = "in")
print(gg2_1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_CLEC_signalingRole.pdf", width = 8, height = 8)
print(gg2_1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_MHC-I_signalingRole.png", width = 8, height = 8, res = 600, units = "in")
print(gg2_2)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_MHC-I_signalingRole.pdf", width = 8, height = 8)
print(gg2_2)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_MIF_signalingRole.png", width = 8, height = 8, res = 600, units = "in")
print(gg2_3)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_MIF_signalingRole.pdf", width = 8, height = 8)
print(gg2_3)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_CD99_signalingRole.png", width = 8, height = 8, res = 600, units = "in")
print(gg2_4)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_CD99_signalingRole.pdf", width = 8, height = 8)
print(gg2_4)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_ADGRE_signalingRole.png", width = 8, height = 8, res = 600, units = "in")
print(gg2_5)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_ADGRE_signalingRole.pdf", width = 8, height = 8)
print(gg2_5)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_TGFb_signalingRole.png", width = 8, height = 8, res = 600, units = "in")
print(gg2_6)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_TGFb_signalingRole.pdf", width = 8, height = 8)
print(gg2_6)
dev.off()


ht1 <- netAnalysis_signalingRole_heatmap(Overall.cellchat, pattern = "outgoing",width=10,height = 10,color.use = colors1)
ht2 <- netAnalysis_signalingRole_heatmap(Overall.cellchat, pattern = "incoming",width=10,height = 10,color.use = colors1)
p = ht1 + ht2
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_all_outgoing.png", width = 10, height = 15, res = 600, units = "in")
print(ht1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_all_outgoing.pdf", width = 10, height = 15)
print(ht1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_all_incoming.png", width = 10, height = 15, res = 600, units = "in")
print(ht2)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_all_incoming.pdf", width = 10, height = 15)
print(ht2)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_all_outin.png", width = 10, height = 15, res = 600, units = "in")
print(p)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_all_outin.pdf", width = 10, height = 15)
print(p)
dev.off()

library(NMF)
library(ggalluvial)

p = netAnalysis_river(Overall.cellchat, pattern = "outgoing",font.size = 1.5)
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_all_outgoing_riverplot.png", width = 8, height = 20, res = 600, units = "in")
print(p)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Overall_all_outgiong_riverplot.pdf", width = 8, height = 20)
print(p)
dev.off()





Control_input <- GetAssayData(Control, layer = 'data')
Control_meta <- Control@meta.data[,c("orig.ident","celltype")]
colnames(Control_meta) <-  c("group","labels")
identical(colnames(Control_input),rownames(Control_meta)) 

Control.cellchat <- createCellChat(object = Control_input, meta = Control_meta, group.by = "labels")

levels(Control.cellchat@idents) 
groupSize <- as.numeric(table(Control.cellchat@idents)) 

CellChatDB <- CellChatDB.human 

Control.cellchat@DB <- CellChatDB


Control.cellchat <- subsetData(Control.cellchat) 
future::plan("multisession", workers = 30) 
Control.cellchat <- identifyOverExpressedGenes(Control.cellchat)
options(future.globals.maxSize = 2000 * 1024^2) 
Control.cellchat <- identifyOverExpressedInteractions(Control.cellchat)

groupSize
current_levels <- levels(Control.cellchat@idents)
print(current_levels)

if (!is.factor(Control.cellchat@idents)) {
  Control.cellchat@idents <- factor(Control.cellchat@idents)
}

Control.cellchat@idents <- droplevels(Control.cellchat@idents)

new_levels <- levels(Control.cellchat@idents)
print(new_levels)

options(future.globals.maxSize = 4000 * 1024^2)  

future::plan("multisession", workers = 32)

Control.cellchat <- computeCommunProb(Control.cellchat, type = "triMean")
dir.create("/thinker/3.tangjiale/shinanxi/CellchatKS")
saveRDS(Control.cellchat,"/thinker/3.tangjiale/shinanxi/CellchatKS/Control.cellchat.rds")
save.image("/thinker/3.tangjiale/shinanxi/CellchatKS/Control.cellchat.Rda")



Control.cellchat <- filterCommunication(Control.cellchat, min.cells = 10)

Control.cellchat <- computeCommunProbPathway(Control.cellchat)

Control.cellchat <- aggregateNet(Control.cellchat)

df.net <- subsetCommunication(Control.cellchat) 



mycolor1_clusters<-c(
  "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
  "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
  "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93")
mycolor2_celltype <-c(
  "#FF34B3","#00F5FF","#BC8F8F","#ADFF2F","#FFFF02",
  "#00CD00","#FF6A6A","#7FFFD4", "#AB82FF")
mycolor3_samples=c(
  "#DC143C","#0000FF","#20B2AA","#FFA500",,"#98FB98","#1E90FF")







table(sceV5$celltype)
colors1 <- c(
  "Transitional activated"="#FF34B3",
  "Cytotoxic effector"="#00F5FF",
  "CSMD1+"="#BC8F8F",
  "NK-like"="#ADFF2F",
  "Activation-regulated"="#FFFF02",
  "Naive-like"="#00CD00",
  "Cytotoxic memory-like"= "#FF6A6A",
  "Terminal-branch" = "#7FFFD4")

groupSize <- as.numeric(table(Control.cellchat@idents)) 
p1 = netVisual_circle(Control.cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Control_Number_of_interactions",color.use = colors1)
p2 = netVisual_circle(Control.cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Control_Interaction_weights_strength",color.use = colors1)
dir.create("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/")
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Number_of_interactions.pdf", width = 12, height = 12)
print(p1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Number_of_interactions.png", width = 12, height = 12, res = 600, units = "in")
print(p1)
dev.off()

pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Interaction_weights_strength.pdf", width = 12, height = 12)
print(p2)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Interaction_weights_strength.png", width = 12, height = 12, res = 600, units = "in")
print(p2)
dev.off()
p3 = pheatmap::pheatmap(Control.cellchat@net$count, border_color = "black", 
                        cluster_cols = F, fontsize = 10, cluster_rows = F,
                        display_numbers = T,number_color="black",number_format = "%.0f")
p4 = pheatmap::pheatmap(Control.cellchat@net$count, border_color = "black", 
                        cluster_cols = T, fontsize = 10, cluster_rows = T,
                        display_numbers = T,number_color="black",number_format = "%.0f")

p5 = pheatmap::pheatmap(Control.cellchat@net$weight, border_color = "black", 
                        cluster_cols = F, fontsize = 10, cluster_rows = F,
                        display_numbers = F,number_color="black",number_format = "%.0f")
p6 = pheatmap::pheatmap(Control.cellchat@net$weight, border_color = "black", 
                        cluster_cols = T, fontsize = 10, cluster_rows = T,
                        display_numbers = F,number_color="black",number_format = "%.0f")
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Number_of_heatmap.png", width = 10, height = 10,dpi = 600,
       plot = p3)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Number_of_heatmap.pdf", width = 10, height = 10, dpi = 300,
       plot = p3)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Number_of_clusterheatmap.png", width = 10, height = 10,dpi = 600,
       plot = p4)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Number_of_clusterheatmap.pdf", width = 10, height = 10, dpi = 300,
       plot = p4)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_of_heatmap.png", width = 10, height = 10,dpi = 600,
       plot = p5)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_of_heatmap.pdf", width = 10, height = 10, dpi = 300,
       plot = p5)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_of_clusterheatmap.png", width = 10, height = 10,dpi = 600,
       plot = p6)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_of_clusterheatmap.pdf", width = 10, height = 10, dpi = 300,
       plot = p6)
table(sceV5$celltype)
celltype_order <- c("Transitional activated","Cytotoxic effector","CSMD1+","NK-like","Activation-regulated","Naive-like","Cytotoxic memory-like","Terminal-branch")
color.use <- colors1

mat <- as.data.frame(Control.cellchat@net$count)
mat <- mat[celltype_order,]
mat <- mat[,celltype_order] %>% as.matrix()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_number_single_interactions.pdf", width = 12, height = 12)
par(mfrow = c(3,3), xpd=TRUE,mar = c(1.5, 1.5, 1.5, 1.5))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, 
                   weight.scale = T, arrow.size=0.05,
                   arrow.width=1, edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   color.use = color.use)
}
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_number_single_interactions.png", width = 12, height = 12, res = 600, units = "in")
par(mfrow = c(3,3), xpd=TRUE,mar = c(1.5, 1.5, 1.5, 1.5))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, 
                   weight.scale = T, arrow.size=0.05,
                   arrow.width=1, edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   color.use = color.use)
}
dev.off()


mat <- as.data.frame(Control.cellchat@net$weight)
mat <- mat[celltype_order,]
mat <- mat[,celltype_order] %>% as.matrix()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_single_interactions.pdf", width = 12, height = 12)
par(mfrow = c(3,3), xpd=TRUE,mar = c(1.5, 1.5, 1.5, 1.5))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, 
                   weight.scale = T, arrow.size=0.05,
                   arrow.width=1, edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   color.use = color.use)
}
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_single_interactions.png", width = 12, height = 12, res = 600, units = "in")
par(mfrow = c(3,3), xpd=TRUE,mar = c(1.5, 1.5, 1.5, 1.5))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, 
                   weight.scale = T, arrow.size=0.05,
                   arrow.width=1, edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   color.use = color.use)
}
dev.off()

Control.cellchat@netP$pathways


p1_1 = netVisual_aggregate(Control.cellchat, signaling = "CLEC",color.use = colors1)
p1_2 = netVisual_aggregate(Control.cellchat, signaling = "MHC-I",color.use = colors1)
p1_3 = netVisual_aggregate(Control.cellchat, signaling = "MIF",color.use = colors1)
p1_4 = netVisual_aggregate(Control.cellchat, signaling = "CD99",color.use = colors1)
p1_5 = netVisual_aggregate(Control.cellchat, signaling = "ADGRE",color.use = colors1)
p1_6 = netVisual_aggregate(Control.cellchat, signaling = "TGFb",color.use = colors1)

pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_CLEC_interactions.pdf", width = 12, height = 12)
print(p1_1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_CLEC_interactions.png", width = 12, height = 12, res = 600, units = "in")
print(p1_1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_MHC-I_interactions.pdf", width = 12, height = 12)
print(p1_2)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_MHC-I_interactions.png", width = 12, height = 12, res = 600, units = "in")
print(p1_2)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_MIF_interactions.pdf", width = 12, height = 12)
print(p1_3)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_MIF_interactions.png", width = 12, height = 12, res = 600, units = "in")
print(p1_3)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_CD99_interactions.pdf", width = 12, height = 12)
print(p1_4)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_CD99_interactions.png", width = 12, height = 12, res = 600, units = "in")
print(p1_4)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_ADGRE_interactions.pdf", width = 12, height = 12)
print(p1_5)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_ADGRE_interactions.png", width = 12, height = 12, res = 600, units = "in")
print(p1_5)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_TGFb_interactions.pdf", width = 12, height = 12)
print(p1_6)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_TGFb_interactions.png", width = 12, height = 12, res = 600, units = "in")
print(p1_6)
dev.off()
p2_1 = netVisual_aggregate(Control.cellchat, signaling = "CLEC", layout = "chord",color.use = colors1)
p2_2 = netVisual_aggregate(Control.cellchat, signaling = "MHC-I", layout = "chord",color.use = colors1)
p2_3 = netVisual_aggregate(Control.cellchat, signaling = "MIF", layout = "chord",color.use = colors1)
p2_4 = netVisual_aggregate(Control.cellchat, signaling = "CD99", layout = "chord",color.use = colors1)
p2_5 = netVisual_aggregate(Control.cellchat, signaling = "ADGRE", layout = "chord",color.use = colors1)
p2_6 = netVisual_aggregate(Control.cellchat, signaling = "TGFb", layout = "chord",color.use = colors1)

pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_CLEC_cirplot.pdf", width = 12, height = 12)
print(p2_1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_CLEC_cirplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_MHC-I_cirplot.pdf", width = 12, height = 12)
print(p2_2)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_MHC-I_cirplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_2)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_MIF_cirplot.pdf", width = 12, height = 12)
print(p2_3)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_MIF_cirplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_3)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_CD99_cirplot.pdf", width = 12, height = 12)
print(p2_4)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_CD99_cirplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_4)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_ADGRE_cirplot.pdf", width = 12, height = 12)
print(p2_5)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_ADGRE_cirplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_5)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_TGFb_cirplot.pdf", width = 12, height = 12)
print(p2_6)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_TGFb_cirplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_6)
dev.off()
p3_1 = netVisual_heatmap(Control.cellchat, signaling = "CLEC",color.heatmap = "Reds",color.use = colors1)
p3_2 = netVisual_heatmap(Control.cellchat, signaling = "MHC-I",color.heatmap = "Reds",color.use = colors1)
p3_3 = netVisual_heatmap(Control.cellchat, signaling = "MIF",color.heatmap = "Reds",color.use = colors1)
p3_4 = netVisual_heatmap(Control.cellchat, signaling = "CD99",color.heatmap = "Reds",color.use = colors1)
p3_5 = netVisual_heatmap(Control.cellchat, signaling = "ADGRE",color.heatmap = "Reds",color.use = colors1)
p3_6 = netVisual_heatmap(Control.cellchat, signaling = "TGFb",color.heatmap = "Reds",color.use = colors1)
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_CLEC_heatmap.pdf", width = 12, height = 12)
print(p3_1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_CLEC_heatmap.png", width = 12, height = 12, res = 600, units = "in")
print(p3_1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_MHC-I_heatmap.pdf", width = 12, height = 12)
print(p3_2)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_MHC-I_heatmap.png", width = 12, height = 12, res = 600, units = "in")
print(p3_2)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_MIF_heatmap.pdf", width = 12, height = 12)
print(p3_3)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_MIF_heatmap.png", width = 12, height = 12, res = 600, units = "in")
print(p3_3)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_CD99_heatmap.pdf", width = 12, height = 12)
print(p3_4)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_CD99_heatmap.png", width = 12, height = 12, res = 600, units = "in")
print(p3_4)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_ADGRE_heatmap.pdf", width = 12, height = 12)
print(p3_5)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_ADGRE_heatmap.png", width = 12, height = 12, res = 600, units = "in")
print(p3_5)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_TGFb_heatmap.pdf", width = 12, height = 12)
print(p3_6)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_TGFb_heatmap.png", width = 12, height = 12, res = 600, units = "in")
print(p3_6)
dev.off()
p4 = netVisual_heatmap(Control.cellchat,color.heatmap = "Reds",color.use = colors1)
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_all_heatmap.pdf", width = 12, height = 12)
print(p4)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_weight_all_heatmap.png", width = 12, height = 12, res = 600, units = "in")
print(p4)
dev.off()


p2_1 = netAnalysis_contribution(Control.cellchat, signaling = "CLEC")
p2_2 = netAnalysis_contribution(Control.cellchat, signaling = "MHC-I")
p2_3 = netAnalysis_contribution(Control.cellchat, signaling = "MIF")
p2_4 = netAnalysis_contribution(Control.cellchat, signaling = "CD99")
p2_5 = netAnalysis_contribution(Control.cellchat, signaling = "ADGRE")
p2_6 = netAnalysis_contribution(Control.cellchat, signaling = "TGFb")

pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_CLEC_dotplot.pdf", width = 12, height = 12)
print(p2_1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_CLEC_dotplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_MHC-I_dotplot.pdf", width = 12, height = 12)
print(p2_2)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_MHC-I_dotplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_2)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_MIF_dotplot.pdf", width = 12, height = 12)
print(p2_3)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_MIF_dotplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_3)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_CD99_dotplot.pdf", width = 12, height = 12)
print(p2_4)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_CD99_dotplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_4)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_ADGRE_dotplot.pdf", width = 12, height = 12)
print(p2_5)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_ADGRE_dotplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_5)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_TGFb_dotplot.pdf", width = 12, height = 12)
print(p2_6)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_TGFb_dotplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_6)
dev.off()
pairLR1 <- extractEnrichedLR(Control.cellchat, signaling = "CLEC", geneLR.return = FALSE)
pairLR2 <- extractEnrichedLR(Control.cellchat, signaling = "MHC-I", geneLR.return = FALSE)
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_CLEC2D_KLRB1_interactions.pdf", width = 12, height = 12)
netVisual_individual(Control.cellchat, signaling = "CLEC", pairLR.use = "CLEC2D_KLRB1", layout = "circle",color.use = colors1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_CLEC2D_KLRB1_interactions.png", width = 12, height = 12, res = 600, units = "in")
netVisual_individual(Control.cellchat, signaling = "CLEC", pairLR.use = "CLEC2D_KLRB1", layout = "circle",color.use = colors1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_CLEC2C_KLRB1_interactions.pdf", width = 12, height = 12)
netVisual_individual(Control.cellchat, signaling = "CLEC", pairLR.use = "CLEC2C_KLRB1", layout = "circle",color.use = colors1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_CLEC2C_KLRB1_interactions.png", width = 12, height = 12, res = 600, units = "in")
netVisual_individual(Control.cellchat, signaling = "CLEC", pairLR.use = "CLEC2C_KLRB1", layout = "circle",color.use = colors1)
dev.off()




pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_HLA-E_KLRK1_interactions.pdf", width = 12, height = 12)
netVisual_individual(Control.cellchat, signaling = "MHC-I", pairLR.use = "HLA-E_KLRK1", layout = "circle",color.use = colors1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_HLA-E_KLRK1_interactions.png", width = 12, height = 12, res = 600, units = "in")
netVisual_individual(Control.cellchat, signaling = "MHC-I", pairLR.use = "HLA-E_KLRK1", layout = "circle",color.use = colors1)
dev.off()


p = netVisual_bubble(Control.cellchat, remove.isolate = FALSE)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_all_sender_receiver.png", width = 15, height = 10,dpi = 600,
       plot = p)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_all_sender_receiver.pdf", width = 15, height = 10, dpi = 300,
       plot = p)



p2 = netVisual_bubble(Control.cellchat, remove.isolate = FALSE,sources.use = "Cytotoxic memory-like")
p3 = netVisual_bubble(Control.cellchat, remove.isolate = FALSE,targets.use = "Cytotoxic memory-like")
p4 = netVisual_bubble(Control.cellchat, remove.isolate = FALSE,sources.use = "Activation-regulated")
p5 = netVisual_bubble(Control.cellchat, remove.isolate = FALSE,targets.use = "Activation-regulated")
p6 = netVisual_bubble(Control.cellchat, remove.isolate = FALSE,sources.use = "CSMD1+")
p7 = netVisual_bubble(Control.cellchat, remove.isolate = FALSE,targets.use = "CSMD1+")
p8 = netVisual_bubble(Control.cellchat, remove.isolate = FALSE,sources.use = "Transitional activated")
p9 = netVisual_bubble(Control.cellchat, remove.isolate = FALSE,targets.use = "Transitional activated")
p10 = netVisual_bubble(Control.cellchat, remove.isolate = FALSE,sources.use = "Naive-like")
p11 = netVisual_bubble(Control.cellchat, remove.isolate = FALSE,targets.use = "Naive-like")
p12 = netVisual_bubble(Control.cellchat, remove.isolate = FALSE,sources.use = "Cytotoxic effector")
p13 = netVisual_bubble(Control.cellchat, remove.isolate = FALSE,targets.use = "Cytotoxic effector")
p14 = netVisual_bubble(Control.cellchat, remove.isolate = FALSE,sources.use = "NK-like")
p15 = netVisual_bubble(Control.cellchat, remove.isolate = FALSE,targets.use = "NK-like")
p16 = netVisual_bubble(Control.cellchat, remove.isolate = FALSE,sources.use = "Terminal-branch")
p17 = netVisual_bubble(Control.cellchat, remove.isolate = FALSE,targets.use = "Terminal-branch")
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Cytotoxic memory-like_sender.png", width = 10, height = 10,dpi = 600,
       plot = p2)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Cytotoxic memory-like_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p2)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Cytotoxic memory-like_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p3)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Cytotoxic memory-like_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p3)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Activation-regulated_sender.png", width = 10, height = 10,dpi = 600,
       plot = p4)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Activation-regulated_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p4)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Activation-regulated_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p5)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Activation-regulated_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p5)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_CSMD1+_sender.png", width = 10, height = 10,dpi = 600,
       plot = p6)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_CSMD1+_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p6)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_CSMD1+_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p7)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_CSMD1+_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p7)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Transitional activated_sender.png", width = 10, height = 10,dpi = 600,
       plot = p8)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Transitional activated_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p8)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Transitional activated_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p9)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Transitional activated_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p9)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Naive-like_sender.png", width = 10, height = 10,dpi = 600,
       plot = p10)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Naive-like_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p10)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Naive-like_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p11)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Naive-like_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p11)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Cytotoxic effector_sender.png", width = 10, height = 10,dpi = 600,
       plot = p12)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Cytotoxic effector_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p12)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Cytotoxic effector_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p13)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Cytotoxic effector_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p13)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_NK-like_sender.png", width = 10, height = 10,dpi = 600,
       plot = p14)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_NK-like_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p14)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_NK-like_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p15)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_NK-like_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p15)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Terminal-branch_sender.png", width = 10, height = 10,dpi = 600,
       plot = p16)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Terminal-branchs_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p16)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Terminal-branch_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p17)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_Terminal-branch_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p17)


Control.cellchat <- netAnalysis_computeCentrality(Control.cellchat, slot.name = "netP")
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_CLEC_charactor.png", width = 10, height = 8, res = 600, units = "in")
netAnalysis_signalingRole_network(Control.cellchat, signaling = "CLEC", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_CLEC_charactor.pdf", width = 10, height = 8)
netAnalysis_signalingRole_network(Control.cellchat, signaling = "CLEC", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_MHC-I_charactor.png", width = 10, height = 8, res = 600, units = "in")
netAnalysis_signalingRole_network(Control.cellchat, signaling = "MHC-I", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_MHC-I_charactor.pdf", width = 10, height = 8)
netAnalysis_signalingRole_network(Control.cellchat, signaling = "MHC-I", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_MIF_charactor.png", width = 10, height = 8, res = 600, units = "in")
netAnalysis_signalingRole_network(Control.cellchat, signaling = "MIF", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_MIF_charactor.pdf", width = 10, height = 8)
netAnalysis_signalingRole_network(Control.cellchat, signaling = "MIF", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_CD99_charactor.png", width = 10, height = 8, res = 600, units = "in")
netAnalysis_signalingRole_network(Control.cellchat, signaling = "CD99", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_CD99_charactor.pdf", width = 10, height = 8)
netAnalysis_signalingRole_network(Control.cellchat, signaling = "CD99", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_ADGRE_charactor.png", width = 10, height = 8, res = 600, units = "in")
netAnalysis_signalingRole_network(Control.cellchat, signaling = "ADGRE", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_ADGRE_charactor.pdf", width = 10, height = 8)
netAnalysis_signalingRole_network(Control.cellchat, signaling = "ADGRE", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_TGFb_charactor.png", width = 10, height = 8, res = 600, units = "in")
netAnalysis_signalingRole_network(Control.cellchat, signaling = "TGFb", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_TGFb_charactor.pdf", width = 10, height = 8)
netAnalysis_signalingRole_network(Control.cellchat, signaling = "TGFb", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
gg1 <- netAnalysis_signalingRole_scatter(Control.cellchat,color.use = colors1)
gg2_1 <- netAnalysis_signalingRole_scatter(Control.cellchat, signaling = "CLEC",color.use = colors1)
gg2_2 <- netAnalysis_signalingRole_scatter(Control.cellchat, signaling = "MHC-I",color.use = colors1)
gg2_3 <- netAnalysis_signalingRole_scatter(Control.cellchat, signaling = "MIF",color.use = colors1)
gg2_4 <- netAnalysis_signalingRole_scatter(Control.cellchat, signaling = "CD99",color.use = colors1)
gg2_5 <- netAnalysis_signalingRole_scatter(Control.cellchat, signaling = "ADGRE",color.use = colors1)
gg2_6 <- netAnalysis_signalingRole_scatter(Control.cellchat, signaling = "TGFb",color.use = colors1)
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_all_signalingRole.png", width = 8, height = 8, res = 600, units = "in")
print(gg1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_all_signalingRole.pdf", width = 8, height = 8)
print(gg1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_CLEC_signalingRole.png", width = 8, height = 8, res = 600, units = "in")
print(gg2_1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_CLEC_signalingRole.pdf", width = 8, height = 8)
print(gg2_1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_MHC-I_signalingRole.png", width = 8, height = 8, res = 600, units = "in")
print(gg2_2)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_MHC-I_signalingRole.pdf", width = 8, height = 8)
print(gg2_2)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_MIF_signalingRole.png", width = 8, height = 8, res = 600, units = "in")
print(gg2_3)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_MIF_signalingRole.pdf", width = 8, height = 8)
print(gg2_3)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_CD99_signalingRole.png", width = 8, height = 8, res = 600, units = "in")
print(gg2_4)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_CD99_signalingRole.pdf", width = 8, height = 8)
print(gg2_4)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_ADGRE_signalingRole.png", width = 8, height = 8, res = 600, units = "in")
print(gg2_5)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_ADGRE_signalingRole.pdf", width = 8, height = 8)
print(gg2_5)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_TGFb_signalingRole.png", width = 8, height = 8, res = 600, units = "in")
print(gg2_6)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_TGFb_signalingRole.pdf", width = 8, height = 8)
print(gg2_6)
dev.off()


ht1 <- netAnalysis_signalingRole_heatmap(Control.cellchat, pattern = "outgoing",width=10,height = 10,color.use = colors1)
ht2 <- netAnalysis_signalingRole_heatmap(Control.cellchat, pattern = "incoming",width=10,height = 10,color.use = colors1)
p = ht1 + ht2
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_all_outgoing.png", width = 10, height = 15, res = 600, units = "in")
print(ht1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_all_outgoing.pdf", width = 10, height = 15)
print(ht1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_all_incoming.png", width = 10, height = 15, res = 600, units = "in")
print(ht2)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_all_incoming.pdf", width = 10, height = 15)
print(ht2)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_all_outin.png", width = 10, height = 15, res = 600, units = "in")
print(p)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_all_outin.pdf", width = 10, height = 15)
print(p)
dev.off()

library(NMF)
library(ggalluvial)

p = netAnalysis_river(Control.cellchat, pattern = "outgoing",font.size = 1.5)
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_all_outgoing_riverplot.png", width = 8, height = 20, res = 600, units = "in")
print(p)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Control_all_outgiong_riverplot.pdf", width = 8, height = 20)
print(p)
dev.off()





Psoriasis_input <- GetAssayData(Psoriasis, layer = 'data')
Psoriasis_meta <- Psoriasis@meta.data[,c("orig.ident","celltype")]
colnames(Psoriasis_meta) <-  c("group","labels")
identical(colnames(Psoriasis_input),rownames(Psoriasis_meta)) 

Psoriasis.cellchat <- createCellChat(object = Psoriasis_input, meta = Psoriasis_meta, group.by = "labels")

levels(Psoriasis.cellchat@idents) 
groupSize <- as.numeric(table(Psoriasis.cellchat@idents)) 

CellChatDB <- CellChatDB.human 

Psoriasis.cellchat@DB <- CellChatDB


Psoriasis.cellchat <- subsetData(Psoriasis.cellchat) 
future::plan("multisession", workers = 30) 
Psoriasis.cellchat <- identifyOverExpressedGenes(Psoriasis.cellchat)
options(future.globals.maxSize = 2000 * 1024^2) 
Psoriasis.cellchat <- identifyOverExpressedInteractions(Psoriasis.cellchat)

groupSize
current_levels <- levels(Psoriasis.cellchat@idents)
print(current_levels)

if (!is.factor(Psoriasis.cellchat@idents)) {
  Psoriasis.cellchat@idents <- factor(Psoriasis.cellchat@idents)
}

Psoriasis.cellchat@idents <- droplevels(Psoriasis.cellchat@idents)

new_levels <- levels(Psoriasis.cellchat@idents)
print(new_levels)

options(future.globals.maxSize = 4000 * 1024^2)  

future::plan("multisession", workers = 32)

Psoriasis.cellchat <- computeCommunProb(Psoriasis.cellchat, type = "triMean")
dir.create("/thinker/3.tangjiale/shinanxi/CellchatKS")
saveRDS(Psoriasis.cellchat,"/thinker/3.tangjiale/shinanxi/CellchatKS/Psoriasis.cellchat.rds")
save.image("/thinker/3.tangjiale/shinanxi/CellchatKS/Psoriasis.cellchat.Rda")



Psoriasis.cellchat <- filterCommunication(Psoriasis.cellchat, min.cells = 10)

Psoriasis.cellchat <- computeCommunProbPathway(Psoriasis.cellchat)

Psoriasis.cellchat <- aggregateNet(Psoriasis.cellchat)

df.net <- subsetCommunication(Psoriasis.cellchat) 



mycolor1_clusters<-c(
  "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
  "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
  "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93")
mycolor2_celltype <-c(
  "#FF34B3","#00F5FF","#BC8F8F","#ADFF2F","#FFFF02",
  "#00CD00","#FF6A6A","#7FFFD4", "#AB82FF")
mycolor3_samples=c(
  "#DC143C","#0000FF","#20B2AA","#FFA500","#98FB98","#1E90FF")







table(sceV5$celltype)
colors1 <- c(
  "Transitional activated"="#FF34B3",
  "Cytotoxic effector"="#00F5FF",
  "CSMD1+"="#BC8F8F",
  "NK-like"="#ADFF2F",
  "Activation-regulated"="#FFFF02",
  "Naive-like"="#00CD00",
  "Cytotoxic memory-like"= "#FF6A6A",
  "Terminal-branch" = "#7FFFD4")

groupSize <- as.numeric(table(Psoriasis.cellchat@idents)) 
p1 = netVisual_circle(Psoriasis.cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Psoriasis_Number_of_interactions",color.use = colors1)
p2 = netVisual_circle(Psoriasis.cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Psoriasis_Interaction_weights_strength",color.use = colors1)
dir.create("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/")
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Number_of_interactions.pdf", width = 12, height = 12)
print(p1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Number_of_interactions.png", width = 12, height = 12, res = 600, units = "in")
print(p1)
dev.off()

pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Interaction_weights_strength.pdf", width = 12, height = 12)
print(p2)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Interaction_weights_strength.png", width = 12, height = 12, res = 600, units = "in")
print(p2)
dev.off()
p3 = pheatmap::pheatmap(Psoriasis.cellchat@net$count, border_color = "black", 
                        cluster_cols = F, fontsize = 10, cluster_rows = F,
                        display_numbers = T,number_color="black",number_format = "%.0f")
p4 = pheatmap::pheatmap(Psoriasis.cellchat@net$count, border_color = "black", 
                        cluster_cols = T, fontsize = 10, cluster_rows = T,
                        display_numbers = T,number_color="black",number_format = "%.0f")

p5 = pheatmap::pheatmap(Psoriasis.cellchat@net$weight, border_color = "black", 
                        cluster_cols = F, fontsize = 10, cluster_rows = F,
                        display_numbers = F,number_color="black",number_format = "%.0f")
p6 = pheatmap::pheatmap(Psoriasis.cellchat@net$weight, border_color = "black", 
                        cluster_cols = T, fontsize = 10, cluster_rows = T,
                        display_numbers = F,number_color="black",number_format = "%.0f")
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Number_of_heatmap.png", width = 10, height = 10,dpi = 600,
       plot = p3)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Number_of_heatmap.pdf", width = 10, height = 10, dpi = 300,
       plot = p3)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Number_of_clusterheatmap.png", width = 10, height = 10,dpi = 600,
       plot = p4)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Number_of_clusterheatmap.pdf", width = 10, height = 10, dpi = 300,
       plot = p4)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_of_heatmap.png", width = 10, height = 10,dpi = 600,
       plot = p5)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_of_heatmap.pdf", width = 10, height = 10, dpi = 300,
       plot = p5)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_of_clusterheatmap.png", width = 10, height = 10,dpi = 600,
       plot = p6)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_of_clusterheatmap.pdf", width = 10, height = 10, dpi = 300,
       plot = p6)
table(sceV5$celltype)
celltype_order <- c("Transitional activated","Cytotoxic effector","CSMD1+","NK-like","Activation-regulated","Naive-like","Cytotoxic memory-like","Terminal-branch")
color.use <- colors1

mat <- as.data.frame(Psoriasis.cellchat@net$count)
mat <- mat[celltype_order,]
mat <- mat[,celltype_order] %>% as.matrix()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_number_single_interactions.pdf", width = 12, height = 12)
par(mfrow = c(3,3), xpd=TRUE,mar = c(1.5, 1.5, 1.5, 1.5))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, 
                   weight.scale = T, arrow.size=0.05,
                   arrow.width=1, edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   color.use = color.use)
}
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_number_single_interactions.png", width = 12, height = 12, res = 600, units = "in")
par(mfrow = c(3,3), xpd=TRUE,mar = c(1.5, 1.5, 1.5, 1.5))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, 
                   weight.scale = T, arrow.size=0.05,
                   arrow.width=1, edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   color.use = color.use)
}
dev.off()


mat <- as.data.frame(Psoriasis.cellchat@net$weight)
mat <- mat[celltype_order,]
mat <- mat[,celltype_order] %>% as.matrix()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_single_interactions.pdf", width = 12, height = 12)
par(mfrow = c(3,3), xpd=TRUE,mar = c(1.5, 1.5, 1.5, 1.5))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, 
                   weight.scale = T, arrow.size=0.05,
                   arrow.width=1, edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   color.use = color.use)
}
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_single_interactions.png", width = 12, height = 12, res = 600, units = "in")
par(mfrow = c(3,3), xpd=TRUE,mar = c(1.5, 1.5, 1.5, 1.5))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, 
                   weight.scale = T, arrow.size=0.05,
                   arrow.width=1, edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   color.use = color.use)
}
dev.off()

Psoriasis.cellchat@netP$pathways


p1_1 = netVisual_aggregate(Psoriasis.cellchat, signaling = "CLEC",color.use = colors1)
p1_2 = netVisual_aggregate(Psoriasis.cellchat, signaling = "MHC-I",color.use = colors1)
p1_3 = netVisual_aggregate(Psoriasis.cellchat, signaling = "MIF",color.use = colors1)
p1_4 = netVisual_aggregate(Psoriasis.cellchat, signaling = "CD99",color.use = colors1)
p1_5 = netVisual_aggregate(Psoriasis.cellchat, signaling = "ADGRE",color.use = colors1)
p1_6 = netVisual_aggregate(Psoriasis.cellchat, signaling = "TGFb",color.use = colors1)

pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_CLEC_interactions.pdf", width = 12, height = 12)
print(p1_1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_CLEC_interactions.png", width = 12, height = 12, res = 600, units = "in")
print(p1_1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_MHC-I_interactions.pdf", width = 12, height = 12)
print(p1_2)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_MHC-I_interactions.png", width = 12, height = 12, res = 600, units = "in")
print(p1_2)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_MIF_interactions.pdf", width = 12, height = 12)
print(p1_3)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_MIF_interactions.png", width = 12, height = 12, res = 600, units = "in")
print(p1_3)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_CD99_interactions.pdf", width = 12, height = 12)
print(p1_4)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_CD99_interactions.png", width = 12, height = 12, res = 600, units = "in")
print(p1_4)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_ADGRE_interactions.pdf", width = 12, height = 12)
print(p1_5)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_ADGRE_interactions.png", width = 12, height = 12, res = 600, units = "in")
print(p1_5)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_TGFb_interactions.pdf", width = 12, height = 12)
print(p1_6)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_TGFb_interactions.png", width = 12, height = 12, res = 600, units = "in")
print(p1_6)
dev.off()
p2_1 = netVisual_aggregate(Psoriasis.cellchat, signaling = "CLEC", layout = "chord",color.use = colors1)
p2_2 = netVisual_aggregate(Psoriasis.cellchat, signaling = "MHC-I", layout = "chord",color.use = colors1)
p2_3 = netVisual_aggregate(Psoriasis.cellchat, signaling = "MIF", layout = "chord",color.use = colors1)
p2_4 = netVisual_aggregate(Psoriasis.cellchat, signaling = "CD99", layout = "chord",color.use = colors1)
p2_5 = netVisual_aggregate(Psoriasis.cellchat, signaling = "ADGRE", layout = "chord",color.use = colors1)
p2_6 = netVisual_aggregate(Psoriasis.cellchat, signaling = "TGFb", layout = "chord",color.use = colors1)

pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_CLEC_cirplot.pdf", width = 12, height = 12)
print(p2_1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_CLEC_cirplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_MHC-I_cirplot.pdf", width = 12, height = 12)
print(p2_2)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_MHC-I_cirplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_2)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_MIF_cirplot.pdf", width = 12, height = 12)
print(p2_3)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_MIF_cirplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_3)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_CD99_cirplot.pdf", width = 12, height = 12)
print(p2_4)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_CD99_cirplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_4)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_ADGRE_cirplot.pdf", width = 12, height = 12)
print(p2_5)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_ADGRE_cirplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_5)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_TGFb_cirplot.pdf", width = 12, height = 12)
print(p2_6)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_TGFb_cirplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_6)
dev.off()
p3_1 = netVisual_heatmap(Psoriasis.cellchat, signaling = "CLEC",color.heatmap = "Reds",color.use = colors1)
p3_2 = netVisual_heatmap(Psoriasis.cellchat, signaling = "MHC-I",color.heatmap = "Reds",color.use = colors1)
p3_3 = netVisual_heatmap(Psoriasis.cellchat, signaling = "MIF",color.heatmap = "Reds",color.use = colors1)
p3_4 = netVisual_heatmap(Psoriasis.cellchat, signaling = "CD99",color.heatmap = "Reds",color.use = colors1)
p3_5 = netVisual_heatmap(Psoriasis.cellchat, signaling = "ADGRE",color.heatmap = "Reds",color.use = colors1)
p3_6 = netVisual_heatmap(Psoriasis.cellchat, signaling = "TGFb",color.heatmap = "Reds",color.use = colors1)
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_CLEC_heatmap.pdf", width = 12, height = 12)
print(p3_1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_CLEC_heatmap.png", width = 12, height = 12, res = 600, units = "in")
print(p3_1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_MHC-I_heatmap.pdf", width = 12, height = 12)
print(p3_2)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_MHC-I_heatmap.png", width = 12, height = 12, res = 600, units = "in")
print(p3_2)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_MIF_heatmap.pdf", width = 12, height = 12)
print(p3_3)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_MIF_heatmap.png", width = 12, height = 12, res = 600, units = "in")
print(p3_3)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_CD99_heatmap.pdf", width = 12, height = 12)
print(p3_4)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_CD99_heatmap.png", width = 12, height = 12, res = 600, units = "in")
print(p3_4)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_ADGRE_heatmap.pdf", width = 12, height = 12)
print(p3_5)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_ADGRE_heatmap.png", width = 12, height = 12, res = 600, units = "in")
print(p3_5)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_TGFb_heatmap.pdf", width = 12, height = 12)
print(p3_6)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_TGFb_heatmap.png", width = 12, height = 12, res = 600, units = "in")
print(p3_6)
dev.off()
p4 = netVisual_heatmap(Psoriasis.cellchat,color.heatmap = "Reds",color.use = colors1)
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_all_heatmap.pdf", width = 12, height = 12)
print(p4)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_weight_all_heatmap.png", width = 12, height = 12, res = 600, units = "in")
print(p4)
dev.off()


p2_1 = netAnalysis_contribution(Psoriasis.cellchat, signaling = "CLEC")
p2_2 = netAnalysis_contribution(Psoriasis.cellchat, signaling = "MHC-I")
p2_3 = netAnalysis_contribution(Psoriasis.cellchat, signaling = "MIF")
p2_4 = netAnalysis_contribution(Psoriasis.cellchat, signaling = "CD99")
p2_5 = netAnalysis_contribution(Psoriasis.cellchat, signaling = "ADGRE")
p2_6 = netAnalysis_contribution(Psoriasis.cellchat, signaling = "TGFb")

pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_CLEC_dotplot.pdf", width = 12, height = 12)
print(p2_1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_CLEC_dotplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_MHC-I_dotplot.pdf", width = 12, height = 12)
print(p2_2)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_MHC-I_dotplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_2)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_MIF_dotplot.pdf", width = 12, height = 12)
print(p2_3)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_MIF_dotplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_3)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_CD99_dotplot.pdf", width = 12, height = 12)
print(p2_4)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_CD99_dotplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_4)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_ADGRE_dotplot.pdf", width = 12, height = 12)
print(p2_5)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_ADGRE_dotplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_5)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_TGFb_dotplot.pdf", width = 12, height = 12)
print(p2_6)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_TGFb_dotplot.png", width = 12, height = 12, res = 600, units = "in")
print(p2_6)
dev.off()
pairLR1 <- extractEnrichedLR(Psoriasis.cellchat, signaling = "CLEC", geneLR.return = FALSE)
pairLR2 <- extractEnrichedLR(Psoriasis.cellchat, signaling = "MHC-I", geneLR.return = FALSE)
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_CLEC2D_KLRB1_interactions.pdf", width = 12, height = 12)
netVisual_individual(Psoriasis.cellchat, signaling = "CLEC", pairLR.use = "CLEC2D_KLRB1", layout = "circle",color.use = colors1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_CLEC2D_KLRB1_interactions.png", width = 12, height = 12, res = 600, units = "in")
netVisual_individual(Psoriasis.cellchat, signaling = "CLEC", pairLR.use = "CLEC2D_KLRB1", layout = "circle",color.use = colors1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_CLEC2C_KLRB1_interactions.pdf", width = 12, height = 12)
netVisual_individual(Psoriasis.cellchat, signaling = "CLEC", pairLR.use = "CLEC2C_KLRB1", layout = "circle",color.use = colors1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_CLEC2C_KLRB1_interactions.png", width = 12, height = 12, res = 600, units = "in")
netVisual_individual(Psoriasis.cellchat, signaling = "CLEC", pairLR.use = "CLEC2C_KLRB1", layout = "circle",color.use = colors1)
dev.off()




pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_HLA-E_KLRK1_interactions.pdf", width = 12, height = 12)
netVisual_individual(Psoriasis.cellchat, signaling = "MHC-I", pairLR.use = "HLA-E_KLRK1", layout = "circle",color.use = colors1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_HLA-E_KLRK1_interactions.png", width = 12, height = 12, res = 600, units = "in")
netVisual_individual(Psoriasis.cellchat, signaling = "MHC-I", pairLR.use = "HLA-E_KLRK1", layout = "circle",color.use = colors1)
dev.off()


p = netVisual_bubble(Psoriasis.cellchat, remove.isolate = FALSE)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_all_sender_receiver.png", width = 15, height = 10,dpi = 600,
       plot = p)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_all_sender_receiver.pdf", width = 15, height = 10, dpi = 300,
       plot = p)



p2 = netVisual_bubble(Psoriasis.cellchat, remove.isolate = FALSE,sources.use = "Cytotoxic memory-like")
p3 = netVisual_bubble(Psoriasis.cellchat, remove.isolate = FALSE,targets.use = "Cytotoxic memory-like")
p4 = netVisual_bubble(Psoriasis.cellchat, remove.isolate = FALSE,sources.use = "Activation-regulated")
p5 = netVisual_bubble(Psoriasis.cellchat, remove.isolate = FALSE,targets.use = "Activation-regulated")
p6 = netVisual_bubble(Psoriasis.cellchat, remove.isolate = FALSE,sources.use = "CSMD1+")
p7 = netVisual_bubble(Psoriasis.cellchat, remove.isolate = FALSE,targets.use = "CSMD1+")
p8 = netVisual_bubble(Psoriasis.cellchat, remove.isolate = FALSE,sources.use = "Transitional activated")
p9 = netVisual_bubble(Psoriasis.cellchat, remove.isolate = FALSE,targets.use = "Transitional activated")
p10 = netVisual_bubble(Psoriasis.cellchat, remove.isolate = FALSE,sources.use = "Naive-like")
p11 = netVisual_bubble(Psoriasis.cellchat, remove.isolate = FALSE,targets.use = "Naive-like")
p12 = netVisual_bubble(Psoriasis.cellchat, remove.isolate = FALSE,sources.use = "Cytotoxic effector")
p13 = netVisual_bubble(Psoriasis.cellchat, remove.isolate = FALSE,targets.use = "Cytotoxic effector")
p14 = netVisual_bubble(Psoriasis.cellchat, remove.isolate = FALSE,sources.use = "NK-like")
p15 = netVisual_bubble(Psoriasis.cellchat, remove.isolate = FALSE,targets.use = "NK-like")
p16 = netVisual_bubble(Psoriasis.cellchat, remove.isolate = FALSE,sources.use = "Terminal-branch")
p17 = netVisual_bubble(Psoriasis.cellchat, remove.isolate = FALSE,targets.use = "Terminal-branch")
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Cytotoxic memory-like_sender.png", width = 10, height = 10,dpi = 600,
       plot = p2)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Cytotoxic memory-like_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p2)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Cytotoxic memory-like_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p3)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Cytotoxic memory-like_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p3)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Activation-regulated_sender.png", width = 10, height = 10,dpi = 600,
       plot = p4)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Activation-regulated_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p4)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Activation-regulated_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p5)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Activation-regulated_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p5)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_CSMD1+_sender.png", width = 10, height = 10,dpi = 600,
       plot = p6)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_CSMD1+_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p6)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_CSMD1+_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p7)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_CSMD1+_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p7)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Transitional activated_sender.png", width = 10, height = 10,dpi = 600,
       plot = p8)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Transitional activated_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p8)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Transitional activated_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p9)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Transitional activated_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p9)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Naive-like_sender.png", width = 10, height = 10,dpi = 600,
       plot = p10)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Naive-like_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p10)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Naive-like_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p11)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Naive-like_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p11)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Cytotoxic effector_sender.png", width = 10, height = 10,dpi = 600,
       plot = p12)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Cytotoxic effector_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p12)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Cytotoxic effector_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p13)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Cytotoxic effector_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p13)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_NK-like_sender.png", width = 10, height = 10,dpi = 600,
       plot = p14)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_NK-like_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p14)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_NK-like_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p15)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_NK-like_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p15)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Terminal-branch_sender.png", width = 10, height = 10,dpi = 600,
       plot = p16)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Terminal-branchs_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p16)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Terminal-branch_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p17)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_Terminal-branch_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p17)


Psoriasis.cellchat <- netAnalysis_computeCentrality(Psoriasis.cellchat, slot.name = "netP")
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_CLEC_charactor.png", width = 10, height = 8, res = 600, units = "in")
netAnalysis_signalingRole_network(Psoriasis.cellchat, signaling = "CLEC", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_CLEC_charactor.pdf", width = 10, height = 8)
netAnalysis_signalingRole_network(Psoriasis.cellchat, signaling = "CLEC", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_MHC-I_charactor.png", width = 10, height = 8, res = 600, units = "in")
netAnalysis_signalingRole_network(Psoriasis.cellchat, signaling = "MHC-I", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_MHC-I_charactor.pdf", width = 10, height = 8)
netAnalysis_signalingRole_network(Psoriasis.cellchat, signaling = "MHC-I", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_MIF_charactor.png", width = 10, height = 8, res = 600, units = "in")
netAnalysis_signalingRole_network(Psoriasis.cellchat, signaling = "MIF", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_MIF_charactor.pdf", width = 10, height = 8)
netAnalysis_signalingRole_network(Psoriasis.cellchat, signaling = "MIF", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_CD99_charactor.png", width = 10, height = 8, res = 600, units = "in")
netAnalysis_signalingRole_network(Psoriasis.cellchat, signaling = "CD99", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_CD99_charactor.pdf", width = 10, height = 8)
netAnalysis_signalingRole_network(Psoriasis.cellchat, signaling = "CD99", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_ADGRE_charactor.png", width = 10, height = 8, res = 600, units = "in")
netAnalysis_signalingRole_network(Psoriasis.cellchat, signaling = "ADGRE", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_ADGRE_charactor.pdf", width = 10, height = 8)
netAnalysis_signalingRole_network(Psoriasis.cellchat, signaling = "ADGRE", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_TGFb_charactor.png", width = 10, height = 8, res = 600, units = "in")
netAnalysis_signalingRole_network(Psoriasis.cellchat, signaling = "TGFb", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_TGFb_charactor.pdf", width = 10, height = 8)
netAnalysis_signalingRole_network(Psoriasis.cellchat, signaling = "TGFb", width = 8, height = 2.5, font.size = 10,color.use = colors1)
dev.off()
gg1 <- netAnalysis_signalingRole_scatter(Psoriasis.cellchat,color.use = colors1)
gg2_1 <- netAnalysis_signalingRole_scatter(Psoriasis.cellchat, signaling = "CLEC",color.use = colors1)
gg2_2 <- netAnalysis_signalingRole_scatter(Psoriasis.cellchat, signaling = "MHC-I",color.use = colors1)
gg2_3 <- netAnalysis_signalingRole_scatter(Psoriasis.cellchat, signaling = "MIF",color.use = colors1)
gg2_4 <- netAnalysis_signalingRole_scatter(Psoriasis.cellchat, signaling = "CD99",color.use = colors1)
gg2_5 <- netAnalysis_signalingRole_scatter(Psoriasis.cellchat, signaling = "ADGRE",color.use = colors1)
gg2_6 <- netAnalysis_signalingRole_scatter(Psoriasis.cellchat, signaling = "TGFb",color.use = colors1)
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_all_signalingRole.png", width = 8, height = 8, res = 600, units = "in")
print(gg1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_all_signalingRole.pdf", width = 8, height = 8)
print(gg1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_CLEC_signalingRole.png", width = 8, height = 8, res = 600, units = "in")
print(gg2_1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_CLEC_signalingRole.pdf", width = 8, height = 8)
print(gg2_1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_MHC-I_signalingRole.png", width = 8, height = 8, res = 600, units = "in")
print(gg2_2)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_MHC-I_signalingRole.pdf", width = 8, height = 8)
print(gg2_2)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_MIF_signalingRole.png", width = 8, height = 8, res = 600, units = "in")
print(gg2_3)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_MIF_signalingRole.pdf", width = 8, height = 8)
print(gg2_3)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_CD99_signalingRole.png", width = 8, height = 8, res = 600, units = "in")
print(gg2_4)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_CD99_signalingRole.pdf", width = 8, height = 8)
print(gg2_4)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_ADGRE_signalingRole.png", width = 8, height = 8, res = 600, units = "in")
print(gg2_5)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_ADGRE_signalingRole.pdf", width = 8, height = 8)
print(gg2_5)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_TGFb_signalingRole.png", width = 8, height = 8, res = 600, units = "in")
print(gg2_6)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_TGFb_signalingRole.pdf", width = 8, height = 8)
print(gg2_6)
dev.off()


ht1 <- netAnalysis_signalingRole_heatmap(Psoriasis.cellchat, pattern = "outgoing",width=10,height = 10,color.use = colors1)
ht2 <- netAnalysis_signalingRole_heatmap(Psoriasis.cellchat, pattern = "incoming",width=10,height = 10,color.use = colors1)
p = ht1 + ht2
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_all_outgoing.png", width = 10, height = 15, res = 600, units = "in")
print(ht1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_all_outgoing.pdf", width = 10, height = 15)
print(ht1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_all_incoming.png", width = 10, height = 15, res = 600, units = "in")
print(ht2)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_all_incoming.pdf", width = 10, height = 15)
print(ht2)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_all_outin.png", width = 10, height = 15, res = 600, units = "in")
print(p)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_all_outin.pdf", width = 10, height = 15)
print(p)
dev.off()

library(NMF)
library(ggalluvial)

p = netAnalysis_river(Psoriasis.cellchat, pattern = "outgoing",font.size = 1.5)
png("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_all_outgoing_riverplot.png", width = 8, height = 20, res = 600, units = "in")
print(p)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/plot/Psoriasis_all_outgiong_riverplot.pdf", width = 8, height = 20)
print(p)
dev.off()
```

###########11.CellChat Multiple Group Comparative Analysis####################
```r
Control.cellchat = readRDS("/thinker/3.tangjiale/shinanxi/CellchatKS/Control.cellchat.rds")
Control.cellchat <- filterCommunication(Control.cellchat, min.cells = 10)
Control.cellchat <- computeCommunProbPathway(Control.cellchat)
Control.cellchat <- aggregateNet(Control.cellchat)

Psoriasis.cellchat = readRDS("/thinker/3.tangjiale/shinanxi/CellchatKS/Psoriasis.cellchat.rds")
Psoriasis.cellchat <- filterCommunication(Psoriasis.cellchat, min.cells = 10)
Psoriasis.cellchat <- computeCommunProbPathway(Psoriasis.cellchat)
Psoriasis.cellchat <- aggregateNet(Psoriasis.cellchat)

print(levels(Control.cellchat@idents))
print(levels(Psoriasis.cellchat@idents))
common_celltypes <- intersect(levels(Control.cellchat@idents), levels(Psoriasis.cellchat@idents))
Control.cellchat <- subsetCellChat(Control.cellchat, idents.use = common_celltypes)
Psoriasis.cellchat <- subsetCellChat(Psoriasis.cellchat, idents.use = common_celltypes)
Control.cellchat@idents <- factor(Control.cellchat@idents, levels = common_celltypes)
Psoriasis.cellchat@idents <- factor(Psoriasis.cellchat@idents, levels = common_celltypes)
object.list <- list(Control = Control.cellchat, Psoriasis = Psoriasis.cellchat)
i = 1
object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
object.list[[i+1]] <- netAnalysis_computeCentrality(object.list[[i+1]])
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

colors1 <- c(
  "Transitional activated"="#FF34B3",
  "Cytotoxic effector"="#00F5FF",
  "CSMD1+"="#BC8F8F",
  "NK-like"="#ADFF2F",
  "Activation-regulated"="#FFFF02",
  "Naive-like"="#00CD00",
  "Cytotoxic memory-like"= "#FF6A6A",
  "Terminal-branch" = "#7FFFD4")

color_Control_Psoriasis = c("#FFE4B5","#800080")

dir.create("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot")
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2),color.use = color_Control_Psoriasis)
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight",color.use = color_Control_Psoriasis)
p = gg1 + gg2
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_number_histogram.png", width = 10, height = 14,dpi = 600,
       plot = gg1)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_number_histogram.pdf", width = 10, height = 14, dpi = 300,
       plot = gg1)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_weight_histogram.png", width = 10, height = 14,dpi = 600,
       plot = gg2)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_weight_histogram.pdf", width = 10, height = 14, dpi = 300,
       plot = gg2)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_nw_histogram.png", width = 10, height = 14,dpi = 600,
       plot = p)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_nw_histogram.pdf", width = 10, height = 14, dpi = 300,
       plot = p)
coloredge = c("#E63863","#2166ac")
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/number_interaction_Psoriasis_compared_to_Control_redis_higher.pdf", width = 10, height = 10)
netVisual_diffInteraction(cellchat, weight.scale = T,color.use = colors1,color.edge = coloredge)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/number_interaction_Psoriasis_compared_to_Control_redis_higher.png", width = 10, height = 10, res = 600, units = "in")
netVisual_diffInteraction(cellchat, weight.scale = T,color.use = colors1,,color.edge = coloredge)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/weight_interaction_Psoriasis_compared_to_Control_redis_higher.pdf", width = 10, height = 10)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",color.use = colors1,color.edge = coloredge)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/weight_interaction_Psoriasis_compared_to_Control_redis_higher.png", width = 10, height = 10, res = 600, units = "in")
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",color.use = colors1,,color.edge = coloredge)
dev.off()

gg1 <- netVisual_heatmap(cellchat,color.use = colors1)
gg2 <- netVisual_heatmap(cellchat, measure = "weight",color.use = colors1)
p = gg1 + gg2
png("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/number_heatmap_Psoriasis_compared_to_Control_redis_higher.png", width = 10, height = 10, res = 600, units = "in")
print(gg1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/number_heatmap_Psoriasis_compared_to_Control_redis_higher.pdf", width = 10, height = 10)
print(gg1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/weight_heatmap_Psoriasis_compared_to_Control_redis_higher.png", width = 10, height = 10, res = 600, units = "in")
print(gg2)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/weight_heatmap_Psoriasis_compared_to_Control_redis_higher.pdf", width = 10, height = 10)
print(gg2)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/nw_heatmap_Psoriasis_compared_to_Control_redis_higher.png", width = 20, height = 10, res = 600, units = "in")
print(p)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/nw_heatmap_Psoriasis_compared_to_Control_redis_higher.pdf", width = 20, height = 10)
print(p)
dev.off()

gg3_1 <- netVisual_heatmap(cellchat,color.use = colors1,signaling = "CLEC",title.name = "Differential number of interactions(CLEC)")
gg4_1 <- netVisual_heatmap(cellchat,color.use = colors1,signaling = "CLEC",measure = "weight",title.name = "Differential weight of interactions(CLEC)")
gg3_2 <- netVisual_heatmap(cellchat,color.use = colors1,signaling = "MHC-I",title.name = "Differential number of interactions(MHC-I)")
gg4_2 <- netVisual_heatmap(cellchat,color.use = colors1,signaling = "MHC-I",measure = "weight",title.name = "Differential weight of interactions(MHC-I)")

png("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/CLEC_number_heatmap_Psoriasis_compared_to_Control_redis_higher.png", width = 10, height = 10, res = 600, units = "in")
print(gg3_1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/CLEC_number_heatmap_Psoriasis_compared_to_Control_redis_higher.pdf", width = 10, height = 10)
print(gg3_1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/CLEC_weight_heatmap_Psoriasis_compared_to_Control_redis_higher.png", width = 10, height = 10, res = 600, units = "in")
print(gg4_1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/CLEC_weight_heatmap_Psoriasis_compared_to_Control_redis_higher.pdf", width = 10, height = 10)
print(gg4_1)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/MHC-I_number_heatmap_Psoriasis_compared_to_Control_redis_higher.png", width = 10, height = 10, res = 600, units = "in")
print(gg3_2)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/MHC-I_number_heatmap_Psoriasis_compared_to_Control_redis_higher.pdf", width = 10, height = 10)
print(gg3_2)
dev.off()
png("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/MHC-I_weight_heatmap_Psoriasis_compared_to_Control_redis_higher.png", width = 10, height = 10, res = 600, units = "in")
print(gg4_2)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/MHC-I_weight_heatmap_Psoriasis_compared_to_Control_redis_higher.pdf", width = 10, height = 10)
print(gg4_2)
dev.off()

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link))
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]],color.use = colors1,title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
p = patchwork::wrap_plots(plots = gg)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_signalingRole.png", width = 20, height = 10,dpi = 600,
       plot = p)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_signalingRole.pdf", width = 20, height = 10, dpi = 300,
       plot = p)

p1 = netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Transitional activated",
                                          color.use = c("grey10", "#FFE4B5", "#800080"),
                                          comparison = c(1, 2))
p2 = netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Cytotoxic effector",
                                          color.use = c("grey10", "#FFE4B5", "#800080"),
                                          comparison = c(1, 2))
p3 = netAnalysis_signalingChanges_scatter(cellchat, idents.use = "CSMD1+",
                                          color.use = c("grey10", "#FFE4B5", "#800080"),
                                          comparison = c(1, 2))
p4 = netAnalysis_signalingChanges_scatter(cellchat, idents.use = "NK-like",
                                          color.use = c("grey10", "#FFE4B5", "#800080"),
                                          comparison = c(1, 2))
p5 = netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Activation-regulated",
                                          color.use = c("grey10", "#FFE4B5", "#800080"),
                                          comparison = c(1, 2))
p6 = netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Naive-like",
                                          color.use = c("grey10", "#FFE4B5", "#800080"),
                                          comparison = c(1, 2))
p7 = netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Cytotoxic memory-like",
                                          color.use = c("grey10", "#FFE4B5", "#800080"),
                                          comparison = c(1, 2))
p8 = netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Terminal-branch",
                                          color.use = c("grey10", "#FFE4B5", "#800080"),
                                          comparison = c(1, 2))
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/from_Control_to_Psoriasis_Transitional activated_signalingChanges.png", width = 6, height = 4.5,dpi = 600,
       plot = p1)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/from_Control_to_Psoriasis_Transitional activated_signalingChanges.pdf", width = 6, height = 4.5, dpi = 300,
       plot = p1)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/from_Control_to_Psoriasis_Cytotoxic effector_signalingChanges.png", width = 6, height = 4.5,dpi = 600,
       plot = p2)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/from_Control_to_Psoriasis_Cytotoxic effector_signalingChanges.pdf", width = 6, height = 4.5, dpi = 300,
       plot = p2)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/from_Control_to_Psoriasis_CSMD1+_signalingChanges.png", width = 6, height = 4.5,dpi = 600,
       plot = p3)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/from_Control_to_Psoriasis_CSMD1+_signalingChanges.pdf", width = 6, height = 4.5, dpi = 300,
       plot = p3)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/from_Control_to_Psoriasis_NK-like_signalingChanges.png", width = 6, height = 4.5,dpi = 600,
       plot = p4)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/from_Control_to_Psoriasis_NK-like_signalingChanges.pdf", width = 6, height = 4.5, dpi = 300,
       plot = p4)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/from_Control_to_Psoriasis_Activation-regulated_signalingChanges.png", width = 6, height = 4.5,dpi = 600,
       plot = p5)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/from_Control_to_Psoriasis_Activation-regulated_signalingChanges.pdf", width = 6, height = 4.5, dpi = 300,
       plot = p5)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/from_Control_to_Psoriasis_Naive-like_signalingChanges.png", width = 6, height = 4.5,dpi = 600,
       plot = p6)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/from_Control_to_Psoriasis_Naive-like_signalingChanges.pdf", width = 6, height = 4.5, dpi = 300,
       plot = p6)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/from_Control_to_Psoriasis_Cytotoxic memory-like_signalingChanges.png", width = 6, height = 4.5,dpi = 600,
       plot = p7)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/from_Control_to_Psoriasis_Cytotoxic memory-like_signalingChanges.pdf", width = 6, height = 4.5, dpi = 300,
       plot = p7)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/from_Control_to_Psoriasis_Terminal-branch_signalingChanges.png", width = 6, height = 4.5,dpi = 600,
       plot = p8)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/from_Control_to_Psoriasis_Terminal-branch_signalingChanges.pdf", width = 6, height = 4.5, dpi = 300,
       plot = p8)
p1 = plotGeneExpression(cellchat, signaling = "MHC-I", split.by = "datasets", colors.ggplot = T,color.use = c("#FFE4B5", "#800080"))
p6 = plotGeneExpression(cellchat, signaling = "CLEC", split.by = "datasets", colors.ggplot = T,color.use = c("#FFE4B5", "#800080"))
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/from_Control_to_Psoriasis_MHC-I_vlnplot.png", width = 10, height = 10,dpi = 600,
       plot = p1)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/from_Control_to_Psoriasis_MHC-I_vlnplot.pdf", width = 10, height = 10, dpi = 300,
       plot = p1)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/from_Control_to_Psoriasis_CLEC_vlnplot.png", width = 10, height = 10,dpi = 600,
       plot = p6)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/from_Control_to_Psoriasis_CLEC_vlnplot.pdf", width = 10, height = 10, dpi = 300,
       plot = p6)

library(uwot)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional",umap.method = "uwot")
cellchat <- netClustering(cellchat, type = "functional")
p = netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5,color.use = mycolor1_clusters)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_embeddingPairwise.png", width = 10, height = 10,dpi = 600,
       plot = p)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_embeddingPairwise.pdf", width = 10, height = 10, dpi = 300,
       plot = p)

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"))
p = gg1 + gg2
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ALL_ControlPsoriasis_rankNet.png", width = 10, height = 10,dpi = 600,
       plot = p)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ALL_ControlPsoriasis_rankNet.pdf", width = 10, height = 10, dpi = 300,
       plot = p)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ALL_ControlPsoriasis_rankNet_1.png", width = 10, height = 10,dpi = 600,
       plot = gg1)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ALL_ControlPsoriasis_rankNet_1.pdf", width = 10, height = 10, dpi = 300,
       plot = gg1)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),signaling = c("MHC-I","CLEC"))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),signaling = c("MHC-I","CLEC"))
p = gg1 + gg2
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/signal_ControlPsoriasis_rankNet.png", width = 10, height = 10,dpi = 600,
       plot = p)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/signal_ControlPsoriasis_rankNet.pdf", width = 10, height = 10, dpi = 300,
       plot = p)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/signal_ControlPsoriasis_rankNet_1.png", width = 10, height = 10,dpi = 600,
       plot = gg1)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/signal_ControlPsoriasis_rankNet_1.pdf", width = 10, height = 10, dpi = 300,
       plot = gg1)
gg1_1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),sources.use = "Transitional activated")
gg2_1 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),sources.use = "Transitional activated")
p1 = gg1_1 + gg2_1
gg1_2 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),sources.use = "Cytotoxic effector")
gg2_2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),sources.use = "Cytotoxic effector")
p2 = gg1_2 + gg2_2
gg1_3 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),sources.use = "CSMD1+")
gg2_3 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),sources.use = "CSMD1+")
p3 = gg1_3 + gg2_3
gg1_4 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),sources.use = "NK-like")
gg2_4 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),sources.use = "NK-like")
p4 = gg1_4 + gg2_4
gg1_5 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),sources.use = "Activation-regulated")
gg2_5 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),sources.use = "Activation-regulated")
p5 = gg1_5 + gg2_5
gg1_6 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),sources.use = "Naive-like")
gg2_6 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),sources.use = "Naive-like")
p6 = gg1_6 + gg2_6
gg1_7 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),sources.use = "Cytotoxic memory-like")
gg2_7 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),sources.use = "Cytotoxic memory-like")
p7 = gg1_7 + gg2_7
gg1_8 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),sources.use = "Terminal-branch")
gg2_8 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),sources.use = "Terminal-branch")
p8 = gg1_8 + gg2_8
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_Transitional activated_ControlPsoriasis_rankNet.png", width = 10, height = 10,dpi = 600,
       plot = p1)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_Transitional activated_ControlPsoriasis_rankNet.pdf", width = 10, height = 10, dpi = 300,
       plot = p1)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_Transitional activated_ControlPsoriasis_rankNet_1.png", width = 10, height = 10,dpi = 600,
       plot = gg1_1)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_Transitional activated_ControlPsoriasis_rankNet_1.pdf", width = 10, height = 10, dpi = 300,
       plot = gg1_1)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_Cytotoxic effector_ControlPsoriasis_rankNet.png", width = 10, height = 10,dpi = 600,
       plot = p2)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_Cytotoxic effector_ControlPsoriasis_rankNet.pdf", width = 10, height = 10, dpi = 300,
       plot = p2)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_Cytotoxic effector_ControlPsoriasis_rankNet_1.png", width = 10, height = 10,dpi = 600,
       plot = gg1_2)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_Cytotoxic effector_ControlPsoriasis_rankNet_1.pdf", width = 10, height = 10, dpi = 300,
       plot = gg1_2)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_CSMD1+_ControlPsoriasis_rankNet.png", width = 10, height = 10,dpi = 600,
       plot = p3)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_CSMD1+_rankNet.pdf", width = 10, height = 10, dpi = 300,
       plot = p3)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_CSMD1+_rankNet_1.png", width = 10, height = 10,dpi = 600,
       plot = gg1_3)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_CSMD1+_rankNet_1.pdf", width = 10, height = 10, dpi = 300,
       plot = gg1_3)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_NK-like_rankNet.png", width = 10, height = 10,dpi = 600,
       plot = p4)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_NK-like_ControlPsoriasis_rankNet.pdf", width = 10, height = 10, dpi = 300,
       plot = p4)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_NK-like_ControlPsoriasis_rankNet_1.png", width = 10, height = 10,dpi = 600,
       plot = gg1_4)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_NK-like_ControlPsoriasis_rankNet_1.pdf", width = 10, height = 10, dpi = 300,
       plot = gg1_4)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_Activation-regulated_ControlPsoriasis_rankNet.png", width = 10, height = 10,dpi = 600,
       plot = p5)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_Activation-regulated_ControlPsoriasis_rankNet.pdf", width = 10, height = 10, dpi = 300,
       plot = p5)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_Activation-regulated_ControlPsoriasis_rankNet_1.png", width = 10, height = 10,dpi = 600,
       plot = gg1_5)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_Activation-regulated_ControlPsoriasis_rankNet_1.pdf", width = 10, height = 10, dpi = 300,
       plot = gg1_5)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_Naive-like_ControlPsoriasis_rankNet.png", width = 10, height = 10,dpi = 600,
       plot = p6)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_Naive-like_ControlPsoriasis_rankNet.pdf", width = 10, height = 10, dpi = 300,
       plot = p6)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_Naive-like_ControlPsoriasis_rankNet_1.png", width = 10, height = 10,dpi = 600,
       plot = gg1_6)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_Naive-like_ControlPsoriasis_rankNet_1.pdf", width = 10, height = 10, dpi = 300,
       plot = gg1_6)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_Cytotoxic memory-like_ControlPsoriasis_rankNet.png", width = 10, height = 10,dpi = 600,
       plot = p7)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_Cytotoxic memory-like_ControlPsoriasis_rankNet.pdf", width = 10, height = 10, dpi = 300,
       plot = p7)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_Cytotoxic memory-like_ControlPsoriasis_rankNet_1.png", width = 10, height = 10,dpi = 600,
       plot = gg1_7)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_Cytotoxic memory-like_ControlPsoriasis_rankNet_1.pdf", width = 10, height = 10, dpi = 300,
       plot = gg1_7)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_Terminal-branch_ControlPsoriasis_rankNet.png", width = 10, height = 10,dpi = 600,
       plot = p8)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_Terminal-branch_ControlPsoriasis_rankNet.pdf", width = 10, height = 10, dpi = 300,
       plot = p8)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_Terminal-branch_ControlPsoriasis_rankNet_1.png", width = 10, height = 10,dpi = 600,
       plot = gg1_8)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/source_Terminal-branch_ControlPsoriasis_rankNet_1.pdf", width = 10, height = 10, dpi = 300,
       plot = gg1_8)

gg1_1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),targets.use = "Transitional activated")
gg2_1 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),targets.use = "Transitional activated")
p1 = gg1_1 + gg2_1
gg1_2 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),targets.use = "Cytotoxic effector")
gg2_2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),targets.use = "Cytotoxic effector")
p2 = gg1_2 + gg2_2
gg1_3 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),targets.use = "CSMD1+")
gg2_3 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),targets.use = "CSMD1+")
p3 = gg1_3 + gg2_3
gg1_4 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),targets.use = "NK-like")
gg2_4 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),targets.use = "NK-like")
p4 = gg1_4 + gg2_4
gg1_5 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),targets.use = "Activation-regulated")
gg2_5 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),targets.use = "Activation-regulated")
p5 = gg1_5 + gg2_5
gg1_6 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),targets.use = "Naive-like")
gg2_6 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),targets.use = "Naive-like")
p6 = gg1_6 + gg2_6
gg1_7 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),targets.use = "Cytotoxic memory-like")
gg2_7 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),targets.use = "Cytotoxic memory-like")
p7 = gg1_7 + gg2_7
gg1_8 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),targets.use = "Terminal-branch")
gg2_8 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use = c("#FFE4B5", "#800080"),targets.use = "Terminal-branch")
p8 = gg1_8 + gg2_8
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_Transitional activated_ControlPsoriasis_rankNet.png", width = 10, height = 10,dpi = 600,
       plot = p1)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_Transitional activated_ControlPsoriasis_rankNet.pdf", width = 10, height = 10, dpi = 300,
       plot = p1)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_Transitional activated_ControlPsoriasis_rankNet_1.png", width = 10, height = 10,dpi = 600,
       plot = gg1_1)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_Transitional activated_ControlPsoriasis_rankNet_1.pdf", width = 10, height = 10, dpi = 300,
       plot = gg1_1)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_Cytotoxic effector_ControlPsoriasis_rankNet.png", width = 10, height = 10,dpi = 600,
       plot = p2)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_Cytotoxic effector_ControlPsoriasis_rankNet.pdf", width = 10, height = 10, dpi = 300,
       plot = p2)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_Cytotoxic effector_ControlPsoriasis_rankNet_1.png", width = 10, height = 10,dpi = 600,
       plot = gg1_2)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_Cytotoxic effector_ControlPsoriasis_rankNet_1.pdf", width = 10, height = 10, dpi = 300,
       plot = gg1_2)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_CSMD1+_ControlPsoriasis_rankNet.png", width = 10, height = 10,dpi = 600,
       plot = p3)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_CSMD1+_rankNet.pdf", width = 10, height = 10, dpi = 300,
       plot = p3)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_CSMD1+_rankNet_1.png", width = 10, height = 10,dpi = 600,
       plot = gg1_3)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_CSMD1+_rankNet_1.pdf", width = 10, height = 10, dpi = 300,
       plot = gg1_3)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_NK-like_rankNet.png", width = 10, height = 10,dpi = 600,
       plot = p4)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_NK-like_ControlPsoriasis_rankNet.pdf", width = 10, height = 10, dpi = 300,
       plot = p4)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_NK-like_ControlPsoriasis_rankNet_1.png", width = 10, height = 10,dpi = 600,
       plot = gg1_4)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_NK-like_ControlPsoriasis_rankNet_1.pdf", width = 10, height = 10, dpi = 300,
       plot = gg1_4)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_Activation-regulated_ControlPsoriasis_rankNet.png", width = 10, height = 10,dpi = 600,
       plot = p5)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_Activation-regulated_ControlPsoriasis_rankNet.pdf", width = 10, height = 10, dpi = 300,
       plot = p5)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_Activation-regulated_ControlPsoriasis_rankNet_1.png", width = 10, height = 10,dpi = 600,
       plot = gg1_5)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_Activation-regulated_ControlPsoriasis_rankNet_1.pdf", width = 10, height = 10, dpi = 300,
       plot = gg1_5)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_Naive-like_ControlPsoriasis_rankNet.png", width = 10, height = 10,dpi = 600,
       plot = p6)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_Naive-like_ControlPsoriasis_rankNet.pdf", width = 10, height = 10, dpi = 300,
       plot = p6)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_Naive-like_ControlPsoriasis_rankNet_1.png", width = 10, height = 10,dpi = 600,
       plot = gg1_6)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_Naive-like_ControlPsoriasis_rankNet_1.pdf", width = 10, height = 10, dpi = 300,
       plot = gg1_6)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_Cytotoxic memory-like_ControlPsoriasis_rankNet.png", width = 10, height = 10,dpi = 600,
       plot = p7)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_Cytotoxic memory-like_ControlPsoriasis_rankNet.pdf", width = 10, height = 10, dpi = 300,
       plot = p7)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_Cytotoxic memory-like_ControlPsoriasis_rankNet_1.png", width = 10, height = 10,dpi = 600,
       plot = gg1_7)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_Cytotoxic memory-like_ControlPsoriasis_rankNet_1.pdf", width = 10, height = 10, dpi = 300,
       plot = gg1_7)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_Terminal-branch_ControlPsoriasis_rankNet.png", width = 10, height = 10,dpi = 600,
       plot = p8)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_Terminal-branch_ControlPsoriasis_rankNet.pdf", width = 10, height = 10, dpi = 300,
       plot = p8)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_Terminal-branch_ControlPsoriasis_rankNet_1.png", width = 10, height = 10,dpi = 600,
       plot = gg1_8)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/target_Terminal-branch_ControlPsoriasis_rankNet_1.pdf", width = 10, height = 10, dpi = 300,
       plot = gg1_8)

library(ComplexHeatmap)
i = 1

pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
object.list[[i+1]] <- netAnalysis_computeCentrality(object.list[[i+1]])
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 20,color.use = colors1)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 20,color.use = colors1)
p1 = draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
object.list[[i+1]] <- netAnalysis_computeCentrality(object.list[[i+1]])
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 20,color.use = colors1)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 20,color.use = colors1)
p2 = draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

png("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_outgoing.png", width = 10, height = 10, res = 600, units = "in")
print(p1)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_outgoing.pdf", width = 10, height = 10)
print(p1)
dev.off()

png("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_incoming.png", width = 10, height = 10, res = 600, units = "in")
print(p2)
dev.off()
pdf("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_incoming.pdf", width = 10, height = 10)
print(p2)
dev.off()

p = netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45,remove.isolate = F,color.text = c("#FFE4B5", "#800080"))
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_all_sender_receiver.png", width = 20, height = 20,dpi = 600,
       plot = p)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_all_sender_receiver.pdf", width = 20, height = 45, dpi = 300,
       plot = p)
p1 = netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45,sources.use = "Transitional activated",color.text = c("#FFE4B5", "#800080"))
p2 = netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45,targets.use = "Transitional activated",color.text = c("#FFE4B5", "#800080"))
p3 = netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45,sources.use = "Cytotoxic effector",color.text = c("#FFE4B5", "#800080"))
p4 = netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45,targets.use = "Cytotoxic effector",color.text = c("#FFE4B5", "#800080"))
p5 = netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45,sources.use = "CSMD1+",color.text = c("#FFE4B5", "#800080"))
p6 = netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45,targets.use = "CSMD1+",color.text = c("#FFE4B5", "#800080"))
p7 = netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45,sources.use = "NK-like",color.text = c("#FFE4B5", "#800080"))
p8 = netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45,targets.use = "NK-like",color.text = c("#FFE4B5", "#800080"))
p9 = netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45,sources.use = "Activation-regulated",color.text = c("#FFE4B5", "#800080"))
p10 = netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45,targets.use = "Activation-regulated",color.text = c("#FFE4B5", "#800080"))
p11 = netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45,sources.use = "Naive-like",color.text = c("#FFE4B5", "#800080"))
p12 = netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45,targets.use = "Naive-like",color.text = c("#FFE4B5", "#800080"))
p13 = netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45,sources.use = "Cytotoxic memory-like",color.text = c("#FFE4B5", "#800080"))
p14 = netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45,targets.use = "Cytotoxic memory-like",color.text = c("#FFE4B5", "#800080"))
p15 = netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45,sources.use = "Terminal-branch",color.text = c("#FFE4B5", "#800080"))
p16 = netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45,targets.use = "Terminal-branch",color.text = c("#FFE4B5", "#800080"))

ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_Transitional activated_sender.png", width = 10, height = 10,dpi = 600,
       plot = p1)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_Transitional activated_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p1)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_Transitional activated_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p2)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_Transitional activated_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p2)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_Cytotoxic effector_sender.png", width = 10, height = 10,dpi = 600,
       plot = p3)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_Cytotoxic effector_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p3)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_Cytotoxic effector_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p4)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_Cytotoxic effector_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p4)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_CSMD1+_sender.png", width = 10, height = 10,dpi = 600,
       plot = p5)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_CSMD1+_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p5)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_CSMD1+_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p6)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_CSMD1+_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p6)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_NK-like_sender.png", width = 10, height = 10,dpi = 600,
       plot = p7)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_NK-like_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p7)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_NK-like_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p8)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_NK-like_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p8)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_Activation-regulated_sender.png", width = 10, height = 10,dpi = 600,
       plot = p9)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_Activation-regulated_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p9)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_Activation-regulated_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p10)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_Activation-regulated_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p10)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_Naive-like_sender.png", width = 10, height = 10,dpi = 600,
       plot = p11)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_Naive-like_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p11)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_Naive-like_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p12)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_Naive-like_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p12)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_Cytotoxic memory-like_sender.png", width = 10, height = 10,dpi = 600,
       plot = p13)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_Cytotoxic memory-like_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p13)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_Cytotoxic memory-like_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p14)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_Cytotoxic memory-like_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p14)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_Terminal-branch_sender.png", width = 10, height = 10,dpi = 600,
       plot = p15)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_Terminal-branch_sender.pdf", width = 10, height = 10, dpi = 300,
       plot = p15)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_Terminal-branch_receiver.png", width = 10, height = 10,dpi = 600,
       plot = p16)
ggsave("/thinker/3.tangjiale/shinanxi/CellchatKS/compare_plot/ControlPsoriasis_Terminal-branch_receiver.pdf", width = 10, height = 10, dpi = 300,
       plot = p16)
```






##########12.CellChat analysis of CSMD1 subgroups using public databases################
#########(1)Setting up the environment################
setwd("/thinker/5.chensihan/shinanxi_cellchat/26.1.24/")
.libPaths(c("/home/02chensihan/R/x86_64-pc-linux-gnu-library/4.3" , .libPaths()))
.libPaths(c("/opt/R/4.3.3/lib/R/library" , .libPaths()))
.libPaths(c("/thinker/5.chensihan/packages_new/", .libPaths()))
.libPaths(c("/thinker/5.chensihan/SeuratV4/", .libPaths()))
library(Seurat)
.libPaths(c("/thinker/5.chensihan/packages_new/", .libPaths()))
library(usethis)
library(devtools)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(batchelor)
#remotes::install_github('satijalab/seurat-wrappers@community-vignette')
library(SeuratWrappers)
library(magrittr)
library(tidyverse)
library(clusterProfiler)
library(GO.db)
library(org.Hs.eg.db)
library(DOSE)
library(DoubletFinder)
library(SingleR)
library(celldex)
# BiocManager::install("DOSE")

library(Rcpp)
library(harmony)
library(pheatmap)
#install.packages("scCustomize")
library(scCustomize)
library(RColorBrewer)
.libPaths(c("/thinker/5.chensihan/SeuratV4/", .libPaths()))
library(Seurat)
.libPaths(c("/thinker/5.chensihan/packages_new/", .libPaths()))
# .libPaths()
library(igraph)
# remotes::install_github("sqjin/CellChat", dependencies = TRUE)
library(CellChat)

#rm(list = ls())
#pacman::p_unload(pacman::p_loaded(), character.only = TRUE)
#suppressWarnings(library(Seurat))

seu <- readRDS("../rds/Overall_shinanxi_New.rds")

##################(2)CSMD1 subgroups and PBMCs from public databases###########################
setwd("/thinker/5.chensihan/shinanxi_cellchat/26.1.24/")

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(ggplot2)
  library(patchwork)
  library(future)
  library(dplyr)
  library(harmony)
  library(NMF)
})

future::plan("multisession", workers = 4) 
options(stringsAsFactors = FALSE)

print("Loading user data...")
seu_user <- readRDS("../rds/Overall_shinanxi_New.rds")

seu_user$Dataset <- "My_PBMC"
seu_user$final_annotation <- as.character(seu_user$celltype)

print(table(seu_user$final_annotation))

sample_folders <- list.dirs(path = "./public_data/PBMC/", full.names = FALSE, recursive = FALSE)
sample_folders <- grep("^PBMC", sample_folders, value = TRUE)

seurat_list <- list()

for (sample in sample_folders) {
  message("Processing: ", sample)
  
  matrix_path   <- file.path("./public_data/PBMC",sample, paste0(sample, ".matrix.mtx.gz"))
  barcodes_path <- file.path("./public_data/PBMC",sample, paste0(sample, ".barcodes.tsv.gz"))
  features_path <- file.path("./public_data/PBMC",sample, paste0(sample, ".features.tsv.gz"))
  
  if (file.exists(matrix_path) & file.exists(barcodes_path) & file.exists(features_path)) {
    
    counts <- ReadMtx(
      mtx = matrix_path,
      cells = barcodes_path,
      features = features_path
    )
    
    seurat_obj <- CreateSeuratObject(
      counts = counts, 
      project = sample, 
      min.cells = 3, 
      min.features = 200
    )
    
    seurat_list[[sample]] <- seurat_obj
    
  } else {
    warning(paste("Skipping sample", sample, ": file path incomplete"))
  }
}

if (length(seurat_list) > 0) {
  seu_PBMC <- merge(
    x = seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = names(seurat_list),
    project = "PBMC_Merged"
  )
  
  print(seu_PBMC)
} else {
  stop("No data successfully read, please check path.")
}

meta_PBMC <- read.table(
  gzfile("/thinker/5.chensihan/shinanxi_cellchat/26.1.24/public_data/PBMC/GSE194315_CellMetadata-PSA_TotalCiteseq_20220103.tsv.gz"),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)
Idents(seu_PBMC) <- "meta.data"
head(meta_PBMC)
rownames(meta_PBMC) <- meta_PBMC$CellName
rownames(meta_PBMC) <- paste0(rownames(meta_PBMC), "-1")
seu_PBMC <- AddMetaData(seu_PBMC, metadata = meta_PBMC)
seu_PBMC <- seu_PBMC[, !is.na(seu_PBMC$CellType)]
seu_PBMC$celltype <- seu_PBMC$CellType

seu_PBMC$group <- NA

seu_PBMC$group[seu_PBMC$Status == "Healthy"] <- "Control"

seu_PBMC$group[seu_PBMC$Status %in% c("PSO")] <- "Psoriasis"

seu_PBMC <- seu_PBMC[,!is.na(seu_PBMC$group)]

print("Current Status column categories:")
table(seu_PBMC$Status)

group_disease <- "PSO"
group_control <- "Healthy"

samples_disease_all <- unique(seu_PBMC@meta.data$Sample[seu_PBMC$Status == group_disease])
samples_control_all <- unique(seu_PBMC@meta.data$Sample[seu_PBMC$Status == group_control])

print(paste(group_disease, "total samples:", length(samples_disease_all)))
print(paste(group_control, "total samples:", length(samples_control_all)))

set.seed(123) 

target_disease <- if(length(samples_disease_all) >= 3) sample(samples_disease_all, 3) else samples_disease_all
target_control <- if(length(samples_control_all) >= 3) sample(samples_control_all, 3) else samples_control_all

print("Selected disease samples:")
print(target_disease)
print("Selected healthy samples:")
print(target_control)

target_samples <- c(target_disease, target_control)

seu_sub <- subset(seu_PBMC, subset = Sample %in% target_samples)

print("Subsampled object information:")
print(seu_sub)
table(seu_sub$Status, seu_sub$Sample)

seu_PBMC <- seu_sub

print("Converted group statistics:")
table(seu_PBMC$group)

print("Correspondence between original Status and new group:")
table(seu_PBMC$Status, seu_PBMC$group)
saveRDS(seu_PBMC,"./seu_PBMC.rds")

seu_combined <- merge(seu_user, y = seu_PBMC, add.cell.ids = c("gdT", "PBMC"))
table(seu_combined$celltype, useNA = "ifany")
saveRDS(seu_combined,"./result_combined_PBMC.rds")
library(loupeR)
loupeR::setup()
create_loupe_from_seurat(seu_PBMC)

setwd("/thinker/5.chensihan/shinanxi_cellchat/26.1.24/PBMC_2.14/")
save_plot_auto_2.5 <- function(plot_code_or_obj, filename_base, width = 18, height = 9.5, type = "both") {
  if(!dir.exists("cellchat")) dir.create("cellchat", recursive = TRUE)
  f_png <- paste0("cellchat/", filename_base, ".png")
  f_pdf <- paste0("cellchat/", filename_base, ".pdf")
  
  draw_it <- function(obj) {
    if (inherits(obj, "grob")) {
      grid::grid.draw(obj)
    } else if (is.expression(obj)) {
      eval(obj)
    } else {
      print(obj)
    }
  }
  
  if(type %in% c("both", "png")) {
    png(f_png, width = width, height = height, units = "in", res = 600)
    grid::grid.newpage() 
    tryCatch({
      draw_it(plot_code_or_obj)
    }, error = function(e) message(paste("Error saving PNG:", filename_base, e$message)))
    dev.off()
  }
  
  if(type %in% c("both", "pdf")) {
    pdf(f_pdf, width = width, height = height)
    grid::grid.newpage()
    tryCatch({
      draw_it(plot_code_or_obj)
    }, error = function(e) message(paste("Error saving PDF:", filename_base, e$message)))
    dev.off()
  }
  
  message(paste("Saved successfully:", filename_base))
}

group.cellType <- union(levels(cellchat_list[[1]]@idents), levels(cellchat_list[[2]]@idents))
for (i in 1:length(cellchat_list)) {
  cellchat_list[[i]] <- liftCellChat(cellchat_list[[i]], group.new = group.cellType)
}
cellchat_merged <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))
compare_vec <- c(2, 1)


setwd("/thinker/5.chensihan/shinanxi_cellchat/26.1.24/")
.libPaths(c("/home/02chensihan/R/x86_64-pc-linux-gnu-library/4.3" , .libPaths()))
.libPaths(c("/opt/R/4.3.3/lib/R/library" , .libPaths()))
.libPaths(c("/thinker/5.chensihan/packages_new/", .libPaths()))
.libPaths(c("/thinker/5.chensihan/SeuratV4/", .libPaths()))
library(Seurat)
.libPaths(c("/thinker/5.chensihan/packages_new/", .libPaths()))
library(usethis)
library(devtools)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(batchelor)
library(SeuratWrappers)
library(magrittr)
library(tidyverse)
library(clusterProfiler)
library(GO.db)
library(org.Hs.eg.db)
library(DOSE)
library(DoubletFinder)
library(SingleR)
library(celldex)
library(ggalluvial)
library(Rcpp)
library(harmony)
library(pheatmap)
library(scCustomize)
library(RColorBrewer)
library(circlize)
.libPaths(c("/thinker/5.chensihan/SeuratV4/", .libPaths()))
library(Seurat)
.libPaths(c("/thinker/5.chensihan/packages_new/", .libPaths()))
library(igraph)
library(CellChat)
suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(ggplot2)
  library(patchwork)
  library(future)
  library(dplyr)
  library(harmony)
  library(NMF)
  library(CellChat)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(dplyr)
  library(pheatmap)
  library(grid)
})
mycolor2 <- c(
  "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F","#FF6A6A","#7FFFD4", 
  "#AB82FF","#90EE90","#00CD00","#008B8B","#6495ED","#FFC1C1","#CD5C5C","#8B008B",
  "#FF3030", "#7CFC00","#000000","#708090"
)

cell_types <- levels(cellchat_merged@meta$celltype)

print(paste("Number of cell types:", length(cell_types)))
print(paste("Number of provided colors:", length(mycolor2)))

all_pathways_1 <- cellchat_list[[1]]@netP$pathways
keep_pathways_1 <- all_pathways_1[!grepl("HLA|MHC", all_pathways_1, ignore.case = TRUE)]

all_pathways_2 <- cellchat_list[[2]]@netP$pathways
keep_pathways_2 <- all_pathways_2[!grepl("HLA|MHC", all_pathways_2, ignore.case = TRUE)]

print(head(cellchat_list[[1]]@netP$pathways))
print(head(cellchat_list[[2]]@netP$pathways))

p1 <- cellchat_list[[1]]@netP$pathways
keep_pathways_1 <- p1[!grepl("HLA|MHC", p1, ignore.case = TRUE)]

p2 <- cellchat_list[[2]]@netP$pathways
keep_pathways_2 <- p2[!grepl("HLA|MHC", p2, ignore.case = TRUE)]

mat_control <- cellchat_list[[1]]@net$count
mat_disease <- cellchat_list[[2]]@net$count

library(pheatmap)
library(grid)
library(gridExtra)

global_max <- max(max(mat_control, na.rm = T), max(mat_disease, na.rm = T))
common_breaks <- seq(0, global_max, length.out = 100)
common_color <- colorRampPalette(c("white", "#E41A1C"))(100)

if(!dir.exists("cellchat")) dir.create("cellchat")

p3_ctrl <- pheatmap::pheatmap(cellchat_list[[1]]@net$count, border_color = "black",
                              cluster_cols = F, cluster_rows = F, fontsize = 10,
                              display_numbers = T, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Control", silent = TRUE)$gtable

p3_case <- pheatmap::pheatmap(mat_disease, border_color = "black",
                              cluster_cols = F, cluster_rows = F, fontsize = 10,
                              display_numbers = T, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks, 
                              main = "Psoriasis", silent = TRUE)$gtable

p3_combined <- arrangeGrob(p3_ctrl, p3_case, ncol = 2, top = textGrob("Interaction Count (Raw)", gp=gpar(fontsize=15, font=2)))
save_plot_auto_2.5(p3_combined, "CellChat_Heatmap_NoCluster_WithNum")

p4_ctrl <- pheatmap::pheatmap(mat_control, border_color = "black",
                              cluster_cols = T, cluster_rows = T, fontsize = 10,
                              display_numbers = T, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Control", silent = TRUE)$gtable

p4_case <- pheatmap::pheatmap(mat_disease, border_color = "black",
                              cluster_cols = T, cluster_rows = T, fontsize = 10,
                              display_numbers = T, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Psoriasis", silent = TRUE)$gtable

p4_combined <- arrangeGrob(p4_ctrl, p4_case, ncol = 2, top = textGrob("Interaction Count (Clustered)", gp=gpar(fontsize=15, font=2)))
save_plot_auto_2.5(p4_combined, "CellChat_Heatmap_Clustered_WithNum")

p5_ctrl <- pheatmap::pheatmap(mat_control, border_color = "black",
                              cluster_cols = F, cluster_rows = F, fontsize = 10,
                              display_numbers = F, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Control", silent = TRUE)$gtable

p5_case <- pheatmap::pheatmap(mat_disease, border_color = "black",
                              cluster_cols = F, cluster_rows = F, fontsize = 10,
                              display_numbers = F, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Psoriasis", silent = TRUE)$gtable

p5_combined <- arrangeGrob(p5_ctrl, p5_case, ncol = 2, top = textGrob("Interaction Count (Raw, No Num)", gp=gpar(fontsize=15, font=2)))
save_plot_auto_2.5(p5_combined, "CellChat_Heatmap_NoCluster_NoNum")

p6_ctrl <- pheatmap::pheatmap(mat_control, border_color = "black",
                              cluster_cols = T, cluster_rows = T, fontsize = 10,
                              display_numbers = F, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Control", silent = TRUE)$gtable

p6_case <- pheatmap::pheatmap(mat_disease, border_color = "black",
                              cluster_cols = T, cluster_rows = T, fontsize = 10,
                              display_numbers = F, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Psoriasis", silent = TRUE)$gtable

p6_combined <- arrangeGrob(p6_ctrl, p6_case, ncol = 2, top = textGrob("Interaction Count (Clustered, No Num)", gp=gpar(fontsize=15, font=2)))
save_plot_auto_2.5(p6_combined, "CellChat_Heatmap_Clustered_NoNum")



mat <- cellchat_merged@net$count
cell_types <- levels(cellchat_merged@idents$joint)
groupSize <- as.numeric(table(cellchat_merged@idents$joint))

code_single_circle <- expression({
  par(mfrow = c(3,3), xpd=TRUE, mar = c(1,1,1,1))
  for (i in 1:length(cell_types)) {
    netVisual_circle(cellchat_list[["Psoriasis"]]@net$count, vertex.weight = groupSize, weight.scale = T, 
                     edge.weight.max = max(cellchat_list[["Psoriasis"]]@net$count), 
                     title.name = paste("Pso -", cell_types[i]), color.use = mycolor2, sources.use = cell_types[i])
  }
})
save_plot_auto_2.5(code_single_circle, "03_Single_CellType_Source_Circle_Psoriasis")
code_single_circle <- expression({
  par(mfrow = c(3,3), xpd=TRUE, mar = c(1,1,1,1))
  for (i in 1:length(cell_types)) {
    netVisual_circle(cellchat_list[["Control"]]@net$count, vertex.weight = groupSize, weight.scale = T, 
                     edge.weight.max = max(cellchat_list[["Control"]]@net$count), 
                     title.name = paste("Pso -", cell_types[i]), color.use = mycolor2, sources.use = cell_types[i])
  }
})
save_plot_auto_2.5(code_single_circle, "03_Single_CellType_Source_Circle_Control")

groupSize <- as.numeric(table(cellchat_list[[1]]@idents))

p1 <- expression({
  par(mfrow = c(1,2), xpd=TRUE)
  
  netVisual_circle(cellchat_list[[1]]@net$count, vertex.weight = groupSize, weight.scale = T,
                   label.edge= F,                    arrow.width = 1.2,
                   arrow.size = 0.24,color.use = mycolor2,
                   title.name = "Number of interactions in Cont")
  netVisual_circle(cellchat_list[[1]]@net$weight, vertex.weight = groupSize, weight.scale = T, 
                   label.edge= F, arrow.width = 1.2, color.use = mycolor2,
                   arrow.size = 0.24 ,title.name = "Interaction weights/strength in Cont")
  
  groupSize <- as.numeric(table(cellchat_list[[2]]@idents))})
save_plot_auto_2.5(p1, "03_Total_CellType_Source_Circle_Control")

p2 <- expression({
  par(mfrow = c(1,2), xpd=TRUE)
  
  netVisual_circle(cellchat_list[[2]]@net$count, vertex.weight = groupSize, weight.scale = T,
                   label.edge= F, 
                   arrow.width = 1.2, color.use = mycolor2,
                   arrow.size = 0.24,
                   title.name = "Number of interactions in Pso")
  netVisual_circle(cellchat_list[[2]]@net$weight, vertex.weight = groupSize, weight.scale = T, 
                   label.edge= F, 
                   arrow.width = 1.2, color.use = mycolor2,
                   arrow.size = 0.24,
                   title.name = "Interaction weights/strength in Pso")})
save_plot_auto_2.5(p2, "03_Total_CellType_Source_Circle_Psoriasis")

cellchat_list[[1]]@netP$pathways

pp1 <- rankNet(cellchat_merged, mode = "comparison", stacked = T, do.stat = TRUE)
pp2 <- rankNet(cellchat_merged, mode = "comparison", stacked = F, do.stat = TRUE)
save_plot_auto_2.5(pp1 + pp2, "03_rankNet")

target_groups   <- c("Control", "Psoriasis")
chat_ctrl <- cellchat_list[[1]] 
chat_pso  <- cellchat_list[[2]]

pathway_name <- "CLEC"
obj_use <- chat_ctrl
prefix <- "Control_CLEC"

p1 = netVisual_aggregate(obj_use, signaling = pathway_name, color.use = mycolor2)
save_plot_auto_2.5(p1, paste0(prefix, "_circle"))

p3 = netVisual_heatmap(obj_use, signaling = pathway_name, color.heatmap = "Reds", color.use = mycolor2_named)
save_plot_auto_2.5(p3, paste0(prefix, "_heatmap"))

p4_expr = expression({
  netVisual_chord_cell(obj_use, signaling = pathway_name, group = group.cellType, 
                       title.name = paste0(pathway_name, " (Control) - Immune vs Other"), color.use = mycolor2_named)
})
save_plot_auto_2.5(p4_expr, paste0(prefix, "_grouped_chord"))

p5 = netAnalysis_contribution(obj_use, signaling = pathway_name)
save_plot_auto_2.5(p5, paste0(prefix, "_contribution"))

pairLR_use <- extractEnrichedLR(obj_use, signaling = pathway_name, geneLR.return = FALSE)[1]
p6_expr = expression({
  netVisual_individual(obj_use, signaling = pathway_name, pairLR.use = pairLR_use, layout = "circle", color.use = mycolor2)
})
save_plot_auto_2.5(p6_expr, paste0(prefix, "_", pairLR_use, "_interaction"))

state_list <- list(
  "Control"   = cellchat_list[[1]],
  "Psoriasis" = cellchat_list[[2]]
)

target_pathways <- c("MIF")

for (current_state in names(state_list)) {
  
  obj_use <- state_list[[current_state]]
  
  for (pathway_name in target_pathways) {
    
    if (!(pathway_name %in% obj_use@netP$pathways)) {
      message(paste("Skip:", pathway_name, "not found in", current_state))
      next
    }
    
    prefix <- paste0(current_state, "_", pathway_name)
    
    message(paste("Processing:", prefix))
    
    p1 = netVisual_aggregate(obj_use, signaling = pathway_name, color.use = mycolor2)
    save_plot_auto_2.5(p1, paste0(prefix, "_circle"))
    
    p3 = netVisual_heatmap(obj_use, signaling = pathway_name, color.heatmap = "Reds", color.use = mycolor2_named)
    save_plot_auto_2.5(p3, paste0(prefix, "_heatmap"))
    
    group.cellType <- c(
      "Transitional activated"    = "gdT",
      "Cytotoxic effector"  = "gdT",
      "CSMD1+"   = "gdT",
      "NK-like"   = "gdT",
      "Activation-regulated"         = "gdT",
      "Naive-like"   = "gdT",
      "Cytotoxic memory-like" = "gdT",
      "Terminal-branch"         = "gdT",
      "gdT"               = "gdT",
      
      "KC-S.Corneum"             = "Other",
      "KC-S.Granulosum"             = "Other",
      "KC-S.Spinosum"           = "Other",
      "KC-S.Basale"           = "Other",
      "KC-Wound/Prolif"              = "Other",
      "CD4_T_cell"          = "Other",
      "Treg"     = "Other",
      "Cytotoxic_NK_T"              = "Other",
      "Myeloid"              = "Other",
      "Stroma "              = "Other"
    )
    
    color_grouped <- c("gdT" = "#FF34B3", "Other" = "grey80")
    p4_expr = expression({
      netVisual_chord_cell(obj_use, signaling = pathway_name, group = group.cellType,
                           title.name = paste0(pathway_name, " (", current_state, ") - gdT vs Other"), color.use = mycolor2_named)
    })
    save_plot_auto_2.5(p4_expr, paste0(prefix, "_grouped_chord"))
    
    p5 = netAnalysis_contribution(obj_use, signaling = pathway_name)
    save_plot_auto_2.5(p5, paste0(prefix, "_contribution"))
    
    tryCatch({
      pairLR_use <- extractEnrichedLR(obj_use, signaling = pathway_name, geneLR.return = FALSE)[1]
      
      if (!is.na(pairLR_use)) {
        p6_expr = expression({
          netVisual_individual(obj_use, signaling = pathway_name, pairLR.use = pairLR_use, layout = "circle", color.use = mycolor2)
        })
        save_plot_auto_2.5(p6_expr, paste0(prefix, "_", pairLR_use, "_interaction"))
      }
    }, error = function(e) { message(paste("Error plotting individual LR for", prefix)) })
    
    p7 = netVisual_bubble(obj_use, remove.isolate = FALSE,signaling = pathway_name)
    save_plot_auto_2.5(p7, paste0(prefix, "_sender_receiver"))
    
    p8 <- expression(netAnalysis_signalingRole_network(obj_use, signaling = pathway_name, width = 8, height = 2.5, font.size = 10,color.use = mycolor2_named))
    save_plot_auto_2.5(p8, paste0(prefix, "_charactor"))
    
    p9 <- netAnalysis_signalingRole_scatter(obj_use, signaling = pathway_name,color.use = mycolor2)
    save_plot_auto_2.5(p9, paste0(prefix, "_signalingRole_scatte"))
  }
}

CSMD1_chord <- expression(netVisual_chord_gene(cellchat_list[[1]], sources.use = "CSMD1+", lab.cex = 0.5,legend.pos.y =50))
save_plot_auto_2.5(CSMD1_chord,"CSMD1_chord")



for (i in 1:length(cellchat_list)) {
  cellchat_list[[i]] <- netAnalysis_computeCentrality(cellchat_list[[i]], slot.name = "netP")
}
ht_out_control <- netAnalysis_signalingRole_heatmap(cellchat_list[[1]], pattern = "outgoing", width=12, height=14, color.use = mycolor2_named,signaling = keep_pathways_1)
ht_in_control <- netAnalysis_signalingRole_heatmap(cellchat_list[[1]], pattern = "incoming", width=12, height=14, color.use = mycolor2_named,signaling = keep_pathways_1)
save_plot_auto_2.5(ht_out_control + ht_in_control, "SignalingRole_Heatmap_control")

ht_out_PSO <- netAnalysis_signalingRole_heatmap(cellchat_list[[2]], pattern = "outgoing", width=12, height=14, color.use = mycolor2_named,signaling = keep_pathways_2)
ht_in_PSO <- netAnalysis_signalingRole_heatmap(cellchat_list[[2]], pattern = "incoming", width=12, height=14, color.use = mycolor2_named,signaling = keep_pathways_2)
save_plot_auto_2.5(ht_out_PSO + ht_in_PSO, "SignalingRole_Heatmapr_PSO")

p1 <- netVisual_aggregate(
  cellchat_list[[1]], 
  signaling = keep_pathways_1, 
  layout = "chord",
  signaling.name = "Control",
  title.name = "NULL" 
)
save_plot_auto_2.5(p1, "_Cont_circle")

p2 <- netVisual_aggregate(
  cellchat_list[[2]],
  signaling = keep_pathways_2, 
  color.use = mycolor2, 
  layout = "chord", 
  signaling.name = "Psoriasis",
  title.name = "NULL"
)
save_plot_auto_2.5(p2, "_Pso_circle")

gg1 <- netAnalysis_signalingRole_scatter(cellchat_list[[1]],color.use = mycolor2,signaling = keep_pathways_1)
save_plot_auto_2.5(gg1,"netAnalysis_Cont")
gg2 <- netAnalysis_signalingRole_scatter(cellchat_list[[2]],color.use = mycolor2,signaling = keep_pathways_2)
save_plot_auto_2.5(gg2,"netAnalysis_Pso")

bubble_Cont = netVisual_bubble(cellchat_list[[1]], remove.isolate = FALSE,signaling = keep_pathways_1 )
save_plot_auto_2.5(bubble_Cont,"bubble_Cont", width = 30)
bubble_Pso = netVisual_bubble(cellchat_list[[2]], remove.isolate = FALSE,signaling = keep_pathways_2)
save_plot_auto_2.5(bubble_Pso,"bubble_Pso", width = 30)

NMF_control <- selectK(cellchat_list[[1]], pattern = "outgoing")
NMF_PSO <-  selectK(cellchat_list[[2]], pattern = "outgoing")
save_plot_auto_2.5(NMF_control + NMF_PSO, "17_NMF_PSO")

k_c = 5
k_p = 6
cellchat_list[[1]] <- identifyCommunicationPatterns(cellchat_list[[1]], pattern = "outgoing", k = 5)
cellchat_list[[2]] <- identifyCommunicationPatterns(cellchat_list[[2]], pattern = "outgoing", k = 6)

png("./cellchat/Pattern_River_Outgoing_Cont.png", width = 8, height = 8,res = 600, units = "in")
netAnalysis_river(cellchat_list[[1]], pattern = "outgoing")
dev.off()
pdf("./cellchat/Pattern_River_Outgoing_Cont.pdf", width = 8, height = 8)
netAnalysis_river(cellchat_list[[1]], pattern = "outgoing")
dev.off()

png("./cellchat/Pattern_Dot_Outgoing_Cont.png", width = 8, height = 8,res = 600, units = "in")
netAnalysis_dot(cellchat_list[[1]], pattern = "outgoing")
dev.off()
pdf("./cellchat/Pattern_Dot_Outgoing_Cont.pdf", width = 8, height = 8)
netAnalysis_dot(cellchat_list[[1]], pattern = "outgoing")
dev.off()

png("./cellchat/Pattern_River_Outgoing_Pso.png", width = 8, height = 8,res = 600, units = "in")
netAnalysis_river(cellchat_list[[2]], pattern = "outgoing")
dev.off()
pdf("./cellchat/Pattern_River_Outgoing_Pso.pdf", width = 8, height = 8)
netAnalysis_river(cellchat_list[[2]], pattern = "outgoing")
dev.off()

png("./cellchat/Pattern_River_Outgoing_Pso.png", width = 8, height = 8,res = 600, units = "in")
netAnalysis_river(cellchat_list[[2]], pattern = "outgoing")
dev.off()
pdf("./cellchat/Pattern_River_Outgoing_Pso.pdf", width = 8, height = 8)
netAnalysis_river(cellchat_list[[2]], pattern = "outgoing")
dev.off()

plot_cellchat_LR <- function(cellchat_obj, 
                             celltype_inter, 
                             celltype_color,
                             ligand_col = "#FB8072",
                             receptor_col = "#1988B0", 
                             text_size,
                             legend_in_plot = T,
                             Group,
                             top_n) {
  
  library(tidyverse)
  library(circlize)
  library(ggsci)
  library(igraph)
  library(gtools)
  library(ComplexHeatmap)
  library(stringr)
  
  net.df <- subsetCommunication(cellchat_obj)
  net.df$source <- str_remove_all(net.df$source, "[^[:alnum:][:space:]]")
  net.df$cell_type_pair <- paste0(net.df$source, '_', net.df$target)
  
  net.df_filtered <- net.df %>%
    filter(source %in% celltype_inter & target %in% celltype_inter) %>%
    group_by(cell_type_pair) %>%
    slice_max(prob, n = top_n, with_ties = FALSE) %>%
    ungroup()
  
  if (nrow(net.df_filtered) == 0) {
    message("No interactions found after filtering")
    return()
  }
  
  node_v1 <- vector()
  node_v2 <- vector()
  weight_v <- vector()
  
  for (i in 1:nrow(net.df_filtered)) {
    cci <- strsplit(net.df_filtered$interaction_name_2[i], split = ' ')[[1]]
    ct.pair <- strsplit(net.df_filtered$cell_type_pair[i], split = '_')[[1]]
    prefix1 <- ct.pair[1]
    prefix2 <- ct.pair[2]
    
    node1 <- paste0(prefix1, '_', cci[1], '_L')  
    
    id2 <- cci[length(cci)]
    weight <- net.df_filtered$prob[i]
    
    if (length(grep('\\+', id2)) > 0) {
      array <- strsplit(id2, split = '\\+')[[1]]
      node2 <- paste0(prefix2, '_', substr(array[1], 2, nchar(array[1])), '_R')
      node3 <- paste0(prefix2, '_', substr(array[2], 1, nchar(array[2]) - 1), '_R')
      node_v1 <- c(node_v1, node1, node1)
      node_v2 <- c(node_v2, node2, node3)
      weight_v <- c(weight_v, weight, weight)
    } else {
      node2 <- paste0(prefix2, '_', id2, '_R')
      node_v1 <- c(node_v1, node1)
      node_v2 <- c(node_v2, node2)
      weight_v <- c(weight_v, weight)
    }
  }
  
  if (length(node_v1) <= 0) { return() }
  
  g_df <- data.frame(node_v1, node_v2, weight_v)
  g <- graph.data.frame(g_df, directed = T)
  E(g)$weight <- g_df[[3]]
  adj <- get.adjacency(g, attr = 'weight')
  
  graph_adj <- as.data.frame(as.matrix(adj))
  Genes <- names(V(g))
  
  arrays <- strsplit(Genes, split = '_')
  ID <- sapply(arrays, function(x) {
    cell_gene <- paste(x[-length(x)], collapse = "_")
    return(strsplit(cell_gene, split = '_')[[1]][1])
  })
  
  category <- data.frame(Genes = Genes, ID = ID)
  graph_module <- category
  graph_module$LR <- ifelse(grepl("_L$", graph_module$Genes), "ligand", "receptor")
  
  if(is.null(celltype_inter)==F){
    graph_module <- graph_module[graph_module$ID %in% celltype_inter,]
    str <- paste0("^(", paste(celltype_inter, collapse = "|"), ")")
    str_cells <- grepl(str, rownames(graph_adj))
    str_cells1 <- rownames(graph_adj)[str_cells]
    str_cells2 <- colnames(graph_adj)[str_cells]
    
    graph_adj <- graph_adj[str_cells1,]
    graph_adj <- graph_adj[,str_cells2]
  }
  
  g <- graph.adjacency(as.matrix(graph_adj), weighted = T)
  
  LR_color <- data.frame(LR = c("ligand", "receptor"),
                         color1 = c(ligand_col, receptor_col))
  graph_module <- left_join(graph_module, LR_color, by = "LR")
  
  if (is.vector(celltype_color) && !is.null(names(celltype_color))) {
    celltype_color <- data.frame(ID = names(celltype_color), 
                                 color2 = as.character(celltype_color),
                                 stringsAsFactors = FALSE)
  } else {
    celltype_color <- as.data.frame(celltype_color)
    colnames(celltype_color) <- c("ID", 'color2')
  }
  
  graph_module <- left_join(graph_module, celltype_color, by = "ID")
  
  raw_edges <- as.data.frame(cbind(get.edgelist(g), E(g)$weight)) %>%
    mutate(
      V1 = gsub('\\.', '-', V1),
      V2 = gsub('\\.', '-', V2),
      V3 = as.numeric(V3),
      V4 = 1
    )
  edges <- raw_edges %>% arrange(V3)
  
  nodes <- unique(c(edges$V1, edges$V2))
  sectors <- sort(unique(c(raw_edges$V1, raw_edges$V2)))
  
  col_fun = colorRamp2(range(edges$V3), c("#FFFDE7", "#013220"))
  
  circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(-0.15, 0.2))
  circos.initialize(sectors, xlim = c(0, 1))
  circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05, bg.border = NA)
  
  circos.track(
    track.index = 1,
    panel.fun = function(x, y) {
      sector.name = get.cell.meta.data("sector.index")
      xlim = get.cell.meta.data("xlim")
      
      display_name <- gsub("_(L|R)$", "", sector.name)
      
      node_text_color <- graph_module %>%
        dplyr::filter(Genes == sector.name) %>%
        pull(color2) %>%
        as.character()
      
      if(length(node_text_color) == 0 || is.na(node_text_color)) node_text_color <- "black"
      
      node_LR_color <- graph_module %>%
        dplyr::filter(Genes == sector.name) %>%
        pull(color1) %>%
        as.character()
      
      circos.rect(
        xlim[1], 0, xlim[2], 1,
        col = node_LR_color,
        border = NA
      )
      
      circos.text(
        mean(xlim), 2,
        display_name, 
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        col = node_text_color, 
        cex = text_size
      )
    },
    bg.border = NA
  ) 
  
  for (i in seq_len(nrow(edges))) {
    link <- edges[i,]
    circos.link(link[[1]], c(0, 1), link[[2]], c(0, 1),
                col = col_fun(link[[3]]), border = NA)
  }
  
  title(paste0("Cell interactions in ", Group))
  
  if(legend_in_plot == T){
    lgd <- Legend(title = "Score", col_fun = col_fun, direction = "horizontal", border = "black")
    grid.draw(lgd)
    
    legend(1.1, 0.5, pch = 15, legend = c("Ligand", "Receptor"), bty = "n",
           col = c(ligand_col, receptor_col), cex = 1, pt.cex = 3, border = "black") 
    
    legend(1.1, 0, pch = 20, legend = celltype_color$ID, bty = "n",
           col = celltype_color$color2, cex = 1, pt.cex = 3, border = "black") 
  }
}

pdf("./cellchat_new_2.5/Sender_recever_LR_Pso.pdf", width = 20, height = 20)
p2 <- 
  plot_cellchat_LR_v2(
    cellchat_obj = cellchat_list[[2]],
    celltype_inter = target_celltypes,
    celltype_color = mycolor2_named,
    ligand_col = "#FB8072",
    receptor_col = "#1988B0",
    text_size = 0.8,
    Group = "Psoriasis",
    legend_in_plot = TRUE,
    top_n = 1
  )
dev.off()
png("./cellchat_new_2.5/Sender_recever_LR_Pso.png", width = 20, height = 20,res = 600, units = "in")
p2 <- 
  plot_cellchat_LR_v2(
    cellchat_obj = cellchat_list[[2]],
    celltype_inter = target_celltypes,
    celltype_color = mycolor2_named,
    ligand_col = "#FB8072",
    receptor_col = "#1988B0",
    text_size = 0.8,
    Group = "Psoriasis",
    legend_in_plot = TRUE,
    top_n = 1
  )
dev.off()

pdf("./cellchat_new_2.5/Sender_recever_LR_Cont.pdf", width = 20, height = 20)
p1 <- 
  plot_cellchat_LR_v2(
    cellchat_obj = cellchat_list[[1]],
    celltype_inter = target_celltypes,
    celltype_color = mycolor2_named,
    ligand_col = "#FB8072",
    receptor_col = "#1988B0",
    text_size = 0.8,
    Group = "Psoriasis",
    legend_in_plot = TRUE,
    top_n = 1
  )
dev.off()
png("./cellchat_new_2.5/Sender_recever_LR_Cont.png", width = 20, height = 20,res = 600, units = "in")
p1 <- 
  plot_cellchat_LR_v2(
    cellchat_obj = cellchat_list[[1]],
    celltype_inter = target_celltypes,
    celltype_color = mycolor2_named,
    ligand_col = "#FB8072",
    receptor_col = "#1988B0",
    text_size = 0.8,
    Group = "Psoriasis",
    legend_in_plot = TRUE,
    top_n = 1
  )
dev.off()

gg1 <- compareInteractions(cellchat_merged, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat_merged, show.legend = F, group = c(1,2), measure = "weight")
save_plot_auto_2.5(gg1 + gg2, "5_Compare_Interactions_Barplot")

p_circle_diff_expr <- expression({
  par(mfrow = c(1,2), xpd=TRUE)
  
  netVisual_diffInteraction(cellchat_merged, weight.scale = T, comparison = c(1,2),
                            margin = 0.1,
                            arrow.width = 1,
                            color.use = mycolor2, 
                            arrow.size = 0.2)
  
  netVisual_diffInteraction(cellchat_merged, weight.scale = T, measure = "weight", comparison = c(1,2),
                            margin = 0.1,
                            arrow.width = 1,
                            color.use = mycolor2, 
                            arrow.size = 0.2)
  
  title("Red: Increased in Psoriasis | Blue: Decreased", outer = TRUE, line = -1)
})
save_plot_auto_2.5(p_circle_diff_expr, "6_Diff_Interaction_CirclePlot")

ht1 <- netVisual_heatmap(cellchat_merged, measure = "count", comparison = c(1,2), color.use = mycolor2_named)
ht2 <- netVisual_heatmap(cellchat_merged, measure = "weight", comparison = c(1,2), color.use = mycolor2_named)
save_plot_auto_2.5(expression({
  print(ht1 + ht2)
}), "7_Diff_Interaction_Heatmap")

if ("celltype" %in% colnames(cellchat_merged@meta)) {
  new_idents <- as.factor(cellchat_merged@meta$celltype)
  cellchat_merged@idents <- new_idents
  print("Successfully updated CellChat identity! Current identity levels:")
  print(levels(cellchat_merged@idents))
} else {
  stop("Error: cannot find 'celltype' column in cellchat_merged@meta!")
}

all_idents <- levels(cellchat_merged@idents)

gdt_cells <- grep("γδT|gdT", all_idents, value = TRUE)

PBMC_2.8_cells <- setdiff(all_idents, gdt_cells)

message(paste("Senders (gdT):", paste(gdt_cells, collapse = ", ")))
message(paste("Receivers (PBMC):", paste(PBMC_2.8_cells, collapse = ", ")))

target_patterns <- PBMC_2.8_cells

regex_pattern <- paste(target_patterns, collapse = "|")
keep_pathways_merged <- union(keep_pathways_1, keep_pathways_2)
if(length(gdt_cells) > 0 && length(PBMC_2.8_cells) > 0) {
  
  p_bubble <- netVisual_bubble(cellchat_merged, 
                               sources.use = gdt_cells, 
                               targets.use = PBMC_2.8_cells, 
                               comparison = c(1,2),
                               angle.x = 45, 
                               signaling = keep_pathways_merged,
                               remove.isolate = FALSE,
                               title.name = "Diff Interactions: gdT -> PBMC_2.8 Cells")
  
  save_plot_auto_2.5(p_bubble, "8_Diff_Bubble_gdT_Sender_PBMC_2.8_Receiver")
  
} else {
  message("Specified gdT or PBMC_2.8 cells not found, please check name matching rules")
}

if(length(gdt_cells) > 0 && length(PBMC_2.8_cells) > 0) {
  
  p_bubble_rev <- netVisual_bubble(cellchat_merged, 
                                   sources.use = PBMC_2.8_cells, 
                                   targets.use = gdt_cells, 
                                   comparison = c(1,2), 
                                   signaling = keep_pathways_merged,
                                   angle.x = 45, 
                                   remove.isolate = FALSE,
                                   title.name = "Diff Interactions: PBMC_2.8 Cells -> gdT")
  
  save_plot_auto_2.5(p_bubble_rev, "8_Diff_Bubble_PBMC_2.8_Sender_gdT_Receiver")
}

saveRDS(cellchat_list, file = "./cellchat_list_PBMC_2.14.rds")




seu_combined_foucus <- subset(seu_combined, 
                              celltype == "CSMD1+" | !grepl("γδT", celltype))

seu_combined_foucus$celltype <- droplevels(seu_combined_foucus$celltype)

table(seu_combined_foucus$celltype)

print("Dataset merge completed.")

save_plot_auto_2.5 <- function(plot_code_or_obj, filename_base, width = 18, height = 9.5, type = "both") {
  if(!dir.exists("CSMD1_cellchat")) dir.create("CSMD1_cellchat", recursive = TRUE)
  f_png <- paste0("CSMD1_cellchat/", filename_base, ".png")
  f_pdf <- paste0("CSMD1_cellchat/", filename_base, ".pdf")
  
  draw_it <- function(obj) {
    if (inherits(obj, "grob")) {
      grid::grid.draw(obj)
    } else if (is.expression(obj)) {
      eval(obj)
    } else {
      print(obj)
    }
  }
  
  if(type %in% c("both", "png")) {
    png(f_png, width = width, height = height, units = "in", res = 600)
    grid::grid.newpage() 
    tryCatch({
      draw_it(plot_code_or_obj)
    }, error = function(e) message(paste("Error saving PNG:", filename_base, e$message)))
    dev.off()
  }
  
  if(type %in% c("both", "pdf")) {
    pdf(f_pdf, width = width, height = height)
    grid::grid.newpage()
    tryCatch({
      draw_it(plot_code_or_obj)
    }, error = function(e) message(paste("Error saving PDF:", filename_base, e$message)))
    dev.off()
  }
  
  message(paste("Saved successfully:", filename_base))
}

run_cellchat_pipeline <- function(seurat_obj, group_name) {
  message(paste0("Processing group: ", group_name))
  
  data.input <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  meta <- seurat_obj@meta.data
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
  
  CellChatDB <- CellChatDB.human
  cellchat@DB <- subsetDB(CellChatDB, search = "Secreted Signaling")
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  return(cellchat)
}

seurat_list <- SplitObject(seu_combined_foucus, split.by = "group")
Idents(seu_combined_foucus) <- "celltype"

cellchat_list <- list()

for (grp in names(seurat_list)) {
  
  seurat_obj <- seurat_list[[grp]]
  
  seurat_obj$celltype <- droplevels(as.factor(seurat_obj$celltype))
  
  Idents(seurat_obj) <- seurat_obj$celltype
  
  print(paste0(">>> Processing group: ", grp))
  print(paste("    Total cells:", ncol(seurat_obj)))
  print(paste("    Number of cell types:", length(levels(Idents(seurat_obj)))))
  
  if(ncol(seurat_obj) > 50 && length(levels(Idents(seurat_obj))) >= 2) {
    
    cellchat_list[[grp]] <- run_cellchat_pipeline(seurat_obj, grp)
    
  } else {
    print(paste("    [Skip] Too few cells or single type"))
  }
}

cellchat_list <- list()
for (grp in names(seurat_list)) {
  if(ncol(seurat_list[[grp]]) > 50) {
    cellchat_list[[grp]] <- run_cellchat_pipeline(seurat_list[[grp]], grp)
  }
}

group.cellType <- union(levels(cellchat_list[[1]]@idents), levels(cellchat_list[[2]]@idents))
for (i in 1:length(cellchat_list)) {
  cellchat_list[[i]] <- liftCellChat(cellchat_list[[i]], group.new = group.cellType)
}
cellchat_merged <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))
compare_vec <- c(2, 1)

message("--- Re-merging CellChat object ---")
cellchat_merged <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))

print("Merged cell type order:")
print(levels(cellchat_merged@idents$joint))

all_pathways_1 <- cellchat_list[[1]]@netP$pathways
keep_pathways_1 <- all_pathways_1[!grepl("HLA|MHC", all_pathways_1, ignore.case = TRUE)]

all_pathways_2 <- cellchat_list[[2]]@netP$pathways
keep_pathways_2 <- all_pathways_2[!grepl("HLA|MHC", all_pathways_2, ignore.case = TRUE)]

print(head(cellchat_list[[1]]@netP$pathways))
print(head(cellchat_list[[2]]@netP$pathways))

p1 <- cellchat_list[[1]]@netP$pathways
keep_pathways_1 <- p1[!grepl("HLA|MHC", p1, ignore.case = TRUE)]

p2 <- cellchat_list[[2]]@netP$pathways
keep_pathways_2 <- p2[!grepl("HLA|MHC", p2, ignore.case = TRUE)]

mat_control <- cellchat_list[[1]]@net$count
mat_disease <- cellchat_list[[2]]@net$count

library(pheatmap)
library(grid)
library(gridExtra)

global_max <- max(max(mat_control, na.rm = T), max(mat_disease, na.rm = T))
common_breaks <- seq(0, global_max, length.out = 100)
common_color <- colorRampPalette(c("white", "#E41A1C"))(100)

if(!dir.exists("cellchat")) dir.create("cellchat")

p3_ctrl <- pheatmap::pheatmap(cellchat_list[[1]]@net$count, border_color = "black",
                              cluster_cols = F, cluster_rows = F, fontsize = 10,
                              display_numbers = T, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Control", silent = TRUE)$gtable

p3_case <- pheatmap::pheatmap(mat_disease, border_color = "black",
                              cluster_cols = F, cluster_rows = F, fontsize = 10,
                              display_numbers = T, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks, 
                              main = "Psoriasis", silent = TRUE)$gtable

p3_combined <- arrangeGrob(p3_ctrl, p3_case, ncol = 2, top = textGrob("Interaction Count (Raw)", gp=gpar(fontsize=15, font=2)))
save_plot_auto_2.5(p3_combined, "CellChat_Heatmap_NoCluster_WithNum")

p4_ctrl <- pheatmap::pheatmap(mat_control, border_color = "black",
                              cluster_cols = T, cluster_rows = T, fontsize = 10,
                              display_numbers = T, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Control", silent = TRUE)$gtable

p4_case <- pheatmap::pheatmap(mat_disease, border_color = "black",
                              cluster_cols = T, cluster_rows = T, fontsize = 10,
                              display_numbers = T, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Psoriasis", silent = TRUE)$gtable

p4_combined <- arrangeGrob(p4_ctrl, p4_case, ncol = 2, top = textGrob("Interaction Count (Clustered)", gp=gpar(fontsize=15, font=2)))
save_plot_auto_2.5(p4_combined, "CellChat_Heatmap_Clustered_WithNum")

p5_ctrl <- pheatmap::pheatmap(mat_control, border_color = "black",
                              cluster_cols = F, cluster_rows = F, fontsize = 10,
                              display_numbers = F, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Control", silent = TRUE)$gtable

p5_case <- pheatmap::pheatmap(mat_disease, border_color = "black",
                              cluster_cols = F, cluster_rows = F, fontsize = 10,
                              display_numbers = F, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Psoriasis", silent = TRUE)$gtable

p5_combined <- arrangeGrob(p5_ctrl, p5_case, ncol = 2, top = textGrob("Interaction Count (Raw, No Num)", gp=gpar(fontsize=15, font=2)))
save_plot_auto_2.5(p5_combined, "CellChat_Heatmap_NoCluster_NoNum")

p6_ctrl <- pheatmap::pheatmap(mat_control, border_color = "black",
                              cluster_cols = T, cluster_rows = T, fontsize = 10,
                              display_numbers = F, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Control", silent = TRUE)$gtable

p6_case <- pheatmap::pheatmap(mat_disease, border_color = "black",
                              cluster_cols = T, cluster_rows = T, fontsize = 10,
                              display_numbers = F, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Psoriasis", silent = TRUE)$gtable

p6_combined <- arrangeGrob(p6_ctrl, p6_case, ncol = 2, top = textGrob("Interaction Count (Clustered, No Num)", gp=gpar(fontsize=15, font=2)))
save_plot_auto_2.5(p6_combined, "CellChat_Heatmap_Clustered_NoNum")

mat <- cellchat_merged@net$count
cell_types <- levels(cellchat_merged@idents$joint)
groupSize <- as.numeric(table(cellchat_merged@idents$joint))

code_single_circle <- expression({
  par(mfrow = c(3,3), xpd=TRUE, mar = c(1,1,1,1))
  for (i in 1:length(cell_types)) {
    netVisual_circle(cellchat_list[["Psoriasis"]]@net$count, vertex.weight = groupSize, weight.scale = T, 
                     edge.weight.max = max(cellchat_list[["Psoriasis"]]@net$count), 
                     title.name = paste("Pso -", cell_types[i]), color.use = mycolor2_foucus, sources.use = cell_types[i])
  }
})
save_plot_auto_2.5(code_single_circle, "03_Single_CellType_Source_Circle_Psoriasis")
code_single_circle <- expression({
  par(mfrow = c(3,3), xpd=TRUE, mar = c(1,1,1,1))
  for (i in 1:length(cell_types)) {
    netVisual_circle(cellchat_list[["Control"]]@net$count, vertex.weight = groupSize, weight.scale = T, 
                     edge.weight.max = max(cellchat_list[["Control"]]@net$count), 
                     title.name = paste("Pso -", cell_types[i]), color.use = mycolor2_foucus, sources.use = cell_types[i])
  }
})
save_plot_auto_2.5(code_single_circle, "03_Single_CellType_Source_Circle_Control")

groupSize <- as.numeric(table(cellchat_list[[1]]@idents))

p1 <- expression({
  par(mfrow = c(1,2), xpd=TRUE)
  
  netVisual_circle(cellchat_list[[1]]@net$count, vertex.weight = groupSize, weight.scale = T,
                   label.edge= F,                    arrow.width = 1.2,
                   arrow.size = 0.24, color.use = mycolor2_foucus,
                   title.name = "Number of interactions in Cont")
  netVisual_circle(cellchat_list[[1]]@net$weight, vertex.weight = groupSize, weight.scale = T, 
                   label.edge= F, arrow.width = 1.2, color.use = mycolor2_foucus,
                   arrow.size = 0.24 ,title.name = "Interaction weights/strength in Cont")
  
  groupSize <- as.numeric(table(cellchat_list[[1]]@idents))
})
save_plot_auto_2.5(p1, "03_Total_CellType_Source_Circle_Control")

p2 <- expression({
  par(mfrow = c(1,2), xpd=TRUE)
  
  netVisual_circle(cellchat_list[[2]]@net$count, vertex.weight = groupSize, weight.scale = T,
                   label.edge= F, 
                   arrow.width = 1.2,
                   arrow.size = 0.24, color.use = mycolor2_foucus,
                   title.name = "Number of interactions in Pso")
  netVisual_circle(cellchat_list[[2]]@net$weight, vertex.weight = groupSize, weight.scale = T, 
                   label.edge= F, 
                   arrow.width = 1.2,
                   arrow.size = 0.24, color.use = mycolor2_foucus,
                   title.name = "Interaction weights/strength in Pso")
})
save_plot_auto_2.5(p2, "03_Total_CellType_Source_Circle_Psoriasis")

cellchat_list[[1]]@netP$pathways

pp1 <- rankNet(cellchat_merged, mode = "comparison", stacked = T, do.stat = TRUE)
pp2 <- rankNet(cellchat_merged, mode = "comparison", stacked = F, do.stat = TRUE)
save_plot_auto_2.5(pp1 + pp2, "03_rankNet")

CSMD1_chord <- expression(netVisual_chord_gene(cellchat_list[[1]], sources.use = "CSMD1+", lab.cex = 0.5,legend.pos.y =50))
save_plot_auto_2.5(CSMD1_chord,"CSMD1_chord")

for (i in 1:length(cellchat_list)) {
  cellchat_list[[i]] <- netAnalysis_computeCentrality(cellchat_list[[i]], slot.name = "netP")
}
ht_out_control <- netAnalysis_signalingRole_heatmap(cellchat_list[[1]], pattern = "outgoing", width=12, height=14, color.use = mycolor2_named_foucus,signaling = keep_pathways_1)
ht_in_control <- netAnalysis_signalingRole_heatmap(cellchat_list[[1]], pattern = "incoming", width=12, height=14, color.use = mycolor2_named_foucus,signaling = keep_pathways_1)
save_plot_auto_2.5(ht_out_control + ht_in_control, "SignalingRole_Heatmap_control")

ht_out_PSO <- netAnalysis_signalingRole_heatmap(cellchat_list[[2]], pattern = "outgoing", width=12, height=14, color.use = mycolor2_named_foucus,signaling = keep_pathways_2)
ht_in_PSO <- netAnalysis_signalingRole_heatmap(cellchat_list[[2]], pattern = "incoming", width=12, height=14, color.use = mycolor2_named_foucus,signaling = keep_pathways_2)
save_plot_auto_2.5(ht_out_PSO + ht_in_PSO, "SignalingRole_Heatmapr_PSO")

p1 <- netVisual_aggregate(
  cellchat_list[[1]],
  signaling = keep_pathways_1,
  layout = "chord",
  signaling.name = "Control",
  title.name = "NULL"
)
save_plot_auto_2.5(p1, "_Cont_circle")

p2 <- netVisual_aggregate(
  cellchat_list[[2]],
  signaling = keep_pathways_2,
  color.use = mycolor2_foucus,
  layout = "chord",
  signaling.name = "Psoriasis",
  title.name = "NULL"
)
save_plot_auto_2.5(p2, "_Pso_circle")

gg1 <- netAnalysis_signalingRole_scatter(cellchat_list[[1]],color.use = mycolor2_foucus,signaling = keep_pathways_1)
save_plot_auto_2.5(gg1,"netAnalysis_Cont")
gg2 <- netAnalysis_signalingRole_scatter(cellchat_list[[2]],color.use = mycolor2_foucus,signaling = keep_pathways_2)
save_plot_auto_2.5(gg2,"netAnalysis_Pso")

bubble_Cont = netVisual_bubble(cellchat_list[[1]], remove.isolate = FALSE,signaling = keep_pathways_1 )
save_plot_auto_2.5(bubble_Cont,"bubble_Cont", width = 30)
bubble_Pso = netVisual_bubble(cellchat_list[[2]], remove.isolate = FALSE,signaling = keep_pathways_2)
save_plot_auto_2.5(bubble_Pso,"bubble_Pso", width = 30)

gg1 <- compareInteractions(cellchat_merged, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat_merged, show.legend = F, group = c(1,2), measure = "weight")
save_plot_auto_2.5(gg1 + gg2, "5_Compare_Interactions_Barplot")

p_circle_diff_expr <- expression({
  par(mfrow = c(1,2), xpd=TRUE)
  
  netVisual_diffInteraction(cellchat_merged, weight.scale = T, comparison = c(1,2),
                            margin = 0.1,
                            arrow.width = 1,
                            arrow.size = 0.2)
  
  netVisual_diffInteraction(cellchat_merged, weight.scale = T, measure = "weight", comparison = c(1,2),
                            margin = 0.1,
                            arrow.width = 1,
                            arrow.size = 0.2)
  
  title("Red: Increased in Psoriasis | Blue: Decreased", outer = TRUE, line = -1)
})
save_plot_auto_2.5(p_circle_diff_expr, "6_Diff_Interaction_CirclePlot")

ht1 <- netVisual_heatmap(cellchat_merged, measure = "count", comparison = c(1,2), color.use = mycolor2_named_foucus)
ht2 <- netVisual_heatmap(cellchat_merged, measure = "weight", comparison = c(1,2), color.use = mycolor2_named_foucus)
save_plot_auto_2.5(expression({
  print(ht1 + ht2)
}), "7_Diff_Interaction_Heatmap")

if ("celltype" %in% colnames(cellchat_merged@meta)) {
  new_idents <- as.factor(cellchat_merged@meta$celltype)
  cellchat_merged@idents <- new_idents
  print("Successfully updated CellChat identity! Current identity levels:")
  print(levels(cellchat_merged@idents))
} else {
  stop("Error: cannot find 'celltype' column in cellchat_merged@meta!")
}

all_idents <- levels(cellchat_merged@idents)

gdt_cells <- grep("γδT|gdT", all_idents, value = TRUE)

PBMC_2.8_cells <- setdiff(all_idents, gdt_cells)

message(paste("Senders (gdT):", paste(gdt_cells, collapse = ", ")))
message(paste("Receivers (PBMC):", paste(PBMC_2.8_cells, collapse = ", ")))

target_patterns <- PBMC_2.8_cells

if(length(gdt_cells) > 0 && length(PBMC_2.8_cells) > 0) {
  
  p_bubble <- netVisual_bubble(cellchat_merged, 
                               sources.use = gdt_cells, 
                               targets.use = PBMC_2.8_cells, 
                               comparison = c(1,2),
                               signaling = keep_pathways_merged,
                               angle.x = 45, 
                               remove.isolate = FALSE,
                               title.name = "Diff Interactions: gdT -> PBMC_2.8 Cells")
  
  save_plot_auto_2.5(p_bubble, "8_Diff_Bubble_gdT_Sender_PBMC_2.8_Receiver")
  
} else {
  message("Specified gdT or PBMC_2.8 cells not found, please check name matching rules")
}

keep_pathways_merged <- union(keep_pathways_1, keep_pathways_2)

if(length(gdt_cells) > 0 && length(PBMC_2.8_cells) > 0) {
  p_bubble_rev <- netVisual_bubble(cellchat_merged, 
                                   sources.use = PBMC_2.8_cells, 
                                   targets.use = gdt_cells, 
                                   comparison = c(1,2), 
                                   signaling = keep_pathways_merged,
                                   angle.x = 45, 
                                   remove.isolate = FALSE,
                                   title.name = "Diff Interactions: PBMC_2.8 Cells -> gdT")
  
  save_plot_auto_2.5(p_bubble_rev, "8_Diff_Bubble_PBMC_2.8_Sender_gdT_Receiver")
}

saveRDS(cellchat_list, file = "./cellchat_list_PBMC_foucus_2.14.rds")





###############(3)CSMD1 subgroups and Skins from public databases############
setwd("/thinker/5.chensihan/shinanxi_cellchat/26.1.24/")

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(ggplot2)
  library(patchwork)
  library(future)
  library(dplyr)
  library(harmony)
})

future::plan("multisession", workers = 4)
options(stringsAsFactors = FALSE)

print("Loading user PBMC data...")
seu_user <- readRDS("../rds/Overall_shinanxi_New.rds")

seu_user$Dataset <- "My_PBMC"
seu_user$final_annotation <- as.character(seu_user$celltype)

print(table(seu_user$final_annotation))


load_geo_data <- function(path, project_name, dataset_label) {
  
  read_one_10x <- function(dir){
    counts <- tryCatch(
      Read10X(data.dir = dir),
      error = function(e) {
        ReadMtx(
          mtx      = file.path(dir, "matrix.mtx.gz"),
          cells    = file.path(dir, "barcodes.tsv.gz"),
          features = file.path(dir, "features.tsv.gz")
        )
      }
    )
    seu <- CreateSeuratObject(counts = counts, project = project_name)
    return(seu)
  }
  
  qc_and_norm <- function(seu){
    seu$Dataset <- dataset_label
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
    seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
    
    seu <- NormalizeData(seu)
    seu <- FindVariableFeatures(seu)
    seu <- ScaleData(seu)
    return(seu)
  }
  
  if (dir.exists(path)) {
    
    has_mtx <- file.exists(file.path(path, "matrix.mtx")) || file.exists(file.path(path, "matrix.mtx.gz"))
    if (has_mtx) {
      obj <- read_one_10x(path)
      obj$sample <- basename(path)
      obj <- qc_and_norm(obj)
      return(obj)
    }
    
    sample_dirs <- list.dirs(path, recursive = FALSE, full.names = TRUE)
    sample_dirs <- sample_dirs[
      file.exists(file.path(sample_dirs, "matrix.mtx.gz")) | file.exists(file.path(sample_dirs, "matrix.mtx"))
    ]
    
    if (length(sample_dirs) == 0) stop(paste("No 10X subdirectory found:", path))
    
    obj_list <- lapply(sample_dirs, function(d){
      seu <- read_one_10x(d)
      smp <- basename(d)
      seu$sample <- smp
      seu$cell_id <- ifelse(grepl("control", smp, ignore.case = TRUE), "Control", "Psoriasis")
      seu <- qc_and_norm(seu)
      seu
    })
    names(obj_list) <- basename(sample_dirs)
    
    obj <- merge(
      x = obj_list[[1]],
      y = obj_list[-1],
      add.cell.ids = names(obj_list),
      project = project_name
    )
    return(obj)
    
  } else if (file.exists(path) && grepl("\\.rds$", path, ignore.case = TRUE)) {
    
    obj <- readRDS(path)
    if (!"percent.mt" %in% colnames(obj@meta.data)) {
      obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    }
    obj$Dataset <- dataset_label
    obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj)
    obj <- ScaleData(obj)
    return(obj)
    
  } else {
    stop(paste("Path does not exist or format not supported:", path))
  }
}

seu_skin <- load_geo_data("/thinker/5.chensihan/shinanxi_cellchat/26.1.24/public_data/Skin", "Skin", "Public_Skin")



seu_skin <- readRDS("../result_skin/seu_skin.rds")

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

save_plot_auto <- function(plot_code_or_obj, filename_base, width = 18, height = 9.5, type = "both") {
  
  f_png <- paste0(filename_base, ".png")
  f_pdf <- paste0( filename_base, ".pdf")
  
  if(type %in% c("both", "png")) {
    png(f_png, width = width, height = height, units = "in", res = 600)
    tryCatch({
      if (is.ggplot(plot_code_or_obj) || inherits(plot_code_or_obj, "HeatmapList")) {
        print(plot_code_or_obj)
      } else {
        eval(plot_code_or_obj)
      }
    }, error = function(e) message(paste("Error saving PNG:", filename_base, e)))
    dev.off()
  }
  
  if(type %in% c("both", "pdf")) {
    pdf(f_pdf, width = width, height = height)
    tryCatch({
      if (is.ggplot(plot_code_or_obj) || inherits(plot_code_or_obj, "HeatmapList")) {
        print(plot_code_or_obj)
      } else {
        eval(plot_code_or_obj)
      }
    }, error = function(e) message(paste("Error saving PDF:", filename_base, e)))
    dev.off()
  }
  
  message(paste("Saved:", filename_base))
}

keep_samples <- c("Control01", "Control02", "Control03", "Control04", "Control05", 
                  "Psoriasis01", "Psoriasis02", "Psoriasis03", "Psoriasis04", 
                  "Psoriasis05", "Psoriasis06", "Psoriasis07", "Psoriasis08", 
                  "Psoriasis09", "Psoriasis10", "Psoriasis11", "Psoriasis12", "Psoriasis13")

seu_skin <- subset(seu_skin, subset = sample %in% grep("F", unique(seu_skin$sample), value = T, invert = T))
seu_skin$sample <- factor(seu_skin$sample)
seu_skin$sample <- droplevels(seu_skin$sample)

obj = SplitObject(seu_skin, split.by = "sample")
obj_rm=list()
doublets_plot = list()
pc.num = 1:20
dir.create("./SingleCell_QC")

RemoveDoublets <-function(
    object,
    doublet.rate,
    pN=0.25,
    pc.num=1:30
){
  sweep.res.list <- paramSweep(object, PCs = pc.num, sct = F)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  homotypic.prop <- modelHomotypic(object$seurat_clusters)
  nExp_poi <- round(doublet.rate*ncol(object)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  seu.scored <- doubletFinder(object, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, 
                              nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
  cname <-colnames(seu.scored[[]])
  DF<-cname[grep('^DF',cname)]
  seu.scored[["doublet"]] <- as.numeric(seu.scored[[DF]]=="Doublet")
  
  seu.removed <- subset(seu.scored, subset = doublet != 1)
  p1 <- DimPlot(seu.scored, group.by = DF)
  res.list <- list("plot"=p1, "obj"=seu.removed)
  return(res.list)
}

for (sample_name in names(obj)) {
  current_seu <- obj[[sample_name]]
  
  current_seu <- NormalizeData(current_seu)
  current_seu <- FindVariableFeatures(current_seu, selection.method = "vst", nfeatures = 2000)
  current_seu <- ScaleData(current_seu)
  current_seu <- RunPCA(current_seu)
  current_seu <- RunUMAP(current_seu, dims = 1:20)
  current_seu <- FindNeighbors(current_seu, dims = pc.num) %>% FindClusters(resolution = 0.3)
  
  cell_count <- ncol(current_seu)
  
  if ( cell_count < 750 ){
    doublet_rate = 0.004
  }else if ( cell_count %in% 750:1499 ) {
    doublet_rate = 0.008
  }else if ( cell_count %in% 1500:2499 ) {
    doublet_rate = 0.016
  }else if ( cell_count %in% 2500:3499 ) {
    doublet_rate = 0.023
  }else if ( cell_count %in% 3500:4499 ) {
    doublet_rate = 0.031
  }else if ( cell_count %in% 4500:5499 ) {
    doublet_rate = 0.039
  }else if ( cell_count %in% 5500:6499 ) {
    doublet_rate = 0.046
  }else if ( cell_count %in% 6500:7499 ) {
    doublet_rate = 0.054
  }else if ( cell_count %in% 7500:8499 ) {
    doublet_rate = 0.061
  }else if ( cell_count %in% 8500:9499 ) {
    doublet_rate = 0.069
  }else if ( cell_count %in% 9500:10499 ) {
    doublet_rate = 0.076
  }else if ( cell_count >= 10500 ) {
    doublet_rate = 0.1
  }
  
  tmp <- RemoveDoublets(
    current_seu, 
    doublet.rate =  doublet_rate,
    pc.num = pc.num
  )
  
  obj_rm[[sample_name]] <- tmp$obj
  doublets_plot[[sample_name]] <- tmp$plot
  
  cat(
    "Sample", sample_name, 
    "cell count:", cell_count, 
    ", doublet rate:", round(doublet_rate * 100, 2), "%\n", 
    sep = ""
  )
}

for (i in 1:18) {
  assign(paste0("p", i), doublets_plot[[i]])
}

dir.create("Doublets")

my_clean_list <- list(
  "Control01" = p1, "Control02" = p2, "Control03" = p3, "Control04" = p4, "Control05" = p5,
  "Psoriasis01" = p6, "Psoriasis02" = p7, "Psoriasis03" = p8, "Psoriasis04" = p9, "Psoriasis05" = p10,
  "Psoriasis06" = p11, "Psoriasis07" = p12, "Psoriasis08" = p13, "Psoriasis09" = p14, "Psoriasis10" = p15,
  "Psoriasis11" = p16, "Psoriasis12" = p17, "Psoriasis13" = p18
)

print("Check Control01 type:")
print(class(my_clean_list[["Control01"]]))

p_con <- wrap_plots(my_clean_list[1:5], ncol = 3) + plot_annotation(title = "Control Group")
ggsave("Doublets/Doublets_Control_Group.png", p_con, width = 16, height = 8, dpi = 600)
ggsave("Doublets/Doublets_Control_Group.pdf", p_con, width = 16, height = 8, dpi = 300)

p_pso_v2 <- wrap_plots(my_clean_list[6:15], ncol = 4) + plot_annotation(title = "Psoriasis V2.0")
ggsave("Doublets/Doublets_Psoriasis_V2.png", p_pso_v2, width = 16, height = 12, dpi = 600)
ggsave("Doublets/Doublets_Psoriasis_V2.pdf", p_pso_v2, width = 16, height = 12, dpi = 300)

p_pso_v3 <- wrap_plots(my_clean_list[16:18], ncol = 3) + plot_annotation(title = "Psoriasis V3.0")
ggsave("Doublets/Doublets_Psoriasis_V3.png", p_pso_v3, width = 16, height = 5, dpi = 600)
ggsave("Doublets/Doublets_Psoriasis_V3.pdf", p_pso_v3, width = 16, height = 5, dpi = 300)

seu_skin <- merge(obj_rm[[1]], y=c(obj_rm[[2]],obj_rm[[3]],obj_rm[[4]],obj_rm[[5]],
                                   obj_rm[[6]],obj_rm[[7]],obj_rm[[8]],obj_rm[[9]],obj_rm[[10]],
                                   obj_rm[[11]],obj_rm[[12]],obj_rm[[13]],obj_rm[[14]],obj_rm[[15]],
                                   obj_rm[[16]],obj_rm[[17]],obj_rm[[18]]))

extracted_samples <- sapply(strsplit(colnames(seu_skin), "_"), `[`, 2)
head(extracted_samples)

seu_skin$sample <- extracted_samples
my_levels <- c(
  "Control01", "Control02", "Control03", "Control04", "Control05",
  "Psoriasis01", "Psoriasis02", "Psoriasis03", "Psoriasis04", 
  "Psoriasis05", "Psoriasis06", "Psoriasis07", "Psoriasis08", 
  "Psoriasis09", "Psoriasis10", "Psoriasis11", "Psoriasis12", "Psoriasis13"
)
seu_skin$sample <- factor(seu_skin$sample, levels = my_levels)

samples <- seu_skin$sample
seu_skin$chemistry <- "V2.0"
v3_samples <- c("Psoriasis11", "Psoriasis12", "Psoriasis13")
seu_skin$chemistry[seu_skin$sample %in% v3_samples] <- "V3.0"
seu_skin$disease <- ifelse(grepl("Psoriasis", seu_skin$sample), "Psoriasis", "Control")
seu_skin$batch_group <- paste(seu_skin$disease, seu_skin$chemistry, sep = "_")
table(seu_skin$batch_group, useNA = "ifany")

Idents(seu_skin) <- 'sample'

rna_genes <- rownames(seu_skin@assays$RNA)

mito_features <- grep("^mt-", rna_genes, value = TRUE, ignore.case = TRUE)
cat("Number of mitochondrial genes found:", length(mito_features), "\n")
if (length(mito_features) > 0) {
  head(mito_features, 20)
  seu_skin[["percent.mt"]] <- PercentageFeatureSet(seu_skin, features = mito_features)
  cat("percent.mt recalculated with detected mitochondrial genes\n")
} else {
  cat("No mt- genes detected, check gene naming\n")
}

HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
seu_skin[["percent.HB"]] <- PercentageFeatureSet(seu_skin, features = HB.genes)
cat("percent.HB recalculated with candidate HB genes\n")

summary(seu_skin@meta.data[, c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB")])

beforeQC_vlnplot = VlnPlot(seu_skin, 
                           features = c("nFeature_RNA", 
                                        "nCount_RNA", 
                                        "percent.mt",
                                        "percent.HB"), 
                           ncol = 4, cols = mycolor1_clusters,
                           pt.size = 0.1)

dir.create("./SingleCell_QC")
ggsave("./SingleCell_QC/BeforeQC_nFeature_nCount_percent.mt_percent.HB_vlnplot.pdf", 
       plot = beforeQC_vlnplot,width = 12, height = 6, dpi=300)
ggsave("./SingleCell_QC/BeforeQC_nFeature_nCount_percent.mt_percent.HB_vlnplot.png",
       plot = beforeQC_vlnplot,width = 12, height = 6, dpi=600)
summary(seu_skin@meta.data[,c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB")])
beforeQC_vlnplot

plot1 <- FeatureScatter(seu_skin, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seu_skin, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(seu_skin, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot <- CombinePlots(plots = list(plot1, plot2, plot3), nrow=1, legend="none")
ggsave("./SingleCell_QC/pearplot_before_qc.pdf", plot = pearplot, width = 12, height = 5) 
ggsave("./SingleCell_QC/pearplot_before_qc.png", plot = pearplot, width = 12, height = 5)

minGene=200
maxGene=6000
pctMT=20
pctHB=3
Idents(seu_skin) <- "sample"
seu_skin = subset(seu_skin, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT & percent.HB < pctHB)

afterQC_vlnplot = VlnPlot(seu_skin, 
                          features = c("nFeature_RNA", 
                                       "nCount_RNA", 
                                       "percent.mt",
                                       "percent.HB"), 
                          ncol = 4, 
                          pt.size = 0.1)
ggsave("./SingleCell_QC/afterQC_nFeature_nCount_percent.mt_percent.HB_vlnplot.pdf", plot = afterQC_vlnplot,width = 12, height = 5)
ggsave("./SingleCell_QC/afterQC_nFeature_nCount_percent.mt_percent.HB_vlnplot.png", plot = afterQC_vlnplot,width = 12, height = 5)

afterQC_vlnplot_clear = VlnPlot(seu_skin, 
                                features = c("nFeature_RNA", 
                                             "nCount_RNA", 
                                             "percent.mt",
                                             "percent.HB"), 
                                ncol = 4, 
                                pt.size = 0)
ggsave("./SingleCell_QC/afterQC_nFeature_nCount_percent.mt_percent.HB_vlnplot_clear.pdf", plot = afterQC_vlnplot_clear,width = 12, height = 5)
ggsave("./SingleCell_QC/afterQC_nFeature_nCount_percent.mt_percent.HB_vlnplot_clear.png", plot = afterQC_vlnplot_clear,width = 12, height = 5)

afterQC_vlnplot

saveRDS(seu_skin,"./SingleCell_QC/After_QC.rds")



mycolor1_clusters <- c(
  "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
  "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
  "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93")

mycolor2 <- c(
  "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F","#FF6A6A","#7FFFD4", "#AB82FF","#90EE90",
  "#00CD00","#008B8B","#6495ED","#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00","#000000","#708090")

mycolor4_samples <- c(
  "#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#228B20",
  "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
  "#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#FFFF02","#E9967A","#4682B4",
  "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

seu_skin <- NormalizeData(seu_skin, normalization.method = "LogNormalize", scale.factor = 10000)
seu_skin <- FindVariableFeatures(seu_skin, selection.method = "vst", nfeatures = 2000)
seu_skin <- ScaleData(seu_skin, features = VariableFeatures(seu_skin)) 
seu_skin <- RunPCA(seu_skin, features = VariableFeatures(object = seu_skin))

Elbow <- ElbowPlot(seu_skin, ndims = 50) 
ggsave("./SingleCell_QC/ElbowPlot.pdf", plot = Elbow,width = 12, height = 5)
ggsave("./SingleCell_QC/ElbowPlot.png", plot = Elbow,width = 12, height = 5)

library(harmony)

Idents(seu_skin) <- "batch_group"
seu_skin <- RunHarmony(seu_skin, group.by.vars = "batch_group", plot_convergence = FALSE)

seu_skin <- RunUMAP(seu_skin, reduction = "harmony", dims = 1:20)
seu_skin <- FindNeighbors(seu_skin, reduction = "harmony", dims = 1:20)
seu_skin <- FindClusters(seu_skin, resolution = 0.6)

Idents(seu_skin) <- "seurat_clusters"
DimPlot(seu_skin, reduction = "umap", label = TRUE) + NoLegend()

dir.create("./Marker")
all.markers = FindAllMarkers(seu_skin, 
                             min.pct = 0.25, 
                             logfc.threshold = 0.25, 
                             only.pos = TRUE)
head(all.markers)

top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top5 = all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

write.table(all.markers, 
            "Marker/all_Markers_of_each_clusters.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")
write.table(top10, 
            "Marker/top10_Markers_of_each_clusters.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")

seu_skin <- ScaleData(seu_skin, features = row.names(seu_skin))

heatmap_plot1 = DoHeatmap(object = seu_skin, label = F, 
                          features = as.character(unique(top10$gene)),   
                          group.by = "seurat_clusters",  
                          assay = "RNA",  
                          group.colors = mycolor1_clusters) +
  scale_fill_gradientn(colors = c("navy","white","firebrick3"))+
  theme(axis.text.y = element_text(size = 4))

heatmap_plot2 = DoHeatmap(object = seu_skin, label = F, 
                          features = as.character(unique(top10$gene)),   
                          group.by = "sample",  
                          assay = "RNA",  
                          group.colors = mycolor4_samples)+
  scale_fill_gradientn(colors = c("navy","white","firebrick3"))+
  theme(axis.text.y = element_text(size = 4))

heatmap_plot3 = DoHeatmap(object = seu_skin, label = F, 
                          features = as.character(unique(top5$gene)),   
                          group.by = "seurat_clusters",  
                          assay = "RNA",  
                          group.colors = mycolor1_clusters) +
  scale_fill_gradientn(colors = c("navy","white","firebrick3"))+
  theme(axis.text.y = element_text(size = 4))

heatmap_plot4 = DoHeatmap(object = seu_skin, label = F, 
                          features = as.character(unique(top5$gene)),   
                          group.by = "sample",  
                          assay = "RNA",  
                          group.colors = mycolor4_samples)+
  scale_fill_gradientn(colors = c("navy","white","firebrick3"))+
  theme(axis.text.y = element_text(size = 4))

ggsave("./Marker/top10_marker_of_each_cluster_heatmap.pdf", width = 8, height = 7,dpi = 300,
       plot = heatmap_plot1)
ggsave("./Marker/top10_marker_of_each_cluster_heatmap.png", width = 8, height = 7, dpi = 600,
       plot = heatmap_plot1)
ggsave("./Marker/top10_marker_of_each_sample_heatmap.pdf", width = 6, height = 8,dpi = 300,
       plot = heatmap_plot2)
ggsave("./Marker/top10_marker_of_each_sample_heatmap.png", width = 6, height = 8, dpi = 600,
       plot = heatmap_plot2)
ggsave("./Marker/top5_marker_of_each_cluster_heatmap.pdf", width = 8, height = 7,dpi = 300,
       plot = heatmap_plot3)
ggsave("./Marker/top5_marker_of_each_cluster_heatmap.png", width = 8, height = 7, dpi = 600,
       plot = heatmap_plot3)
ggsave("./Marker/top5_marker_of_each_sample_heatmap.pdf", width = 6, height = 8,dpi = 300,
       plot = heatmap_plot4)
ggsave("./Marker/top5_marker_of_each_sample_heatmap.png", width = 6, height = 8, dpi = 600,
       plot = heatmap_plot4)

dir.create("./Diffexp", showWarnings = FALSE, recursive = TRUE)
library(loupeR)

library(stringr)

seqs <- str_extract(colnames(seu_skin), "[ACGT]{10,}")
head(seqs)

new_names <- paste0(seqs, "-", as.numeric(as.factor(seu_skin$sample)))

if(any(duplicated(new_names))){
  stop("Error: duplicate barcodes generated!")
} else {
  message("Barcodes unique, renaming...")
}

seu_skin <- RenameCells(seu_skin, new.names = new_names)
loupeR::setup()
loupeR::create_loupe_from_seurat(seu_skin, output_name = "Skin_article")

markers <- c(
  "KLRB1", "GNLY", "CD3D", "TRAC", "IL7R", "GZMK", "CD8A", "CD8B",
  "FOXP3", "CTLA4", "TIGIT", "LAMP3", "CCR7", "CD40",
  "HLA-DRA", "HLA-DQB1", "CD14", "LYZ", "CD163", "C1QA",
  "DCT", "MLANA", "PMEL",
  "COL1A1", "DCN", "LUM", "CCL21", "TFPI", "PECAM1",
  "SPRR2G", "LCE3D", "LCE2C", "CDSN", "SPRR2A",
  "KRT1", "KRT10", "FABP5", "KRT14", "KRT5", "KRT15", "KRT6A"
)
markers <- unique(markers)

p2 <- DotPlot(seu_skin, features = markers, group.by = "seurat_clusters") + 
  RotatedAxis() +
  scale_color_gradientn(colours = c("white", "orange", "red"))
save_plot_auto (p2,"clusters_Dotplot")

p3 <- VlnPlot(seu_skin, features = markers, stack = TRUE, flip = TRUE) + NoLegend()
save_plot_auto (p3,"clusters_VlnPlot",width = 8.5)

new_cluster_ids <- c(
  "6"  = "CD4_T_cell",
  "1"  = "KC-S.Granulosum",
  "2"  = "KC-S.Corneum",
  "3"  = "Mature_DC",
  "4"  = "KC-S.Spinosum",
  "5"  = "KC-Wound/Prolif",
  "0"  = "Treg",
  "7"  = "CD8_T_cell",
  "8"  = "Macrophage",
  "9"  = "KC-S.Basale",
  "10" = "KC-S.Corneum",
  "11" = "Endothelial",
  "12" = "CD161_T_cell",
  "13" = "Melanocyte",
  "14" = "NK_cell",
  "15" = "KC-S.Corneum",
  "16" = "Mature_DC",
  "17" = "Mature_DC",
  "18" = "Mature_DC",
  "19" = "Fibroblast",
  "20" = "Mature_DC"
)

seu_skin <- RenameIdents(seu_skin, new_cluster_ids)
seu_skin$celltype <- Idents(seu_skin)

cell_type_levels <- c(
  "NK_cell", "CD161_T_cell", "CD4_T_cell", "CD8_T_cell", "Treg",
  "Mature_DC", "Macrophage", "Melanocyte", "Fibroblast", "Endothelial",
  "KC-S.Corneum", "KC-S.Granulosum", "KC-S.Spinosum", "KC-S.Basale", "KC-Wound/Prolif"
)
existing_cells <- unique(seu_skin@meta.data$celltype)
final_levels <- cell_type_levels[cell_type_levels %in% existing_cells]
seu_skin$celltype <- factor(seu_skin$celltype, levels = final_levels)
Idents(seu_skin) <- seu_skin$celltype

p2 <- DotPlot(seu_skin, features = markers, group.by = "celltype") + 
  RotatedAxis() + 
  scale_color_gradientn(colours = c("white", "orange", "red")) +
  labs(x = "Markers", y = "Cell Type")
save_plot_auto(p2,"annotion_dotplot")
p3 = DimPlot(seu_skin, reduction = "umap", label = TRUE) + NoLegend()
save_plot_auto(p3,"annotted_umap",width = 15, height = 18 )

saveRDS(seu_skin,"./seu_skin_2.8.rds")



seu_user <- readRDS("../../rds/Overall_shinanxi_New.rds")

gdt_order <- c(
  "Transitional activated",
  "Activation-regulated",
  "CSMD1+",
  "Cytotoxic effector",
  "Naive-like",
  "Terminal-branch",
  "NK-like",
  "Cytotoxic memory-like"
)

mycolor2 <- c(
  "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F","#FF6A6A","#7FFFD4", 
  "#AB82FF","#90EE90","#00CD00","#008B8B","#6495ED","#FFC1C1","#CD5C5C","#8B008B",
  "#FF3030", "#7CFC00","#000000","#708090"
)

if (is.null(seu_skin$celltype)) {
  seu_skin$celltype <- Idents(seu_skin)
}

Idents(seu_skin) <- seu_skin$seurat_clusters
seu_skin <- RenameIdents(seu_skin, new_cluster_ids)
seu_skin$celltype <- Idents(seu_skin)

print("seu_skin annotation mapping complete.")

seu_skin$group <- seu_skin$cell_id
seu_combined <- merge(seu_skin, y = seu_user, add.cell.ids = c("skin", "gdT"))

print("Data merge complete.")

seu_combined$celltype_old <- seu_combined$celltype

current_types <- as.character(seu_combined$celltype)

new_types <- case_when(
  current_types %in% gdt_order ~ current_types,
  current_types == "KC-S.Corneum"   ~ "KC-S.Corneum",
  current_types == "KC-S.Spinosum"   ~ "KC-S.Spinosum",
  current_types == "KC-S.Basale"     ~ "KC-S.Basale",
  current_types == "Inflammatory_KC" ~ "Inflammatory_KC",
  current_types == "CD4_T_cell" ~ "CD4_T_cell",
  current_types == "Treg"       ~ "Treg",
  current_types %in% c("CD8_T_cell", "NK_cell", "CD161_T_cell") ~ "Cytotoxic_NK_T",
  current_types %in% c("Mature_DC", "Macrophage","Semimature_DC") ~ "Myeloid",
  current_types %in% c( "Endothelial_Cell", "Melanocyte") ~ "Stroma",
  TRUE ~ current_types
)

seu_combined$celltype <- new_types

final_levels <- c(gdt_order, 
                  "KC-S.Corneum",  "KC-S.Spinosum", "KC-S.Basale", "Inflammatory_KC",
                  "CD4_T_cell", "Treg", "Cytotoxic_NK_T", "Myeloid", "Stroma")

existing_final <- unique(seu_combined$celltype)
final_levels <- final_levels[final_levels %in% existing_final]

seu_combined$celltype <- factor(seu_combined$celltype, levels = final_levels)
Idents(seu_combined) <- seu_combined$celltype

print(paste0("CellChat preparation complete. Number of cell subtypes: ", length(final_levels)))
print(table(seu_combined$celltype))

new_cell_types <- levels(seu_combined$celltype)
n_types <- length(new_cell_types)

if(length(mycolor2) < n_types){
  warning(paste0("Warning: insufficient colors! Need ", n_types, " but only ", length(mycolor2), " provided."))
}

mycolor2_named <- mycolor2[1:n_types]
names(mycolor2_named) <- new_cell_types



run_cellchat_pipeline <- function(seurat_obj, group_name) {
  message(paste0("Processing group: ", group_name))
  
  data.input <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  meta <- seurat_obj@meta.data
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
  
  CellChatDB <- CellChatDB.human
  cellchat@DB <- subsetDB(CellChatDB, search = "Secreted Signaling")
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  return(cellchat)
}

seurat_list <- SplitObject(seu_combined, split.by = "group")
Idents(seu_combined) <- "celltype"

cellchat_list <- list()

for (grp in names(seurat_list)) {
  
  seurat_obj <- seurat_list[[grp]]
  
  seurat_obj$celltype <- droplevels(as.factor(seurat_obj$celltype))
  
  Idents(seurat_obj) <- seurat_obj$celltype
  
  print(paste0(">>> Processing group: ", grp))
  print(paste("    Total cells:", ncol(seurat_obj)))
  print(paste("    Number of cell types:", length(levels(Idents(seurat_obj)))))
  
  if(ncol(seurat_obj) > 50 && length(levels(Idents(seurat_obj))) >= 2) {
    
    cellchat_list[[grp]] <- run_cellchat_pipeline(seurat_obj, grp)
    
  } else {
    print(paste("    [Skipped] Too few cells or single type"))
  }
}

cellchat_list <- list()
for (grp in names(seurat_list)) {
  if(ncol(seurat_list[[grp]]) > 50) {
    cellchat_list[[grp]] <- run_cellchat_pipeline(seurat_list[[grp]], grp)
  }
}

orig.seu_skin <- readRDS("./SingleCell_QC/After_QC.rds")
orig.seu_skin <- NormalizeData(orig.seu_skin, normalization.method = "LogNormalize", scale.factor = 10000)
orig.seu_skin <- FindVariableFeatures(orig.seu_skin, selection.method = "vst", nfeatures = 2000)
orig.seu_skin <- ScaleData(orig.seu_skin, features = VariableFeatures(orig.seu_skin)) 
orig.seu_skin <- RunPCA(orig.seu_skin, features = VariableFeatures(object = orig.seu_skin))

library(harmony)
library(SeuratWrappers)
library(batchelor)

seu_skin_20_MNN_old <- readRDS("./20/harmony/with_doublet_skin.rds")

Idents(seu_skin_20_MNN_old) <- "batch_group"

scRNAlist <- SplitObject(seu_skin_20_MNN_old, split.by = "batch_group")
common_genes <- Reduce(intersect, lapply(scRNAlist, rownames))
message(paste("Common genes count:", length(common_genes)))

scRNAlist <- lapply(scRNAlist, function(x) {
  x <- x[common_genes, ]
  x <- DietSeurat(x, counts = TRUE, data = TRUE, scale.data = FALSE) 
  return(x)
})

features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 2000)
seu_skin_20_MNN_old <- RunFastMNN(object.list = scRNAlist, features = features)

features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 2000)
seu_skin_20_MNN_old <- RunFastMNN(object.list = scRNAlist)
seu_skin_20_MNN_old <- FindVariableFeatures(seu_skin_20_MNN_old)
seu_skin_20_MNN_old <- RunUMAP(seu_skin_20_MNN_old, reduction = "mnn", dims = 1:20)
seu_skin_20_MNN_old <- FindNeighbors(seu_skin_20_MNN_old, reduction = "mnn", dims = 1:20)
seu_skin_20_MNN_old <- FindClusters(seu_skin_20_MNN_old,resolution = 0.85)

Idents(seu_skin_20_MNN_old) <- "seurat_clusters"
p1 = DimPlot(seu_skin_20_MNN_old, reduction = "umap", label = TRUE) + NoLegend()
save_plot_auto(p1,"./20/MNN_old/UMAP_without_annotion", width = 13, height = 13)

dir.create("./20/MNN_old/Diffexp", showWarnings = FALSE, recursive = TRUE)
library(loupeR)
library(stringr)

seqs <- str_extract(colnames(seu_skin_20_MNN_old), "[ACGT]{10,}")
head(seqs)

new_names <- paste0(seqs, "-", as.numeric(as.factor(seu_skin_20_MNN_old$sample)))

if(any(duplicated(new_names))){
  stop("Error: duplicate barcodes generated!")
} else {
  message("Barcodes unique, renaming...")
}

seu_skin_20_MNN_old <- RenameCells(seu_skin_20_MNN_old, new.names = new_names)
loupeR::setup()
loupeR::create_loupe_from_seurat(seu_skin_20_MNN_old, output_name = "./20/MNN_old/Skin_20_MNN")

markers <- c(
  "KLRB1", "GNLY", "CD3D", "TRAC", "TRBC1", "GZMH", "GZMK", "CD8A", "CD8B",
  "TIGIT", "IL2RA", "FOXP3", "CTLA4",
  "LAMP3", "LY75", "CIITA", "CD40", "HLA-DQA1", "HLA-DQB1",
  "HLA-DRB1", "HLA-DRA", "HLA-DRB5", "LYZ", "CD14", "CD163",
  "DCT", "TYRP1", "MLANA",
  "SPRR2G", "LCE3D", "CDSN", "FABP5", "KRT10", "KRT1", "KRT14", "KRT5"
)
markers <- unique(markers)

p2 <- DotPlot(seu_skin_20_MNN_old, features = markers, group.by = "seurat_clusters") + 
  RotatedAxis() +
  scale_color_gradientn(colours = c("white", "orange", "red"))
save_plot_auto (p2,"./20/MNN_old/clusters_Dotplot")

p3 <- VlnPlot(seu_skin_20_MNN_old, features = markers, stack = TRUE, flip = TRUE) + NoLegend()
save_plot_auto (p3,"./20/MNN_old/clusters_VlnPlot",width = 8.5)

Idents(seu_skin_20_MNN_old) = "seurat_clusters"
dir.create("./20/MNN_old/Marker")
all.markers = FindAllMarkers(seu_skin_20_MNN_old, 
                             min.pct = 0.25, 
                             logfc.threshold = 0.25, 
                             only.pos = TRUE)
head(all.markers)

top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top5 = all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

write.table(all.markers, 
            "./20/MNN_old/Marker/all_Markers_of_each_clusters.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")
write.table(top10, 
            "./20/MNN_old/Marker/top10_Markers_of_each_clusters.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")

seu_skin_20_MNN_old <- ScaleData(seu_skin_20_MNN_old, features = row.names(seu_skin_20_MNN_old))

heatmap_plot1 = DoHeatmap(object = seu_skin_20_MNN_old, label = F, 
                          features = as.character(unique(top10$gene)),   
                          group.by = "seurat_clusters",  
                          assay = "RNA",  
                          group.colors = mycolor1_clusters) +
  scale_fill_gradientn(colors = c("navy","white","firebrick3"))+
  theme(axis.text.y = element_text(size = 4))

heatmap_plot2 = DoHeatmap(object = seu_skin_20_MNN_old, label = F, 
                          features = as.character(unique(top10$gene)),   
                          group.by = "sample",  
                          assay = "RNA",  
                          group.colors = mycolor4_samples)+
  scale_fill_gradientn(colors = c("navy","white","firebrick3"))+
  theme(axis.text.y = element_text(size = 4))

heatmap_plot3 = DoHeatmap(object = seu_skin_20_MNN_old, label = F, 
                          features = as.character(unique(top5$gene)),   
                          group.by = "seurat_clusters",  
                          assay = "RNA",  
                          group.colors = mycolor1_clusters) +
  scale_fill_gradientn(colors = c("navy","white","firebrick3"))+
  theme(axis.text.y = element_text(size = 4))

heatmap_plot4 = DoHeatmap(object = seu_skin_20_MNN_old, label = F, 
                          features = as.character(unique(top5$gene)),   
                          group.by = "sample",  
                          assay = "RNA",  
                          group.colors = mycolor4_samples)+
  scale_fill_gradientn(colors = c("navy","white","firebrick3"))+
  theme(axis.text.y = element_text(size = 4))

ggsave("./20/MNN_old/Marker/top10_marker_of_each_cluster_heatmap.pdf", width = 8, height = 7,dpi = 300,
       plot = heatmap_plot1)
ggsave("./20/MNN_old/Marker/top10_marker_of_each_cluster_heatmap.png", width = 8, height = 7, dpi = 600,
       plot = heatmap_plot1)
ggsave("./20/MNN_old/Marker/top10_marker_of_each_sample_heatmap.pdf", width = 6, height = 8,dpi = 300,
       plot = heatmap_plot2)
ggsave("./20/MNN_old/Marker/top10_marker_of_each_sample_heatmap.png", width = 6, height = 8, dpi = 600,
       plot = heatmap_plot2)
ggsave("./20/MNN_old/Marker/top5_marker_of_each_cluster_heatmap.pdf", width = 8, height = 7,dpi = 300,
       plot = heatmap_plot3)
ggsave("./20/MNN_old/Marker/top5_marker_of_each_cluster_heatmap.png", width = 8, height = 7, dpi = 600,
       plot = heatmap_plot3)
ggsave("./20/MNN_old/Marker/top5_marker_of_each_sample_heatmap.pdf", width = 6, height = 8,dpi = 300,
       plot = heatmap_plot4)
ggsave("./20/MNN_old/Marker/top5_marker_of_each_sample_heatmap.png", width = 6, height = 8, dpi = 600,
       plot = heatmap_plot4)

new_annotations_supplement <- c(
  "16" = "NK_cell",
  "11" = "CD161_T_cell",
  "18" = "CD8_T_cell",
  "3"  = "CD4_T_cell",
  "0"  = "Treg",
  "17" = "Treg",
  "20" = "Mature_DC",
  "9"  = "Mature_DC",
  "5"  = "Mature_DC",
  "8"  = "Semimature_DC",
  "22" = "Macrophage",
  "15" = "Endothelial_Cell",
  "12" = "Melanocyte",
  "21" = "Melanocyte",
  "7"  = "Inflammatory_KC",
  "6"  = "KC-S.Basale",
  "10" = "KC-S.Spinosum",
  "1"  = "KC-S.Corneum",
  "2"  = "KC-S.Corneum",
  "4"  = "KC-S.Corneum",
  "13" = "KC-S.Corneum",
  "14" = "KC-S.Corneum",
  "19" = "Doublet",
  "23" = "Remove"
)

current_ids <- levels(seu_skin_20_MNN_old)
missing_ids <- setdiff(current_ids, names(new_annotations_supplement))
if(length(missing_ids) > 0) {
  message("Warning: following clusters unassigned, marking as Unidentified: ", paste(missing_ids, collapse=", "))
  temp_assignment <- rep("Unidentified", length(missing_ids))
  names(temp_assignment) <- missing_ids
  new_annotations_supplement <- c(new_annotations_supplement, temp_assignment)
}

Idents(seu_skin_20_MNN_old) <- "seurat_clusters"
seu_skin_20_MNN_old <- RenameIdents(seu_skin_20_MNN_old, new_annotations_supplement)

seu_skin_final <- subset(seu_skin_20_MNN_old, idents = c("Remove", "Doublet"), invert = TRUE)
seu_skin_final$celltype <- NULL

final_levels <- c(
  "NK_cell", "CD161_T_cell", "CD8_T_cell", "CD4_T_cell", "Treg",
  "Mature_DC", "Semimature_DC", "Macrophage",
  "Melanocyte", "Endothelial_Cell", 
  "KC-S.Corneum", "Inflammatory_KC", "KC-S.Spinosum", "KC-S.Basale"
)
seu_skin_final$celltype <- Idents(seu_skin_final)
final_levels <- final_levels[final_levels %in% unique(Idents(seu_skin_final))]
Idents(seu_skin_final) <- factor(Idents(seu_skin_final), levels = final_levels)
seu_skin_final$celltype <- droplevels(factor(seu_skin_final$celltype))
Idents(seu_skin_final) <- seu_skin_final$celltype

print(levels(Idents(seu_skin_final)))
seu_skin_final$celltype <- factor(seu_skin_final$celltype, levels = final_levels)
Idents(seu_skin_final) <- seu_skin_final$celltype

DefaultAssay(seu_skin_final) <- "RNA"
features.supplement <- c(
  "KLRB1", "GNLY", "CD3D", "TRAC", "TRBC1", "GZMH", "GZMK", "CD8A", "CD8B",
  "TIGIT", "IL2RA", "FOXP3", "CTLA4",
  "LAMP3", "LY75", "CIITA", "CD40", "HLA-DQA1", "HLA-DQB1",
  "HLA-DRB1", "HLA-DRA", "HLA-DRB5", "LYZ", "CD14", "CD163",
  "DCT", "TYRP1", "MLANA", "CCL21", "TFF3",
  "SPRR2G", "LCE3D", "CCL20", "CXCL8", "FABP5", "KRT10", "KRT1", "KRT14", "KRT5"
)

p4 <- DotPlot(seu_skin_final, features = features.supplement, group.by = "celltype") + 
  RotatedAxis() +
  scale_color_gradientn(colours = c("white", "orange", "red"))
save_plot_auto (p4,"./20/MNN_old/celltype_Dotplot")

saveRDS(seu_skin_final, "./Skin_Final_Annotated_Supplemented.rds")



seu_skin <- readRDS( "./Skin_Final_Annotated_Supplemented.rds")

seu_user <- readRDS("../../rds/Overall_shinanxi_New.rds")

gdt_order <- c(
  "Transitional activated",
  "Activation-regulated",
  "CSMD1+",
  "Cytotoxic effector",
  "Naive-like",
  "Terminal-branch",
  "NK-like",
  "Cytotoxic memory-like"
)

mycolor2 <- c(
  "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F","#FF6A6A","#7FFFD4", 
  "#AB82FF","#90EE90","#00CD00","#008B8B","#6495ED","#FFC1C1","#CD5C5C","#8B008B",
  "#FF3030", "#7CFC00","#000000","#708090"
)

if (is.null(seu_skin$celltype)) {
  seu_skin$celltype <- Idents(seu_skin)
}

print("seu_skin annotation mapping complete.")

seu_skin$group <- seu_skin$cell_id
seu_combined <- merge(seu_skin, y = seu_user, add.cell.ids = c("skin", "gdT"))

print("Data merge complete.")

seu_combined$celltype_old <- seu_combined$celltype

current_types <- as.character(seu_combined$celltype)

new_types <- case_when(
  current_types %in% gdt_order ~ current_types,
  current_types == "KC-S.Corneum"   ~ "KC-S.Corneum",
  current_types == "KC-S.Spinosum"   ~ "KC-S.Spinosum",
  current_types == "KC-S.Basale"     ~ "KC-S.Basale",
  current_types == "Inflammatory_KC" ~ "Inflammatory_KC",
  current_types == "CD4_T_cell" ~ "CD4_T_cell",
  current_types == "Treg"       ~ "Treg",
  current_types %in% c("CD8_T_cell", "NK_cell", "CD161_T_cell") ~ "Cytotoxic_NK_T",
  current_types %in% c("Mature_DC", "Macrophage","Semimature_DC") ~ "Myeloid",
  current_types %in% c( "Endothelial_Cell", "Melanocyte") ~ "Stroma",
  TRUE ~ current_types
)

seu_combined$celltype <- new_types

final_levels <- c(gdt_order, 
                  "KC-S.Corneum",  "KC-S.Spinosum", "KC-S.Basale", "Inflammatory_KC",
                  "CD4_T_cell", "Treg", "Cytotoxic_NK_T", "Myeloid", "Stroma")

existing_final <- unique(seu_combined$celltype)
final_levels <- final_levels[final_levels %in% existing_final]

seu_combined$celltype <- factor(seu_combined$celltype, levels = final_levels)
Idents(seu_combined) <- seu_combined$celltype

print(paste0("CellChat preparation complete. Number of cell subtypes: ", length(final_levels)))
print(table(seu_combined$celltype))

new_cell_types <- levels(seu_combined$celltype)
n_types <- length(new_cell_types)

if(length(mycolor2) < n_types){
  warning(paste0("Warning: insufficient colors! Need ", n_types, " but only ", length(mycolor2), " provided."))
}

mycolor2_named <- mycolor2[1:n_types]
names(mycolor2_named) <- new_cell_types

run_cellchat_pipeline <- function(seurat_obj, group_name) {
  message(paste0("Processing group: ", group_name))
  
  data.input <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  meta <- seurat_obj@meta.data
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
  
  CellChatDB <- CellChatDB.human
  cellchat@DB <- subsetDB(CellChatDB, search = "Secreted Signaling")
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  return(cellchat)
}

seurat_list <- SplitObject(seu_combined, split.by = "group")
Idents(seu_combined) <- "celltype"

cellchat_list <- list()

for (grp in names(seurat_list)) {
  
  seurat_obj <- seurat_list[[grp]]
  
  seurat_obj$celltype <- droplevels(as.factor(seurat_obj$celltype))
  
  Idents(seurat_obj) <- seurat_obj$celltype
  
  print(paste0(">>> Processing group: ", grp))
  print(paste("    Total cells:", ncol(seurat_obj)))
  print(paste("    Number of cell types:", length(levels(Idents(seurat_obj)))))
  
  if(ncol(seurat_obj) > 50 && length(levels(Idents(seurat_obj))) >= 2) {
    
    cellchat_list[[grp]] <- run_cellchat_pipeline(seurat_obj, grp)
    
  } else {
    print(paste("    [Skipped] Too few cells or single type"))
  }
}

cellchat_list <- list()
for (grp in names(seurat_list)) {
  if(ncol(seurat_list[[grp]]) > 50) {
    cellchat_list[[grp]] <- run_cellchat_pipeline(seurat_list[[grp]], grp)
  }
}



save_plot_auto_2.5 <- function(plot_code_or_obj, filename_base, width = 18, height = 9.5, type = "both") {
  if(!dir.exists("cellchat")) dir.create("cellchat", recursive = TRUE)
  f_png <- paste0("cellchat/", filename_base, ".png")
  f_pdf <- paste0("cellchat/", filename_base, ".pdf")
  
  draw_it <- function(obj) {
    if (inherits(obj, "grob")) {
      grid::grid.draw(obj)
    } else if (is.expression(obj)) {
      eval(obj)
    } else {
      print(obj)
    }
  }
  
  if(type %in% c("both", "png")) {
    png(f_png, width = width, height = height, units = "in", res = 600)
    grid::grid.newpage() 
    tryCatch({
      draw_it(plot_code_or_obj)
    }, error = function(e) message(paste("Error saving PNG:", filename_base, e$message)))
    dev.off()
  }
  
  if(type %in% c("both", "pdf")) {
    pdf(f_pdf, width = width, height = height)
    grid::grid.newpage()
    tryCatch({
      draw_it(plot_code_or_obj)
    }, error = function(e) message(paste("Error saving PDF:", filename_base, e$message)))
    dev.off()
  }
  
  message(paste("Saved successfully:", filename_base))
}

group.cellType <- union(levels(cellchat_list[[1]]@idents), levels(cellchat_list[[2]]@idents))
for (i in 1:length(cellchat_list)) {
  cellchat_list[[i]] <- liftCellChat(cellchat_list[[i]], group.new = group.cellType)
}
cellchat_merged <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))
compare_vec <- c(2, 1)

message("--- Merging CellChat objects ---")
cellchat_merged <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))

print("Merged cell type order:")
print(levels(cellchat_merged@idents$joint))

mycolor2 <- c(
  "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F","#FF6A6A","#7FFFD4", 
  "#AB82FF","#90EE90","#00CD00","#008B8B","#6495ED","#FFC1C1","#CD5C5C","#8B008B",
  "#FF3030", "#7CFC00","#000000","#708090"
)

cell_types <- levels(cellchat_merged@meta$celltype)
print(paste("Number of cell types:", length(cell_types)))
print(paste("Provided colors:", length(mycolor2)))

mycolor2_named <- mycolor2[1:length(cell_types)]
names(mycolor2_named) <- cell_types

all_pathways_1 <- cellchat_list[[1]]@netP$pathways
keep_pathways_1 <- all_pathways_1[!grepl("HLA|MHC", all_pathways_1, ignore.case = TRUE)]

all_pathways_2 <- cellchat_list[[2]]@netP$pathways
keep_pathways_2 <- all_pathways_2[!grepl("HLA|MHC", all_pathways_2, ignore.case = TRUE)]

print(head(cellchat_list[[1]]@netP$pathways))
print(head(cellchat_list[[2]]@netP$pathways))

p1 <- cellchat_list[[1]]@netP$pathways
keep_pathways_1 <- p1[!grepl("HLA|MHC", p1, ignore.case = TRUE)]

p2 <- cellchat_list[[2]]@netP$pathways
keep_pathways_2 <- p2[!grepl("HLA|MHC", p2, ignore.case = TRUE)]

mat_control <- cellchat_list[[1]]@net$count
mat_disease <- cellchat_list[[2]]@net$count

library(pheatmap)
library(grid)
library(gridExtra)

global_max <- max(max(mat_control, na.rm = T), max(mat_disease, na.rm = T))
common_breaks <- seq(0, global_max, length.out = 100)
common_color <- colorRampPalette(c("white", "#E41A1C"))(100)

if(!dir.exists("cellchat")) dir.create("cellchat")

p3_ctrl <- pheatmap::pheatmap(cellchat_list[[1]]@net$count, border_color = "black",
                              cluster_cols = F, cluster_rows = F, fontsize = 10,
                              display_numbers = T, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Control", silent = TRUE)$gtable

p3_case <- pheatmap::pheatmap(mat_disease, border_color = "black",
                              cluster_cols = F, cluster_rows = F, fontsize = 10,
                              display_numbers = T, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Psoriasis", silent = TRUE)$gtable

p3_combined <- arrangeGrob(p3_ctrl, p3_case, ncol = 2, top = textGrob("Interaction Count (Raw)", gp=gpar(fontsize=15, font=2)))
save_plot_auto_2.5(p3_combined, "CellChat_Heatmap_NoCluster_WithNum")

p4_ctrl <- pheatmap::pheatmap(mat_control, border_color = "black",
                              cluster_cols = T, cluster_rows = T, fontsize = 10,
                              display_numbers = T, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Control", silent = TRUE)$gtable

p4_case <- pheatmap::pheatmap(mat_disease, border_color = "black",
                              cluster_cols = T, cluster_rows = T, fontsize = 10,
                              display_numbers = T, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Psoriasis", silent = TRUE)$gtable

p4_combined <- arrangeGrob(p4_ctrl, p4_case, ncol = 2, top = textGrob("Interaction Count (Clustered)", gp=gpar(fontsize=15, font=2)))
save_plot_auto_2.5(p4_combined, "CellChat_Heatmap_Clustered_WithNum")

p5_ctrl <- pheatmap::pheatmap(mat_control, border_color = "black",
                              cluster_cols = F, cluster_rows = F, fontsize = 10,
                              display_numbers = F, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Control", silent = TRUE)$gtable

p5_case <- pheatmap::pheatmap(mat_disease, border_color = "black",
                              cluster_cols = F, cluster_rows = F, fontsize = 10,
                              display_numbers = F, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Psoriasis", silent = TRUE)$gtable

p5_combined <- arrangeGrob(p5_ctrl, p5_case, ncol = 2, top = textGrob("Interaction Count (Raw, No Num)", gp=gpar(fontsize=15, font=2)))
save_plot_auto_2.5(p5_combined, "CellChat_Heatmap_NoCluster_NoNum")

p6_ctrl <- pheatmap::pheatmap(mat_control, border_color = "black",
                              cluster_cols = T, cluster_rows = T, fontsize = 10,
                              display_numbers = F, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Control", silent = TRUE)$gtable

p6_case <- pheatmap::pheatmap(mat_disease, border_color = "black",
                              cluster_cols = T, cluster_rows = T, fontsize = 10,
                              display_numbers = F, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Psoriasis", silent = TRUE)$gtable

p6_combined <- arrangeGrob(p6_ctrl, p6_case, ncol = 2, top = textGrob("Interaction Count (Clustered, No Num)", gp=gpar(fontsize=15, font=2)))
save_plot_auto_2.5(p6_combined, "CellChat_Heatmap_Clustered_NoNum")

mat <- cellchat_merged@net$count
cell_types <- levels(cellchat_merged@idents$joint)
groupSize <- as.numeric(table(cellchat_merged@idents$joint))

code_single_circle <- expression({
  par(mfrow = c(3,3), xpd=TRUE, mar = c(1,1,1,1))
  for (i in 1:length(cell_types)) {
    netVisual_circle(cellchat_list[["Psoriasis"]]@net$count, vertex.weight = groupSize, weight.scale = T, 
                     edge.weight.max = max(cellchat_list[["Psoriasis"]]@net$count), 
                     title.name = paste("Pso -", cell_types[i]), color.use = mycolor2, sources.use = cell_types[i])
  }
})
save_plot_auto_2.5(code_single_circle, "03_Single_CellType_Source_Circle_Psoriasis")

code_single_circle <- expression({
  par(mfrow = c(3,3), xpd=TRUE, mar = c(1,1,1,1))
  for (i in 1:length(cell_types)) {
    netVisual_circle(cellchat_list[["Control"]]@net$count, vertex.weight = groupSize, weight.scale = T, 
                     edge.weight.max = max(cellchat_list[["Control"]]@net$count), 
                     title.name = paste("Pso -", cell_types[i]), color.use = mycolor2, sources.use = cell_types[i])
  }
})
save_plot_auto_2.5(code_single_circle, "03_Single_CellType_Source_Circle_Control")

groupSize <- as.numeric(table(cellchat_list[[1]]@idents))

p1 <- expression({
  par(mfrow = c(1,2), xpd=TRUE)
  
  netVisual_circle(cellchat_list[[1]]@net$count, vertex.weight = groupSize, weight.scale = T,
                   label.edge= F,                    arrow.width = 1.2,
                   arrow.size = 0.24, color.use = mycolor2,
                   title.name = "Number of interactions in Cont")
  netVisual_circle(cellchat_list[[1]]@net$weight, vertex.weight = groupSize, weight.scale = T, 
                   label.edge= F, arrow.width = 1.2,  color.use = mycolor2,
                   arrow.size = 0.24 ,title.name = "Interaction weights/strength in Cont")
  
  groupSize <- as.numeric(table(cellchat_list[[2]]@idents))})
save_plot_auto_2.5(p1, "03_Total_CellType_Source_Circle_Control")

p2 <- expression({
  par(mfrow = c(1,2), xpd=TRUE)
  
  netVisual_circle(cellchat_list[[2]]@net$count, vertex.weight = groupSize, weight.scale = T,
                   label.edge= F, 
                   arrow.width = 1.2,
                   arrow.size = 0.24, color.use = mycolor2,
                   title.name = "Number of interactions in Pso")
  netVisual_circle(cellchat_list[[2]]@net$weight, vertex.weight = groupSize, weight.scale = T, 
                   label.edge= F, 
                   arrow.width = 1.2,
                   arrow.size = 0.24, color.use = mycolor2,
                   title.name = "Interaction weights/strength in Pso")})
save_plot_auto_2.5(p2, "03_Total_CellType_Source_Circle_Psoriasis")

cellchat_list[[1]]@netP$pathways

pp1 <- rankNet(cellchat_merged, mode = "comparison", stacked = T, do.stat = TRUE)
pp2 <- rankNet(cellchat_merged, mode = "comparison", stacked = F, do.stat = TRUE)
save_plot_auto_2.5(pp1 + pp2, "03_rankNet")

target_groups   <- c("Control", "Psoriasis")
chat_ctrl <- cellchat_list[[1]] 
chat_pso  <- cellchat_list[[2]]

pathway_name <- "CLEC"
obj_use <- chat_ctrl
prefix <- "Control_CLEC"

p1 = netVisual_aggregate(obj_use, signaling = pathway_name, color.use = mycolor2)
save_plot_auto_2.5(p1, paste0(prefix, "_circle"))

p3 = netVisual_heatmap(obj_use, signaling = pathway_name, color.heatmap = "Reds", color.use = mycolor2_named)
save_plot_auto_2.5(p3, paste0(prefix, "_heatmap"))

p4_expr = expression({
  netVisual_chord_cell(obj_use, signaling = pathway_name, group = group.cellType, 
                       title.name = paste0(pathway_name, " (Control) - Immune vs Other"), color.use = mycolor2_named)
})
save_plot_auto_2.5(p4_expr, paste0(prefix, "_grouped_chord"))

p5 = netAnalysis_contribution(obj_use, signaling = pathway_name)
save_plot_auto_2.5(p5, paste0(prefix, "_contribution"))

pairLR_use <- extractEnrichedLR(obj_use, signaling = pathway_name, geneLR.return = FALSE)[1]
p6_expr = expression({
  netVisual_individual(obj_use, signaling = pathway_name, pairLR.use = pairLR_use, layout = "circle", color.use = mycolor2)
})
save_plot_auto_2.5(p6_expr, paste0(prefix, "_", pairLR_use, "_interaction"))

state_list <- list(
  "Control"   = cellchat_list[[1]],
  "Psoriasis" = cellchat_list[[2]]
)

target_pathways <- c("MIF")

for (current_state in names(state_list)) {
  
  obj_use <- state_list[[current_state]]
  
  for (pathway_name in target_pathways) {
    
    if (!(pathway_name %in% obj_use@netP$pathways)) {
      message(paste("Skip:", pathway_name, "not found in", current_state))
      next
    }
    
    prefix <- paste0(current_state, "_", pathway_name)
    
    message(paste("Processing:", prefix))
    
    p1 = netVisual_aggregate(obj_use, signaling = pathway_name, color.use = mycolor2)
    save_plot_auto_2.5(p1, paste0(prefix, "_circle"))
    
    p3 = netVisual_heatmap(obj_use, signaling = pathway_name, color.heatmap = "Reds", color.use = mycolor2_named)
    save_plot_auto_2.5(p3, paste0(prefix, "_heatmap"))
    
    group.cellType <- c(
      "Transitional activated"    = "gdT",
      "Cytotoxic effector"  = "gdT",
      "CSMD1+"   = "gdT",
      "NK-like"   = "gdT",
      "Activation-regulated"         = "gdT",
      "Naive-like"   = "gdT",
      "Cytotoxic memory-like" = "gdT",
      "Terminal-branch"         = "gdT",
      "gdT"               = "gdT",
      
      "KC-S.Corneum"             = "Other",
      "KC-S.Granulosum"             = "Other",
      "KC-S.Spinosum"           = "Other",
      "KC-S.Basale"           = "Other",
      "KC-Wound/Prolif"              = "Other",
      "CD4_T_cell"          = "Other",
      "Treg"     = "Other",
      "Cytotoxic_NK_T"              = "Other",
      "Myeloid"              = "Other",
      "Stroma "              = "Other"
    )
    
    color_grouped <- c("gdT" = "#FF34B3", "Other" = "grey80")
    p4_expr = expression({
      netVisual_chord_cell(obj_use, signaling = pathway_name, group = group.cellType,
                           title.name = paste0(pathway_name, " (", current_state, ") - gdT vs Other"), color.use = mycolor2_named)
    })
    save_plot_auto_2.5(p4_expr, paste0(prefix, "_grouped_chord"))
    
    p5 = netAnalysis_contribution(obj_use, signaling = pathway_name)
    save_plot_auto_2.5(p5, paste0(prefix, "_contribution"))
    
    tryCatch({
      pairLR_use <- extractEnrichedLR(obj_use, signaling = pathway_name, geneLR.return = FALSE)[1]
      
      if (!is.na(pairLR_use)) {
        p6_expr = expression({
          netVisual_individual(obj_use, signaling = pathway_name, pairLR.use = pairLR_use, layout = "circle", color.use = mycolor2)
        })
        save_plot_auto_2.5(p6_expr, paste0(prefix, "_", pairLR_use, "_interaction"))
      }
    }, error = function(e) { message(paste("Error plotting individual LR for", prefix)) })
    
    p7 = netVisual_bubble(obj_use, remove.isolate = FALSE,signaling = pathway_name)
    save_plot_auto_2.5(p7, paste0(prefix, "_sender_receiver"))
    
    p8 <- expression(netAnalysis_signalingRole_network(obj_use, signaling = pathway_name, width = 8, height = 2.5, font.size = 10,color.use = mycolor2_named))
    save_plot_auto_2.5(p8, paste0(prefix, "_charactor"))
    
    p9 <- netAnalysis_signalingRole_scatter(obj_use, signaling = pathway_name,color.use = mycolor2)
    save_plot_auto_2.5(p9, paste0(prefix, "_signalingRole_scatte"))
  }
}

CSMD1_chord <- expression(netVisual_chord_gene(cellchat_list[[1]], sources.use = "CSMD1+", lab.cex = 0.5,legend.pos.y =50))
save_plot_auto_2.5(CSMD1_chord,"CSMD1_chord")

for (i in 1:length(cellchat_list)) {
  cellchat_list[[i]] <- netAnalysis_computeCentrality(cellchat_list[[i]], slot.name = "netP")
}
ht_out_control <- netAnalysis_signalingRole_heatmap(cellchat_list[[1]], pattern = "outgoing", width=12, height=14, color.use = mycolor2_named,signaling = keep_pathways_1)
ht_in_control <- netAnalysis_signalingRole_heatmap(cellchat_list[[1]], pattern = "incoming", width=12, height=14, color.use = mycolor2_named,signaling = keep_pathways_1)
save_plot_auto_2.5(ht_out_control + ht_in_control, "SignalingRole_Heatmap_control")

ht_out_PSO <- netAnalysis_signalingRole_heatmap(cellchat_list[[2]], pattern = "outgoing", width=12, height=14, color.use = mycolor2_named,signaling = keep_pathways_2)
ht_in_PSO <- netAnalysis_signalingRole_heatmap(cellchat_list[[2]], pattern = "incoming", width=12, height=14, color.use = mycolor2_named,signaling = keep_pathways_2)
save_plot_auto_2.5(ht_out_PSO + ht_in_PSO, "SignalingRole_Heatmapr_PSO")

p1 <- netVisual_aggregate(
  cellchat_list[[1]], 
  signaling = keep_pathways_1, 
  layout = "chord",
  signaling.name = "Control",
  title.name = "NULL" 
)
save_plot_auto_2.5(p1, "_Cont_circle")

p2 <- netVisual_aggregate(
  cellchat_list[[2]],
  signaling = keep_pathways_2, 
  color.use = mycolor2, 
  layout = "chord", 
  signaling.name = "Psoriasis",
  title.name = "NULL"
)
save_plot_auto_2.5(p2, "_Pso_circle")

gg1 <- netAnalysis_signalingRole_scatter(cellchat_list[[1]],color.use = mycolor2,signaling = keep_pathways_1)
save_plot_auto_2.5(gg1,"netAnalysis_Cont")
gg2 <- netAnalysis_signalingRole_scatter(cellchat_list[[2]],color.use = mycolor2,signaling = keep_pathways_2)
save_plot_auto_2.5(gg2,"netAnalysis_Pso")

bubble_Cont = netVisual_bubble(cellchat_list[[1]], remove.isolate = FALSE,signaling = keep_pathways_1 )
save_plot_auto_2.5(bubble_Cont,"bubble_Cont", width = 30)
bubble_Pso = netVisual_bubble(cellchat_list[[2]], remove.isolate = FALSE,signaling = keep_pathways_2)
save_plot_auto_2.5(bubble_Pso,"bubble_Pso", width = 30)

NMF_control <- selectK(cellchat_list[[1]], pattern = "outgoing")
NMF_PSO <-  selectK(cellchat_list[[2]], pattern = "outgoing")
save_plot_auto_2.5(NMF_control + NMF_PSO, "17_NMF_PSO")

k_c = 5
k_p = 2
cellchat_list[[1]] <- identifyCommunicationPatterns(cellchat_list[[1]], pattern = "outgoing", k = 5)
cellchat_list[[2]] <- identifyCommunicationPatterns(cellchat_list[[2]], pattern = "outgoing", k = 2)

png("./cellchat/Pattern_River_Outgoing_Cont.png", width = 8, height = 8,res = 600, units = "in")
netAnalysis_river(cellchat_list[[1]], pattern = "outgoing")
dev.off()
pdf("./cellchat/Pattern_River_Outgoing_Cont.pdf", width = 8, height = 8)
netAnalysis_river(cellchat_list[[1]], pattern = "outgoing")
dev.off()

png("./cellchat/Pattern_Dot_Outgoing_Cont.png", width = 8, height = 8,res = 600, units = "in")
netAnalysis_dot(cellchat_list[[1]], pattern = "outgoing")
dev.off()
pdf("./cellchat/Pattern_Dot_Outgoing_Cont.pdf", width = 8, height = 8)
netAnalysis_dot(cellchat_list[[1]], pattern = "outgoing")
dev.off()

png("./cellchat/Pattern_River_Outgoing_Pso.png", width = 8, height = 8,res = 600, units = "in")
netAnalysis_river(cellchat_list[[2]], pattern = "outgoing")
dev.off()
pdf("./cellchat/Pattern_River_Outgoing_Pso.pdf", width = 8, height = 8)
netAnalysis_river(cellchat_list[[2]], pattern = "outgoing")
dev.off()

png("./cellchat/Pattern_River_Outgoing_Pso.png", width = 8, height = 8,res = 600, units = "in")
netAnalysis_river(cellchat_list[[2]], pattern = "outgoing")
dev.off()
pdf("./cellchat/Pattern_River_Outgoing_Pso.pdf", width = 8, height = 8)
netAnalysis_river(cellchat_list[[2]], pattern = "outgoing")
dev.off()

plot_cellchat_LR <- function(cellchat_obj, 
                             celltype_inter, 
                             celltype_color,
                             ligand_col = "#FB8072",
                             receptor_col = "#1988B0", 
                             text_size,
                             legend_in_plot = T,
                             Group,
                             top_n) {
  
  library(tidyverse)
  library(circlize)
  library(ggsci)
  library(igraph)
  library(gtools)
  library(ComplexHeatmap)
  library(stringr)
  
  net.df <- subsetCommunication(cellchat_obj)
  net.df$source <- str_remove_all(net.df$source, "[^[:alnum:][:space:]]")
  net.df$cell_type_pair <- paste0(net.df$source, '_', net.df$target)
  
  net.df_filtered <- net.df %>%
    filter(source %in% celltype_inter & target %in% celltype_inter) %>%
    group_by(cell_type_pair) %>%
    slice_max(prob, n = top_n, with_ties = FALSE) %>%
    ungroup()
  
  if (nrow(net.df_filtered) == 0) {
    message("No interactions found after filtering")
    return()
  }
  
  node_v1 <- vector()
  node_v2 <- vector()
  weight_v <- vector()
  
  for (i in 1:nrow(net.df_filtered)) {
    cci <- strsplit(net.df_filtered$interaction_name_2[i], split = ' ')[[1]]
    ct.pair <- strsplit(net.df_filtered$cell_type_pair[i], split = '_')[[1]]
    prefix1 <- ct.pair[1]
    prefix2 <- ct.pair[2]
    
    node1 <- paste0(prefix1, '_', cci[1], '_L')  
    
    id2 <- cci[length(cci)]
    weight <- net.df_filtered$prob[i]
    
    if (length(grep('\\+', id2)) > 0) {
      array <- strsplit(id2, split = '\\+')[[1]]
      node2 <- paste0(prefix2, '_', substr(array[1], 2, nchar(array[1])), '_R')
      node3 <- paste0(prefix2, '_', substr(array[2], 1, nchar(array[2]) - 1), '_R')
      node_v1 <- c(node_v1, node1, node1)
      node_v2 <- c(node_v2, node2, node3)
      weight_v <- c(weight_v, weight, weight)
    } else {
      node2 <- paste0(prefix2, '_', id2, '_R')
      node_v1 <- c(node_v1, node1)
      node_v2 <- c(node_v2, node2)
      weight_v <- c(weight_v, weight)
    }
  }
  
  if (length(node_v1) <= 0) { return() }
  
  g_df <- data.frame(node_v1, node_v2, weight_v)
  g <- graph.data.frame(g_df, directed = T)
  E(g)$weight <- g_df[[3]]
  adj <- get.adjacency(g, attr = 'weight')
  
  graph_adj <- as.data.frame(as.matrix(adj))
  Genes <- names(V(g))
  
  arrays <- strsplit(Genes, split = '_')
  ID <- sapply(arrays, function(x) {
    cell_gene <- paste(x[-length(x)], collapse = "_")
    return(strsplit(cell_gene, split = '_')[[1]][1])
  })
  
  category <- data.frame(Genes = Genes, ID = ID)
  graph_module <- category
  graph_module$LR <- ifelse(grepl("_L$", graph_module$Genes), "ligand", "receptor")
  
  if(is.null(celltype_inter)==F){
    graph_module <- graph_module[graph_module$ID %in% celltype_inter,]
    str <- paste0("^(", paste(celltype_inter, collapse = "|"), ")")
    str_cells <- grepl(str, rownames(graph_adj))
    str_cells1 <- rownames(graph_adj)[str_cells]
    str_cells2 <- colnames(graph_adj)[str_cells]
    
    graph_adj <- graph_adj[str_cells1,]
    graph_adj <- graph_adj[,str_cells2]
  }
  
  g <- graph.adjacency(as.matrix(graph_adj), weighted = T)
  
  LR_color <- data.frame(LR = c("ligand", "receptor"),
                         color1 = c(ligand_col, receptor_col))
  graph_module <- left_join(graph_module, LR_color, by = "LR")
  
  if (is.vector(celltype_color) && !is.null(names(celltype_color))) {
    celltype_color <- data.frame(ID = names(celltype_color), 
                                 color2 = as.character(celltype_color),
                                 stringsAsFactors = FALSE)
  } else {
    celltype_color <- as.data.frame(celltype_color)
    colnames(celltype_color) <- c("ID", 'color2')
  }
  
  graph_module <- left_join(graph_module, celltype_color, by = "ID")
  
  raw_edges <- as.data.frame(cbind(get.edgelist(g), E(g)$weight)) %>%
    mutate(
      V1 = gsub('\\.', '-', V1),
      V2 = gsub('\\.', '-', V2),
      V3 = as.numeric(V3),
      V4 = 1
    )
  edges <- raw_edges %>% arrange(V3)
  
  nodes <- unique(c(edges$V1, edges$V2))
  sectors <- sort(unique(c(raw_edges$V1, raw_edges$V2)))
  
  col_fun = colorRamp2(range(edges$V3), c("#FFFDE7", "#013220"))
  
  circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(-0.15, 0.2))
  circos.initialize(sectors, xlim = c(0, 1))
  circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05, bg.border = NA)
  
  circos.track(
    track.index = 1,
    panel.fun = function(x, y) {
      sector.name = get.cell.meta.data("sector.index")
      xlim = get.cell.meta.data("xlim")
      
      display_name <- gsub("_(L|R)$", "", sector.name)
      
      node_text_color <- graph_module %>%
        dplyr::filter(Genes == sector.name) %>%
        pull(color2) %>%
        as.character()
      
      if(length(node_text_color) == 0 || is.na(node_text_color)) node_text_color <- "black"
      
      node_LR_color <- graph_module %>%
        dplyr::filter(Genes == sector.name) %>%
        pull(color1) %>%
        as.character()
      
      circos.rect(
        xlim[1], 0, xlim[2], 1,
        col = node_LR_color,
        border = NA
      )
      
      circos.text(
        mean(xlim), 2,
        display_name, 
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        col = node_text_color, 
        cex = text_size
      )
    },
    bg.border = NA
  ) 
  
  for (i in seq_len(nrow(edges))) {
    link <- edges[i,]
    circos.link(link[[1]], c(0, 1), link[[2]], c(0, 1),
                col = col_fun(link[[3]]), border = NA)
  }
  
  title(paste0("Cell interactions in ", Group))
  
  if(legend_in_plot == T){
    lgd <- Legend(title = "Score", col_fun = col_fun, direction = "horizontal", border = "black")
    grid.draw(lgd)
    
    legend(1.1, 0.5, pch = 15, legend = c("Ligand", "Receptor"), bty = "n",
           col = c(ligand_col, receptor_col), cex = 1, pt.cex = 3, border = "black") 
    
    legend(1.1, 0, pch = 20, legend = celltype_color$ID, bty = "n",
           col = celltype_color$color2, cex = 1, pt.cex = 3, border = "black") 
  }
}

target_celltypes <- c(
  "Transitional activated",   "Activation-regulated",        "CSMD1+", 
  "Cytotoxic effector", "Naive-like",  "Terminal-branch", 
  "NK-like",  "Cytotoxic memory-like",
  "CD4 T",            "CD8 T",            "Myeloid", 
  "B cells",          "Treg",             "NK cells", 
  "gdT",              "Innate-like T",    "HSPC"
)

pdf("./cellchat_new_2.5/Sender_recever_LR_Pso.pdf", width = 20, height = 20)
p2 <- 
  plot_cellchat_LR_v2(
    cellchat_obj = cellchat_list[[2]],
    celltype_inter = target_celltypes,
    celltype_color = mycolor2_named,
    ligand_col = "#FB8072",
    receptor_col = "#1988B0",
    text_size = 0.8,
    Group = "Psoriasis",
    legend_in_plot = TRUE,
    top_n = 1
  )
dev.off()
png("./cellchat_new_2.5/Sender_recever_LR_Pso.png", width = 20, height = 20,res = 600, units = "in")
p2 <- 
  plot_cellchat_LR_v2(
    cellchat_obj = cellchat_list[[2]],
    celltype_inter = target_celltypes,
    celltype_color = mycolor2_named,
    ligand_col = "#FB8072",
    receptor_col = "#1988B0",
    text_size = 0.8,
    Group = "Psoriasis",
    legend_in_plot = TRUE,
    top_n = 1
  )
dev.off()

pdf("./cellchat_new_2.5/Sender_recever_LR_Cont.pdf", width = 20, height = 20)
p1 <- 
  plot_cellchat_LR_v2(
    cellchat_obj = cellchat_list[[1]],
    celltype_inter = target_celltypes,
    celltype_color = mycolor2_named,
    ligand_col = "#FB8072",
    receptor_col = "#1988B0",
    text_size = 0.8,
    Group = "Psoriasis",
    legend_in_plot = TRUE,
    top_n = 1
  )
dev.off()
png("./cellchat_new_2.5/Sender_recever_LR_Cont.png", width = 20, height = 20,res = 600, units = "in")
p1 <- 
  plot_cellchat_LR_v2(
    cellchat_obj = cellchat_list[[1]],
    celltype_inter = target_celltypes,
    celltype_color = mycolor2_named,
    ligand_col = "#FB8072",
    receptor_col = "#1988B0",
    text_size = 0.8,
    Group = "Psoriasis",
    legend_in_plot = TRUE,
    top_n = 1
  )
dev.off()

gg1 <- compareInteractions(cellchat_merged, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat_merged, show.legend = F, group = c(1,2), measure = "weight")
save_plot_auto_2.5(gg1 + gg2, "5_Compare_Interactions_Barplot")

p_circle_diff_expr <- expression({
  par(mfrow = c(1,2), xpd=TRUE)
  
  netVisual_diffInteraction(cellchat_merged, weight.scale = T, comparison = c(1,2),
                            margin = 0.1,
                            arrow.width = 1, color.use = mycolor2_foucus,
                            arrow.size = 0.2)
  
  netVisual_diffInteraction(cellchat_merged, weight.scale = T, measure = "weight", comparison = c(1,2),
                            margin = 0.1,
                            arrow.width = 1, color.use = mycolor2_foucus,
                            arrow.size = 0.2)
  
  title("Red: Increased in Psoriasis | Blue: Decreased", outer = TRUE, line = -1)
})
save_plot_auto_2.5(p_circle_diff_expr, "6_Diff_Interaction_CirclePlot")

ht1 <- netVisual_heatmap(cellchat_merged, measure = "count", comparison = c(1,2), color.use = mycolor2_named)
ht2 <- netVisual_heatmap(cellchat_merged, measure = "weight", comparison = c(1,2), color.use = mycolor2_named)
save_plot_auto_2.5(expression({
  print(ht1 + ht2)
}), "7_Diff_Interaction_Heatmap")

if ("celltype" %in% colnames(cellchat_merged@meta)) {
  new_idents <- as.factor(cellchat_merged@meta$celltype)
  cellchat_merged@idents <- new_idents
  print("Successfully updated CellChat identity! Current identity levels:")
  print(levels(cellchat_merged@idents))
} else {
  stop("Error: 'celltype' column not found in cellchat_merged@meta!")
}

all_idents <- levels(cellchat_merged@idents)

gdt_cells <- grep("γδT|gdT", all_idents, value = TRUE)
skin_2.8_cells <- grep("KC-|Stroma|Myeloid|CD4_T_cell|Treg|Cytotoxic_NK_T", all_idents, value = TRUE)

message(paste("Senders (gdT):", paste(gdt_cells, collapse = ", ")))
message(paste("Receivers (skin):", paste(skin_2.8_cells, collapse = ", ")))

target_patterns <- c(
  "KC-",
  "Stroma",
  "Myeloid",
  "CD4_T_cell",
  "Treg",
  "Cytotoxic_NK_T"
)

regex_pattern <- paste(target_patterns, collapse = "|")

keep_pathways_merged <- union(keep_pathways_1, keep_pathways_2)

if(length(gdt_cells) > 0 && length(skin_2.8_cells) > 0) {
  
  p_bubble <- netVisual_bubble(cellchat_merged, 
                               sources.use = gdt_cells, 
                               targets.use = skin_2.8_cells, 
                               comparison = c(1,2),
                               angle.x = 45, 
                               signaling = keep_pathways_merged,
                               remove.isolate = FALSE,
                               title.name = "Diff Interactions: gdT -> skin_2.8 Cells")
  
  save_plot_auto_2.5(p_bubble, "8_Diff_Bubble_gdT_Sender_skin_2.8_Receiver")
  
} else {
  message("gdT or skin_2.8 cells not found, check naming")
}

if(length(gdt_cells) > 0 && length(skin_2.8_cells) > 0) {
  
  p_bubble_rev <- netVisual_bubble(cellchat_merged, 
                                   sources.use = skin_2.8_cells, 
                                   targets.use = gdt_cells, 
                                   comparison = c(1,2), 
                                   signaling = keep_pathways_merged,
                                   angle.x = 45, 
                                   remove.isolate = FALSE,
                                   title.name = "Diff Interactions: skin_2.8 Cells -> gdT")
  
  save_plot_auto_2.5(p_bubble_rev, "8_Diff_Bubble_skin_2.8_Sender_gdT_Receiver")
}

saveRDS(cellchat_list, file = "./cellchat_list_skin_2.14.rds")



all_gdt_names <- grep("γδT", names(mycolor2_named), value = TRUE)
target_gdt <- "CSMD1+"
remove_gdt_names <- setdiff(all_gdt_names, target_gdt)
mycolor2_named_foucus <- mycolor2_named[ !names(mycolor2_named) %in% remove_gdt_names ]
mycolor2_foucus <- unname(mycolor2_named_foucus)

print(names(mycolor2_named_foucus))

seu_user <- readRDS("../../rds/Overall_shinanxi_New.rds") 
cells_to_keep <- names(mycolor2_named)
seu_user_foucus <- subset(seu_user, subset = celltype %in% "CSMD1+")
seu_user_foucus$celltype <- droplevels(seu_user_foucus$celltype)
print(table(seu_user_foucus$celltype))

seu_skin$group <- seu_skin$cell_id
seu_combined_foucus <- merge(seu_skin, y = seu_user_foucus, add.cell.ids = c("skin", "gdT"))

print("Data merge complete.")

seu_combined_foucus$celltype_old <- seu_combined_foucus$celltype

current_types <- as.character(seu_combined_foucus$celltype)

new_types <- case_when(
  current_types %in% gdt_order ~ current_types,
  current_types == "KC-S.Corneum"   ~ "KC-S.Corneum",
  current_types == "KC-S.Spinosum"   ~ "KC-S.Spinosum",
  current_types == "KC-S.Basale"     ~ "KC-S.Basale",
  current_types == "Inflammatory_KC" ~ "Inflammatory_KC",
  current_types == "CD4_T_cell" ~ "CD4_T_cell",
  current_types == "Treg"       ~ "Treg",
  current_types %in% c("CD8_T_cell", "NK_cell", "CD161_T_cell") ~ "Cytotoxic_NK_T",
  current_types %in% c("Mature_DC", "Macrophage","Semimature_DC") ~ "Myeloid",
  current_types %in% c( "Endothelial_Cell", "Melanocyte") ~ "Stroma",
  TRUE ~ current_types
)

seu_combined_foucus$celltype <- new_types

final_levels <- c("CSMD1+", 
                  "KC-S.Corneum",  "KC-S.Spinosum", "KC-S.Basale", "Inflammatory_KC",
                  "CD4_T_cell", "Treg", "Cytotoxic_NK_T", "Myeloid", "Stroma")

existing_final <- unique(seu_combined_foucus$celltype)
final_levels <- final_levels[final_levels %in% existing_final]

seu_combined_foucus$celltype <- factor(seu_combined_foucus$celltype, levels = final_levels)
Idents(seu_combined_foucus) <- seu_combined_foucus$celltype

print(paste0("CellChat preparation complete. Number of cell subtypes: ", length(final_levels)))
print(table(seu_combined_foucus$celltype))

run_cellchat_pipeline <- function(seurat_obj, group_name) {
  message(paste0("Processing group: ", group_name))
  
  data.input <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  meta <- seurat_obj@meta.data
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
  
  CellChatDB <- CellChatDB.human
  cellchat@DB <- subsetDB(CellChatDB, search = "Secreted Signaling")
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  return(cellchat)
}

seurat_list <- SplitObject(seu_combined_foucus, split.by = "group")
Idents(seu_combined_foucus) <- "celltype"

cellchat_list <- list()

for (grp in names(seurat_list)) {
  
  seurat_obj <- seurat_list[[grp]]
  
  seurat_obj$celltype <- droplevels(as.factor(seurat_obj$celltype))
  
  Idents(seurat_obj) <- seurat_obj$celltype
  
  print(paste0(">>> Processing group: ", grp))
  print(paste("    Total cells:", ncol(seurat_obj)))
  print(paste("    Number of cell types:", length(levels(Idents(seurat_obj)))))
  
  if(ncol(seurat_obj) > 50 && length(levels(Idents(seurat_obj))) >= 2) {
    
    cellchat_list[[grp]] <- run_cellchat_pipeline(seurat_obj, grp)
    
  } else {
    print(paste("    [Skipped] Too few cells or single type"))
  }
}

cellchat_list <- list()
for (grp in names(seurat_list)) {
  if(ncol(seurat_list[[grp]]) > 50) {
    cellchat_list[[grp]] <- run_cellchat_pipeline(seurat_list[[grp]], grp)
  }
}

save_plot_auto_2.5 <- function(plot_code_or_obj, filename_base, width = 18, height = 9.5, type = "both") {
  if(!dir.exists("cellchat_CSMD1")) dir.create("cellchat_CSMD1", recursive = TRUE)
  f_png <- paste0("cellchat_CSMD1/", filename_base, ".png")
  f_pdf <- paste0("cellchat_CSMD1/", filename_base, ".pdf")
  
  draw_it <- function(obj) {
    if (inherits(obj, "grob")) {
      grid::grid.draw(obj)
    } else if (is.expression(obj)) {
      eval(obj)
    } else {
      print(obj)
    }
  }
  
  if(type %in% c("both", "png")) {
    png(f_png, width = width, height = height, units = "in", res = 600)
    grid::grid.newpage() 
    tryCatch({
      draw_it(plot_code_or_obj)
    }, error = function(e) message(paste("Error saving PNG:", filename_base, e$message)))
    dev.off()
  }
  
  if(type %in% c("both", "pdf")) {
    pdf(f_pdf, width = width, height = height)
    grid::grid.newpage()
    tryCatch({
      draw_it(plot_code_or_obj)
    }, error = function(e) message(paste("Error saving PDF:", filename_base, e$message)))
    dev.off()
  }
  
  message(paste("Saved successfully:", filename_base))
}

group.cellType <- union(levels(cellchat_list[[1]]@idents), levels(cellchat_list[[2]]@idents))
for (i in 1:length(cellchat_list)) {
  cellchat_list[[i]] <- liftCellChat(cellchat_list[[i]], group.new = group.cellType)
}
cellchat_merged <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))
compare_vec <- c(2, 1)

message("--- Merging CellChat objects ---")
cellchat_merged <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))

print("Merged cell type order:")
print(levels(cellchat_merged@idents$joint))

cell_types <- levels(cellchat_merged@meta$celltype)
print(paste("Number of cell types:", length(cell_types)))
print(paste("Provided colors:", length(mycolor2)))

mycolor2_named_foucus <- mycolor2[1:length(cell_types)]
names(mycolor2_named_foucus) <- cell_types

all_pathways_1 <- cellchat_list[[1]]@netP$pathways
keep_pathways_1 <- all_pathways_1[!grepl("HLA|MHC", all_pathways_1, ignore.case = TRUE)]

all_pathways_2 <- cellchat_list[[2]]@netP$pathways
keep_pathways_2 <- all_pathways_2[!grepl("HLA|MHC", all_pathways_2, ignore.case = TRUE)]

print(head(cellchat_list[[1]]@netP$pathways))
print(head(cellchat_list[[2]]@netP$pathways))

p1 <- cellchat_list[[1]]@netP$pathways
keep_pathways_1 <- p1[!grepl("HLA|MHC", p1, ignore.case = TRUE)]

p2 <- cellchat_list[[2]]@netP$pathways
keep_pathways_2 <- p2[!grepl("HLA|MHC", p2, ignore.case = TRUE)]

mat_control <- cellchat_list[[1]]@net$count
mat_disease <- cellchat_list[[2]]@net$count

library(pheatmap)
library(grid)
library(gridExtra)

global_max <- max(max(mat_control, na.rm = T), max(mat_disease, na.rm = T))
common_breaks <- seq(0, global_max, length.out = 100)
common_color <- colorRampPalette(c("white", "#E41A1C"))(100)

if(!dir.exists("cellchat_CSMD1")) dir.create("cellchat_CSMD1")

p3_ctrl <- pheatmap::pheatmap(cellchat_list[[1]]@net$count, border_color = "black",
                              cluster_cols = F, cluster_rows = F, fontsize = 10,
                              display_numbers = T, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Control", silent = TRUE)$gtable

p3_case <- pheatmap::pheatmap(mat_disease, border_color = "black",
                              cluster_cols = F, cluster_rows = F, fontsize = 10,
                              display_numbers = T, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Psoriasis", silent = TRUE)$gtable

p3_combined <- arrangeGrob(p3_ctrl, p3_case, ncol = 2, top = textGrob("Interaction Count (Raw)", gp=gpar(fontsize=15, font=2)))
save_plot_auto_2.5(p3_combined, "CellChat_Heatmap_NoCluster_WithNum")

p4_ctrl <- pheatmap::pheatmap(mat_control, border_color = "black",
                              cluster_cols = T, cluster_rows = T, fontsize = 10,
                              display_numbers = T, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Control", silent = TRUE)$gtable

p4_case <- pheatmap::pheatmap(mat_disease, border_color = "black",
                              cluster_cols = T, cluster_rows = T, fontsize = 10,
                              display_numbers = T, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Psoriasis", silent = TRUE)$gtable

p4_combined <- arrangeGrob(p4_ctrl, p4_case, ncol = 2, top = textGrob("Interaction Count (Clustered)", gp=gpar(fontsize=15, font=2)))
save_plot_auto_2.5(p4_combined, "CellChat_Heatmap_Clustered_WithNum")

p5_ctrl <- pheatmap::pheatmap(mat_control, border_color = "black",
                              cluster_cols = F, cluster_rows = F, fontsize = 10,
                              display_numbers = F, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Control", silent = TRUE)$gtable

p5_case <- pheatmap::pheatmap(mat_disease, border_color = "black",
                              cluster_cols = F, cluster_rows = F, fontsize = 10,
                              display_numbers = F, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Psoriasis", silent = TRUE)$gtable

p5_combined <- arrangeGrob(p5_ctrl, p5_case, ncol = 2, top = textGrob("Interaction Count (Raw, No Num)", gp=gpar(fontsize=15, font=2)))
save_plot_auto_2.5(p5_combined, "CellChat_Heatmap_NoCluster_NoNum")

p6_ctrl <- pheatmap::pheatmap(mat_control, border_color = "black",
                              cluster_cols = T, cluster_rows = T, fontsize = 10,
                              display_numbers = F, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Control", silent = TRUE)$gtable

p6_case <- pheatmap::pheatmap(mat_disease, border_color = "black",
                              cluster_cols = T, cluster_rows = T, fontsize = 10,
                              display_numbers = F, number_color = "black", number_format = "%.0f",
                              color = common_color, breaks = common_breaks,
                              main = "Psoriasis", silent = TRUE)$gtable

p6_combined <- arrangeGrob(p6_ctrl, p6_case, ncol = 2, top = textGrob("Interaction Count (Clustered, No Num)", gp=gpar(fontsize=15, font=2)))
save_plot_auto_2.5(p6_combined, "CellChat_Heatmap_Clustered_NoNum")

mat <- cellchat_merged@net$count
cell_types <- levels(cellchat_merged@idents$joint)
groupSize <- as.numeric(table(cellchat_merged@idents$joint))

code_single_circle <- expression({
  par(mfrow = c(3,3), xpd=TRUE, mar = c(1,1,1,1))
  for (i in 1:length(cell_types)) {
    netVisual_circle(cellchat_list[["Psoriasis"]]@net$count, vertex.weight = groupSize, weight.scale = T, 
                     edge.weight.max = max(cellchat_list[["Psoriasis"]]@net$count), 
                     title.name = paste("Pso -", cell_types[i]), color.use = mycolor2_foucus, sources.use = cell_types[i])
  }
})
save_plot_auto_2.5(code_single_circle, "03_Single_CellType_Source_Circle_Psoriasis")

code_single_circle <- expression({
  par(mfrow = c(3,3), xpd=TRUE, mar = c(1,1,1,1))
  for (i in 1:length(cell_types)) {
    netVisual_circle(cellchat_list[["Control"]]@net$count, vertex.weight = groupSize, weight.scale = T, 
                     edge.weight.max = max(cellchat_list[["Control"]]@net$count), 
                     title.name = paste("Pso -", cell_types[i]), color.use = mycolor2_foucus, sources.use = cell_types[i])
  }
})
save_plot_auto_2.5(code_single_circle, "03_Single_CellType_Source_Circle_Control")

groupSize <- as.numeric(table(cellchat_list[[1]]@idents))

p1 <- expression({
  par(mfrow = c(1,2), xpd=TRUE)
  
  netVisual_circle(cellchat_list[[1]]@net$count, vertex.weight = groupSize, weight.scale = T,
                   label.edge= F,                    arrow.width = 1.2,
                   arrow.size = 0.24, color.use = mycolor2_foucus,
                   title.name = "Number of interactions in Cont")
  netVisual_circle(cellchat_list[[1]]@net$weight, vertex.weight = groupSize, weight.scale = T, 
                   label.edge= F, arrow.width = 1.2,  color.use = mycolor2_foucus,
                   arrow.size = 0.24 ,title.name = "Interaction weights/strength in Cont")
  
  groupSize <- as.numeric(table(cellchat_list[[2]]@idents))})
save_plot_auto_2.5(p1, "03_Total_CellType_Source_Circle_Control")

p2 <- expression({
  par(mfrow = c(1,2), xpd=TRUE)
  
  netVisual_circle(cellchat_list[[2]]@net$count, vertex.weight = groupSize, weight.scale = T,
                   label.edge= F, 
                   arrow.width = 1.2,
                   arrow.size = 0.24, color.use = mycolor2_foucus,
                   title.name = "Number of interactions in Pso")
  netVisual_circle(cellchat_list[[2]]@net$weight, vertex.weight = groupSize, weight.scale = T, 
                   label.edge= F, 
                   arrow.width = 1.2,
                   arrow.size = 0.24, color.use = mycolor2_foucus,
                   title.name = "Interaction weights/strength in Pso")})
save_plot_auto_2.5(p2, "03_Total_CellType_Source_Circle_Psoriasis")

cellchat_list[[1]]@netP$pathways

pp1 <- rankNet(cellchat_merged, mode = "comparison", stacked = T, do.stat = TRUE)
pp2 <- rankNet(cellchat_merged, mode = "comparison", stacked = F, do.stat = TRUE)
save_plot_auto_2.5(pp1 + pp2, "03_rankNet")

target_groups   <- c("Control", "Psoriasis")
chat_ctrl <- cellchat_list[[1]] 
chat_pso  <- cellchat_list[[2]]

pathway_name <- "CLEC"
obj_use <- chat_ctrl
prefix <- "Control_CLEC"

p1 = netVisual_aggregate(obj_use, signaling = pathway_name, color.use = mycolor2_foucus)
save_plot_auto_2.5(p1, paste0(prefix, "_circle"))

p3 = netVisual_heatmap(obj_use, signaling = pathway_name, color.heatmap = "Reds", color.use = mycolor2_named_foucus)
save_plot_auto_2.5(p3, paste0(prefix, "_heatmap"))

p4_expr = expression({
  netVisual_chord_cell(obj_use, signaling = pathway_name, group = group.cellType, 
                       title.name = paste0(pathway_name, " (Control) - Immune vs Other"), color.use = mycolor2_named_foucus)
})
save_plot_auto_2.5(p4_expr, paste0(prefix, "_grouped_chord"))

p5 = netAnalysis_contribution(obj_use, signaling = pathway_name)
save_plot_auto_2.5(p5, paste0(prefix, "_contribution"))

pairLR_use <- extractEnrichedLR(obj_use, signaling = pathway_name, geneLR.return = FALSE)[1]
p6_expr = expression({
  netVisual_individual(obj_use, signaling = pathway_name, pairLR.use = pairLR_use, layout = "circle", color.use = mycolor2_foucus)
})
save_plot_auto_2.5(p6_expr, paste0(prefix, "_", pairLR_use, "_interaction"))

state_list <- list(
  "Control"   = cellchat_list[[1]],
  "Psoriasis" = cellchat_list[[2]]
)

target_pathways <- c("MIF")

for (current_state in names(state_list)) {
  
  obj_use <- state_list[[current_state]]
  
  for (pathway_name in target_pathways) {
    
    if (!(pathway_name %in% obj_use@netP$pathways)) {
      message(paste("Skip:", pathway_name, "not found in", current_state))
      next
    }
    
    prefix <- paste0(current_state, "_", pathway_name)
    
    message(paste("Processing:", prefix))
    
    p1 = netVisual_aggregate(obj_use, signaling = pathway_name, color.use = mycolor2_foucus)
    save_plot_auto_2.5(p1, paste0(prefix, "_circle"))
    
    p3 = netVisual_heatmap(obj_use, signaling = pathway_name, color.heatmap = "Reds", color.use = mycolor2_named_foucus)
    save_plot_auto_2.5(p3, paste0(prefix, "_heatmap"))
    
    group.cellType <- c(
      "Transitional activated"    = "gdT",
      "Cytotoxic effector"  = "gdT",
      "CSMD1+"   = "gdT",
      "NK-like"   = "gdT",
      "Activation-regulated"         = "gdT",
      "Naive-like"   = "gdT",
      "Cytotoxic memory-like" = "gdT",
      "Terminal-branch"         = "gdT",
      "gdT"               = "gdT",
      
      "KC-S.Corneum"             = "Other",
      "KC-S.Granulosum"             = "Other",
      "KC-S.Spinosum"           = "Other",
      "KC-S.Basale"           = "Other",
      "KC-Wound/Prolif"              = "Other",
      "CD4_T_cell"          = "Other",
      "Treg"     = "Other",
      "Cytotoxic_NK_T"              = "Other",
      "Myeloid"              = "Other",
      "Stroma "              = "Other"
    )
    
    color_grouped <- c("gdT" = "#FF34B3", "Other" = "grey80")
    p4_expr = expression({
      netVisual_chord_cell(obj_use, signaling = pathway_name, group = group.cellType,
                           title.name = paste0(pathway_name, " (", current_state, ") - gdT vs Other"), color.use = mycolor2_named_foucus)
    })
    save_plot_auto_2.5(p4_expr, paste0(prefix, "_grouped_chord"))
    
    p5 = netAnalysis_contribution(obj_use, signaling = pathway_name)
    save_plot_auto_2.5(p5, paste0(prefix, "_contribution"))
    
    tryCatch({
      pairLR_use <- extractEnrichedLR(obj_use, signaling = pathway_name, geneLR.return = FALSE)[1]
      
      if (!is.na(pairLR_use)) {
        p6_expr = expression({
          netVisual_individual(obj_use, signaling = pathway_name, pairLR.use = pairLR_use, layout = "circle", color.use = mycolor2_foucus)
        })
        save_plot_auto_2.5(p6_expr, paste0(prefix, "_", pairLR_use, "_interaction"))
      }
    }, error = function(e) { message(paste("Error plotting individual LR for", prefix)) })
    
    p7 = netVisual_bubble(obj_use, remove.isolate = FALSE,signaling = pathway_name)
    save_plot_auto_2.5(p7, paste0(prefix, "_sender_receiver"))
    
    p8 <- expression(netAnalysis_signalingRole_network(obj_use, signaling = pathway_name, width = 8, height = 2.5, font.size = 10,color.use = mycolor2_named_foucus))
    save_plot_auto_2.5(p8, paste0(prefix, "_charactor"))
    
    p9 <- netAnalysis_signalingRole_scatter(obj_use, signaling = pathway_name,color.use = mycolor2_foucus)
    save_plot_auto_2.5(p9, paste0(prefix, "_signalingRole_scatte"))
  }
}

CSMD1_chord <- expression(netVisual_chord_gene(cellchat_list[[1]], sources.use = "CSMD1+", lab.cex = 0.5,legend.pos.y =50))
save_plot_auto_2.5(CSMD1_chord,"CSMD1_chord")

for (i in 1:length(cellchat_list)) {
  cellchat_list[[i]] <- netAnalysis_computeCentrality(cellchat_list[[i]], slot.name = "netP")
}
ht_out_control <- netAnalysis_signalingRole_heatmap(cellchat_list[[1]], pattern = "outgoing", width=12, height=14, color.use = mycolor2_named_foucus,signaling = keep_pathways_1)
ht_in_control <- netAnalysis_signalingRole_heatmap(cellchat_list[[1]], pattern = "incoming", width=12, height=14, color.use = mycolor2_named_foucus,signaling = keep_pathways_1)
save_plot_auto_2.5(ht_out_control + ht_in_control, "SignalingRole_Heatmap_control")

ht_out_PSO <- netAnalysis_signalingRole_heatmap(cellchat_list[[2]], pattern = "outgoing", width=12, height=14, color.use = mycolor2_named_foucus,signaling = keep_pathways_2)
ht_in_PSO <- netAnalysis_signalingRole_heatmap(cellchat_list[[2]], pattern = "incoming", width=12, height=14, color.use = mycolor2_named_foucus,signaling = keep_pathways_2)
save_plot_auto_2.5(ht_out_PSO + ht_in_PSO, "SignalingRole_Heatmapr_PSO")

p1 <- netVisual_aggregate(
  cellchat_list[[1]], 
  signaling = keep_pathways_1, 
  layout = "chord",
  signaling.name = "Control",
  title.name = "NULL" 
)
save_plot_auto_2.5(p1, "_Cont_circle")

p2 <- netVisual_aggregate(
  cellchat_list[[2]],
  signaling = keep_pathways_2, 
  color.use = mycolor2, 
  layout = "chord", 
  signaling.name = "Psoriasis",
  title.name = "NULL"
)
save_plot_auto_2.5(p2, "_Pso_circle")

gg1 <- netAnalysis_signalingRole_scatter(cellchat_list[[1]],color.use = mycolor2_foucus,signaling = keep_pathways_1)
save_plot_auto_2.5(gg1,"netAnalysis_Cont")
gg2 <- netAnalysis_signalingRole_scatter(cellchat_list[[2]],color.use = mycolor2_foucus,signaling = keep_pathways_2)
save_plot_auto_2.5(gg2,"netAnalysis_Pso")

bubble_Cont = netVisual_bubble(cellchat_list[[1]], remove.isolate = FALSE,signaling = keep_pathways_1 )
save_plot_auto_2.5(bubble_Cont,"bubble_Cont", width = 30)
bubble_Pso = netVisual_bubble(cellchat_list[[2]], remove.isolate = FALSE,signaling = keep_pathways_2)
save_plot_auto_2.5(bubble_Pso,"bubble_Pso", width = 30)

NMF_control <- selectK(cellchat_list[[1]], pattern = "outgoing")
NMF_PSO <-  selectK(cellchat_list[[2]], pattern = "outgoing")
save_plot_auto_2.5(NMF_control + NMF_PSO, "17_NMF_PSO")

k_c = 3
k_p = 5
cellchat_list[[1]] <- identifyCommunicationPatterns(cellchat_list[[1]], pattern = "outgoing", k = 3)
cellchat_list[[2]] <- identifyCommunicationPatterns(cellchat_list[[2]], pattern = "outgoing", k = 5)

png("./cellchat_CSMD1/Pattern_River_Outgoing_Cont.png", width = 8, height = 8,res = 600, units = "in")
netAnalysis_river(cellchat_list[[1]], pattern = "outgoing")
dev.off()
pdf("./cellchat_CSMD1/Pattern_River_Outgoing_Cont.pdf", width = 8, height = 8)
netAnalysis_river(cellchat_list[[1]], pattern = "outgoing")
dev.off()

png("./cellchat_CSMD1/Pattern_Dot_Outgoing_Cont.png", width = 8, height = 8,res = 600, units = "in")
netAnalysis_dot(cellchat_list[[1]], pattern = "outgoing")
dev.off()
pdf("./cellchat_CSMD1/Pattern_Dot_Outgoing_Cont.pdf", width = 8, height = 8)
netAnalysis_dot(cellchat_list[[1]], pattern = "outgoing")
dev.off()

png("./cellchat_CSMD1/Pattern_River_Outgoing_Pso.png", width = 8, height = 8,res = 600, units = "in")
netAnalysis_river(cellchat_list[[2]], pattern = "outgoing")
dev.off()
pdf("./cellchat_CSMD1/Pattern_River_Outgoing_Pso.pdf", width = 8, height = 8)
netAnalysis_river(cellchat_list[[2]], pattern = "outgoing")
dev.off()

png("./cellchat_CSMD1/Pattern_River_Outgoing_Pso.png", width = 8, height = 8,res = 600, units = "in")
netAnalysis_river(cellchat_list[[2]], pattern = "outgoing")
dev.off()
pdf("./cellchat_CSMD1/Pattern_River_Outgoing_Pso.pdf", width = 8, height = 8)
netAnalysis_river(cellchat_list[[2]], pattern = "outgoing")
dev.off()

pdf("./cellchat_new_2.5/Sender_recever_LR_Pso.pdf", width = 20, height = 20)
p2 <- 
  plot_cellchat_LR_v2(
    cellchat_obj = cellchat_list[[2]],
    celltype_inter = target_celltypes,
    celltype_color = mycolor2_named_foucus,
    ligand_col = "#FB8072",
    receptor_col = "#1988B0",
    text_size = 0.8,
    Group = "Psoriasis",
    legend_in_plot = TRUE,
    top_n = 1
  )
dev.off()
png("./cellchat_new_2.5/Sender_recever_LR_Pso.png", width = 20, height = 20,res = 600, units = "in")
p2 <- 
  plot_cellchat_LR_v2(
    cellchat_obj = cellchat_list[[2]],
    celltype_inter = target_celltypes,
    celltype_color = mycolor2_named_foucus,
    ligand_col = "#FB8072",
    receptor_col = "#1988B0",
    text_size = 0.8,
    Group = "Psoriasis",
    legend_in_plot = TRUE,
    top_n = 1
  )
dev.off()

pdf("./cellchat_new_2.5/Sender_recever_LR_Cont.pdf", width = 20, height = 20)
p1 <- 
  plot_cellchat_LR_v2(
    cellchat_obj = cellchat_list[[1]],
    celltype_inter = target_celltypes,
    celltype_color = mycolor2_named_foucus,
    ligand_col = "#FB8072",
    receptor_col = "#1988B0",
    text_size = 0.8,
    Group = "Psoriasis",
    legend_in_plot = TRUE,
    top_n = 1
  )
dev.off()
png("./cellchat_new_2.5/Sender_recever_LR_Cont.png", width = 20, height = 20,res = 600, units = "in")
p1 <- 
  plot_cellchat_LR_v2(
    cellchat_obj = cellchat_list[[1]],
    celltype_inter = target_celltypes,
    celltype_color = mycolor2_named_foucus,
    ligand_col = "#FB8072",
    receptor_col = "#1988B0",
    text_size = 0.8,
    Group = "Psoriasis",
    legend_in_plot = TRUE,
    top_n = 1
  )
dev.off()

gg1 <- compareInteractions(cellchat_merged, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat_merged, show.legend = F, group = c(1,2), measure = "weight")
save_plot_auto_2.5(gg1 + gg2, "5_Compare_Interactions_Barplot")

p_circle_diff_expr <- expression({
  par(mfrow = c(1,2), xpd=TRUE)
  
  netVisual_diffInteraction(cellchat_merged, weight.scale = T, comparison = c(1,2),
                            margin = 0.1,
                            arrow.width = 1, color.use = mycolor2_foucus,
                            arrow.size = 0.2)
  
  netVisual_diffInteraction(cellchat_merged, weight.scale = T, measure = "weight", comparison = c(1,2),
                            margin = 0.1,
                            arrow.width = 1, color.use = mycolor2_foucus,
                            arrow.size = 0.2)
  
  title("Red: Increased in Psoriasis | Blue: Decreased", outer = TRUE, line = -1)
})
save_plot_auto_2.5(p_circle_diff_expr, "6_Diff_Interaction_CirclePlot")

ht1 <- netVisual_heatmap(cellchat_merged, measure = "count", comparison = c(1,2), color.use = mycolor2_named_foucus)
ht2 <- netVisual_heatmap(cellchat_merged, measure = "weight", comparison = c(1,2), color.use = mycolor2_named_foucus)
save_plot_auto_2.5(expression({
  print(ht1 + ht2)
}), "7_Diff_Interaction_Heatmap")

if ("celltype" %in% colnames(cellchat_merged@meta)) {
  new_idents <- as.factor(cellchat_merged@meta$celltype)
  cellchat_merged@idents <- new_idents
  print("Successfully updated CellChat identity! Current identity levels:")
  print(levels(cellchat_merged@idents))
} else {
  stop("Error: 'celltype' column not found in cellchat_merged@meta!")
}

all_idents <- levels(cellchat_merged@idents)

gdt_cells <- grep("γδT|gdT", all_idents, value = TRUE)
skin_2.8_cells <- grep("KC-|Stroma|Myeloid|CD4_T_cell|Treg|Cytotoxic_NK_T", all_idents, value = TRUE)

message(paste("Senders (gdT):", paste(gdt_cells, collapse = ", ")))
message(paste("Receivers (skin):", paste(skin_2.8_cells, collapse = ", ")))

target_patterns <- c(
  "KC-",
  "Stroma",
  "Myeloid",
  "CD4_T_cell",
  "Treg",
  "Cytotoxic_NK_T"
)

regex_pattern <- paste(target_patterns, collapse = "|")

keep_pathways_merged <- union(keep_pathways_1, keep_pathways_2)

if(length(gdt_cells) > 0 && length(skin_2.8_cells) > 0) {
  
  p_bubble <- netVisual_bubble(cellchat_merged, 
                               sources.use = gdt_cells, 
                               targets.use = skin_2.8_cells, 
                               comparison = c(1,2),
                               angle.x = 45, 
                               signaling = keep_pathways_merged,
                               remove.isolate = FALSE,
                               title.name = "Diff Interactions: gdT -> skin_2.8 Cells")
  
  save_plot_auto_2.5(p_bubble, "8_Diff_Bubble_gdT_Sender_skin_2.8_Receiver")
  
} else {
  message("gdT or skin_2.8 cells not found, check naming")
}

if(length(gdt_cells) > 0 && length(skin_2.8_cells) > 0) {
  
  p_bubble_rev <- netVisual_bubble(cellchat_merged, 
                                   sources.use = skin_2.8_cells, 
                                   targets.use = gdt_cells, 
                                   comparison = c(1,2), 
                                   signaling = keep_pathways_merged,
                                   angle.x = 45, 
                                   remove.isolate = FALSE,
                                   title.name = "Diff Interactions: skin_2.8 Cells -> gdT")
  
  save_plot_auto_2.5(p_bubble_rev, "8_Diff_Bubble_skin_2.8_Sender_gdT_Receiver")
}

saveRDS(cellchat_list, file = "./cellchat_list_skin_foucus2.14.rds")




#######################13.Monocle2#########
```r
seurat_ob = readRDS("/thinker/3.tangjiale/shinanxi/monocle2/subset_ob.rds")
cds <- as.CellDataSet(seurat_ob)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds = detectGenes(cds,min_expr = 1)
expressed_genes = row.names(subset(fData(cds),num_cells_expressed>10))
clustering_DEGs = differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr ="~celltype",cores = 1)
class(clustering_DEGs)
write.csv(clustering_DEGs,"/thinker/3.tangjiale/shinanxi/monocle2_OE/diff_test_res.csv")
write.table(clustering_DEGs,"/thinker/3.tangjiale/shinanxi/monocle2_OE/diff_test_res.xls")
featureData(cds)@data[rownames(clustering_DEGs),"pval"]=clustering_DEGs$pval
featureData(cds)@data[rownames(clustering_DEGs),"qval"]=clustering_DEGs$qval
ordering_genes <- row.names (subset(clustering_DEGs, qval < 0.01))
gbm_cds = setOrderingFilter(cds,ordering_genes = ordering_genes)
plot_ordering_genes(gbm_cds)
save(gbm_cds, file = "/thinker/3.tangjiale/shinanxi/monocle2_OE/gbm_cds.RData")
gbm_cds = reduceDimension(gbm_cds,max_components = 2,verbose = T)
saveRDS(gbm_cds,file = "/thinker/3.tangjiale/shinanxi/monocle2_OE/celltype_gbm_cds.rds")
gbm_cds = readRDS("/thinker/3.tangjiale/shinanxi/monocle2_OE/celltype_gbm_cds.rds")
gbm_cds = orderCells(gbm_cds,reverse = F)
gbm_cds <- orderCells(cds_DGT, root_state = 2)
saveRDS(gbm_cds,file = "/thinker/3.tangjiale/shinanxi/monocle2_OE/subcell_cds_DGT_Pseudotime.rds")
cds_DGT = gbm_cds 
p1 = plot_cell_trajectory(cds_DGT, cell_size = 1, color_by = "subcelltype",show_branch_points = F) + facet_wrap(~subcelltype, nrow = 2) + scale_color_manual(values = mycolor2_celltype)
p2 = plot_cell_trajectory(cds_DGT, cell_size = 1, color_by = "subcelltype",show_branch_points = F) + scale_color_manual(values = mycolor2_celltype)
p5 = plot_cell_trajectory(cds_DGT, cell_size = 1, color_by = "orig.ident",show_branch_points = F) + facet_wrap(~orig.ident, nrow = 2) + scale_color_manual(values = mycolor3_samples)
p6 = plot_cell_trajectory(cds_DGT, cell_size = 1, color_by = "orig.ident",show_branch_points = F) + scale_color_manual(values = mycolor3_samples)
p7 = plot_cell_trajectory(cds_DGT, cell_size = 1, color_by = "State",show_branch_points = F) + facet_wrap(~State, nrow = 2) 
p8 = plot_cell_trajectory(cds_DGT, cell_size = 1, color_by = "State",show_branch_points = F)
p9 = plot_cell_trajectory(cds_DGT, color_by = "Pseudotime",show_branch_points = F)+scale_colour_viridis_c(option = "inferno",direction = -1)
p10 = plot_cell_trajectory(cds_DGT, cell_size = 1, color_by = "group",show_branch_points = F) + facet_wrap(~group, nrow = 2) + scale_color_manual(values = mycolor3_samples)
p11 = plot_cell_trajectory(cds_DGT, cell_size = 1, color_by = "group",show_branch_points = F) + scale_color_manual(values = mycolor3_samples)
dir.create("/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/")
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/subcell_plottrajectory_celltype.png', p1, width=8, height=8, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/subcell_plottrajectory_celltype.pdf', p1, width=8, height=8, dpi=300)
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/plottrajectory_celltype.png', p2, width=6.5, height=7, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/plottrajectory_celltype.pdf', p2, width=6.5, height=7, dpi=300)
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/subcell_plottrajectory_sample.png', p5, width=6.5, height=7, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/subcell_plottrajectory_sample.pdf', p5, width=6.5, height=7, dpi=300)
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/plottrajectory_sample.png', p6, width=6.5, height=7, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/plottrajectory_sample.pdf', p6, width=6.5, height=7, dpi=300)
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/subcell_plottrajectory_State.png', p7, width=10, height=7.5, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/subcell_plottrajectory_State.pdf', p7, width=10, height=7.5, dpi=300)
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/plottrajectory_State.png', p8, width=6.5, height=7, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/plottrajectory_State.pdf', p8, width=6.5, height=7, dpi=300)
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/plottrajectory_Pseudotime.png', p9, width=6.5, height=7, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/plottrajectory_Pseudotime.pdf', p9, width=6.5, height=7, dpi=300)
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/subcell_plottrajectory_group.png', p10, width=6.5, height=14, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/subcell_plottrajectory_group.pdf', p10, width=6.5, height=14, dpi=300)
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/plottrajectory_group.png', p11, width=6.5, height=7, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/plottrajectory_group.pdf', p11, width=6.5, height=7, dpi=300)
plot_cell_trajectory(cds_DGT, markers = "KAZN",use_color_gradient=T,cell_size = 1,cell_link_size = 1.5)
exprData = cds_DGT@assayData$exprs
exprData = LogNormalize(exprData)
cds_DGT$KAZN = exprData["KAZN",]
p10 <- plot_cell_trajectory(cds_DGT, color_by = "KAZN",cell_size=1.5)+ scale_color_gradient(low = "transparent",high = "red",limits= range(cds_DGT$KAZN))
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/plottrajectory_KAZN_state4.png', p10, width=10, height=7.5, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/plottrajectory_KAZN_state4.pdf', p10, width=10, height=7.5, dpi=300)
plot_cell_trajectory(cds_DGT, markers = "MICOS10",use_color_gradient=T,cell_size = 1,cell_link_size = 1.5)
exprData = cds_DGT@assayData$exprs
exprData = LogNormalize(exprData)
cds_DGT$MICOS10 = exprData["MICOS10",]
p10 <- plot_cell_trajectory(cds_DGT, color_by = "MICOS10",cell_size=1.5)+ scale_color_gradient(low = "transparent",high = "red",limits= range(cds_DGT$MICOS10))
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/plottrajectory_MICOS10_4.png', p10, width=10, height=7.5, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/plottrajectory_MICOS10_4.pdf', p10, width=10, height=7.5, dpi=300)
plot_cell_trajectory(cds_DGT, markers = "UBXN11",use_color_gradient=T,cell_size = 1,cell_link_size = 1.5)
exprData = cds_DGT@assayData$exprs
exprData = LogNormalize(exprData)
cds_DGT$UBXN11 = exprData["UBXN11",]
p10 <- plot_cell_trajectory(cds_DGT, color_by = "UBXN11",cell_size=1.5)+ scale_color_gradient(low = "transparent",high = "red",limits= range(cds_DGT$UBXN11))
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/plottrajectory_UBXN11_4.png', p10, width=10, height=7.5, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/plottrajectory_UBXN11_4.pdf', p10, width=10, height=7.5, dpi=300)
plot_cell_trajectory(cds_DGT, markers = "UBXN11",use_color_gradient=T,cell_size = 1,cell_link_size = 1.5)
exprData = cds_DGT@assayData$exprs
exprData = LogNormalize(exprData)
cds_DGT$UBXN11 = exprData["UBXN11",]
p10 <- plot_cell_trajectory(cds_DGT, color_by = "UBXN11",cell_size=1.5)+ scale_color_gradient(low = "transparent",high = "red",limits= range(cds_DGT$UBXN11))
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/plottrajectory_UBXN11_4.png', p10, width=10, height=7.5, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/plottrajectory_UBXN11_4.pdf', p10, width=10, height=7.5, dpi=300)
p10 = plot_cell_trajectory(cds_DGT, markers = "AGO4",use_color_gradient=T,cell_size = 1,cell_link_size = 1.5)
exprData = cds_DGT@assayData$exprs
exprData = LogNormalize(exprData)
cds_DGT$AGO4 = exprData["AGO4",]
p10 <- plot_cell_trajectory(cds_DGT, color_by = "AGO4",cell_size=1.5)+ scale_color_gradient(low = "transparent",high = "red",limits= range(cds_DGT$AGO4))
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/plottrajectory_AGO4_3.png', p10, width=10, height=7.5, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/plottrajectory_AGO4_3.pdf', p10, width=10, height=7.5, dpi=300)
cds_DGT_pseudotimegenes <- differentialGeneTest(cds_DGT,fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 1)
cds_DGT_pseudotimegenes_sig <- subset(cds_DGT_pseudotimegenes)
write.csv(cds_DGT_pseudotimegenes_sig,"/thinker/3.tangjiale/shinanxi/monocle2_OE/cds_DGT_pseudotimegenes_sig.csv")
saveRDS(cds_DGT_pseudotimegenes_sig,file = "/thinker/3.tangjiale/shinanxi/monocle2_OE/cds_DGT_pseudotimegenes_sig.rds")
Time_genes <- cds_DGT_pseudotimegenes_sig %>% pull(gene_short_name) %>% as.character()
p <- plot_pseudotime_heatmap(cds_DGT[Time_genes,], 
                             num_cluster = 4, 
                             show_rownames = T, 
                             return_heatmap = T)
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/top15marker_pseudotime_heatmap.png', p, width=7.5, height=7.5, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/monocle2_OE/plot/top15marker_pseudotime_heatmap.pdf', p, width=7.5, height=7.5, dpi=300)
genes <- as.factor(subset(gbm_cds@featureData@data, use_for_ordering == TRUE)$gene_short_name)
to_be_tested <- row.names(subset(fData(gbm_cds), gene_short_name %in% levels(genes)))
to_be_tested = to_be_tested[1:500]
gbm_cds <- gbm_cds[to_be_tested, ]
p <- plot_pseudotime_heatmap(gbm_cds, cluster_rows = T, num_clusters = 4, show_rownames = T, return_heatmap = T)
p = plot_cell_trajectory(cds_DGT, color_by = "Pseudotime",show_branch_points = T)+scale_colour_viridis_c(option = "inferno")
BEAM_res <- BEAM(gbm_cds, branch_point = 1)
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
write.csv(BEAM_res,file = "/thinker/3.tangjiale/shinanxi/monocle2_OE/BEAM_res.csv")
saveRDS(BEAM_res,file = "/thinker/3.tangjiale/shinanxi/monocle2_OE/BEAM_res.rds")
num_clusters = 4
p = plot_genes_branched_heatmap(gbm_cds[row.names(subset(BEAM_res,
                                                         qval < 1e-4)),],
                                branch_point = 1,
                                num_clusters = num_clusters,
                                cores = 1,
                                use_gene_short_name = T,
                                show_rownames = F,return_heatmap = T)
annotation_row <- data.frame(Module = factor(cutree(p$ph_res$tree_row, k = num_clusters)))
head(seurat_ob@reductions$umap@cell.embeddings)
head(cds_DGT@phenoData@data)
seurat_ob@meta.data$Pseudotime <- cds_DGT@phenoData@data$Pseudotime 
head(seurat_ob@meta.data)
saveRDS(seurat_ob,"/thinker/3.tangjiale/shinanxi/monocle2_OE/seurat_ob_pseudotime.rds")
mydata<- FetchData(seurat_ob,vars = c("UMAP_1","UMAP_2","Pseudotime"))
p <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = Pseudotime))+
  geom_point(size = 1)+scale_colour_viridis_c(option = "inferno")
p4 <- p + theme_bw() + theme(panel.border = element_blank(), 
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), 
                             axis.line = element_line(colour = "black"))
p4
ggsave('/thinker/3.tangjiale/LZH/figures/pseudotimeUMAP.png', p4, width=9, height=8, dpi=600)
ggsave('/thinker/3.tangjiale/LZH/figures/pseudotimeUMAP.pdf', p4, width=9, height=8, dpi=300)
```


######14.cytotrace1#################
library(monocle3)
library(Seurat)
library(SeuratObject)
library(Matrix)
library(data.table)
library(CellChat)
BiocManager::install("sva")
devtools::install_local("/thinker/3.tangjiale/02_package/CytoTRACE_0.3.3.tar.gz")
library(CytoTRACE)
sc <- readRDS("/thinker/3.tangjiale/shinanxi/rds/Overall_shinanxi_New.rds")

colors1 <- c(
  "Transitional activated"="#FF34B3",
  "Cytotoxic effector"="#00F5FF",
  "CSMD1+"="#BC8F8F",
  "NK-like"="#ADFF2F",
  "Activation-regulated"="#FFFF02",
  "Naive-like"="#00CD00",
  "Cytotoxic memory-like"= "#FF6A6A",
  "Terminal-branch" = "#7FFFD4")
exp1 <- as.matrix(sc@assays$RNA@counts)
exp1 <- exp1[apply(exp1 > 0,1,sum) >= 5,]
results <- CytoTRACE(exp1,ncores = 1)
phenot <- sc$celltype
phenot <- as.character(phenot)
names(phenot) <- rownames(sc@meta.data)
emb <- sc@reductions[["umap"]]@cell.embeddings
plotCytoTRACE(results, phenotype = phenot, emb = emb, outputDir = '/thinker/3.tangjiale/shinanxi/cytotrace1/color',colors =c(
  "Activation-regulated"="#FFFF02",
  "Transitional activated"="#FF34B3",
  "Cytotoxic effector"="#00F5FF",
  "NK-like"="#ADFF2F",
  "Naive-like"="#00CD00",
  "Cytotoxic memory-like"= "#FF6A6A",
  "Terminal-branch" = "#7FFFD4",
  "CSMD1+"="#BC8F8F"))
plotCytoTRACE(results, phenotype = phenot, emb = emb, outputDir = '/thinker/3.tangjiale/shinanxi/cytotrace1/')
plotCytoGenes(results, numOfGenes = 30, outputDir = '/thinker/3.tangjiale/shinanxi/cytotrace1/')
dim(sc)




#########15.cytotrace2 #######
```r
library(Seurat)
sce1 = readRDS("/thinker/3.tangjiale/shinanxi/rds/Overall_shinanxi_New.rds")
table(sce1$celltype)
DimPlot(sce1, label = T)
sce_sub <- sce1[,sce1$celltype %in% c("Transitional activated","Cytotoxic effector","CSMD1+","NK-like",
                                      "Activation-regulated","Naive-like","Cytotoxic memory-like","Terminal-branch")]
library(CytoTRACE2)

cytotrace2_sce <- cytotrace2(sce_sub,
                             is_seurat = TRUE,
                             slot_type = "counts",
                             species = 'human',ncores = 1,seed = 42)

class(cytotrace2_sce)
dir.create("/thinker/3.tangjiale/shinanxi/cytotrace")
dir.create("/thinker/3.tangjiale/shinanxi/cytotrace/plot")
saveRDS(cytotrace2_sce,"/thinker/3.tangjiale/shinanxi/cytotrace/cytotrace2_sce.rds")

annotation <- data.frame(phenotype = sce_sub@meta.data$celltype) %>% 
  set_rownames(., colnames(sce_sub))

plots1 <- plotData(cytotrace2_result = cytotrace2_sce,
                   annotation = annotation,
                   is_seurat = TRUE)

umap_raw <- as.data.frame(sce_sub@reductions$umap@cell.embeddings)
plot_names <- c("CytoTRACE2_UMAP", "CytoTRACE2_Potency_UMAP",
                "CytoTRACE2_Relative_UMAP", "Phenotype_UMAP")
for (plot_name in plot_names) {
  if (!is.null(plots[[plot_name]][[1]]$data)) {
    plots[[plot_name]][[1]]$data$UMAP_1 <- umap_raw$UMAP_1
    plots[[plot_name]][[1]]$data$UMAP_2 <- umap_raw$UMAP_2
  }
}

class(plots$CytoTRACE2_UMAP[[1]])
class(plots$CytoTRACE2_Potency_UMAP[[1]])
class(plots$CytoTRACE2_Relative_UMAP[[1]])
class(plots$Phenotype_UMAP[[1]])

x_limits <- range(umap_raw$UMAP_1, na.rm = TRUE)
y_limits <- range(umap_raw$UMAP_2, na.rm = TRUE)

plots$CytoTRACE2_UMAP[[1]] <- plots$CytoTRACE2_UMAP[[1]] +
  scale_x_continuous(limits = x_limits) +
  scale_y_continuous(limits = y_limits)
plots$CytoTRACE2_UMAP[[1]] <- plots$CytoTRACE2_UMAP[[1]] +
  coord_cartesian(xlim = x_limits, ylim = y_limits)
plots$CytoTRACE2_UMAP

plots$CytoTRACE2_Potency_UMAP[[1]] <- plots$CytoTRACE2_Potency_UMAP[[1]] +
  scale_x_continuous(limits = x_limits) +
  scale_y_continuous(limits = y_limits)
plots$CytoTRACE2_Potency_UMAP[[1]] <- plots$CytoTRACE2_Potency_UMAP[[1]] +
  coord_cartesian(xlim = x_limits, ylim = y_limits)
plots$CytoTRACE2_Potency_UMAP

plots$CytoTRACE2_Relative_UMAP[[1]] <- plots$CytoTRACE2_Relative_UMAP[[1]] +
  scale_x_continuous(limits = x_limits) +
  scale_y_continuous(limits = y_limits)
plots$CytoTRACE2_Relative_UMAP[[1]] <- plots$CytoTRACE2_Relative_UMAP[[1]] +
  coord_cartesian(xlim = x_limits, ylim = y_limits)
plots$CytoTRACE2_Relative_UMAP

plots$Phenotype_UMAP[[1]] <- plots$Phenotype_UMAP[[1]] +
  scale_x_continuous(limits = x_limits) +
  scale_y_continuous(limits = y_limits)
plots$Phenotype_UMAP[[1]] <- plots$Phenotype_UMAP[[1]] +
  coord_cartesian(xlim = x_limits, ylim = y_limits)
plots$Phenotype_UMAP

p1 = plots[[1]]
p2 = plots[[3]]
p3 = plots[[4]]
p4 = plots[[5]]
dir.create("/thinker/3.tangjiale/shinanxi/cytotrace/plot/")
dir.create("/thinker/3.tangjiale/shinanxi/cytotrace/plot/celltype")
ggsave('/thinker/3.tangjiale/shinanxi/cytotrace/plot/celltype/CytoTRACE2_Potency_UMAP.png', p1, width=10, height=8, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/cytotrace/plot/celltype/CytoTRACE2_Potency_UMAP.pdf', p1, width=10, height=8, dpi=300)
ggsave('/thinker/3.tangjiale/shinanxi/cytotrace/plot/celltype/CytoTRACE2_Relative_UMAP.png', p2, width=10, height=8, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/cytotrace/plot/celltype/CytoTRACE2_Relative_UMAP.pdf', p2, width=10, height=8, dpi=300)
ggsave('/thinker/3.tangjiale/shinanxi/cytotrace/plot/celltype/celltype_UMAP.png', p3, width=10, height=8, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/cytotrace/plot/celltype/celltype_UMAP.pdf', p3, width=10, height=8, dpi=300)
ggsave('/thinker/3.tangjiale/shinanxi/cytotrace/plot/celltype/boxplot.png', p4, width=10, height=8, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/cytotrace/plot/celltype/boxplot.pdf', p4, width=10, height=8, dpi=300)

annotation <- data.frame(phenotype = sce_sub@meta.data$orig.ident) %>%
  set_rownames(., colnames(sce_sub))

plots2 <- plotData(cytotrace2_result = cytotrace2_sce,
                   annotation = annotation,
                   is_seurat = TRUE)

umap_raw <- as.data.frame(sce_sub@reductions$umap@cell.embeddings)
plot_names <- c("CytoTRACE2_UMAP", "CytoTRACE2_Potency_UMAP",
                "CytoTRACE2_Relative_UMAP", "Phenotype_UMAP"
)
for (plot_name in plot_names) {
  if (!is.null(plots[[plot_name]][[1]]$data)) {
    plots[[plot_name]][[1]]$data$UMAP_1 <- umap_raw$UMAP_1
    plots[[plot_name]][[1]]$data$UMAP_2 <- umap_raw$UMAP_2
  }
}

class(plots$CytoTRACE2_UMAP[[1]])
class(plots$CytoTRACE2_Potency_UMAP[[1]])
class(plots$CytoTRACE2_Relative_UMAP[[1]])
class(plots$Phenotype_UMAP[[1]])

x_limits <- range(umap_raw$UMAP_1, na.rm = TRUE)
y_limits <- range(umap_raw$UMAP_2, na.rm = TRUE)

plots$CytoTRACE2_UMAP[[1]] <- plots$CytoTRACE2_UMAP[[1]] +
  scale_x_continuous(limits = x_limits) +
  scale_y_continuous(limits = y_limits)
plots$CytoTRACE2_UMAP[[1]] <- plots$CytoTRACE2_UMAP[[1]] +
  coord_cartesian(xlim = x_limits, ylim = y_limits)
plots$CytoTRACE2_UMAP

plots$CytoTRACE2_Potency_UMAP[[1]] <- plots$CytoTRACE2_Potency_UMAP[[1]] +
  scale_x_continuous(limits = x_limits) +
  scale_y_continuous(limits = y_limits)
plots$CytoTRACE2_Potency_UMAP[[1]] <- plots$CytoTRACE2_Potency_UMAP[[1]] +
  coord_cartesian(xlim = x_limits, ylim = y_limits)
plots$CytoTRACE2_Potency_UMAP

plots$CytoTRACE2_Relative_UMAP[[1]] <- plots$CytoTRACE2_Relative_UMAP[[1]] +
  scale_x_continuous(limits = x_limits) +
  scale_y_continuous(limits = y_limits)
plots$CytoTRACE2_Relative_UMAP[[1]] <- plots$CytoTRACE2_Relative_UMAP[[1]] +
  coord_cartesian(xlim = x_limits, ylim = y_limits)
plots$CytoTRACE2_Relative_UMAP

plots$Phenotype_UMAP[[1]] <- plots$Phenotype_UMAP[[1]] +
  scale_x_continuous(limits = x_limits) +
  scale_y_continuous(limits = y_limits)
plots$Phenotype_UMAP[[1]] <- plots$Phenotype_UMAP[[1]] +
  coord_cartesian(xlim = x_limits, ylim = y_limits)
plots$Phenotype_UMAP

p1 = plots[[1]]
p2 = plots[[3]]
p3 = plots[[4]]
p4 = plots[[5]]
dir.create("/thinker/3.tangjiale/shinanxi/cytotrace/plot/orig.ident")
ggsave('/thinker/3.tangjiale/shinanxi/cytotrace/plot/orig.ident/CytoTRACE2_Potency_UMAP.png', p1, width=10, height=8, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/cytotrace/plot/orig.ident/CytoTRACE2_Potency_UMAP.pdf', p1, width=10, height=8, dpi=300)
ggsave('/thinker/3.tangjiale/shinanxi/cytotrace/plot/orig.ident/CytoTRACE2_Relative_UMAP.png', p2, width=10, height=8, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/cytotrace/plot/orig.ident/CytoTRACE2_Relative_UMAP.pdf', p2, width=10, height=8, dpi=300)
ggsave('/thinker/3.tangjiale/shinanxi/cytotrace/plot/orig.ident/sample_UMAP.png', p3, width=10, height=8, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/cytotrace/plot/orig.ident/sample_UMAP.pdf', p3, width=10, height=8, dpi=300)
ggsave('/thinker/3.tangjiale/shinanxi/cytotrace/plot/orig.ident/boxplot.png', p4, width=10, height=8, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/cytotrace/plot/orig.ident/boxplot.pdf', p4, width=10, height=8, dpi=300)
```




#########16.GeneSwitches#######
```r
list.of.packages <- c("SingleCellExperiment", "Biobase", "fastglm", "ggplot2", "monocle",
                      "plyr", "RColorBrewer", "ggrepel", "ggridges", "gridExtra", "devtools",
                      "mixtools")

## for package "fastglm", "ggplot2", "plyr", "RColorBrewer", "ggrepel", "ggridges", "gridExtra", "mixtools"
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## for package "SingleCellExperiment", "Biobase"
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

devtools::install_github("SGDDNB/GeneSwitches")
library(GeneSwitches)
library(SingleCellExperiment)

scRNA = readRDS("/thinker/3.tangjiale/shinanxi/monocle2/subset_ob.rds")
dim(scRNA)
cds_DGT = readRDS("/thinker/3.tangjiale/shinanxi/monocle2_OE/subcell_cds_DGT_Pseudotime.rds")
dim(cds_DGT)
cardiac_monocle2 = cds_DGT
plot_monocle_State(cardiac_monocle2)


logexpdata <- as.matrix(GetAssayData(scRNA, assay = "RNA", slot = "data"))

monocle_cells <- colnames(cds_DGT)
seurat_cells <- colnames(logexpdata)
length(monocle_cells)
length(monocle_cells)

common_cells <- intersect(monocle_cells, seurat_cells)
length(common_cells)

cds_DGT <- cds_DGT[, common_cells]
logexpdata <- logexpdata[, common_cells]

GeneSwitches::plot_monocle_State(cds_DGT)


p1 = monocle::plot_cell_trajectory(cds_DGT, color_by = "State")


##
## Input log-normalized gene expression, Monocle2 pseudo-time and dimensionality reduction
## Path1 containing cells in states 3,2,1
sce_p1 <- convert_monocle2(monocle2_obj = cardiac_monocle2, 
                           states = c(2,1), expdata = logexpdata)
## Path2 containing cells in states 3,2,5
sce_p2 <- convert_monocle2(monocle2_obj = cardiac_monocle2, 
                           states = c(2,3,7,4,5,6), expdata = logexpdata)

# Trajectory 1
exp_matrix <- assay(sce_p1, "expdata")

genes_to_keep <- rowSums(exp_matrix > 0) >= 10 &
  rowMeans(exp_matrix) >= 0 &
  apply(exp_matrix, 1, var) >= 0.01

sce_p1_filtered <- sce_p1[genes_to_keep, ]

sce_p1_binary <- binarize_exp(sce_p1_filtered, ncores = 1)

sce_p1 <- find_switch_logistic_fastglm(sce_p1_binary , downsample = TRUE, show_warning = FALSE)

sg_allgenes <- filter_switchgenes(sce_p1, allgenes = TRUE, topnum = 15)

sg_gtypes <- filter_switchgenes(sce_p1, allgenes = FALSE, topnum = 20,
                                genelists = gs_genelists, genetype = c("Surface proteins", "TFs"))

sg_vis <- rbind(sg_gtypes, sg_allgenes[setdiff(rownames(sg_allgenes), rownames(sg_gtypes)),])

p1 = plot_timeline_ggplot(sg_vis, timedata = sce_p1$Pseudotime, txtsize = 3)

ggsave('/thinker/3.tangjiale/shinanxi/GeneSwitches/plot/21.png', p1, width=10, height=7.5, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/GeneSwitches/plot/21.pdf', p1, width=10, height=7.5, dpi=300)
```

```r
# Trajectory 2
exp_matrix <- assay(sce_p2, "expdata")

genes_to_keep <- rowSums(exp_matrix > 0) >= 10 &
  rowMeans(exp_matrix) >= 0 &
  apply(exp_matrix, 1, var) >= 0.01

sce_p2_filtered <- sce_p2[genes_to_keep, ]

sce_p2_binary <- binarize_exp(sce_p2_filtered, ncores = 1)

sce_p2 <- find_switch_logistic_fastglm(sce_p2_binary , downsample = TRUE, show_warning = FALSE)

sg_allgenes <- filter_switchgenes(sce_p2, allgenes = TRUE, topnum = 15)

sg_gtypes <- filter_switchgenes(sce_p2, allgenes = FALSE, topnum = 20,
                                genelists = gs_genelists, genetype = c("Surface proteins", "TFs"))

sg_vis <- rbind(sg_gtypes, sg_allgenes[setdiff(rownames(sg_allgenes), rownames(sg_gtypes)),])

p1 = plot_timeline_ggplot(sg_vis, timedata = sce_p2$Pseudotime, txtsize = 3)
ggsave('/thinker/3.tangjiale/shinanxi/GeneSwitches/plot/237456.png', p1, width=10, height=7.5, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/GeneSwitches/plot/237456.pdf', p1, width=10, height=7.5, dpi=300)
```

```r
# Compare switching genes between the two trajectories
sg_p1_FS <- filter_switchgenes(sce_p1, allgenes = TRUE, r2cutoff = 0.25)
sg_p2_FS <- filter_switchgenes(sce_p2, allgenes = TRUE, r2cutoff = 0)

sg_com <- common_genes(sg_p1_FS, sg_p2_FS, r2cutoff = 0,
                       path1name = "path1", path2name = "path2")
common_genes_plot(sg_com, timedata = sce_p1$Pseudotime)

sg_disgs <- distinct_genes(sg_p1_FS, sg_p2_FS, r2cutoff = 0.1,
                           path1name = "path1", path2name = "path2",
                           path1time = sce_p1$Pseudotime, path2time = sce_p2$Pseudotime)
p2 = plot_timeline_ggplot(sg_disgs, timedata = sce_p1$Pseudotime, color_by = "Paths",
                          iffulltml = FALSE, txtsize = 3)
ggsave('/thinker/3.tangjiale/shinanxi/GeneSwitches/plot/branch.png', p2, width=13, height=7.5, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/GeneSwitches/plot/branch.pdf', p2, width=13, height=7.5, dpi=300)

sg_disgs_scale <- distinct_genes(sg_p1_FS, sg_p2_FS, r2cutoff = 0.1,
                                 path1name = "path1", path2name = "path2",
                                 path1time = sce_p1$Pseudotime, path2time = sce_p2$Pseudotime,
                                 scale_timeline = T, bin = 100)
p3 = plot_timeline_ggplot(sg_disgs_scale, timedata = 1:100, color_by = "Paths",
                          iffulltml = FALSE, txtsize = 3)

ggsave('/thinker/3.tangjiale/shinanxi/GeneSwitches/plot/scale_branch.png', p3, width=13, height=7.5, dpi=600)
ggsave('/thinker/3.tangjiale/shinanxi/GeneSwitches/plot/scale_branch.pdf', p3, width=13, height=7.5, dpi=300)
```









#############17.SCENIC###################
#!/usr/bin/env Rscript
Sys.setenv(RETICULATE_PYTHON = "/data/software/conda_envs/scrna_envs/Scenic/bin/python")
# construct the gene regulation network from the scRNA-seq
# expression matrix with TF annotations
# this work is based on the manuscript:
# Aibar et al. (2017) SCENIC: single-cell regulatory network inference and clustering.
# Nature Methods. doi: 10.1038/nmeth.4463.

#' fast parallel block-based cor function for large gene expression matrix adapted from
#' https://rmazing.wordpress.com/2013/02/22/bigcor-large-correlation-matrices-in-r/
#'
#' @param x the gene expression matrix, where sample/cell is row and gene in on the column.
#' @param y NULL (default) or a vector, matrix or data frame with compatible dimensions to x. The default is equivalent to y = x (but more efficient).
#' @param size chunck size for chunck calculation and parallization.
#' @param cores the number of threads used to run
#' @param fun the function used to call.
#' @param method the method for cor, options can be peason, spearman,kendrall. pearson as default.
#' @param verbose
#'
bigcor <- function(
    x,
    y = NULL,
    # fun = c("cor", "cov"),
    size = 2000,
    cores = 8,
    verbose = TRUE,
    ...
) {
  tictoc::tic()
  # if (fun == "cor") FUN <- cor else FUN <- cov
  # if (fun == "cor") STR <- "Correlation" else STR <- "Covariance"
  if (!is.null(y) & NROW(x) != NROW(y)) stop("'x' and 'y' must have compatible dimensions!")
  
  NCOL <- ncol(x)
  if (!is.null(y)) YCOL <- NCOL(y)
  
  ## calculate remainder, largest 'size'-divisible integer and block size
  REST <- NCOL %% size
  LARGE <- NCOL - REST
  NBLOCKS <- NCOL %/% size
  
  ## preallocate square matrix of dimension
  ## ncol(x) in 'ff' single format
  if (is.null(y)) {
    resMAT <- ff::ff(vmode = "double", dim = c(NCOL, NCOL))
  }else{
    resMAT <- ff::ff(vmode = "double", dim = c(NCOL, YCOL))
  }
  
  ## split column numbers into 'nblocks' groups + remaining block
  GROUP <- rep(1:NBLOCKS, each = size)
  if (REST > 0){
    GROUP <- c(GROUP, rep(NBLOCKS + 1, REST))
  }
  SPLIT <- split(1:NCOL, GROUP)
  
  ## create all unique combinations of blocks
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)
  if (!is.null(y)) COMBS <- cbind(1:length(SPLIT), rep(1, length(SPLIT)))
  
  require(doMC)
  ncore = min(future::availableCores(), cores)
  doMC::registerDoMC(cores = ncore)
  ## iterate through each block combination, calculate correlation matrix
  ## between blocks and store them in the preallocated matrix on both
  ## symmetric sides of the diagonal
  results <- foreach(i = 1:nrow(COMBS)) %dopar% {
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    ## if y = NULL
    if (is.null(y)) {
      if (verbose) message("bigcor: ", sprintf("#%d:Block %s and Block %s (%s x %s) ... ",
                                               i, COMB[1], COMB[2], length(G1),  length(G2)))
      flush.console()
      RES<- do.call("cor", list(x = x[, G1], y = x[, G2], ... ))
      # RES <- FUN(x[, G1], x[, G2], ...)
      resMAT[G1, G2] <- RES
      resMAT[G2, G1] <- t(RES)
    } else {## if y = smaller matrix or vector
      if (verbose) message("bigcor: ", sprintf("#%d:Block %s and 'y' (%s x %s) ... ",
                                               i, COMB[1], length(G1),  YCOL))
      flush.console()
      RES<- do.call("cor", list(x = x[, G1], y = y, ... ))
      # RES <- FUN(x[, G1], y, ...)
      resMAT[G1, ] <- RES
    }
  }
  
  if ( is.null(y) ){
    resMAT <- resMAT[1:ncol(x),1:ncol(x)]
    colnames(resMAT) <- colnames(x)
    rownames(resMAT) <- colnames(x)
  }else{
    resMAT <- resMAT[1:ncol(x),1:ncol(y)]
    colnames(resMAT) <- colnames(x)
    rownames(resMAT) <- colnames(y)
  }
  tictoc::toc()
  return(resMAT)
}

write.gmt <- function(geneSet=kegg2symbol_list,gmt_file='kegg2symbol.gmt'){
  sink( gmt_file )
  for (i in 1:length(geneSet)){
    cat(names(geneSet)[i])
    #cat('\tNA\t')
    cat('\t')
    cat(paste(geneSet[[i]],collapse = '\t'))
    cat('\n')
  }
  sink()
}

FilterGenes <- function (object, min.value=1, min.cells = 0, filter.genes = NULL ) {
  genes.use <- rownames(object)
  
  if (min.cells > 0) {
    num.cells <- Matrix::rowSums( GetAssayData(object,slot="counts") > min.value)
    genes.use <- names(num.cells[which(num.cells >= min.cells)])
    object = subset( object, features = genes.use)
    # object = SetAssayData(object, new.data = GetAssayData(object, slot ="data")[genes.use,])
    # object@data <- object@data[genes.use, ] # Seurat V2.x
  }
  if (!is.null(filter.genes)) {
    filter.genes = CaseMatch(search = filter.genes, match = rownames(seurat_ob))
    genes.use <- setdiff(genes.use, filter.genes) #keep genes not in filter.genes
    object = subset( object, features = genes.use)
    # object = SetAssayData(object, new.data = GetAssayData(object, slot ="data")[genes.use,])
    # object[["RNA"]]@data = object[["RNA"]]@data[genes.use,]
    # object@data <- object@data[genes.use, ] #seurat V2.x
  }
  object <- LogSeuratCommand(object)
  return(object)
}

suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("SCENIC"))
suppressPackageStartupMessages(library("reticulate"))
suppressPackageStartupMessages(library("glue"))
suppressPackageStartupMessages(library("tibble"))
suppressPackageStartupMessages(library("dplyr"))
#suppressPackageStartupMessages(library("OESingleCell"))
#=command line parameters setting
option_list = list(
  make_option( c("--input", "-i" ), type = "character",
               help = "The filtered exprssion matrix in prefered format."),
  make_option( c("--informat", "-f" ), type = "character", default = "seurat",
               help = "The indication of type of input expression matrix, the possible type can be:
                            seurat: the seurat object from the clustering results."),
  make_option( c("--cisTargetdb", "-d"), type = "character",
               help = "[REQUIRED]the glob of downloaded official motif annotation database for specified spieces.
                The files are usually suffixed with .feather, hg19-500bp-upstream-7species.mc9nr.feather as example."),
  make_option( c("--species", "-s"), type = "character",
               help = "[REQUIRED]the spieces abstraction, the current options can be 'hgnc' for human and 'mgi' for mouse."),
  make_option( c("--minCell4gene","-x" ),type="double", default = 0.01,metavar = "minimium proportion",
               help="the minimium cell number one gene detected.If the value is less than 1,
                it will be interpreted as a proportion."),
  make_option( c("--tfs"), type = "character",
               help = "the TF list of interested species in file"),
  make_option( c("--coexMethod"), type = "character",
               help = "the co-expression method for finding regulons. The current supported methods are
                w001,w005,top50,top50perTarget,top10perTarget,top5perTarget."),
  make_option( c("--ncores", "-j"), type = "integer", default = 10,
               help = "[OPTIONAL]the CPUs used to run this job, the more the better for project with more than 10k cells."),
  make_option(c("--downsample", "-e"),type = "character", default = "30000",
              help = "the downsample number of cells "),
  make_option( c("--hvg","-v"),type="logical", default = F,
               help="use the high variable features to do analysis", metavar="high variable features"),
  make_option( c("--extended"),type="logical", default = F,
               help="whether to use the extended regulons for calculation and visualization"),
  make_option( c("--output","-o"),type="character", default = "./",
               help="the output directory.", metavar="outputdir")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# setting the output directory
if ( is.null(opt$output) ){
  print("NO output directory specified,the current directory will be used!")
  output_dir = getwd()
}else{
  if ( file.exists(opt$output) ){
    output_dir = opt$output
  }else{
    output_dir = opt$output
    dir.create(output_dir,recursive = T)
  }
}
output_dir = normalizePath(output_dir )

if ( !is.null(opt$informat) & !is.null(opt$input) ){
  if ( opt$informat == "seurat" ){# the input is a seurat object which may contain more than one sample
    seurat_ob = readRDS( opt$input )###i.e. Overall_shinanxi_New.rds
    # in case of seurat pbject derived from different version, update it to be consensus with the
    # current version in this script
    if ( seurat_ob@version < 3){
      seurat_ob = UpdateSeuratObject(seurat_ob)
    }
  }
}


if ("RNA" %in% names(seurat_ob@assays)) {
  DefaultAssay(seurat_ob)="RNA"
} else if ("Spatial" %in% names(seurat_ob@assays)) {
  DefaultAssay(seurat_ob)="Spatial"
} else {
  stop("No counts data available")
}

if ( opt$extended == "F" ){
  non_extended = FALSE
}else{
  non_extended = TRUE
}

if ( is.null( opt$downsample) ){
  downsample = 30000
}else{
  downsample = opt$downsample
}
if (ncol(seurat_ob) > 40000) {
  # library(sampling) need install,instead with manual function
  ratio <- as.numeric(downsample) / ncol(seurat_ob)
  metadata_temp <- as.data.frame(seurat_ob@meta.data)
  # strata(metadata_temp,stratanames="clusters",ratio,description=FALSE)
  cells_sample <- c()
  for (i in unique(seurat_ob$clusters)) {
    cells_temp <- rownames(metadata_temp)[which(metadata_temp$clusters == i)]
    set.seed(2024)
    cells_temp_sample <- sample(cells_temp, ceiling(length(cells_temp) * ratio), replace = FALSE, prob = NULL)
    cells_sample <- append(cells_sample, cells_temp_sample)
  }
  seurat_ob <- subset(seurat_ob, cells = cells_sample)
  saveRDS(seurat_ob,"sub_cells.rds")
}
print(dim(seurat_ob))


if ( is.null( opt$species ) ){
  stop("NO motif annotation for your specified species.")
}else{
  # human = "hgnc", mouse = "mgi"
  org = tolower(opt$species)
}
scenicOptions <- initializeScenic(org=org, dbDir= opt$cisTargetdb, nCores=opt$ncores)
#scenicOptions@settings$seed <- 123

# do primary filtering of genes if it was not carried out before.
# filter genes by using the minimium cell number one gene is detected
# this step shoud be run after cell filtering
if ( opt$minCell4gene < 1 ){ #the parameter is a percentage
  minCell4gene = round(opt$minCell4gene * ncol(seurat_ob))
}else{ #the parameter is a integer
  minCell4gene = opt$minCell4gene
}

# determine the final expression matrix according to the usage of variable features
# filte genes
# seurat_ob = FilterGenes(seurat_ob, min.cells = minCell4gene, filter.genes = NULL )
if ( opt$hvg ){
  exprMat = GetAssayData(seurat_ob, slot = "counts")[VariableFeatures(seurat_ob),]
}else{
  exprMat = GetAssayData(seurat_ob, slot = "counts")
}

genesKept <- geneFiltering(as.matrix(exprMat),
                           scenicOptions=scenicOptions,
                           minCountsPerGene=1,
                           minSamples=minCell4gene)
exprMat_filtered <- exprMat[genesKept, ] # now the final matrix used for GRN inference


# Step1. Inference of co-expression modules

# load the list of known TFs for specified species from the command line file or annotated database
# R version:run Genies3 to do co-expression network analysis
# but for the high performance we change to python version here
# runGenie3(exprMat_filtered, scenicOptions)
# python version for this step using GRNBoost from arboreto
if ( !is.null(opt$tfs) ){
  arb.util = import("arboreto.utils")
  tf_names = arb.util$load_tf_names(opt$tfs)
}else{
  tf_names = getDbTfs( scenicOptions )
}
tf_names = CaseMatch(search=tf_names, match = rownames(seurat_ob)) # all the tf should be in the genes of matrix
arb.algo = import("arboreto.algo")
#adjacencies = arb.algo$grnboost2(as.data.frame(as.matrix(t(exprMat_filtered))), tf_names=tf_names, verbose=T, seed=123L)
adjacencies = arb.algo$grnboost2(as.data.frame(t(as.matrix(exprMat_filtered))), tf_names=tf_names, verbose=T, seed=123L)
colnames(adjacencies) = c( "TF", "Target", "weight" )
saveRDS( adjacencies, file = getIntName(scenicOptions, "genie3ll") )

# correlation analysis to distinguish positive from negative activity for each TF target
# detect the block size according to the number of genes in count matrix
corrMat = bigcor(t(as.matrix(exprMat_filtered)),size = 2000, cores = opt$ncores, method = "spearman")
saveRDS(corrMat, file= getIntName(scenicOptions, "corrMat"))

### Build and score the GRN
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod = opt$coexMethod)

scenicOptions <- initializeScenic(org=org, dbDir= opt$cisTargetdb, nCores=1)
runSCENIC_3_scoreCells(scenicOptions,log2(as.matrix(exprMat_filtered)+1))

regulonAUC  = loadInt(scenicOptions, "aucell_regulonAUC")
if ( non_extended ){
  regulonAUC = regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
}
regulonAUC_mat = getAUC(regulonAUC)
rownames(regulonAUC_mat) = gsub("_", "-", rownames(regulonAUC_mat))
regulonAUC_mat_out = regulonAUC_mat[-grep(pattern="-extended",rownames(regulonAUC_mat)),]
write.table(as.data.frame(regulonAUC_mat_out) %>% tibble::rownames_to_column(var = "regulon"),
            file.path(output_dir,"regulon_activity.xls"),
            sep = "\t", col.names =T, row.names =F)

seurat_ob[["SCENIC"]] = CreateAssayObject(counts = regulonAUC_mat)
seurat_ob = ScaleData(seurat_ob, assay = "SCENIC")
seurat_ob@tools$RunAUCell = regulonAUC
# Tool(seurat_ob) = regulonAUC
# saveRDSMC(seurat_ob,"SCENIC_seurat.rds")

regulonTargetsInfo = loadInt(scenicOptions, "regulonTargetsInfo")
write.table(regulonTargetsInfo,
            file.path(output_dir, "0.1.TF_target_enrichment_annotation.xls"),
            sep = "\t", col.names =T, row.names =F, quote =F)

regulons <- loadInt(scenicOptions, "regulons")
sub_regulons = gsub(" .*","",rownames(regulonAUC_mat_out))
regulons = regulons[sub_regulons]
write.gmt(regulons, gmt_file = file.path(output_dir, "0.2.regulon_annotation.xls"))

if(!file.exists(file.path(output_dir, "Regulon_Explanation_of_Moderator_Analysis.docx"))){
  file.copy("/public/scRNA_works/Documents/Regulon_Explanation_of_Moderator_Analysis.docx",
            file.path(output_dir, "Regulon_Explanation_of_Moderator_Analysis.docx"))
}
if(!file.exists(file.path(output_dir, "0.3.MotifEnrichment_preview.html"))){
  file.copy("./output/Step2_MotifEnrichment_preview.html",
            file.path(output_dir, "0.3.MotifEnrichment_preview.html"))
}







#!/usr/bin/env Rscript
# Title     : CSI analysis of reuglons
# Objective : Processing and visualization of SCENIC results
# Created by: hanmin
# Created on: 2020/5/20

#' Creates a list of unique color values used for plotting
#'
#' @return A named vector of unique hexedecimal color values, either generated from a preselected
#'         vector of 20 unique colors, or from a sequence of colors in hsv colorspace.
#'
#' @param seurat.obj A singular preprocessed Seurat object
#' @param gradient Setting to TRUE will use a sequence of hsv colors instead of 20 unique colors,
#'                 useful for comparisons of more than 20 cell types.
#' @param value The Seurat metadata slot to generate colors for. Defaults to "celltype".
#'
#' @import SingleCellExperiment
#' @import Seurat
#'
#' @seealso \code{\link{as.SingleCellExperimentList}}
#' @seealso \code{\link{ExtractGenes}}
#' @seealso \code{\link{DecoderVariance}}
#' @seealso \code{\link{MeanDecoderVariance}}
#' @seealso \code{\link{GetCharMetadata}}
#'
#' @examples
#' DimPlot(object = seurat.obj,
#'         reduction = "tsne",
#'         cols = SelectColors(seurat.obj),
#'         group.by = "celltype",
#'         label = TRUE,
#'         repel = TRUE)
#'
#' @export

save_ggplots <- function( filename = NULL,
                          plot = NULL,
                          width = 6,
                          height = 8,
                          dpi = 300,
                          to.pdf = TRUE,
                          to.png = FALSE,
                          ...
) {
  if ( to.png ) {
    ggplot2::ggsave(paste(filename, "png", sep = "."),
                    plot,
                    device = "png",
                    width = width,
                    height = height,
                    units = "in",
                    dpi = dpi,
                    ...
    )
  }
  if ( to.pdf ) {
    ggplot2::ggsave(paste(filename, "pdf", sep = "."),
                    plot,
                    device = "pdf",
                    width = width,
                    height = height,
                    units = "in",
                    dpi = dpi,
                    ...
    )
  }
}

SelectColors <- function(
    object = NULL,
    palette = "blindless",
    value = "celltype",
    n = NULL
){
  if ( !is.null(object) ){
    if ( class(object) == "data.frame" ){
      colid <- ifelse( is.null(value), colnames(object)[1], value )
      if (is.factor(object[[colid]])) {
        names= levels(object[[colid]])
      } else {
        names <- unique(object[[colid]])
      }
    }
    if ( is.factor(object) ){
      names <- levels(object)
    }
    if ( is.vector(object) ){
      names <- unique(object)
    }
    n = length(names)
  }else if ( !is.null(n) ) {
    names = NULL
  }
  
  colors2pick = switch(palette,
                       ##ref: http://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
                       ditto = c("#E69F00", "#56B4E9", "#009E73", "#0072B2",
                                 "#D55E00", "#CC79A7", "#666666", "#AD7700", "#1C91D4",
                                 "#007756", "#D5C711", "#005685", "#A04700", "#B14380",
                                 "#4D4D4D", "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71",
                                 "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C"),
                       CustomCol2 = c(
                         "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
                         "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
                         "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
                         "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
                         "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17"),
                       seurat = hcl( h = seq(15, 375, length = n+1), l = 65, c = 100),
                       col50 =  c("#982f29", "#5ddb53", "#8b35d6", "#a9e047", "#4836be",
                                  "#e0dc33", "#d248d5", "#61a338", "#9765e5", "#69df96",
                                  "#7f3095", "#d0d56a", "#371c6b", "#cfa738", "#5066d1",
                                  "#e08930", "#6a8bd3", "#da4f1e", "#83e6d6", "#df4341",
                                  "#6ebad4", "#e34c75", "#50975f", "#d548a4", "#badb97",
                                  "#b377cf", "#899140", "#564d8b", "#ddb67f", "#292344",
                                  "#d0cdb8", "#421b28", "#5eae99", "#a03259", "#406024",
                                  "#e598d7", "#343b20", "#bbb5d9", "#975223", "#576e8b",
                                  "#d97f5e", "#253e44", "#de959b", "#417265", "#712b5b",
                                  "#8c6d30", "#a56c95", "#5f3121", "#8f846e", "#8f5b5c"),
                       paired = brewer.pal(n = n, 'Paired'),
                       colx22 = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4',
                                  '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
                                  '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
                                  '#000075', '#808080', '#4f34ff', '#f340F0'),
                       jet = c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
                               "#FF7F00", "red", "#7F0000" ),
                       tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
                                     "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
                                     "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
                                     "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5"),
                       tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                                           "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                                           "#CDCC5D", "#6DCCDA"),
                       colorblind10 = c("#006BA4", "#FF800E", "#ABABAB", "#595959",
                                        "#5F9ED1", "#C85200", "#898989", "#A2C8EC",
                                        "#FFBC79", "#CFCFCF"),
                       trafficlight = c("#B10318", "#DBA13A", "#309343", "#D82526",
                                        "#FFC156", "#69B764", "#F26C64", "#FFDD71",
                                        "#9FCD99"),
                       purplegray12 = c("#7B66D2", "#A699E8", "#DC5FBD", "#FFC0DA",
                                        "#5F5A41", "#B4B19B", "#995688", "#D898BA",
                                        "#AB6AD5", "#D098EE", "#8B7C6E", "#DBD4C5"),
                       bluered12 = c("#2C69B0", "#B5C8E2", "#F02720", "#FFB6B0", "#AC613C",
                                     "#E9C39B", "#6BA3D6", "#B5DFFD", "#AC8763", "#DDC9B4",
                                     "#BD0A36", "#F4737A"),
                       greenorange12 = c("#32A251", "#ACD98D", "#FF7F0F", "#FFB977",
                                         "#3CB7CC", "#98D9E4", "#B85A0D", "#FFD94A",
                                         "#39737C", "#86B4A9", "#82853B", "#CCC94D"),
                       cyclic = c("#1F83B4", "#1696AC", "#18A188", "#29A03C", "#54A338",
                                  "#82A93F", "#ADB828", "#D8BD35", "#FFBD4C", "#FFB022",
                                  "#FF9C0E", "#FF810E", "#E75727", "#D23E4E", "#C94D8C",
                                  "#C04AA7", "#B446B3", "#9658B1", "#8061B4", "#6F63BB"),
                       CustomCol = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17",
                                     "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
                                     "#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
                                     "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"),
                       blindless = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F",
                                     "#BF5B17", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                                     "#66A61E", "#E6AB02", "#A6761D", "#A6CEE3", "#1F78B4",
                                     "#B3DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
                                     "#CAB2D6", "#6A3D9A", "#B15928", "#FBB4AE", "#B3CDE3", "#BC80BD",
                                     "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC",
                                     "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9",
                                     "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8",
                                     "#984EA3",  "#FFFF33", "#A65628", "#F781BF", "#999999","#FFED6F",
                                     "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
                                     "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
                                     "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9","#666666"),
  )
  
  if ( is.null(n) ){
    colors_use <- colors2pick
  }else{
    colors_use <- colors2pick[1:n]
  }
  if ( !is.null(names) ){ names(colors_use) <- names }
  return(colors_use)
}


my_color <- function(n){
  my_palette=pal_d3("category20")(10)
  return(my_palette[n])
}

#' Calculates CSI values for all regulon pairs
#'
#' @param regulonAUC_mat The AUC  matrix for all regulons as calculated by SCENIC.
#' @keywords SCENIC, regulons, binary activity, kmeans, thresholds
#' @export
#' @examples
calculate_csi <- function(
    regulonAUC_mat = AUCell::getAUC(regulonAUC),
    calc_extended = FALSE,
    verbose = FALSE
){
  if (calc_extended == FALSE){
    regulonAUC_mat <- subset(regulonAUC_mat,!grepl("extended",rownames(regulonAUC_mat)))
  }
  
  pearson_cor <- cor(t(regulonAUC_mat))
  pearson_cor_df <- as.data.frame(pearson_cor)
  pearson_cor_df$regulon_1 <- rownames(pearson_cor_df)
  pearson_cor_long <- pearson_cor_df %>%
    tidyr::gather(regulon_2,pcc,-regulon_1) %>%
    dplyr::mutate("regulon_pair" = paste(regulon_1,regulon_2,sep="_"))
  
  regulon_names <- unique(colnames(pearson_cor))
  num_of_calculations <- length(regulon_names)*length(regulon_names)
  csi_regulons <- data.frame(matrix(nrow=num_of_calculations,ncol = 3))
  colnames(csi_regulons) <- c("regulon_1", "regulon_2", "CSI")
  num_regulons <- length(regulon_names)
  
  f <- 0
  for(reg in regulon_names){
    for(reg2 in regulon_names){
      f <- f + 1
      # fraction_lower <- calc_csi(reg,reg2,pearson_cor)
      test_cor <- pearson_cor[reg,reg2]
      total_n <- ncol(pearson_cor)
      pearson_cor_sub <- subset(pearson_cor,rownames(pearson_cor) == reg | rownames(pearson_cor) == reg2)
      
      # sums <- apply(pearson_cor_sub,MARGIN = 2, FUN = compare_pcc, pcc = test_cor)
      sums <- apply(pearson_cor_sub, 2,function(m) ifelse( length(m[m>test_cor]) == length(m), 0, length(m)) )
      fraction_lower <- length(sums[sums == nrow(pearson_cor_sub)]) / total_n
      csi_regulons[f,] <- c(reg,reg2,fraction_lower)
    }
  }
  csi_regulons$CSI <- as.numeric(csi_regulons$CSI)
  return(csi_regulons)
}

#' Calculate CSI module activity over all cell types
#'
#' @param clusters_df -
#' @param regulonAUC_mat The AUC  matrix for all regulons as calculated by SCENIC.
#' @param metadata -
#' @keywords SCENIC, regulons, CSI activity
#' @export
#' @examples
calc_csi_module_activity <- function(
    clusters_df,
    regulonAUC_mat = AUCell::getAUC(regulonAUC),
    metadata
){
  cell_types<- unique(metadata$cell_type)
  regulons <- unique(clusters_df$regulon)
  
  regulonAUC_mat <- regulonAUC_mat[as.character(regulons),]
  
  csi_activity_matrix_list <- list()
  csi_cluster_activity <- data.frame("csi_module" = c(),
                                     "mean_activity" = c(),
                                     "cell_type" = c())
  
  cell_type_counter <- 0
  groupby_1 =c()   # select group with only one cell
  for(i in names(table(metadata$cell_type)) ){
    if( table(metadata$cell_type)[i]==1 ){
      groupby_1= i
    }
  }
  cell_types <- setdiff(cell_types, groupby_1)
  
  regulon_counter <-
    for(ct in cell_types) {
      cell_type_counter <- cell_type_counter + 1
      if (unname(table(metadata$cell_type)[ct]) == 1){
        cell_type_aucs = regulonAUC_mat[,rownames(subset(metadata,cell_type == ct))]
      } else {
        cell_type_aucs <- rowMeans(regulonAUC_mat[,rownames(subset(metadata,cell_type == ct))])
      }
      cell_type_aucs_df <- data.frame("regulon" = names(cell_type_aucs),
                                      "activtiy"= cell_type_aucs,
                                      "cell_type" = ct)
      csi_activity_matrix_list[[ct]] <- cell_type_aucs_df
    }
  if(!is.null(groupby_1)){
    regulonAUC_mat_1 = regulonAUC_mat[,rownames(subset(metadata,cell_type == groupby_1))]
    csi_activity_matrix_list_1=data.frame("regulon" = names(regulonAUC_mat_1),
                                          "activtiy"= as.vector(regulonAUC_mat_1),
                                          "cell_type" = groupby_1)
    csi_activity_matrix_list[[groupby_1]] <- csi_activity_matrix_list_1
  }
  for(ct in names(csi_activity_matrix_list)){
    for(cluster in unique(clusters_df$csi_module)){
      csi_regulon <- subset(clusters_df,csi_module == cluster)
      csi_regulon_activtiy <- subset(csi_activity_matrix_list[[ct]],regulon %in% csi_regulon$regulon)
      csi_activtiy_mean <- mean(csi_regulon_activtiy$activtiy)
      this_cluster_ct_activity <- data.frame("csi_module" = cluster,
                                             "mean_activity" = csi_activtiy_mean,
                                             "cell_type" = ct)
      csi_cluster_activity <- rbind(csi_cluster_activity,this_cluster_ct_activity)
    }
  }
  
  csi_cluster_activity[is.na(csi_cluster_activity)] <- 0
  
  csi_cluster_activity_wide <- csi_cluster_activity %>%
    spread(cell_type,mean_activity)
  
  rownames(csi_cluster_activity_wide) <- csi_cluster_activity_wide$csi_module
  csi_cluster_activity_wide <- as.matrix(csi_cluster_activity_wide[2:ncol(csi_cluster_activity_wide)])
  
  return(csi_cluster_activity_wide)
}

#' Plots a heatmap for the connection specificty indices for all regulons.
#'
#' @param csi_df Data frame containing CSI values for all pairwise regulons.
#' @param nclust Number of clusters to divide the heatmap into
#' @param font_size_regulons Font size for regulon names.
#' @keywords SCENIC, regulons, CSI
#' @export
#' @examples
#' @import tidyr
#' @import pheatmap
#' @import viridis
plot_csi_modules <- function(
    csi_df,
    nclust = 10,
    row_anno = NULL,
    font_size_regulons = 6,
    vcolors = "viridis"
){
  ## subset csi data frame based on threshold
  csi_test_mat <- csi_df %>% tidyr::spread(regulon_2,CSI)
  
  future_rownames <- csi_test_mat$regulon_1
  csi_test_mat <- as.matrix(csi_test_mat[,2:ncol(csi_test_mat)])
  rownames(csi_test_mat) <- future_rownames
  
  color_map = SelectColors(object=NULL, n = nclust, palette = "ditto")
  names(color_map) <- 1:nclust
  color_use = list()
  #for (i in names(color_map)){
  #    color_use[["csi_module"]][[i]]=color_map[[i]]
  #}
  color_use[[opt$groupby]] <- color_map  ##different with local script 
  
  pheatmap::pheatmap(csi_test_mat,
                     show_colnames = TRUE,
                     #border_color = NA,
                     color = colorRampPalette(continuous_palette[[vcolors]])(10),
                     fontsize_row = font_size_regulons,
                     fontsize_col = font_size_regulons,
                     angle_col = 90,
                     cutree_cols = nclust,
                     cutree_rows = nclust,
                     cluster_cols = TRUE,
                     cluster_rows = TRUE,
                     annotation_row = row_anno,
                     annotation_colors = color_use,
                     treeheight_row = 20,
                     treeheight_col = 20,
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean",
                     width = 4500,
                     height = 4000,
                     fontface = "bold")
}




# package loading
#
suppressWarnings({
  suppressPackageStartupMessages(library("optparse"))
  #suppressPackageStartupMessages(library("OESingleCell"))
  suppressPackageStartupMessages(library("tidyverse"))
  suppressPackageStartupMessages(library("viridis"))
  suppressPackageStartupMessages(library("Seurat"))
})
source("/public/scRNA_works/pipeline/scRNA-seq_further_analysis/Get_colors.R")
#=command line parameters setting
option_list = list(
  make_option( c("--input", "-i"), type = "character",
               help = "the seurat object saved as R object in RDS format." ),
  make_option( c("--auc", "-v"), type = "character", 
               help = "The regulon activity scores results from the SCENIC output,
		 usually the 3.4_regulonAUC.Rds from SCENIC" ),
  make_option( c("--aucformat", "-f" ), type = "character", default = "rds",
               help = "The indication of type of input AUC, the possible type can be:
              rds: the aucellResults object in RDS format,
	      xsv: the delimited table from AUC.
	      [default \"%default\"]"),
  make_option( c("--groupby", "-c"), type = "character",
               help = "the groupping column of cells in the metadata of seurat object." ),
  make_option( c("--nclust", "-n" ), type="integer", default = 5,
               help="the number of csi modules to use [default \"%default\"]"),
  make_option( c("--extended"),type="logical", default = FALSE,
               help="whether to use the extended regulons for calculation and visualization. [default \"%default\"]"),
  make_option( c("--predicate" ), type = "character", default = NULL, 
               help = "[OPTIONAL]The column name in cell metadata used as identity of each cell combined with which_cell."),
  make_option( c("--output", "-o"), type = "character", default = "./",
               help = "the output directory of QC results.[default \"%default\"]", metavar = "outdir" ),
  make_option(c("--groupby_levels"), type = "character", default = NULL, help = "Tcell,Bcell,NK"),
  make_option(c("--vcolors"), type = "character", default = "viridis",help = "3.2 csi_Color_Palette_Selection"),
  make_option(c("--rowcluster"), type = "character", default="T",help="Whether to cluster the rows, boolean"),
  make_option(c("--colcluster"), type = "character", default="T",help="Whether to cluster the rows, boolean")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if ( opt$extended){
  non_extended = TRUE 
}else{
  non_extended = FALSE 
}

# setting the output directory
if ( is.null(opt$output) ){
  print("NO output directory specified,the current directory will be used!")
  output_dir = getwd()
}else{
  if ( file.exists(opt$output) ){
    output_dir = opt$output
  }else{
    output_dir = opt$output
    dir.create(output_dir,recursive = T)
  }
}
output_dir = normalizePath(output_dir )

seurat_ob = readRDS( opt$input)
if (is.null(opt$groupby_levels)) {
  # 确保是因子型
  if (!is.factor(seurat_ob@meta.data[[opt$groupby]])) {
    seurat_ob@meta.data[[opt$groupby]] <- factor(seurat_ob@meta.data[[opt$groupby]])
  }
  group_levels <- levels(seurat_ob@meta.data[[opt$groupby]])
} else {
  group_levels <- strsplit(opt$groupby_levels, ",")[[1]]
}

## subset 
if (!is.null(opt$predicate)) {
  df <- seurat_ob@meta.data
  desired_cells <- subset(df, eval(parse(text = opt$predicate)))
  seurat_ob <- seurat_ob[, rownames(desired_cells)]
}

if ( seurat_ob@version < 3 ){
  seurat_ob = UpdateSeuratObject(seurat_ob)
}

if ( !is.null(opt[["auc"]])){ # using the user supplied regulon activity score matrix
  if ( tolower(opt$aucformat) == "rds" ){
    auc = readRDS(opt$auc)
    seurat_ob=seurat_ob[,colnames(auc)]
    auc = auc[,rownames(seurat_ob@meta.data)]
    auc = auc[which(rowSums(getAUC(auc))!=0),]
    auc = AUCell::getAUC(auc)
  }else if ( tolower(opt$aucformat) == "xsv" ){
    auc = vroom::vroom(opt$auc, col_names = T) %>% as.data.frame()
    rownames(auc) = auc[,1]
    auc = as.matrix(auc[,-1])
  }
}else {
  auc = Tool(seurat_ob, slot = "RunAUCell")
  auc = auc[,rownames(seurat_ob@meta.data)]
  auc = auc[which(rowSums(getAUC(auc))!=0),]
  auc = AUCell::getAUC(auc)
  #stop("NO regulon activity results supplied!")
}

regulons_csi <- calculate_csi(auc,calc_extended = non_extended)

csi_csi_wide <- regulons_csi %>% tidyr::spread(regulon_2,CSI)
future_rownames <- csi_csi_wide$regulon_1
csi_csi_wide <- as.matrix(csi_csi_wide[,2:ncol(csi_csi_wide)])
rownames(csi_csi_wide) <- future_rownames
regulons_hclust <- hclust(dist(csi_csi_wide,method = "euclidean"))

nclust <- opt$nclust
while (TRUE) {
  clusters <- cutree(regulons_hclust, k = nclust)
  cluster_sizes <- table(clusters)
  if (any(cluster_sizes == 1)) {
    nclust <- nclust - 1
    if (nclust <= 1) {
      stop("The number of clusters has been reduced to 1, and it can't be reduced any further.")
    }
  } else {
    break
  }
}
clusters <- cutree(regulons_hclust,k= nclust)
clusters_df <- data.frame("regulon" = names(clusters),"csi_module" = clusters)
write.table(clusters_df, file.path(output_dir, "3.1.csi_module_annotation.xls"),
            sep = "\t", col.names = T, row.names =F, quote =F)

row_anno <- clusters_df
rownames(row_anno) <- NULL
row_anno$csi_module <- as.character(row_anno$csi_module)
plot = plot_csi_modules(regulons_csi,
                        nclust = nclust,
                        row_anno = row_anno %>% column_to_rownames(var="regulon"),
                        font_size_regulons = 9,vcolors = opt$vcolors)
save_ggplots(file.path(output_dir, paste0("3.2.regulons_csi_correlation_heatmap", collapse = "")),
             plot = plot,
             width = length(unique(regulons_csi[,1])) * 0.15+2,
             height = length(unique(regulons_csi[,1])) * 0.15,
             dpi = 1000,
             limitsize = F,bg="white"
)

cellinfo = seurat_ob@meta.data
cellinfo[,"cell_type"] = cellinfo[,opt$groupby]
csi_cluster_activity_wide <- calc_csi_module_activity(clusters_df,auc, cellinfo)
csi_cluster_activity_wide = csi_cluster_activity_wide[,group_levels]
plot = pheatmap::pheatmap(csi_cluster_activity_wide,
                          show_colnames = TRUE,
                          color = colorRampPalette(continuous_palette[[opt$vcolors]])(10),,
                          cellwidth = 24,
                          cellheight = 24,
                          cluster_rows=as.logical(opt$rowcluster),
                          cluster_cols=as.logical(opt$colcluster),
                          clustering_distance_rows = "euclidean",
                          clustering_distance_cols = "euclidean",
                          fontsize = 15,  
                          fontsize_row = 15,
                          fontsize_col = 15,
                          fontface = "bold")
if (dim(csi_cluster_activity_wide)[1]<7) {
  hig <- 6
} else {
  hig <- dim(csi_cluster_activity_wide)[1]*0.8
}
save_ggplots(file.path(output_dir, paste0("3.3.csi_module_activity_heatmap", collapse = "")),
             plot = plot,
             width = dim(csi_cluster_activity_wide)[2]/2+2,
             height = hig,
             dpi = 1000,
             limitsize = F,bg="white"
)
if(!file.exists(file.path(output_dir, "Regulon_Explanation_of_Moderator_Analysis.docx"))){
  file.copy("/public/scRNA_works/Documents/Regulon_Explanation_of_Moderator_Analysis.docx",
            file.path(output_dir, "Regulon_Explanation_of_Moderator_Analysis.docx"))
}
setwd(output_dir)
print("Convert pdf to png...")
system("module load ImageMagick && ls *.pdf | xargs -I {} -P 8 bash -c '/data/software/ImageMagick/ImageMagick-v7.0.8-14/bin/convert  -verbose -density 200 -trim  ${0}  -quality 100  -flatten  ${0%.pdf}.png' {}")


#!/usr/bin/env Rscript
# Objective : Processing and visualization of SCENIC results

#
rm(list=ls())
save_ggplots <- function( filename = NULL,
                          plot = NULL,
                          width = 6,
                          height = 8,
                          dpi = 300,
                          to.pdf = TRUE,
                          to.png = TRUE,
                          ...
) {
  if ( to.png ) {
    ggplot2::ggsave(paste(filename, "png", sep = "."),
                    plot,
                    device = "png",
                    width = width,
                    height = height,
                    units = "in",
                    dpi = dpi,
                    ...
    )
  }
  if ( to.pdf ) {
    ggplot2::ggsave(paste(filename, "pdf", sep = "."),
                    plot,
                    device = "pdf",
                    width = width,
                    height = height,
                    units = "in",
                    dpi = dpi,
                    ...
    )
  }
}

## RAS
RemoveOutlier <- function(
    metric,
    nmads = 5,
    type = c("both", "lower", "higher"),
    log = FALSE,
    subset = NULL,
    batch = NULL,
    min_diff = NA
) {
  if (log) {
    metric <- log10(metric)
  }
  if (any(is.na(metric))) {
    warning("missing values ignored during outlier detection")
  }
  
  if (!is.null(batch)) {
    N <- length(metric)
    if (length(batch) != N) {
      stop("length of 'batch' must equal length of 'metric'")
    }
    
    # Coercing non-NULL subset into a logical vector.
    if (!is.null(subset)) {
      new.subset <- logical(N)
      names(new.subset) <- names(metric)
      new.subset[subset] <- TRUE
      subset <- new.subset
    }
    
    # Computing QC metrics for each batch.
    by.batch <- split(seq_len(N), batch)
    collected <- logical(N)
    all.threshold <- vector("list", length(by.batch))
    for (b in seq_along(by.batch)) {
      bdx <- by.batch[[b]]
      current <- Recall(metric[bdx], nmads = nmads, type = type, log = FALSE, subset = subset[bdx], batch = NULL, min_diff = min_diff)
      all.threshold[[b]] <- attr(x, "thresholds")
      collected[bdx] <- current
    }
    
    all.threshold <- do.call(cbind, all.threshold)
    colnames(all.threshold) <- names(by.batch)
    # return(.store_thresholds(collected, all.threshold, logged=log))
    if ( log ){ val <- 10^all.threshold }
    attr(collected, "thresholds") <- val
    return( collected )
  }
  # Computing median/MAD (possibly based on subset of the data).
  if (!is.null(subset)) {
    submetric <- metric[subset]
    if (length(submetric) == 0L) {
      warning("no observations remaining after subsetting")
    }
  } else {
    submetric <- metric
  }
  cur.med <- median(submetric, na.rm = TRUE)
  cur.mad <- mad(submetric, center = cur.med, na.rm = TRUE)
  
  diff.val <- max(min_diff, nmads * cur.mad, na.rm = TRUE)
  upper.limit <- cur.med + diff.val
  lower.limit <- cur.med - diff.val
  
  type <- match.arg(type)
  if (type == "lower") {
    upper.limit <- Inf
  } else if (type == "higher") {
    lower.limit <- -Inf
  }
  
  kx = metric < lower.limit | upper.limit < metric
  val = c(lower=lower.limit, higher=upper.limit)
  if ( log ){
    val <- 10^val
  }
  attr(kx, "thresholds") <- val
  return( kx )
  
  # .store_thresholds <- function(x, val, logged=FALSE) {
  #     if (logged) val <- 10^val
  #     attr(x, "thresholds") <- val
  #     x
  # }
}
CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}
## RSS
BinaryCount <- function(
    object,
    method = c("kmeans", "aucell" ),
    auc = NULL,
    nCores = 10,
    do.tfidf = FALSE,
    ...
){
  if ( !is.null(auc) ){
    if ( class(auc) == "matrix" ){
      AUC_mat <- auc
    }else if ( class(auc) == "aucellResults"){
      regulonAUC <- auc
      AUC_mat <- AUCell::getAUC(regulonAUC)
    }
  }else{
    regulonAUC <- Tool(object, slot = "RunAUCell")
    if ( is.null(regulonAUC) ){
      stop("NO regulon AUC matrix supplied or found in the object!")
    }else{
      AUC_mat <- AUCell::getAUC(regulonAUC)
    }
  }
  
  binary_mat <- switch (tolower(method),
                        "aucell" = {
                          cells_AUCellThresholds <- AUCell::AUCell_exploreThresholds(regulonAUC,
                                                                                     smallestPopPercent=0.25,
                                                                                     assignCells=TRUE, plotHist=FALSE,
                                                                                     verbose=FALSE, nCores=nCores)
                          # Get cells assigned to each regulon
                          cellsAssigned <- AUCell::getAssignments(cells_AUCellThresholds)
                          # cellsAssigned   <- lapply(cells_AUCellThresholds, function(x) x$assignment)
                          assignmentTable <- reshape2::melt(cellsAssigned, value.name="cell")
                          colnames(assignmentTable)[2] <- "geneSet"
                          assignmentMat <- table(assignmentTable[,"geneSet"], assignmentTable[,"cell"])
                          binary_mat <- tidyr::spread(as.data.frame.table(assignmentMat),key = "Var2","Freq" ) %>%
                            tibble::column_to_rownames(var = "Var1") %>% as.matrix()
                        },
                        "kmeans"= {
                          # Iterate over each regulon in the AUC matrix
                          AUC_df <- as.data.frame(AUC_mat) %>%
                            tibble::rownames_to_column(var = "regulon") %>%
                            tidyr::gather("cell", "auc", -regulon) %>%
                            dplyr::filter( auc > 0) %>%
                            dplyr::group_by( regulon ) %>%
                            dplyr::mutate( cluster = as.factor(kmeans(auc, centers = 2)$cluster)) %>%
                            dplyr::ungroup()
                          kmeans_thresholds <- lapply(split(AUC_df,as.factor(AUC_df$regulon),drop = F), function(df){
                            cluster1_max <- max(subset(df,cluster == 1)$auc)
                            cluster2_max <- max(subset(df,cluster == 2)$auc)
                            
                            if(cluster1_max > cluster2_max){
                              df <- df %>% mutate("cluster" = gsub(2,3,cluster)) %>%
                                mutate("cluster" = gsub(1,2,cluster)) %>%
                                mutate("cluster" = gsub(3,1,cluster))
                            }
                            
                            df <- df %>% arrange(desc(auc))
                            df_sub <- df %>% subset(cluster == 1)
                            auc_thresholds <- df_sub[1,]$auc
                          })
                          binary_mat <- as.data.frame(AUC_mat) %>%
                            tibble::rownames_to_column(var = "regulon") %>%
                            tidyr::gather("cells", "auc", -regulon) %>%
                            dplyr::group_by( regulon ) %>%
                            mutate("values"= if_else(auc >= kmeans_thresholds[regulon],1,0)) %>%
                            dplyr::ungroup() %>% select(-auc) %>%
                            tidyr::spread("cells", "values") %>%
                            tibble::column_to_rownames(var = "regulon") %>% as.matrix()
                        })
  
  object <- subset(object, cells = colnames(binary_mat) )
  row.names(binary_mat) = gsub("_", "-",rownames(binary_mat) )
  object[["SCENIC"]] <- CreateAssayObject(data = binary_mat)
  if ( do.tfidf ){
    object[["SCENIC"]] <- SetAssayData(object[["SCENIC"]], slot = "data",
                                       new.data =as(TF.IDF(binary_mat), 'dgCMatrix') )
  }
  
  object = LogSeuratCommand(object)
  return( object )
}


##' Calculates Regulon specificity score (RSS) from binary regulon activity.
##' Iterate over all cell types and perform jensen shannon divergence test
##' using binary regulon activity and genotype
##'
##' @param object the Seurat object with binarized assay "SCENIC".
##' @param assay the SCENIC assay as default.
##' @param slot "data" as default.
##' @param metadata Dataframe containing metadata about cells. Has to create a column named cell_type that assigns groupings to cells.
##' Can be the meta.data slot from a Seurat object.
##' @param binary_regulons Data frame with binary regulons, where regulons are rows and columns are cells. Can be created from output of binarize_regulons().
##' @param group.by the groupping factor for RSS calculaion.
##' @import  philentropy JSD
##' @keywords SCENIC, regulons, binary activity, kmeans, thresholds
##' @export
##' @examples
RunRSS <- function(
    object,
    assay = "SCENIC",
    slot = "data",
    binary_regulons = NULL,
    metadata = NULL,
    group.by = "cell_type"
){
  if ( !is.null( binary_regulons) ){
    regulons <- rownames(binary_regulons)
  }else{
    binary_regulons <- Seurat::GetAssayData(object, assay = assay, slot = slot)
    regulons <- rownames(binary_regulons)
  }
  
  if ( is.null(metadata) ){
    metadata <- object@meta.data
  }
  cell_types <- unique(metadata[,group.by])
  jsd_matrix_ct <- data.frame("regulon" = c(), "cell_type" = c(), "jsd" = c())
  
  cell_type_counter <- 0
  for(ct in unique(cell_types)) {
    cell_type_counter <- cell_type_counter + 1
    print(paste("Processing cell type:",cell_type_counter,ct,sep=" "))
    for(regulon_no in 1:length(regulons)) {
      regulon <- regulons[regulon_no]
      regulon_vec <- binary_regulons[regulon,]
      regulon_vec_sum <- sum(regulon_vec)
      ## Check that there are cells with binary activity > 0 for this regulon
      if(regulon_vec_sum > 0){
        #progress(regulon_no)
        regulon_norm <- regulon_vec/regulon_vec_sum
        genotype_vec <- metadata[colnames(binary_regulons),]
        genotype_vec <- genotype_vec %>%
          mutate("cell_class" = if_else(get(group.by) == ct,1,0))
        genotype_vec <- genotype_vec$cell_class
        genotype_norm <- genotype_vec/sum(genotype_vec)
        dist_df <- rbind(regulon_norm,genotype_norm)
        ## Calculate the Jensen-Shannon divergence
        jsd_divergence <- suppressMessages(philentropy::JSD(dist_df))
        ## Calculate Jensen-Shannon distance
        rss <- 1-sqrt(jsd_divergence)
        regulon_jsd <- data.frame("regulon" = regulon, "cell_type" = ct, "RSS" = rss[1])
        jsd_matrix_ct <- rbind(jsd_matrix_ct,regulon_jsd)
      }else if(regulon_vec_sum == 0){
        print(paste("Filtered out:",regulon,". No cells with binary activity > 0 identified. Please check your threshold for this regulon!",sep=""))
      }
    }
  }
  
  jsd_matrix_ct <- jsd_matrix_ct %>% dplyr::arrange(desc(RSS))
  jsd_matrix_ct[,group.by] <- jsd_matrix_ct$cell_type
  jsd_matrix_ct <- jsd_matrix_ct %>% dplyr::select(-cell_type)
  return(jsd_matrix_ct)
}

##' Calculates Regulon specificity score (RSS) from binary regulon activity.
##'
##' @param rrs_df Data frame containing RSS scores for all regulons over all cell types. Can be created with calculate_rrs.
##' @param cell_type Cell type for which to plot jsd ranking. Select "all" to plot a facet plot over all cell types.
##' @param ggrepel_force same as the force parameter for geom_text_repel.
##' @param ggrepel_point_padding same as the force parameter for geom_text_repel.
##' @param top_genes Number of top genes to label in plot using ggrepel.
##' @keywords SCENIC, regulons, RRS, cell type classification
##' @export
##' @examples
##'
## Plot JSD enrichment plot for specific cell type
RSSRanking <- function(
    rss_df,
    group.by,
    ggrepel_force = 2,
    ggrepel_point_padding = 0.5,
    top_genes = 3,
    plot_extended = FALSE,
    size = 4
){
  require(ggrepel)
  require(cowplot)
  require(dplyr)
  require(ggplot2)
  
  # 筛选是否显示 extended regulons
  if(plot_extended){
    rss_df <- rss_df %>% subset(grepl("extended", regulon))
  } else {
    rss_df <- rss_df %>% subset(!grepl("extended", regulon))
  }
  
  # 排名
  rss_df_sub <- rss_df %>%
    dplyr::group_by(.data[[group.by]]) %>%
    mutate(rank = order(order(RSS, decreasing = TRUE)))
  
  # 绘图
  rrs_ranking_plot <- ggplot(rss_df_sub, aes(rank, RSS, label = regulon)) +
    geom_point(color = "grey20", size = 2) +
    geom_point(data = subset(rss_df_sub, rank <= top_genes), color = "red", size = 2) +
    geom_text_repel(
      data = subset(rss_df_sub, rank <= top_genes),
      force = ggrepel_force,
      point.padding = ggrepel_point_padding,
      size = size # 标签字体大小，可调
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.border = element_rect(color = "black", size = 1.5),
      panel.grid = element_blank(),
      strip.text = element_text(size = 12, face = "bold"),       # facet标题加粗
      axis.text = element_text(face = "bold", size = 12),        # 坐标轴刻度加粗
      axis.title = element_text(face = "bold", size = 12),       # 坐标轴标题加粗
      panel.spacing = unit(1.5, "lines"),
      plot.margin = margin(10, 40, 10, 10)
    ) +
    labs(x = "Rank", y = "RSS", title = group.by) +
    facet_wrap(as.formula(paste("~", group.by)), ncol = 2, scales = "free_y")+
    scale_y_continuous(
      expand = expansion(mult = c(0.1, 0.2))  # 增加顶部留白（20%）
    )
  rrs_ranking_plot = rrs_ranking_plot + oe_theme(base_size = 12,base_family = "sans")
  
  return(rrs_ranking_plot)
}

#
# package loading
#
suppressWarnings({
  #========import packages
  suppressPackageStartupMessages( library("Seurat") )
  suppressPackageStartupMessages(library("optparse"))
  suppressPackageStartupMessages(library("tidyverse"))
  #suppressPackageStartupMessages( library("OESingleCell") )
  suppressPackageStartupMessages( library("dplyr") )
  suppressPackageStartupMessages( library("SCENIC") )
  suppressPackageStartupMessages( library("ggplot2") )
  suppressPackageStartupMessages( library("circlize") )
  suppressPackageStartupMessages( library("ComplexHeatmap") )
})

source("/public/scRNA_works/pipeline/scRNA-seq_further_analysis/Get_colors.R")
source("/gpfs/oe-scrna/pipeline/scRNA-seq_further_analysis/function/seuratFc.r")
source("/gpfs/oe-scrna/guokaiqi/test/report_plot/oe_theme.R")
#=command line parameters setting
option_list = list(
  make_option( c("--input", "-i"), type = "character",
               help = "the seurat object saved as R object in RDS format." ),
  make_option( c("--auc", "-v"), type = "character",  default = NULL, 
               help = "The regulon activity scores results from the SCENIC output." ),
  make_option( c("--aucformat", "-f" ), type = "character", default = "rds",
               help = "The indication of type for input AUC, the choices can be:
                            rds: the aucellResults object in RDS format,
                                 usually the 3.4_regulonAUC.Rds from SCENIC;
                            xsv: the delimited table from AUC."),
  make_option( c("--topGenes", "-t"), type = "integer", default = NULL,
               help = "Number of top genes to label in plot of rrs ranking." ),
  make_option( c("--groupby", "-c"), type = "character",
               help = "the groupping column of cells in the metadata of seurat object." ),
  make_option( c("--binmethod", "-m"), type = "character", default = "aucell",
               help = "the binary methods used to binarize the regulon activity matrix element to 0/1.
            Options can be aucell and kmeans. Note that 'aucell' is only avaiable for aucellResults object" ),
  make_option( c("--ncores", "-j" ), type="integer", default = 10,
               help="the number of CPUs used to improve the performace."),
  make_option( c("--threshold", "-s" ), type="double",
               help="subset the regulon according to the threshold of RSS"),
  make_option( c("--extended"),type="logical", default = FALSE,
               help="whether to use the extended regulons for calculation and visualization"),
  make_option( c("--predicate" ), type = "character", default = NULL, 
               help = "[OPTIONAL]The column name in cell metadata used as identity of each cell combined with which_cell."),
  make_option( c("--output", "-o"), type = "character", default = "./",
               help = "the output directory of QC results.", metavar = "outdir" ),
  make_option( c("--use_color_anno" ), type = "logical",  default = TRUE,
               help = "[Optional]Whether to use the color information commented in the RDS, default is to use it, if not, the color will be automatically re-annotated."),
  make_option( c("--color_file" ), type = "character",  default = NULL,
               help = "Optional, enter a tab-separated file where the first column's header is the metadata column name, the first column contains the elements of that column, and the second column contains the corresponding colors."),
  make_option( c("--palette" ), type = "character",  default = "customecol2",
               help = "[Optional] Fill in if needed, specify the name of the discrete color palette in Get_colors.R."),
  make_option(c("--ccolors"), type = "character", default = "rdwibl",help = "Continuous Color Palette Selection"),
  make_option(c("--rowcluster"), type = "character", default="T",help="Whether to cluster the rows, boolean"),
  make_option(c("--colcluster"), type = "character", default="T",help="Whether to cluster the rows, boolean"),
  make_option(c("--groupby_levels"), type = "character", default = NULL, help = "Tcell,Bcell,NK")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if ( is.null(opt$topGenes) ){
  topGenes = 3
}else{
  topGenes = opt$topGenes
}

if ( opt$extended ){
  non_extended = TRUE 
}else{
  non_extended = FALSE 
}

# setting the output directory
if ( is.null(opt$output) ){
  print("NO output directory specified,the current directory will be used!")
  output_dir = getwd()
}else{
  if ( file.exists(opt$output) ){
    output_dir = opt$output
  }else{
    output_dir = opt$output
    dir.create(output_dir,recursive = T)
  }
}
output_dir = normalizePath(output_dir )

seurat_ob = readRDS( opt$input)
#
if (is.null(opt$groupby_levels)) {
 
  if (!is.factor(seurat_ob@meta.data[[opt$groupby]])) {
    seurat_ob@meta.data[[opt$groupby]] <- factor(seurat_ob@meta.data[[opt$groupby]])
  }
  group_levels <- levels(seurat_ob@meta.data[[opt$groupby]])
} else {
  group_levels <- strsplit(opt$groupby_levels, ",")[[1]]
}
if ( ! opt$use_color_anno ){
  seurat_ob@meta.data = seurat_ob@meta.data[ ,!grepl(paste0("^",opt$groupby,"_col$" ), colnames(seurat_ob@meta.data))]
}
if ( !is.null(opt$color_file)){
  color_file = read.delim(opt$color_file, sep="\t", header = T)
  meta_anno = color_anno(seurat_ob@meta.data, color_file)
} else {
  meta_anno = seurat_ob@meta.data
}
color_use = get_colors(meta_anno, opt$groupby, opt$palette)
seurat_ob = AddMetaData( seurat_ob, metadata = color_use[["object_meta"]])
color_map = color_use[["new_celltype_pal"]]
color_map = na.omit(color_map)
color_use = list()
color_use[[opt$groupby]] <- color_map

##
## subset 
if (!is.null(opt$predicate)) {
  df <- seurat_ob@meta.data
  desired_cells <- subset(df, eval(parse(text = opt$predicate)))
  seurat_ob <- seurat_ob[, rownames(desired_cells)]
}

if ( seurat_ob@version < 3 ){
  seurat_ob = UpdateSeuratObject(seurat_ob)
}

if ( length(levels(seurat_ob@meta.data[[opt$groupby]])) != length(unique(seurat_ob@meta.data[[opt$groupby]])) ) {
  seurat_ob@meta.data[[opt$groupby]] = factor(seurat_ob@meta.data[[opt$groupby]], levels= sort(unique(seurat_ob@meta.data[[opt$groupby]])) )
}

if ( !is.null(opt[["auc"]])){ # using the user supplied regulon activity score matrix
  if ( tolower(opt$aucformat) == "rds" ){
    regulonAUC = readRDS(opt$auc)
    seurat_ob = seurat_ob[,colnames(regulonAUC)]
    regulonAUC = regulonAUC[,rownames(seurat_ob@meta.data)]
    regulonAUC = regulonAUC[which(rowSums(AUCell::getAUC(regulonAUC))!=0),]
  }else if ( tolower(opt$aucformat) == "xsv" ){
    regulonAUC = vroom::vroom(opt$auc, col_names = T) %>% as.data.frame()
    rownames(regulonAUC) = regulonAUC[,1]
    regulonAUC = as.matrix(regulonAUC[,-1])
  }
}else { # use the AUC result integrated in the Seurat object by running RunAUCell
  regulonAUC = Tool(seurat_ob, slot = "RunAUCell")
  regulonAUC = regulonAUC[,rownames(seurat_ob@meta.data)]
  regulonAUC = regulonAUC[which(rowSums(AUCell::getAUC(regulonAUC))!=0),]
}


#
# RAS
#
# matrix prep
regulonAUC= regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonAUC_mat = AUCell::getAUC(regulonAUC)
rownames(regulonAUC_mat) = gsub("_", "-", rownames(regulonAUC_mat))
regulonAUC_mat_out = regulonAUC_mat[grep(pattern="-extended",rownames(regulonAUC_mat),invert=T),]
## groupby 
cellInfo = seurat_ob@meta.data
col_anno = as.data.frame(seurat_ob@meta.data) %>% rownames_to_column(var="barcodes")
if ( dim(cellInfo)[1] > 65536 ){
  col_anno = col_anno[,c("barcodes",opt$groupby)] %>% group_by(get(opt$groupby)) %>% sample_frac(65536/length(col_anno$barcodes))
}else{
  col_anno = col_anno[,c("barcodes",opt$groupby)]
}

col_anno = col_anno %>% arrange(get(opt$groupby)) %>% column_to_rownames(var="barcodes")
col_anno[[opt$groupby]]<- factor(col_anno[[opt$groupby]],levels = group_levels)
annotation_colors <- list(
  setNames(getColorOrder(seurat_ob, opt$groupby), 
           group_levels)
)
names(annotation_colors) <- opt$groupby
regulonAUC_plotdata = regulonAUC_mat_out[,rownames(col_anno)]
col_order <- order(col_anno[[opt$groupby]])   # 或者按 levels 排序
col_anno <- col_anno[col_order, ]
#col_list <- setNames(annotation_colors[[opt$groupby]], levels(col_anno[[opt$groupby]]))
col_list <- setNames(annotation_colors[[opt$groupby]], levels(col_anno))
annot_df <- data.frame(
  group = col_anno
)
names(annot_df) <- opt$groupby
annot_colors <- list()
annot_colors[[opt$groupby]] <- annotation_colors[[opt$groupby]]
col_ha <- HeatmapAnnotation(
  df = annot_df,
  col = annot_colors,
  annotation_name_gp = gpar(fontsize = 12),
  annotation_legend_param = list(
    labels_gp = gpar(fontsize = 12),
    title_gp = gpar(fontsize = 12)
  )
)
regulonAUC_scaled <- t(scale(t(regulonAUC_plotdata)))

ht <- Heatmap(
  regulonAUC_scaled,
  name = "AUC",
  col = colorRamp2(
    breaks = c(seq(-2.5, 0, length.out = 100), seq(0, 2.5, length.out = 100)),
    colors = colorRampPalette(continuous_palette[[opt$ccolors]])(200)
  ),show_row_dend = FALSE,
  top_annotation = col_ha,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 10), 
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 12),
    labels_gp = gpar(fontsize = 12)
  )
)
max_length <- max(nchar(group_levels))
width <- if (max_length < 15) {
  8 
} else {
  10 
}
pdf(file.path(output_dir, "1.1.regulon_activity_heatmap_groupby_cells.pdf"),
    width = width, height = 6)
draw(ht)
dev.off()
png(file.path(output_dir, "1.1.regulon_activity_heatmap_groupby_cells.png"),
    width = width, height = 6, units = "in", res = 300)
draw(ht)
dev.off()
detach("package:ComplexHeatmap", unload = TRUE)


# to remove outliers before mean calculation
groupby_1 =c()   # select group with only one cell
for(i in names(table(cellInfo[[opt$groupby]])) ){
  if( table(cellInfo[[opt$groupby]])[i]==1 ){
    groupby_1 = c(groupby_1,i)
  }
}
if(is.null(groupby_1)){
  regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo[[opt$groupby]]),function(cells){
    rowmean = apply( regulonAUC_mat_out[,cells], 1, function(x){
      mean(x[!RemoveOutlier( x, nmads = 3, type = "higher")])
    })
  } )
}else{
  cellInfo_rm1 = cellInfo[-which(cellInfo[[opt$groupby]] %in% groupby_1),]
  cellInfo_rm1[[opt$groupby]] = droplevels(cellInfo_rm1[[opt$groupby]])
  regulonActivity_byCellType0 <- sapply(split(rownames(cellInfo_rm1), cellInfo_rm1[[opt$groupby]]),function(cells){
    rowmean = apply( regulonAUC_mat_out[,cells], 1, function(x){
      mean(x[!RemoveOutlier( x, nmads = 3, type = "higher")])
    })
  } )
  regulonAUC_mat_out_1 = regulonAUC_mat_out[,rownames(cellInfo[which(cellInfo[[opt$groupby]] %in% groupby_1),] )]
  regulonActivity_byCellType <- cbind(regulonActivity_byCellType0, regulonAUC_mat_out_1)
  colnames(regulonActivity_byCellType) <- c(colnames(regulonActivity_byCellType0), groupby_1)
}


# caculate the mean activity of each regulon in each group and plot heatmap
regulonActivity_byCellType_processed <- t(scale(t(regulonActivity_byCellType), center = T, scale=ifelse(length(unique(cellInfo[[opt$groupby]])) ==  2,F,T)))
regulonActivity_byCellType_processed <- na.omit(regulonActivity_byCellType_processed)
#regulonActivity_byCellType_processed <- regulonActivity_byCellType_processed[which(rowSums(regulonActivity_byCellType_processed) != 0),]
df = as.data.frame(regulonActivity_byCellType_processed) %>% tibble::rownames_to_column(var = "regulon")
write.table(df,file.path(output_dir, "1.2.centered_regulon_activity_groupby_design.xls"),
            sep = "\t", col.names =T, row.names =F, quote =F)

# regulonActivity_byCellType_processed <- regulonActivity_byCellType_processed[which(rowSums(abs(regulonActivity_byCellType_processed)) != 0),]
regulonActivity_byCellType_processed <- regulonActivity_byCellType_processed[which(rowSums(abs(regulonActivity_byCellType_processed)) > max(abs(regulonActivity_byCellType_processed))/4),]
if (dim(regulonActivity_byCellType_processed)[1]<11) {
  hig <-5
} else {
  hig <- dim(regulonActivity_byCellType_processed)[1]*0.45
}
regulonActivity_byCellType_processed = regulonActivity_byCellType_processed[,group_levels]
plot = pheatmap::pheatmap( regulonActivity_byCellType_processed,
                           # scale = "row",
                           cellwidth = 18,
                           cellheight = 18,
                           color=colorRampPalette(continuous_palette[[opt$ccolors]])(299),,
                           # annotation_col = factor(cellInfo[[opt$groupby]]),
                           angle_col = 45,
                           treeheight_col=20, 
                           treeheight_row=20,
                           border_color=NA,
                           fontsize = 15,  
                           fontsize_row = 13,
                           fontsize_col = 15,
                           cluster_rows=as.logical(opt$rowcluster),
                           cluster_cols=as.logical(opt$colcluster),
                           fontface = "bold")
save_ggplots(file.path(output_dir, paste0("1.3.regulon_activity_heatmap", collapse = "")),
             plot = plot,
             width = dim(regulonActivity_byCellType_processed)[2]/2+5,
             height = hig,
             dpi = 1000,
             limitsize = F,bg="white"
)

#
# RSS
#
seurat_ob = BinaryCount(seurat_ob, method = opt$binmethod, auc = regulonAUC, nCores = opt$ncores)
#if (!is.null(opt[["auc"]])) saveRDSMC(seurat_ob, opt$input)
rss_df = RunRSS(seurat_ob, group.by = opt$groupby)
rss_df_out = rss_df %>% subset(!grepl("extended",regulon))
write.table(rss_df_out, file.path(output_dir, "2.1.regulon_RSS_annotation.xls"),
            sep = "\t", col.names = T, row.names =F, quote =F)
rss_df[[opt$groupby]] = factor(rss_df[[opt$groupby]], levels = group_levels)
gg_rss_rank = RSSRanking(rss_df, group.by = opt$groupby, top_genes = topGenes, plot_extended = non_extended )
save_ggplots(file.path(output_dir, paste0("2.2.RSS_ranking_plot", collapse = "")),
             plot = gg_rss_rank,
             width = 7,
             height = ceiling(length(unique(rss_df[,opt$groupby]))/2) * 4,
             dpi = 1000,
             limitsize = F,bg="white"
)

rss_df_wide <- rss_df_out %>% tidyr::spread_( opt$groupby,"RSS")
rownames(rss_df_wide) <- rss_df_wide$regulon
rss_df_wide <- rss_df_wide[,2:ncol(rss_df_wide)]
## Subset all regulons that don't have at least an RSS of 0.4 for one cell type
if ( is.null(opt$threshold) ){
  rss_threshold = 0
}else{
  rss_threshold = as.numeric(opt$threshold)
}

rss_df_wide_specific <- rss_df_wide[apply(rss_df_wide,MARGIN = 1 ,FUN = function(x) any(x > rss_threshold)),]
rss_df_wide_specific = rss_df_wide_specific[,group_levels]
plot = pheatmap::pheatmap(rss_df_wide_specific,
                          cellwidth = 18,
                          cellheight = 18,
                          color=colorRampPalette(continuous_palette[[opt$ccolors]])(299),
                          angle_col = 45,
                          treeheight_col=20,
                          treeheight_row=20, 
                          border_color=NA,
                          fontsize = 15,  
                          fontsize_row = 15,
                          fontsize_col = 15,
                          cluster_rows=as.logical(opt$rowcluster),
                          cluster_cols=as.logical(opt$colcluster),
                          fontface = "bold")

save_ggplots(file.path(output_dir, paste0("2.3.RSS_heatmap", collapse = "")),
             plot = plot,
             width = dim(rss_df_wide_specific)[2]/2+5,
             height = dim(rss_df_wide_specific)[1]*0.6,
             dpi = 1000,
             limitsize = F,bg="white"
)
setwd(output_dir)
print("Convert pdf to png...")
system("module load ImageMagick && ls *.pdf | xargs -I {} -P 8 bash -c '/data/software/ImageMagick/ImageMagick-v7.0.8-14/bin/convert  -verbose -density 200 -trim  ${0}  -quality 100  -flatten  ${0%.pdf}.png' {}")









###18.MEBOCOST#######
###########(1)Overall##############
#########The following code runs on a Linux server
# MEBOCOST predicts cell-cell communication events where metabolites (e.g., lipids) are secreted by sender cells and sensed by sensor proteins on receiver cells.
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
# Initialize shell (choose according to your shell type)
source ~/anaconda3/bin/activate  # or ~/miniconda3/bin/activate
conda init bash                  # Permanent effect: auto-configure ~/.bashrc

#conda env list
#conda env remove --name mebocost

conda create -n mebocost python=3.12
conda activate mebocost          # Enter the environment
# Create a new directory, download mebocost package, then enter the downloaded folder
mkdir metabocost_analysis
git clone https://github.com/kaifuchenlab/MEBOCOST.git
cd MEBOCOST
# Install dependencies
pip install -r requirements.txt
# Install mebocost
python -m pip install .
# mebocost v1.0.4
mkdir /thinker/3.tangjiale/shinanxi/metabocost_analysis
cd  /thinker/3.tangjiale/shinanxi/metabocost_analysis
git clone https://github.com/kaifuchenlab/MEBOCOST.git
cd /thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/
  # Install dependencies
  pip install -r requirements.txt
# Install mebocost
python -m pip install .

# mebocost v1.0.4

python3.12

##########In Python
import os,sys
import scanpy as sc
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

from mebocost import mebocost



######In R, save expression matrix and cell annotation
library(Seurat)
library(dplyr)
# Expression matrix
exp <-  GetAssayData(seurat_obj, slot = 'data') %>% as.data.frame() 
exp <- exp[rowSums(exp) > 0,]  # Remove genes with zero expression across all cells
# Can save as csv or txt (both saving and Python reading are too slow, not recommended)
#write.csv(exp,file = "exp.txt",quote = F,sep = "\t",row.names = T,col.names = T)
# Save as feather file
library(arrow)
exp$row_id <- rownames(exp) 
rownames(exp) <- NULL           
write_feather(exp, "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Overall/exp.feather")
# Cell annotation
cell_anno <- seurat_obj@meta.data %>% select("celltype")
write.csv(cell_anno, file = "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Overall/cell_anno.csv",row.names = T,col.names = T)

####In Python
exp_mat = pd.read_feather("/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Overall/exp.feather")
## If error, maybe pyarrow is not installed. Use pip in Linux (recommended)
#pip install pyarrow
##### If pandas and pyarrow extensions conflict, restart the Python3 interface in the mebocost environment.
exp_mat.set_index("row_id", inplace=True)  # Restore row names as index
exp_mat
# Cell annotation, row names are cells
cell_ann = pd.read_csv('/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Overall/cell_anno.csv',index_col=0)
cell_ann
mebo_obj = mebocost.create_obj(
  adata = None,
  ## "Celltype" must be a column name of adata.obs table, otherwise, change it according to your data.
  group_col = 'celltype',
  condition_col = None,
  met_est = 'mebocost',
  # make sure mebocost.conf file in the same folder of this notebook, otherwise, provide an absolute path.
  config_path = './mebocost.conf', 
  exp_mat=exp_mat,
  cell_ann=cell_ann,
  ## make sure you set the right species
  species='human',
  met_pred=None,
  met_enzyme=None,
  met_sensor=None,
  met_ann=None,
  scFEA_ann=None,
  compass_met_ann=None,
  compass_rxn_ann=None,
  cutoff_exp='auto', ## automated cutoff to exclude lowly ranked 25% sensors across all cells
  cutoff_met='auto', ## automated cutoff to exclude lowly ranked 25% metabolites across all cells
  cutoff_prop=0.15, ## at least 15% of cells should be expressed the sensor or present the metabolite in the cell group (specified by group_col)
  sensor_type='All',
  thread=8
)
# This construction also works, but takes a long time
mebo_obj._load_config_()
mebo_obj.estimator()

## check the aggregated enzyme expression for metabolites
met_mat = pd.DataFrame(mebo_obj.met_mat.toarray(),
                       index = mebo_obj.met_mat_indexer,
                       columns = mebo_obj.met_mat_columns)
## print head
met_mat.head()

# This step takes some time
commu_res = mebo_obj.infer_commu(
  n_shuffle=1000,
  seed=12345, 
  Return=True, 
  thread=None,
  save_permuation=False,
  min_cell_number = 10,
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05
)

print('Number of mCCC detected by enzyme and sensor co-expression: ', commu_res.shape[0])
# You can see the number of communications
commu_res
# Result can be saved as pk file
# ### save 
mebocost.save_obj(obj = mebo_obj, path = '/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Overall/demo_HNSC_200cell_commu.pk')
## re-load the previous object if needed
mebo_obj = mebocost.load_obj('/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Overall/demo_HNSC_200cell_commu.pk')

print('sensor_exp cutoff: %s'%mebo_obj.cutoff_exp)
print('metabolite_agg_enzyme cutoff: %s'%mebo_obj.cutoff_met)

#### The above steps have already completed the mebocost analysis. The result is a data frame of cell-cell metabolic communications. The author's tutorial also includes an integration with COMPASS, using flux constraints to filter mCCCs. Here, COMPASS is run on the average gene expression of each cell type.
### Running COMPASS for each cell type by the average gene expression
### output average gene expression
#avg_exp = sc.get.aggregate(adata, by = 'celltype', func='mean')
#avg_exp = pd.DataFrame(avg_exp.layers['mean'], index = avg_exp.obs_names, columns = avg_exp.var_names).T
# ## do un log since COMPASS will take log in the algorithm
#avg_exp = avg_exp.apply(lambda col: np.exp(col)-1)
#avg_exp

#avg_exp.to_csv('avg_exp_mat.tsv', sep = '\t')  # Save file for COMPASS analysis
# Because we take the average, the number of columns is small, equivalent to only 10 cells, so it runs quickly!


####In R:
library(Seurat)
library(tidyverse)

# Read RDS file (assuming Seurat object)
seurat_obj <- readRDS('/thinker/3.tangjiale/shinanxi/rds/Overall_shinanxi_New.rds')

# 1. Compute average expression per cell type (log space)
avg_log <- AverageExpression(seurat_obj, 
                             group.by = "celltype",
                             assays = "RNA",        # Specify the assay to use
                             slot = "data")$RNA     # Default uses data slot (log1p transformed values)

# 2. Inverse log transformation: expm1(x) = exp(x)-1
avg_linear <- expm1(avg_log)

# 3. Save as TSV (rows: genes, columns: cell types)
write.table(avg_linear, 
            file = "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Overall/avg_exp_mat.tsv",
            sep = "\t",
            col.names = NA,
            quote = FALSE)

# Run COMPASS pipeline in another terminal window:
conda activate compass
## The number of processes (10) should ideally be greater than the number of samples (9)
cd /thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Overall/
  compass --data avg_exp_mat.tsv --num-processes 10 --species homo_sapiens --output-dir ./compass_res --temp-dir ./compass_res_tmp --calc-metabolites --lambda 0
### The core issue of AttributeError: 'DataFrame' object has no attribute 'iteritems' is that the pandas DataFrame method .iteritems() is deprecated and removed in newer pandas versions.
# Alternatively use conda (if conda can resolve dependency conflicts)
# Execute under the activated compass environment: conda install pandas=1.4.0 (if conda can resolve dependency conflicts)



###Execute in Python###
#python3
updated_res = mebo_obj._ConstrainCompassFlux_(compass_folder='/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Overall/compass_res/', efflux_cut='auto', influx_cut='auto', inplace=False)
## update to the object
mebo_obj.commu_res = updated_res
print('Number of mCCC detected by further flux constraints', updated_res.shape[0])
# You can see that after COMPASS constraint, the original 330 communications are reduced to 259

##########Visualize results
### ## filter by FDR less than 0.05, save significant results
commu_res = mebo_obj.commu_res.copy()
commu_res = commu_res[commu_res['permutation_test_fdr']<=0.05]
## write to tsv file
commu_res.to_csv('communication_result.tsv', sep = '\t', index = None)
#commu_res = pd.read_csv('communication_result.tsv', sep='\t')

###Bar plot##
## sender and receiver event number
p1 = mebo_obj.eventnum_bar(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  xorder=[],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  comm_score_col='Commu_Score',
  comm_score_cutoff = 0,
  cutoff_prop = 0.15,
  figsize=(6,4.5),
  save=None,
  show_plot=True,
  show_num = True,
  include=['sender-receiver'],
  group_by_cell=True,
  colorcmap='tab20',
  return_fig= True
)
# Custom save (PNG+PDF)
p1.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Overall/plot/barplot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p1.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Overall/plot/barplot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p1)  # Free memory

## circle plot to show communications between cell groups
custom_colors = {
  "Transitional activated":"#FF34B3",
  "Cytotoxic effector":"#00F5FF",
  "CSMD1+":"#BC8F8F",
  "NK-like":"#ADFF2F",
  "Activation-regulated":"#FFFF02",
  "Naive-like":"#00CD00",
  "Cytotoxic memory-like":"#FF6A6A",
  "Terminal-branch" : "#7FFFD4"
}
p2 = mebo_obj.commu_network_plot(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  node_cmap= custom_colors,
  figsize=(6,3),
  line_cmap='bwr',
  line_color_vmin=None,
  line_color_vmax=None,
  linewidth_norm=(0.2, 1),
  linewidth_value_range = None,
  node_size_norm=(50, 200),
  node_value_range = None,
  adjust_text_pos_node=True,
  node_text_hidden = False,
  node_text_font=4,  # Set font size
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=True,
  return_fig=True
)

# Custom save (PNG+PDF)
p2.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Overall/plot/network_plot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p2.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Overall/plot/network_plot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p2)  # Free memory

p3 = mebo_obj.count_dot_plot(
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  cmap='bwr',
  figsize=(10,6),
  save=None,
  dot_size_norm =(20, 200),
  dot_value_range = None,
  dot_color_vmin=None,
  dot_color_vmax=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  dendrogram_cluster=True,
  sender_order=[],
  receiver_order=[],
  return_fig = True
)

# Custom save (PNG+PDF)
p3.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Overall/plot/count_dot_plot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p3.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Overall/plot/count_dot_plot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p3)  # Free memory

## Focus on all cell types, use receiver_focus=[] to include all cell types
p4 = mebo_obj.commu_dotmap(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  figsize='auto',
  cmap='bwr',
  cmap_vmin = None,
  cmap_vmax = None,
  cellpair_order=[],
  met_sensor_order=[],
  dot_size_norm=(10, 150),
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_range = None,
  comm_score_cutoff=0,
  swap_axis = False,
  return_fig = True
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Overall/plot/commu_dotmap.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Overall/plot/commu_dotmap.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p4)  # Free memory

p4 = mebo_obj.commu_dotmap(
  sender_focus=["Transitional activated"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  figsize='auto',
  cmap='bwr',
  cmap_vmin = None,
  cmap_vmax = None,
  cellpair_order=[],
  met_sensor_order=[],
  dot_size_norm=(10, 150),
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_range = None,
  comm_score_cutoff=0,
  swap_axis = False,
  return_fig = True
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Overall/plot/sender_focus_GZMK+IL7R+_commu_dotmap.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Overall/plot/sender_focus_GZMK+IL7R+_commu_dotmap.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p4)

## Focus on all cell types, use receiver_focus=[] to include all cell types
p5 = mebo_obj.FlowPlot(
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  ## set to focus on mCCC of interest
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  remove_unrelevant = False,
  and_or='and',
  node_label_size=4,
  node_alpha=0.6,
  figsize='auto',
  node_cmap='Set1',
  line_cmap='bwr',
  line_cmap_vmin = None,
  line_cmap_vmax = None,
  node_size_norm=(20, 150),
  node_value_range = None,
  linewidth_norm=(0.5, 5),
  linewidth_value_range = None,
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=False,
  return_fig = True
)
p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Overall/plot/FlowPlot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Overall/plot/FlowPlot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p5)  

## Focus on all cell types, use receiver_focus=[] to include all cell types
p5 = mebo_obj.FlowPlot(
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  ## set to focus on mCCC of interest
  sender_focus=["Transitional activated"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  remove_unrelevant = False,
  and_or='and',
  node_label_size=4,
  node_alpha=0.6,
  figsize='auto',
  node_cmap='Set1',
  line_cmap='bwr',
  line_cmap_vmin = None,
  line_cmap_vmax = None,
  node_size_norm=(20, 150),
  node_value_range = None,
  linewidth_norm=(0.5, 5),
  linewidth_value_range = None,
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=False,
  return_fig = True
)
p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Overall/plot/sender_focus_GZMK+IL7R+_FlowPlot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Overall/plot/sender_focus_GZMK+IL7R+_FlowPlot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p5)  

commu_df = mebo_obj.commu_res.copy()
good_met = commu_df[(commu_df['permutation_test_fdr']<=0.05)]['Metabolite_Name'].sort_values().unique()
p6 = mebo_obj.violin_plot(
  sensor_or_met=good_met[:10], ## Show top10 metabolite expression only, or choose as desired
  cell_focus=[],
  cell_order = [],
  row_zscore = False,
  cmap=None,
  vmin=None,
  vmax=None,
  figsize='auto',
  cbar_title='',
  save=None,
  show_plot=True,
  return_fig = True
)
p6.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Overall/plot/Metabolite_Name_violin_plot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p6.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Overall/plot/Metabolite_Name_violin_plot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p6) 

good_sensor = commu_df[(commu_df['permutation_test_fdr']<=0.05)]['Sensor'].sort_values().unique()

p7 = mebo_obj.violin_plot(
  sensor_or_met=good_sensor[:10],  ## only top 10 as example, or choose as desired
  cell_focus=[],
  cell_order = [],
  row_zscore = False,
  cmap=None,
  vmin=None,
  vmax=None,
  figsize='auto',
  cbar_title='',
  save=None,
  show_plot=True,
  return_fig = True
)

p7.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Overall/plot/sensor_violin_plot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p7.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Overall/plot/sensor_violin_plot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p7) 


################(2)Control#################
seurat_obj <- readRDS("/thinker/3.tangjiale/shinanxi/monocle2/subset_ob.rds")
Idents(seurat_obj) <- 'group'
levels(Idents(seurat_obj))
## Check the number of cells per cell type
table(seurat_obj@meta.data$group)

scRNA_Control <- subset(seurat_obj, idents = 'Control')
seurat_obj = scRNA_Control

######In R, save expression matrix and cell annotation
library(Seurat)
library(dplyr)
# Expression matrix
exp <-  GetAssayData(seurat_obj, slot = 'data') %>% as.data.frame() 
exp <- exp[rowSums(exp) > 0,]  # Remove genes with all zero expression
# Can save as csv or txt (both saving and Python reading are too slow, not recommended)
#write.csv(exp,file = "exp.txt",quote = F,sep = "\t",row.names = T,col.names = T)
# Save as feather file
library(arrow)
exp$row_id <- rownames(exp) 
rownames(exp) <- NULL           
write_feather(exp, "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/exp.feather")
# Cell annotation
cell_anno <- seurat_obj@meta.data %>% select("celltype")
write.csv(cell_anno, file = "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/cell_anno.csv",row.names = T,col.names = T)


####In Python
exp_mat = pd.read_feather("/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/exp.feather")
## If error, maybe pyarrow is not installed. Use pip in Linux (recommended)
#pip install pyarrow
##### If pandas and pyarrow extensions conflict, restart the Python3 interface in the mebocost environment.
exp_mat.set_index("row_id", inplace=True)  # Restore row names as index
exp_mat
# Cell annotation, row names are cells
cell_ann = pd.read_csv('/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/cell_anno.csv',index_col=0)
cell_ann
mebo_obj = mebocost.create_obj(
  adata = None,
  ## "Celltype" must be a column name of adata.obs table, otherwise, change it according to your data.
  group_col = 'celltype',
  condition_col = None,
  met_est = 'mebocost',
  # make sure mebocost.conf file in the same folder of this notebook, otherwise, provide an absolute path.
  config_path = './mebocost.conf', 
  exp_mat=exp_mat,
  cell_ann=cell_ann,
  ## make sure you set the right species
  species='human',
  met_pred=None,
  met_enzyme=None,
  met_sensor=None,
  met_ann=None,
  scFEA_ann=None,
  compass_met_ann=None,
  compass_rxn_ann=None,
  cutoff_exp='auto', ## automated cutoff to exclude lowly ranked 25% sensors across all cells
  cutoff_met='auto', ## automated cutoff to exclude lowly ranked 25% metabolites across all cells
  cutoff_prop=0.15, ## at least 15% of cells should be expressed the sensor or present the metabolite in the cell group (specified by group_col)
  sensor_type='All',
  thread=8
)
# This construction also works, but takes a long time
mebo_obj._load_config_()
mebo_obj.estimator()

## check the aggregated enzyme expression for metabolites
met_mat = pd.DataFrame(mebo_obj.met_mat.toarray(),
                       index = mebo_obj.met_mat_indexer,
                       columns = mebo_obj.met_mat_columns)
## print head
met_mat.head()

# This step takes some time
commu_res = mebo_obj.infer_commu(
  n_shuffle=1000,
  seed=12345, 
  Return=True, 
  thread=None,
  save_permuation=False,
  min_cell_number = 10,
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05
)

print('Number of mCCC detected by enzyme and sensor co-expression: ', commu_res.shape[0])
# You can see the number of communications
commu_res
# Result can be saved as pk file
# ### save 
mebocost.save_obj(obj = mebo_obj, path = '/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/demo_HNSC_200cell_commu.pk')
## re-load the previous object if needed
mebo_obj = mebocost.load_obj('/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/demo_HNSC_200cell_commu.pk')

print('sensor_exp cutoff: %s'%mebo_obj.cutoff_exp)
print('metabolite_agg_enzyme cutoff: %s'%mebo_obj.cutoff_met)

#### The above steps have already completed the mebocost analysis. The result is a data frame of cell-cell metabolic communications. The author's tutorial also includes an integration with COMPASS, using flux constraints to filter mCCCs. Here, COMPASS is run on the average gene expression of each cell type.
### Running COMPASS for each cell type by the average gene expression
### output average gene expression
#avg_exp = sc.get.aggregate(adata, by = 'celltype', func='mean')
#avg_exp = pd.DataFrame(avg_exp.layers['mean'], index = avg_exp.obs_names, columns = avg_exp.var_names).T
# ## do un log since COMPASS will take log in the algorithm
#avg_exp = avg_exp.apply(lambda col: np.exp(col)-1)
#avg_exp

#avg_exp.to_csv('avg_exp_mat.tsv', sep = '\t')  # Save file for COMPASS analysis
# Because we take the average, the number of columns is small, equivalent to only 10 cells, so it runs quickly!



####In R:
library(Seurat)
library(tidyverse)

# Read RDS file (assuming Seurat object)
unique(seurat_obj$group)

# 1. Compute average expression per cell type (log space)
avg_log <- AverageExpression(seurat_obj, 
                             group.by = "celltype",
                             assays = "RNA",        # Specify the assay to use
                             slot = "data")$RNA     # Default uses data slot (log1p transformed values)

# 2. Inverse log transformation: expm1(x) = exp(x)-1
avg_linear <- expm1(avg_log)

# 3. Save as TSV (rows: genes, columns: cell types)
write.table(avg_linear, 
            file = "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/avg_exp_mat.tsv",
            sep = "\t",
            col.names = NA,
            quote = FALSE)



# Run COMPASS pipeline in another terminal window:
conda activate compass
## The number of processes (10) should ideally be greater than the number of samples (9)
cd /thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/
  compass --data avg_exp_mat.tsv --num-processes 10 --species homo_sapiens --output-dir ./compass_res --temp-dir ./compass_res_tmp --calc-metabolites --lambda 0
### The core issue of AttributeError: 'DataFrame' object has no attribute 'iteritems' is that the pandas DataFrame method .iteritems() is deprecated and removed in newer pandas versions.
# Alternatively use conda (if conda can resolve dependency conflicts)
# Execute under the activated compass environment: conda install pandas=1.4.0 (if conda can resolve dependency conflicts)


###Execute in Python###
#python3
updated_res = mebo_obj._ConstrainCompassFlux_(compass_folder='/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/compass_res/', efflux_cut='auto', influx_cut='auto', inplace=False)
## update to the object
mebo_obj.commu_res = updated_res
print('Number of mCCC detected by further flux constraints', updated_res.shape[0])
# You can see that after COMPASS constraint, the original 330 communications are reduced to 259

##########Visualize results
### ## filter by FDR less than 0.05, save significant results
commu_res = mebo_obj.commu_res.copy()
commu_res = commu_res[commu_res['permutation_test_fdr']<=0.05]
## write to tsv file
commu_res.to_csv('./Control/communication_result.tsv', sep = '\t', index = None)
#commu_res = pd.read_csv('communication_result.tsv', sep='\t')

###Bar plot##
## sender and receiver event number
p1 = mebo_obj.eventnum_bar(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  xorder=[],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  comm_score_col='Commu_Score',
  comm_score_cutoff = 0,
  cutoff_prop = 0.15,
  figsize=(6,4.5),
  save=None,
  show_plot=True,
  show_num = True,
  include=['sender-receiver'],
  group_by_cell=True,
  colorcmap='tab20',
  return_fig= True
)
# Custom save (PNG+PDF)
p1.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/barplot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p1.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/barplot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p1)  # Free memory

## circle plot to show communications between cell groups
custom_colors = {
  "Transitional activated":"#FF34B3",
  "Cytotoxic effector":"#00F5FF",
  "CSMD1+":"#BC8F8F",
  "NK-like":"#ADFF2F",
  "Activation-regulated":"#FFFF02",
  "Naive-like":"#00CD00",
  "Cytotoxic memory-like":"#FF6A6A",
  "Terminal-branch" : "#7FFFD4"
}
p2 = mebo_obj.commu_network_plot(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  node_cmap= custom_colors,
  figsize=(6,3),
  line_cmap='bwr',
  line_color_vmin=None,
  line_color_vmax=None,
  linewidth_norm=(0.2, 1),
  linewidth_value_range = None,
  node_size_norm=(50, 200),
  node_value_range = None,
  adjust_text_pos_node=True,
  node_text_hidden = False,
  node_text_font=4,  # Set font size
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=True,
  return_fig=True
)

# Custom save (PNG+PDF)
p2.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/network_plot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p2.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/network_plot.png", 
  dpi=1600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p2)  # Free memory

p3 = mebo_obj.count_dot_plot(
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  cmap='bwr',
  figsize=(10,6),
  save=None,
  dot_size_norm =(20, 200),
  dot_value_range = None,
  dot_color_vmin=None,
  dot_color_vmax=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  dendrogram_cluster=True,
  sender_order=[],
  receiver_order=[],
  return_fig = True
)

# Custom save (PNG+PDF)
p3.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/count_dot_plot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p3.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/count_dot_plot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p3)  # Free memory

## Focus on all cell types, use receiver_focus=[] to include all cell types
p4 = mebo_obj.commu_dotmap(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  figsize='auto',
  cmap='bwr',
  cmap_vmin = None,
  cmap_vmax = None,
  cellpair_order=[],
  met_sensor_order=[],
  dot_size_norm=(10, 150),
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_range = None,
  comm_score_cutoff=0,
  swap_axis = False,
  return_fig = True
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/commu_dotmap.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/commu_dotmap.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p4)  # Free memory

p4 = mebo_obj.commu_dotmap(
  sender_focus=["Transitional activated"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  figsize='auto',
  cmap='bwr',
  cmap_vmin = None,
  cmap_vmax = None,
  cellpair_order=[],
  met_sensor_order=[],
  dot_size_norm=(10, 150),
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_range = None,
  comm_score_cutoff=0,
  swap_axis = False,
  return_fig = True
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/sender_focus_GZMK+IL7R+_commu_dotmap.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/sender_focus_GZMK+IL7R+_commu_dotmap.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p4)

p4 = mebo_obj.commu_dotmap(
  sender_focus=["Cytotoxic effector"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  figsize='auto',
  cmap='bwr',
  cmap_vmin = None,
  cmap_vmax = None,
  cellpair_order=[],
  met_sensor_order=[],
  dot_size_norm=(10, 150),
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_range = None,
  comm_score_cutoff=0,
  swap_axis = False,
  return_fig = True
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/sender_focus_ZEB2+FGFBP2+_commu_dotmap.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/sender_focus_ZEB2+FGFBP2+_commu_dotmap.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p4)

## Focus on all cell types, use receiver_focus=[] to include all cell types
p5 = mebo_obj.FlowPlot(
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  ## set to focus on mCCC of interest
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  remove_unrelevant = False,
  and_or='and',
  node_label_size=4,
  node_alpha=0.6,
  figsize='auto',
  node_cmap='Set1',
  line_cmap='bwr',
  line_cmap_vmin = None,
  line_cmap_vmax = None,
  node_size_norm=(20, 150),
  node_value_range = None,
  linewidth_norm=(0.5, 5),
  linewidth_value_range = None,
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=False,
  return_fig = True
)
p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/FlowPlot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/FlowPlot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p5)  

## Focus on all cell types, use receiver_focus=[] to include all cell types
p5 = mebo_obj.FlowPlot(
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  ## set to focus on mCCC of interest
  sender_focus=["Transitional activated"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  remove_unrelevant = False,
  and_or='and',
  node_label_size=4,
  node_alpha=0.6,
  figsize='auto',
  node_cmap='Set1',
  line_cmap='bwr',
  line_cmap_vmin = None,
  line_cmap_vmax = None,
  node_size_norm=(20, 150),
  node_value_range = None,
  linewidth_norm=(0.5, 5),
  linewidth_value_range = None,
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=False,
  return_fig = True
)
p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/sender_focus_GZMK+IL7R+_FlowPlot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/sender_focus_GZMK+IL7R+_FlowPlot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p5)  

p5 = mebo_obj.FlowPlot(
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  ## set to focus on mCCC of interest
  sender_focus=["Cytotoxic effector"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  remove_unrelevant = False,
  and_or='and',
  node_label_size=4,
  node_alpha=0.6,
  figsize='auto',
  node_cmap='Set1',
  line_cmap='bwr',
  line_cmap_vmin = None,
  line_cmap_vmax = None,
  node_size_norm=(20, 150),
  node_value_range = None,
  linewidth_norm=(0.5, 5),
  linewidth_value_range = None,
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=False,
  return_fig = True
)
p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/sender_focus_ZEB2+FGFBP2+_FlowPlot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/sender_focus_ZEB2+FGFBP2+_FlowPlot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p5)  

commu_df = mebo_obj.commu_res.copy()
good_met = commu_df[(commu_df['permutation_test_fdr']<=0.05)]['Metabolite_Name'].sort_values().unique()
p6 = mebo_obj.violin_plot(
  sensor_or_met=good_met[:10], ## Show top10 metabolite expression only, or choose as desired
  cell_focus=[],
  cell_order = [],
  row_zscore = False,
  cmap=None,
  vmin=None,
  vmax=None,
  figsize='auto',
  cbar_title='',
  save=None,
  show_plot=True,
  return_fig = True
)
p6.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/Metabolite_Name_violin_plot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p6.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/Metabolite_Name_violin_plot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p6) 

good_sensor = commu_df[(commu_df['permutation_test_fdr']<=0.05)]['Sensor'].sort_values().unique()

p7 = mebo_obj.violin_plot(
  sensor_or_met=good_sensor[:10],  ## only top 10 as example, or choose as desired
  cell_focus=[],
  cell_order = [],
  row_zscore = False,
  cmap=None,
  vmin=None,
  vmax=None,
  figsize='auto',
  cbar_title='',
  save=None,
  show_plot=True,
  return_fig = True
)

p7.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/sensor_violin_plot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p7.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/sensor_violin_plot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p7) 

#######sender

p2 = mebo_obj.commu_network_plot(
  sender_focus=["CSMD1+"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  node_cmap= custom_colors,
  figsize=(6,3),
  line_cmap='bwr',
  line_color_vmin=None,
  line_color_vmax=None,
  linewidth_norm=(0.2, 1),
  linewidth_value_range = None,
  node_size_norm=(50, 200),
  node_value_range = None,
  adjust_text_pos_node=True,
  node_text_hidden = True,
  node_text_font=4,  # Set font size
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=True,
  return_fig=True
)

p2.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/sender_focus_GZMK+CSMD1+_network_plot.png", 
  dpi=1600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p2)  

p4 = mebo_obj.commu_dotmap(
  sender_focus=["CSMD1+"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  figsize='auto',
  cmap='bwr',
  cmap_vmin = None,
  cmap_vmax = None,
  cellpair_order=[],
  met_sensor_order=[],
  dot_size_norm=(10, 150),
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_range = None,
  comm_score_cutoff=0,
  swap_axis = False,
  return_fig = True
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/sender_focus_GZMK+CSMD1+_commu_dotmap.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/sender_focus_GZMK+CSMD1+_commu_dotmap.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p4)

p5 = mebo_obj.FlowPlot(
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  ## set to focus on mCCC of interest
  sender_focus=["CSMD1+"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  remove_unrelevant = False,
  and_or='and',
  node_label_size=4,
  node_alpha=0.6,
  figsize='auto',
  node_cmap='Set1',
  line_cmap='bwr',
  line_cmap_vmin = None,
  line_cmap_vmax = None,
  node_size_norm=(20, 150),
  node_value_range = None,
  linewidth_norm=(0.5, 5),
  linewidth_value_range = None,
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=False,
  return_fig = True
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/sender_focus_GZMK+CSMD1+_FlowPlot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/sender_focus_GZMK+CSMD1+_FlowPlot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p5)  

p2 = mebo_obj.commu_network_plot(
  sender_focus=["Transitional activated"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  node_cmap= custom_colors,
  figsize=(6,3),
  line_cmap='bwr',
  line_color_vmin=None,
  line_color_vmax=None,
  linewidth_norm=(0.2, 1),
  linewidth_value_range = None,
  node_size_norm=(50, 200),
  node_value_range = None,
  adjust_text_pos_node=True,
  node_text_hidden = True,
  node_text_font=4,  # Set font size
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=True,
  return_fig=True
)

p2.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/sender_focus_GZMK+IL7R+_network_plot.png", 
  dpi=1600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p2)  

p4 = mebo_obj.commu_dotmap(
  sender_focus=["Transitional activated"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  figsize='auto',
  cmap='bwr',
  cmap_vmin = None,
  cmap_vmax = None,
  cellpair_order=[],
  met_sensor_order=[],
  dot_size_norm=(10, 150),
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_range = None,
  comm_score_cutoff=0,
  swap_axis = False,
  return_fig = True
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/sender_focus_GZMK+IL7R+_commu_dotmap1.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/sender_focus_GZMK+IL7R+_commu_dotmap1.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p4)

p5 = mebo_obj.FlowPlot(
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  ## set to focus on mCCC of interest
  sender_focus=["Transitional activated"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  remove_unrelevant = False,
  and_or='and',
  node_label_size=4,
  node_alpha=0.6,
  figsize='auto',
  node_cmap='Set1',
  line_cmap='bwr',
  line_cmap_vmin = None,
  line_cmap_vmax = None,
  node_size_norm=(20, 150),
  node_value_range = None,
  linewidth_norm=(0.5, 5),
  linewidth_value_range = None,
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=False,
  return_fig = True
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/sender_focus_GZMK+IL7R+_FlowPlot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/sender_focus_GZMK+IL7R+_FlowPlot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p5)  

p2 = mebo_obj.commu_network_plot(
  sender_focus=["Cytotoxic effector"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  node_cmap= custom_colors,
  figsize=(6,3),
  line_cmap='bwr',
  line_color_vmin=None,
  line_color_vmax=None,
  linewidth_norm=(0.2, 1),
  linewidth_value_range = None,
  node_size_norm=(50, 200),
  node_value_range = None,
  adjust_text_pos_node=True,
  node_text_hidden = True,
  node_text_font=4,  # Set font size
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=True,
  return_fig=True
)

p2.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/sender_focus_ZEB2+FGFBP2+_network_plot.png", 
  dpi=1600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p2)  

p4 = mebo_obj.commu_dotmap(
  sender_focus=["Cytotoxic effector"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  figsize='auto',
  cmap='bwr',
  cmap_vmin = None,
  cmap_vmax = None,
  cellpair_order=[],
  met_sensor_order=[],
  dot_size_norm=(10, 150),
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_range = None,
  comm_score_cutoff=0,
  swap_axis = False,
  return_fig = True
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/sender_focus_ZEB2+FGFBP2+_commu_dotmap.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/sender_focus_ZEB2+FGFBP2+_commu_dotmap.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p4)

p5 = mebo_obj.FlowPlot(
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  ## set to focus on mCCC of interest
  sender_focus=["Cytotoxic effector"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  remove_unrelevant = False,
  and_or='and',
  node_label_size=4,
  node_alpha=0.6,
  figsize='auto',
  node_cmap='Set1',
  line_cmap='bwr',
  line_cmap_vmin = None,
  line_cmap_vmax = None,
  node_size_norm=(20, 150),
  node_value_range = None,
  linewidth_norm=(0.5, 5),
  linewidth_value_range = None,
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=False,
  return_fig = True
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/sender_focus_ZEB2+FGFBP2+_FlowPlot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/sender_focus_ZEB2+FGFBP2+_FlowPlot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p5)  

##############receiver
p2 = mebo_obj.commu_network_plot(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["NK-like"],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  node_cmap= custom_colors,
  figsize=(6,3),
  line_cmap='bwr',
  line_color_vmin=None,
  line_color_vmax=None,
  linewidth_norm=(0.2, 1),
  linewidth_value_range = None,
  node_size_norm=(50, 200),
  node_value_range = None,
  adjust_text_pos_node=True,
  node_text_hidden = True,
  node_text_font=4,  # Set font size
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=True,
  return_fig=True
)

p2.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/receiver_focus_ZEB2+KLRF1+_network_plot.png", 
  dpi=1600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p2)  

p4 = mebo_obj.commu_dotmap(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["NK-like"],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  figsize='auto',
  cmap='bwr',
  cmap_vmin = None,
  cmap_vmax = None,
  cellpair_order=[],
  met_sensor_order=[],
  dot_size_norm=(10, 150),
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_range = None,
  comm_score_cutoff=0,
  swap_axis = False,
  return_fig = True
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/receiver_focus_ZEB2+KLRF1+_commu_dotmap.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/receiver_focus_ZEB2+KLRF1+_commu_dotmap.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p4)

p5 = mebo_obj.FlowPlot(
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  ## set to focus on mCCC of interest
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["NK-like"],
  remove_unrelevant = False,
  and_or='and',
  node_label_size=4,
  node_alpha=0.6,
  figsize='auto',
  node_cmap='Set1',
  line_cmap='bwr',
  line_cmap_vmin = None,
  line_cmap_vmax = None,
  node_size_norm=(20, 150),
  node_value_range = None,
  linewidth_norm=(0.5, 5),
  linewidth_value_range = None,
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=False,
  return_fig = True
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/receiver_focus_ZEB2+KLRF1+_FlowPlot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/receiver_focus_ZEB2+KLRF1+_FlowPlot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p5)  

p2 = mebo_obj.commu_network_plot(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["Cytotoxic effector"],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  node_cmap= custom_colors,
  figsize=(6,3),
  line_cmap='bwr',
  line_color_vmin=None,
  line_color_vmax=None,
  linewidth_norm=(0.2, 1),
  linewidth_value_range = None,
  node_size_norm=(50, 200),
  node_value_range = None,
  adjust_text_pos_node=True,
  node_text_hidden = True,
  node_text_font=4,  # Set font size
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=True,
  return_fig=True
)

p2.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/receiver_focus_ZEB2+FGFBP2+_network_plot.png", 
  dpi=1600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p2)  

p4 = mebo_obj.commu_dotmap(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["Cytotoxic effector"],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  figsize='auto',
  cmap='bwr',
  cmap_vmin = None,
  cmap_vmax = None,
  cellpair_order=[],
  met_sensor_order=[],
  dot_size_norm=(10, 150),
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_range = None,
  comm_score_cutoff=0,
  swap_axis = False,
  return_fig = True
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/receiver_focus_ZEB2+FGFBP2+_commu_dotmap.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/receiver_focus_ZEB2+FGFBP2+_commu_dotmap.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p4)

p5 = mebo_obj.FlowPlot(
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  ## set to focus on mCCC of interest
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["Cytotoxic effector"],
  remove_unrelevant = False,
  and_or='and',
  node_label_size=4,
  node_alpha=0.6,
  figsize='auto',
  node_cmap='Set1',
  line_cmap='bwr',
  line_cmap_vmin = None,
  line_cmap_vmax = None,
  node_size_norm=(20, 150),
  node_value_range = None,
  linewidth_norm=(0.5, 5),
  linewidth_value_range = None,
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=False,
  return_fig = True
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/receiver_focus_ZEB2+FGFBP2+_FlowPlot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/receiver_focus_ZEB2+FGFBP2+_FlowPlot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p5)  

p2 = mebo_obj.commu_network_plot(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["Terminal-branch"],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  node_cmap= custom_colors,
  figsize=(6,3),
  line_cmap='bwr',
  line_color_vmin=None,
  line_color_vmax=None,
  linewidth_norm=(0.2, 1),
  linewidth_value_range = None,
  node_size_norm=(50, 200),
  node_value_range = None,
  adjust_text_pos_node=True,
  node_text_hidden = True,
  node_text_font=4,  # Set font size
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=True,
  return_fig=True
)

p2.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/receiver_focus_ZEB2+_network_plot.png", 
  dpi=1600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p2)  

p4 = mebo_obj.commu_dotmap(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["Terminal-branch"],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  figsize='auto',
  cmap='bwr',
  cmap_vmin = None,
  cmap_vmax = None,
  cellpair_order=[],
  met_sensor_order=[],
  dot_size_norm=(10, 150),
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_range = None,
  comm_score_cutoff=0,
  swap_axis = False,
  return_fig = True
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/receiver_focus_ZEB2+_commu_dotmap.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/receiver_focus_ZEB2+_commu_dotmap.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p4)

p5 = mebo_obj.FlowPlot(
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  ## set to focus on mCCC of interest
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["Terminal-branch"],
  remove_unrelevant = False,
  and_or='and',
  node_label_size=4,
  node_alpha=0.6,
  figsize='auto',
  node_cmap='Set1',
  line_cmap='bwr',
  line_cmap_vmin = None,
  line_cmap_vmax = None,
  node_size_norm=(20, 150),
  node_value_range = None,
  linewidth_norm=(0.5, 5),
  linewidth_value_range = None,
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=False,
  return_fig = True
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/receiver_focus_ZEB2+_FlowPlot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/receiver_focus_ZEB2+_FlowPlot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p5)  

p2 = mebo_obj.commu_network_plot(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["Naive-like"],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  node_cmap= custom_colors,
  figsize=(6,3),
  line_cmap='bwr',
  line_color_vmin=None,
  line_color_vmax=None,
  linewidth_norm=(0.2, 1),
  linewidth_value_range = None,
  node_size_norm=(50, 200),
  node_value_range = None,
  adjust_text_pos_node=True,
  node_text_hidden = True,
  node_text_font=4,  # Set font size
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=True,
  return_fig=True
)

p2.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/receiver_focus_LEF1+NELL2+_network_plot.png", 
  dpi=1600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p2)  

p4 = mebo_obj.commu_dotmap(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["Naive-like"],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  figsize='auto',
  cmap='bwr',
  cmap_vmin = None,
  cmap_vmax = None,
  cellpair_order=[],
  met_sensor_order=[],
  dot_size_norm=(10, 150),
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_range = None,
  comm_score_cutoff=0,
  swap_axis = False,
  return_fig = True
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/receiver_focus_LEF1+NELL2+_commu_dotmap.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/receiver_focus_LEF1+NELL2+_commu_dotmap.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p4)

p5 = mebo_obj.FlowPlot(
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  ## set to focus on mCCC of interest
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["Naive-like"],
  remove_unrelevant = False,
  and_or='and',
  node_label_size=4,
  node_alpha=0.6,
  figsize='auto',
  node_cmap='Set1',
  line_cmap='bwr',
  line_cmap_vmin = None,
  line_cmap_vmax = None,
  node_size_norm=(20, 150),
  node_value_range = None,
  linewidth_norm=(0.5, 5),
  linewidth_value_range = None,
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=False,
  return_fig = True
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/receiver_focus_LEF1+NELL2+_FlowPlot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/receiver_focus_LEF1+NELL2+_FlowPlot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p5)  

p2 = mebo_obj.commu_network_plot(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["Cytotoxic memory-like"],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  node_cmap= custom_colors,
  figsize=(6,3),
  line_cmap='bwr',
  line_color_vmin=None,
  line_color_vmax=None,
  linewidth_norm=(0.2, 1),
  linewidth_value_range = None,
  node_size_norm=(50, 200),
  node_value_range = None,
  adjust_text_pos_node=True,
  node_text_hidden = True,
  node_text_font=4,  # Set font size
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=True,
  return_fig=True
)

p2.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/receiver_focus_FGFBP2+ITGB1+_network_plot.png", 
  dpi=1600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p2)  

p4 = mebo_obj.commu_dotmap(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["Cytotoxic memory-like"],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  figsize='auto',
  cmap='bwr',
  cmap_vmin = None,
  cmap_vmax = None,
  cellpair_order=[],
  met_sensor_order=[],
  dot_size_norm=(10, 150),
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_range = None,
  comm_score_cutoff=0,
  swap_axis = False,
  return_fig = True
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/receiver_focus_FGFBP2+ITGB1+_commu_dotmap1.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/receiver_focus_FGFBP2+ITGB1+_commu_dotmap1.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p4)

p5 = mebo_obj.FlowPlot(
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  ## set to focus on mCCC of interest
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["Cytotoxic memory-like"],
  remove_unrelevant = False,
  and_or='and',
  node_label_size=4,
  node_alpha=0.6,
  figsize='auto',
  node_cmap='Set1',
  line_cmap='bwr',
  line_cmap_vmin = None,
  line_cmap_vmax = None,
  node_size_norm=(20, 150),
  node_value_range = None,
  linewidth_norm=(0.5, 5),
  linewidth_value_range = None,
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=False,
  return_fig = True
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/receiver_focus_FGFBP2+ITGB1+_FlowPlot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Control/plot/receiver_focus_FGFBP2+ITGB1+_FlowPlot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p5)  



################(3)Psoriasis#################
seurat_obj <- readRDS("/thinker/3.tangjiale/shinanxi/monocle2/subset_ob.rds")
Idents(seurat_obj) <- 'group'
levels(Idents(seurat_obj))
## Check the number of cells per cell type
table(seurat_obj@meta.data$group)

scRNA_Psoriasis <- subset(seurat_obj, idents = 'Psoriasis')
seurat_obj = scRNA_Psoriasis

######In R, save expression matrix and cell annotation
library(Seurat)
library(dplyr)
# Expression matrix
exp <-  GetAssayData(seurat_obj, slot = 'data') %>% as.data.frame() 
exp <- exp[rowSums(exp) > 0,]  # Remove genes with all zero expression
# Can save as csv or txt (both saving and Python reading are too slow, not recommended)
#write.csv(exp,file = "exp.txt",quote = F,sep = "\t",row.names = T,col.names = T)
# Save as feather file
library(arrow)
exp$row_id <- rownames(exp) 
rownames(exp) <- NULL           
write_feather(exp, "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/exp.feather")
# Cell annotation
cell_anno <- seurat_obj@meta.data %>% select("celltype")
write.csv(cell_anno, file = "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/cell_anno.csv",row.names = T,col.names = T)



####In Python
exp_mat = pd.read_feather("/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/exp.feather")
## If error, maybe pyarrow is not installed. Use pip in Linux (recommended)
#pip install pyarrow
##### If pandas and pyarrow extensions conflict, restart the Python3 interface in the mebocost environment.
exp_mat.set_index("row_id", inplace=True)  # Restore row names as index
exp_mat
# Cell annotation, row names are cells
cell_ann = pd.read_csv('/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/cell_anno.csv',index_col=0)
cell_ann
mebo_obj = mebocost.create_obj(
  adata = None,
  ## "Celltype" must be a column name of adata.obs table, otherwise, change it according to your data.
  group_col = 'celltype',
  condition_col = None,
  met_est = 'mebocost',
  # make sure mebocost.conf file in the same folder of this notebook, otherwise, provide an absolute path.
  config_path = './mebocost.conf', 
  exp_mat=exp_mat,
  cell_ann=cell_ann,
  ## make sure you set the right species
  species='human',
  met_pred=None,
  met_enzyme=None,
  met_sensor=None,
  met_ann=None,
  scFEA_ann=None,
  compass_met_ann=None,
  compass_rxn_ann=None,
  cutoff_exp='auto', ## automated cutoff to exclude lowly ranked 25% sensors across all cells
  cutoff_met='auto', ## automated cutoff to exclude lowly ranked 25% metabolites across all cells
  cutoff_prop=0.15, ## at least 15% of cells should be expressed the sensor or present the metabolite in the cell group (specified by group_col)
  sensor_type='All',
  thread=8
)
# This construction also works, but takes a long time
mebo_obj._load_config_()
mebo_obj.estimator()

## check the aggregated enzyme expression for metabolites
met_mat = pd.DataFrame(mebo_obj.met_mat.toarray(),
                       index = mebo_obj.met_mat_indexer,
                       columns = mebo_obj.met_mat_columns)
## print head
met_mat.head()

# This step takes some time
commu_res = mebo_obj.infer_commu(
  n_shuffle=1000,
  seed=12345, 
  Return=True, 
  thread=None,
  save_permuation=False,
  min_cell_number = 10,
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05
)

print('Number of mCCC detected by enzyme and sensor co-expression: ', commu_res.shape[0])
# You can see the number of communications
commu_res
# Result can be saved as pk file
# ### save 
mebocost.save_obj(obj = mebo_obj, path = '/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/demo_HNSC_200cell_commu.pk')
## re-load the previous object if needed
mebo_obj = mebocost.load_obj('/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/demo_HNSC_200cell_commu.pk')

print('sensor_exp cutoff: %s'%mebo_obj.cutoff_exp)
print('metabolite_agg_enzyme cutoff: %s'%mebo_obj.cutoff_met)

#### The above steps have already completed the mebocost analysis. The result is a data frame of cell-cell metabolic communications. The author's tutorial also includes an integration with COMPASS, using flux constraints to filter mCCCs. Here, COMPASS is run on the average gene expression of each cell type.
### Running COMPASS for each cell type by the average gene expression
### output average gene expression
#avg_exp = sc.get.aggregate(adata, by = 'celltype', func='mean')
#avg_exp = pd.DataFrame(avg_exp.layers['mean'], index = avg_exp.obs_names, columns = avg_exp.var_names).T
# ## do un log since COMPASS will take log in the algorithm
#avg_exp = avg_exp.apply(lambda col: np.exp(col)-1)
#avg_exp

#avg_exp.to_csv('avg_exp_mat.tsv', sep = '\t')  # Save file for COMPASS analysis
# Because we take the average, the number of columns is small, equivalent to only 10 cells, so it runs quickly!



####In R:
library(Seurat)
library(tidyverse)

# Read RDS file (assuming Seurat object)
unique(seurat_obj$group)

# 1. Compute average expression per cell type (log space)
avg_log <- AverageExpression(seurat_obj, 
                             group.by = "celltype",
                             assays = "RNA",        # Specify the assay to use
                             slot = "data")$RNA     # Default uses data slot (log1p transformed values)

# 2. Inverse log transformation: expm1(x) = exp(x)-1
avg_linear <- expm1(avg_log)

# 3. Save as TSV (rows: genes, columns: cell types)
write.table(avg_linear, 
            file = "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/avg_exp_mat.tsv",
            sep = "\t",
            col.names = NA,
            quote = FALSE)


# Run COMPASS pipeline in another terminal window:
conda activate compass
## The number of processes (10) should ideally be greater than the number of samples (9)
cd /thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/
  compass --data avg_exp_mat.tsv --num-processes 10 --species homo_sapiens --output-dir ./compass_res --temp-dir ./compass_res_tmp --calc-metabolites --lambda 0
### The core issue of AttributeError: 'DataFrame' object has no attribute 'iteritems' is that the pandas DataFrame method .iteritems() is deprecated and removed in newer pandas versions.
# Alternatively use conda (if conda can resolve dependency conflicts)
# Execute under the activated compass environment: conda install pandas=1.4.0 (if conda can resolve dependency conflicts)


###Execute in Python###
#python3
updated_res = mebo_obj._ConstrainCompassFlux_(compass_folder='/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/compass_res', efflux_cut='auto', influx_cut='auto', inplace=False)
## update to the object
mebo_obj.commu_res = updated_res
print('Number of mCCC detected by further flux constraints', updated_res.shape[0])
# You can see that after COMPASS constraint, the original 330 communications are reduced to 259

##########Visualize results
### ## filter by FDR less than 0.05, save significant results
commu_res = mebo_obj.commu_res.copy()
commu_res = commu_res[commu_res['permutation_test_fdr']<=0.05]
## write to tsv file
commu_res.to_csv('./Psoriasis/communication_result.tsv', sep = '\t', index = None)
#commu_res = pd.read_csv('communication_result.tsv', sep='\t')

###Bar plot##
## sender and receiver event number
p1 = mebo_obj.eventnum_bar(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  xorder=[],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  comm_score_col='Commu_Score',
  comm_score_cutoff = 0,
  cutoff_prop = 0.15,
  figsize=(6,4.5),
  save=None,
  show_plot=True,
  show_num = True,
  include=['sender-receiver'],
  group_by_cell=True,
  colorcmap='tab20',
  return_fig= True
)
# Custom save (PNG+PDF)
p1.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/barplot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p1.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/barplot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p1)  # Free memory

## circle plot to show communications between cell groups
custom_colors = {
  "Transitional activated":"#FF34B3",
  "Cytotoxic effector":"#00F5FF",
  "CSMD1+":"#BC8F8F",
  "NK-like":"#ADFF2F",
  "Activation-regulated":"#FFFF02",
  "Naive-like":"#00CD00",
  "Cytotoxic memory-like":"#FF6A6A",
  "Terminal-branch" : "#7FFFD4"
}
p2 = mebo_obj.commu_network_plot(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  node_cmap= custom_colors,
  figsize=(6,3),
  line_cmap='bwr',
  line_color_vmin=None,
  line_color_vmax=None,
  linewidth_norm=(0.2, 1),
  linewidth_value_range = None,
  node_size_norm=(50, 200),
  node_value_range = None,
  adjust_text_pos_node=True,
  node_text_hidden = False,
  node_text_font=4,  # Set font size
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=True,
  return_fig=True
)

# Custom save (PNG+PDF)
p2.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/network_plot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p2.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/network_plot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p2)  # Free memory

p3 = mebo_obj.count_dot_plot(
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  cmap='bwr',
  figsize=(10,6),
  save=None,
  dot_size_norm =(20, 200),
  dot_value_range = None,
  dot_color_vmin=None,
  dot_color_vmax=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  dendrogram_cluster=True,
  sender_order=[],
  receiver_order=[],
  return_fig = True
)

# Custom save (PNG+PDF)
p3.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/count_dot_plot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p3.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/count_dot_plot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p3)  # Free memory

## Focus on all cell types, use receiver_focus=[] to include all cell types
p4 = mebo_obj.commu_dotmap(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  figsize='auto',
  cmap='bwr',
  cmap_vmin = None,
  cmap_vmax = None,
  cellpair_order=[],
  met_sensor_order=[],
  dot_size_norm=(10, 150),
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_range = None,
  comm_score_cutoff=0,
  swap_axis = False,
  return_fig = True
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/commu_dotmap.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/commu_dotmap.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p4)  # Free memory

p4 = mebo_obj.commu_dotmap(
  sender_focus=["Transitional activated"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  figsize='auto',
  cmap='bwr',
  cmap_vmin = None,
  cmap_vmax = None,
  cellpair_order=[],
  met_sensor_order=[],
  dot_size_norm=(10, 150),
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_range = None,
  comm_score_cutoff=0,
  swap_axis = False,
  return_fig = True
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/sender_focus_GZMK+IL7R+_commu_dotmap.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/sender_focus_GZMK+IL7R+_commu_dotmap.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p4)

p4 = mebo_obj.commu_dotmap(
  sender_focus=["Cytotoxic effector"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  figsize='auto',
  cmap='bwr',
  cmap_vmin = None,
  cmap_vmax = None,
  cellpair_order=[],
  met_sensor_order=[],
  dot_size_norm=(10, 150),
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_range = None,
  comm_score_cutoff=0,
  swap_axis = False,
  return_fig = True
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/sender_focus_ZEB2+FGFBP2+_commu_dotmap.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/sender_focus_ZEB2+FGFBP2+_commu_dotmap.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p4)

## Focus on all cell types, use receiver_focus=[] to include all cell types
p5 = mebo_obj.FlowPlot(
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  ## set to focus on mCCC of interest
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  remove_unrelevant = False,
  and_or='and',
  node_label_size=4,
  node_alpha=0.6,
  figsize='auto',
  node_cmap='Set1',
  line_cmap='bwr',
  line_cmap_vmin = None,
  line_cmap_vmax = None,
  node_size_norm=(20, 150),
  node_value_range = None,
  linewidth_norm=(0.5, 5),
  linewidth_value_range = None,
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=False,
  return_fig = True
)
p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/FlowPlot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/FlowPlot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p5)  

## Focus on all cell types, use receiver_focus=[] to include all cell types
p5 = mebo_obj.FlowPlot(
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  ## set to focus on mCCC of interest
  sender_focus=["Transitional activated"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  remove_unrelevant = False,
  and_or='and',
  node_label_size=4,
  node_alpha=0.6,
  figsize='auto',
  node_cmap='Set1',
  line_cmap='bwr',
  line_cmap_vmin = None,
  line_cmap_vmax = None,
  node_size_norm=(20, 150),
  node_value_range = None,
  linewidth_norm=(0.5, 5),
  linewidth_value_range = None,
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=False,
  return_fig = True
)
p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/sender_focus_GZMK+IL7R+_FlowPlot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/sender_focus_GZMK+IL7R+_FlowPlot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p5)  

p5 = mebo_obj.FlowPlot(
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  ## set to focus on mCCC of interest
  sender_focus=["Cytotoxic effector"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  remove_unrelevant = False,
  and_or='and',
  node_label_size=4,
  node_alpha=0.6,
  figsize='auto',
  node_cmap='Set1',
  line_cmap='bwr',
  line_cmap_vmin = None,
  line_cmap_vmax = None,
  node_size_norm=(20, 150),
  node_value_range = None,
  linewidth_norm=(0.5, 5),
  linewidth_value_range = None,
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=False,
  return_fig = True
)
p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/sender_focus_ZEB2+FGFBP2+_FlowPlot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/sender_focus_ZEB2+FGFBP2+_FlowPlot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p5)  

commu_df = mebo_obj.commu_res.copy()
good_met = commu_df[(commu_df['permutation_test_fdr']<=0.05)]['Metabolite_Name'].sort_values().unique()
p6 = mebo_obj.violin_plot(
  sensor_or_met=good_met[:10], ## Show top10 metabolite expression only, or choose as desired
  cell_focus=[],
  cell_order = [],
  row_zscore = False,
  cmap=None,
  vmin=None,
  vmax=None,
  figsize='auto',
  cbar_title='',
  save=None,
  show_plot=True,
  return_fig = True
)
p6.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/Metabolite_Name_violin_plot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p6.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/Metabolite_Name_violin_plot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p6) 

good_sensor = commu_df[(commu_df['permutation_test_fdr']<=0.05)]['Sensor'].sort_values().unique()

p7 = mebo_obj.violin_plot(
  sensor_or_met=good_sensor[:10],  ## only top 10 as example, or choose as desired
  cell_focus=[],
  cell_order = [],
  row_zscore = False,
  cmap=None,
  vmin=None,
  vmax=None,
  figsize='auto',
  cbar_title='',
  save=None,
  show_plot=True,
  return_fig = True
)

p7.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/sensor_violin_plot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p7.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/sensor_violin_plot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p7) 

#######sender#

p2 = mebo_obj.commu_network_plot(
  sender_focus=["CSMD1+"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  node_cmap= custom_colors,
  figsize=(6,3),
  line_cmap='bwr',
  line_color_vmin=None,
  line_color_vmax=None,
  linewidth_norm=(0.2, 1),
  linewidth_value_range = None,
  node_size_norm=(50, 200),
  node_value_range = None,
  adjust_text_pos_node=True,
  node_text_hidden = True,
  node_text_font=4,  # Set font size
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=True,
  return_fig=True
)

p2.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/sender_focus_GZMK+CSMD1+_network_plot.png", 
  dpi=1600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p2)  

p4 = mebo_obj.commu_dotmap(
  sender_focus=["CSMD1+"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  figsize='auto',
  cmap='bwr',
  cmap_vmin = None,
  cmap_vmax = None,
  cellpair_order=[],
  met_sensor_order=[],
  dot_size_norm=(10, 150),
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_range = None,
  comm_score_cutoff=0,
  swap_axis = False,
  return_fig = True
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/sender_focus_GZMK+CSMD1+_commu_dotmap.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/sender_focus_GZMK+CSMD1+_commu_dotmap.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p4)

p5 = mebo_obj.FlowPlot(
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  ## set to focus on mCCC of interest
  sender_focus=["CSMD1+"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  remove_unrelevant = False,
  and_or='and',
  node_label_size=4,
  node_alpha=0.6,
  figsize='auto',
  node_cmap='Set1',
  line_cmap='bwr',
  line_cmap_vmin = None,
  line_cmap_vmax = None,
  node_size_norm=(20, 150),
  node_value_range = None,
  linewidth_norm=(0.5, 5),
  linewidth_value_range = None,
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=False,
  return_fig = True
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/sender_focus_GZMK+CSMD1+_FlowPlot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/sender_focus_GZMK+CSMD1+_FlowPlot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p5)  

p2 = mebo_obj.commu_network_plot(
  sender_focus=["Transitional activated"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  node_cmap= custom_colors,
  figsize=(6,3),
  line_cmap='bwr',
  line_color_vmin=None,
  line_color_vmax=None,
  linewidth_norm=(0.2, 1),
  linewidth_value_range = None,
  node_size_norm=(50, 200),
  node_value_range = None,
  adjust_text_pos_node=True,
  node_text_hidden = True,
  node_text_font=4,  # Set font size
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=True,
  return_fig=True
)

p2.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/sender_focus_GZMK+IL7R+_network_plot.png", 
  dpi=1600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p2)  

p4 = mebo_obj.commu_dotmap(
  sender_focus=["Transitional activated"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  figsize='auto',
  cmap='bwr',
  cmap_vmin = None,
  cmap_vmax = None,
  cellpair_order=[],
  met_sensor_order=[],
  dot_size_norm=(10, 150),
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_range = None,
  comm_score_cutoff=0,
  swap_axis = False,
  return_fig = True
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/sender_focus_GZMK+IL7R+_commu_dotmap.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/sender_focus_GZMK+IL7R+_commu_dotmap.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p4)

p5 = mebo_obj.FlowPlot(
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  ## set to focus on mCCC of interest
  sender_focus=["Transitional activated"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  remove_unrelevant = False,
  and_or='and',
  node_label_size=4,
  node_alpha=0.6,
  figsize='auto',
  node_cmap='Set1',
  line_cmap='bwr',
  line_cmap_vmin = None,
  line_cmap_vmax = None,
  node_size_norm=(20, 150),
  node_value_range = None,
  linewidth_norm=(0.5, 5),
  linewidth_value_range = None,
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=False,
  return_fig = True
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/sender_focus_GZMK+IL7R+_FlowPlot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/sender_focus_GZMK+IL7R+_FlowPlot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p5)  

p2 = mebo_obj.commu_network_plot(
  sender_focus=["Cytotoxic effector"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  node_cmap= custom_colors,
  figsize=(6,3),
  line_cmap='bwr',
  line_color_vmin=None,
  line_color_vmax=None,
  linewidth_norm=(0.2, 1),
  linewidth_value_range = None,
  node_size_norm=(50, 200),
  node_value_range = None,
  adjust_text_pos_node=True,
  node_text_hidden = True,
  node_text_font=4,  # Set font size
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=True,
  return_fig=True
)

p2.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/sender_focus_ZEB2+FGFBP2+_network_plot.png", 
  dpi=1600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p2)  

p4 = mebo_obj.commu_dotmap(
  sender_focus=["Cytotoxic effector"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  figsize='auto',
  cmap='bwr',
  cmap_vmin = None,
  cmap_vmax = None,
  cellpair_order=[],
  met_sensor_order=[],
  dot_size_norm=(10, 150),
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_range = None,
  comm_score_cutoff=0,
  swap_axis = False,
  return_fig = True
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/sender_focus_ZEB2+FGFBP2+_commu_dotmap.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/sender_focus_ZEB2+FGFBP2+_commu_dotmap.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p4)

p5 = mebo_obj.FlowPlot(
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  ## set to focus on mCCC of interest
  sender_focus=["Cytotoxic effector"],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=[],
  remove_unrelevant = False,
  and_or='and',
  node_label_size=4,
  node_alpha=0.6,
  figsize='auto',
  node_cmap='Set1',
  line_cmap='bwr',
  line_cmap_vmin = None,
  line_cmap_vmax = None,
  node_size_norm=(20, 150),
  node_value_range = None,
  linewidth_norm=(0.5, 5),
  linewidth_value_range = None,
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=False,
  return_fig = True
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/sender_focus_ZEB2+FGFBP2+_FlowPlot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/sender_focus_ZEB2+FGFBP2+_FlowPlot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p5)  

##############receiver
p2 = mebo_obj.commu_network_plot(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["NK-like"],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  node_cmap= custom_colors,
  figsize=(6,3),
  line_cmap='bwr',
  line_color_vmin=None,
  line_color_vmax=None,
  linewidth_norm=(0.2, 1),
  linewidth_value_range = None,
  node_size_norm=(50, 200),
  node_value_range = None,
  adjust_text_pos_node=True,
  node_text_hidden = True,
  node_text_font=4,  # Set font size
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=True,
  return_fig=True
)

p2.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/receiver_focus_ZEB2+KLRF1+_network_plot.png", 
  dpi=1600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p2)  

p4 = mebo_obj.commu_dotmap(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["NK-like"],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  figsize='auto',
  cmap='bwr',
  cmap_vmin = None,
  cmap_vmax = None,
  cellpair_order=[],
  met_sensor_order=[],
  dot_size_norm=(10, 150),
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_range = None,
  comm_score_cutoff=0,
  swap_axis = False,
  return_fig = True
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/receiver_focus_ZEB2+KLRF1+_commu_dotmap1.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/receiver_focus_ZEB2+KLRF1+_commu_dotmap1.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p4)

p5 = mebo_obj.FlowPlot(
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  ## set to focus on mCCC of interest
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["NK-like"],
  remove_unrelevant = False,
  and_or='and',
  node_label_size=4,
  node_alpha=0.6,
  figsize='auto',
  node_cmap='Set1',
  line_cmap='bwr',
  line_cmap_vmin = None,
  line_cmap_vmax = None,
  node_size_norm=(20, 150),
  node_value_range = None,
  linewidth_norm=(0.5, 5),
  linewidth_value_range = None,
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=False,
  return_fig = True
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/receiver_focus_ZEB2+KLRF1+_FlowPlot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/receiver_focus_ZEB2+KLRF1+_FlowPlot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p5)  

p2 = mebo_obj.commu_network_plot(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["Cytotoxic effector"],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  node_cmap= custom_colors,
  figsize=(6,3),
  line_cmap='bwr',
  line_color_vmin=None,
  line_color_vmax=None,
  linewidth_norm=(0.2, 1),
  linewidth_value_range = None,
  node_size_norm=(50, 200),
  node_value_range = None,
  adjust_text_pos_node=True,
  node_text_hidden = True,
  node_text_font=4,  # Set font size
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=True,
  return_fig=True
)

p2.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/receiver_focus_ZEB2+FGFBP2+_network_plot.png", 
  dpi=1600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p2)  

p4 = mebo_obj.commu_dotmap(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["Cytotoxic effector"],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  figsize='auto',
  cmap='bwr',
  cmap_vmin = None,
  cmap_vmax = None,
  cellpair_order=[],
  met_sensor_order=[],
  dot_size_norm=(10, 150),
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_range = None,
  comm_score_cutoff=0,
  swap_axis = False,
  return_fig = True
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/receiver_focus_ZEB2+FGFBP2+_commu_dotmap.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/receiver_focus_ZEB2+FGFBP2+_commu_dotmap.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p4)

p5 = mebo_obj.FlowPlot(
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  ## set to focus on mCCC of interest
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["Cytotoxic effector"],
  remove_unrelevant = False,
  and_or='and',
  node_label_size=4,
  node_alpha=0.6,
  figsize='auto',
  node_cmap='Set1',
  line_cmap='bwr',
  line_cmap_vmin = None,
  line_cmap_vmax = None,
  node_size_norm=(20, 150),
  node_value_range = None,
  linewidth_norm=(0.5, 5),
  linewidth_value_range = None,
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=False,
  return_fig = True
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/receiver_focus_ZEB2+FGFBP2+_FlowPlot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/receiver_focus_ZEB2+FGFBP2+_FlowPlot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p5)  

p2 = mebo_obj.commu_network_plot(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["Terminal-branch"],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  node_cmap= custom_colors,
  figsize=(6,3),
  line_cmap='bwr',
  line_color_vmin=None,
  line_color_vmax=None,
  linewidth_norm=(0.2, 1),
  linewidth_value_range = None,
  node_size_norm=(50, 200),
  node_value_range = None,
  adjust_text_pos_node=True,
  node_text_hidden = True,
  node_text_font=4,  # Set font size
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=True,
  return_fig=True
)

p2.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/receiver_focus_ZEB2+_network_plot.png", 
  dpi=1600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p2)  

p4 = mebo_obj.commu_dotmap(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["Terminal-branch"],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  figsize='auto',
  cmap='bwr',
  cmap_vmin = None,
  cmap_vmax = None,
  cellpair_order=[],
  met_sensor_order=[],
  dot_size_norm=(10, 150),
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_range = None,
  comm_score_cutoff=0,
  swap_axis = False,
  return_fig = True
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/receiver_focus_ZEB2+_commu_dotmap1.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/receiver_focus_ZEB2+_commu_dotmap1.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p4)

p5 = mebo_obj.FlowPlot(
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  ## set to focus on mCCC of interest
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["Terminal-branch"],
  remove_unrelevant = False,
  and_or='and',
  node_label_size=4,
  node_alpha=0.6,
  figsize='auto',
  node_cmap='Set1',
  line_cmap='bwr',
  line_cmap_vmin = None,
  line_cmap_vmax = None,
  node_size_norm=(20, 150),
  node_value_range = None,
  linewidth_norm=(0.5, 5),
  linewidth_value_range = None,
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=False,
  return_fig = True
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/receiver_focus_ZEB2+_FlowPlot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/receiver_focus_ZEB2+_FlowPlot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p5)  

p2 = mebo_obj.commu_network_plot(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["Naive-like"],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  node_cmap= custom_colors,
  figsize=(6,3),
  line_cmap='bwr',
  line_color_vmin=None,
  line_color_vmax=None,
  linewidth_norm=(0.2, 1),
  linewidth_value_range = None,
  node_size_norm=(50, 200),
  node_value_range = None,
  adjust_text_pos_node=True,
  node_text_hidden = True,
  node_text_font=4,  # Set font size
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=True,
  return_fig=True
)

p2.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/receiver_focus_LEF1+NELL2+_network_plot.png", 
  dpi=1600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p2)  

p4 = mebo_obj.commu_dotmap(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["Naive-like"],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  figsize='auto',
  cmap='bwr',
  cmap_vmin = None,
  cmap_vmax = None,
  cellpair_order=[],
  met_sensor_order=[],
  dot_size_norm=(10, 150),
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_range = None,
  comm_score_cutoff=0,
  swap_axis = False,
  return_fig = True
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/receiver_focus_LEF1+NELL2+_commu_dotmap.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/receiver_focus_LEF1+NELL2+_commu_dotmap.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p4)

p5 = mebo_obj.FlowPlot(
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  ## set to focus on mCCC of interest
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["Naive-like"],
  remove_unrelevant = False,
  and_or='and',
  node_label_size=4,
  node_alpha=0.6,
  figsize='auto',
  node_cmap='Set1',
  line_cmap='bwr',
  line_cmap_vmin = None,
  line_cmap_vmax = None,
  node_size_norm=(20, 150),
  node_value_range = None,
  linewidth_norm=(0.5, 5),
  linewidth_value_range = None,
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=False,
  return_fig = True
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/receiver_focus_LEF1+NELL2+_FlowPlot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/receiver_focus_LEF1+NELL2+_FlowPlot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p5)  

p2 = mebo_obj.commu_network_plot(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["Cytotoxic memory-like"],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  node_cmap= custom_colors,
  figsize=(6,3),
  line_cmap='bwr',
  line_color_vmin=None,
  line_color_vmax=None,
  linewidth_norm=(0.2, 1),
  linewidth_value_range = None,
  node_size_norm=(50, 200),
  node_value_range = None,
  adjust_text_pos_node=True,
  node_text_hidden = True,
  node_text_font=4,  # Set font size
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=True,
  return_fig=True
)

p2.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/receiver_focus_FGFBP2+ITGB1+_network_plot.png", 
  dpi=1600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p2)  

p4 = mebo_obj.commu_dotmap(
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["Cytotoxic memory-like"],
  and_or='and',
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  figsize='auto',
  cmap='bwr',
  cmap_vmin = None,
  cmap_vmax = None,
  cellpair_order=[],
  met_sensor_order=[],
  dot_size_norm=(10, 150),
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_range = None,
  comm_score_cutoff=0,
  swap_axis = False,
  return_fig = True
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/receiver_focus_FGFBP2+ITGB1+_commu_dotmap1.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p4.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/receiver_focus_FGFBP2+ITGB1+_commu_dotmap1.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p4)

p5 = mebo_obj.FlowPlot(
  pval_method='permutation_test_fdr',
  pval_cutoff=0.05,
  ## set to focus on mCCC of interest
  sender_focus=[],
  metabolite_focus=[],
  sensor_focus=[],
  receiver_focus=["Cytotoxic memory-like"],
  remove_unrelevant = False,
  and_or='and',
  node_label_size=4,
  node_alpha=0.6,
  figsize='auto',
  node_cmap='Set1',
  line_cmap='bwr',
  line_cmap_vmin = None,
  line_cmap_vmax = None,
  node_size_norm=(20, 150),
  node_value_range = None,
  linewidth_norm=(0.5, 5),
  linewidth_value_range = None,
  save=None,
  show_plot=True,
  comm_score_col='Commu_Score',
  comm_score_cutoff=0,
  text_outline=False,
  return_fig = True
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/receiver_focus_FGFBP2+ITGB1+_FlowPlot.pdf",  
  dpi=600,              # Set DPI
  bbox_inches="tight",  # Keep edges compact
  format="pdf"          # Save as PDF
)

p5.savefig(
  "/thinker/3.tangjiale/shinanxi/metabocost_analysis/MEBOCOST/Psoriasis/plot/receiver_focus_FGFBP2+ITGB1+_FlowPlot.png", 
  dpi=600, 
  bbox_inches="tight",
  format="png"        
)

plt.close(p5)  




##########19.CellOracle#########
######## (1) Convert Seurat object to anndata object
# In R
library(Seurat)
library(sceasy)
# devtools::install_github("cellgeni/sceasy")
library(reticulate)

# Run in linux terminal
conda create -n sceasy_env python=3.9
conda activate sceasy_env
conda install anndata scipy loompy numpy==1.24.4  # fix numpy version to avoid compatibility issues
conda env list

# In R
use_condaenv("/home/ye/envs/sceasy_env")  # specify Python environment
# If Seurat V5, convert Assay type
seurat_obj = readRDS("/thinker/3.tangjiale/shinanxi/monocle2_OE/seurat_ob_pseudotime_20000.rds")
seurat_obj[["RNA"]] <- as(seurat_obj[["RNA"]], "Assay")
# Convert and save as .h5ad
sceasy::convertFormat(seurat_obj, from = "seurat", to = "anndata", outFile = "/thinker/3.tangjiale/shinanxi/celloracle/output_20000.h5ad")


######## Run in linux terminal (linux terminal can install python environment, execute in python environment)
conda activate celloracle_env
cd /thinker/3.tangjiale/shinanxi/celloracle/
  pip install plotly  
python3.10


# In python
# 0.1. Import public libraries
import copy
import glob
import time
import os
import shutil
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from tqdm.auto import tqdm

# 0.2. Import our library
import celloracle as co
from celloracle.applications import Pseudotime_calculator

# 0.3. Plotting parameter settings
# plt.rcParams["font.family"] = "arial"
plt.rcParams["figure.figsize"] = [5,5]
# %config InlineBackend.figure_format = 'retina'
plt.rcParams["savefig.dpi"] = 300
# %matplotlib inline

# 1.1. Load anndata object
adata = sc.read_h5ad("output.h5ad")
# Verify dimensionality reduction coordinate names in Python
print("Available reduction keys:", adata.obsm.keys())
# Available reduction keys: KeysView(AxisArrays with keys: X_harmony, X_mnn, X_pca, X_umap)

# Instantiate pseudotime calculator object
pt = Pseudotime_calculator(adata=adata,
                           obsm_key="X_umap",  # dimensionality reduction data name
                           cluster_column_name="celltype"  # cluster annotation column name
)

# 2. Pseudotime calculation

# 2.1. Add lineage information
# Pseudotime can be calculated for each lineage. We can define these lineages below.

# 2.1.1 Check cluster units
print("Clustering name: ", pt.cluster_column_name)
# Clustering name: subcelltype
print("Cluster list", pt.cluster_list)
# Cluster list ['C1_APOC1', 'C2_NSUN6', 'C3_LTF', 'C4_DEFA3', 'C5_SKAP1', 'C6_MT1X']

# 2.1.2. Define lineages
clusters_in_lineage1 = ["Naive-like","Transitional activated","Activation-regulated","Cytotoxic effector","CSMD1+"]
clusters_in_lineage2 = ["Naive-like","Transitional activated","Activation-regulated","Cytotoxic effector","Cytotoxic memory-like","NK-like","Terminal-branch"]
# Make a dictionary
lineage_dictionary = {"Lineage_1": clusters_in_lineage1,
  "Lineage_2": clusters_in_lineage2}
lineage_dictionary1 = {"Lineage_1": clusters_in_lineage1}
lineage_dictionary2 = {"Lineage_2": clusters_in_lineage2}

# Input lineage information into pseudotime object
pt.set_lineage(lineage_dictionary=lineage_dictionary1)
# Visualize lineage information
pt.plot_lineages()  # plot lineages
plt.savefig("./figures/lineages1_plot.png",
            dpi=600,
            bbox_inches="tight")
plt.savefig("./figures/lineages1_plot.pdf",
            dpi=600,
            bbox_inches="tight")

pt.set_lineage(lineage_dictionary=lineage_dictionary2)
# Visualize lineage information
pt.plot_lineages()  # plot lineages
plt.savefig("./figures/lineages2_plot.png",
            dpi=600,
            bbox_inches="tight")
plt.savefig("./figures/lineages2_plot.pdf",
            dpi=600,
            bbox_inches="tight")

# Use plotly to visualize interactive UMAP to find root cell names
import plotly.express as px
import pandas as pd
# Create UMAP chart
df = pd.DataFrame(adata.obsm["X_umap"], columns=["x", "y"])
df["cluster"] = adata.obs["celltype"]
df["cell_id"] = adata.obs.index
fig = px.scatter(df, x="x", y="y", hover_name="cell_id", color="cluster")

fig.write_html("cell_umap_interactive.html", include_plotlyjs="cdn")
## Open the html file in a browser on Windows

### Select root cells from the html
pt.set_lineage(lineage_dictionary=lineage_dictionary1)
root_cells = {"Lineage_1": "NPBMC1_ATGCGTACGGATGCTTGGTT"}
pt.set_root_cells(root_cells=root_cells)
pt.plot_root_cells()
plt.savefig("./figures/lineages1_root_plot.png",
            dpi=600,
            bbox_inches="tight")
plt.savefig("./figures/lineages1_root_plot.pdf",
            dpi=600,
            bbox_inches="tight")

pt.set_lineage(lineage_dictionary=lineage_dictionary2)
root_cells = {"Lineage_2": "NPBMC1_ATGCGTACGGATGCTTGGTT"}
pt.set_root_cells(root_cells=root_cells)
pt.plot_root_cells()
plt.savefig("./figures/lineages2_root_plot.png",
            dpi=600,
            bbox_inches="tight")
plt.savefig("./figures/lineages2_root_plot.pdf",
            dpi=600,
            bbox_inches="tight")

pt.set_lineage(lineage_dictionary=lineage_dictionary)
root_cells = {"Lineage_1": "NPBMC1_ATGCGTACGGATGCTTGGTT","Lineage_2": "NPBMC1_ATGCGTACGGATGCTTGGTT"}
pt.set_root_cells(root_cells=root_cells)

# 2.3.1. Check diffusion map
# Check diffusion map data.
"X_diffmap" in pt.adata.obsm
# False
# Calculate diffusion map if your anndata object does not have diffusion map data.
sc.pp.neighbors(pt.adata, n_neighbors=30)
sc.tl.diffmap(pt.adata)
"X_diffmap" in pt.adata.obsm
# True

# 2.3.2. Calculate pseudotime
# Calculate pseudotime
pt.get_pseudotime_per_each_lineage()

# Check results
pt.plot_pseudotime(cmap="rainbow")
plt.savefig("./figures/lineages_pseudotime_plot.png",
            dpi=600,
            bbox_inches="tight")
plt.savefig("./figures/lineages_pseudotime_plot.pdf",
            dpi=600,
            bbox_inches="tight")

pt.adata.obs[["Pseudotime"]].head()

# 3. Save data
# Add calculated pseudotime data to the oracle object
adata.obs = pt.adata.obs

# Save updated anndata object
adata.write_h5ad("/thinker/3.tangjiale/shinanxi/celloracle/adata.h5ad")


# In linux terminal
# Initialization
conda activate celloracle_env
cd /thinker/3.tangjiale/shinanxi/celloracle/
  python3.10


import os
import sys

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

import celloracle as co
plt.rcParams["figure.figsize"] = [6,6]
plt.rcParams["savefig.dpi"] = 600

save_folder = "figures"
os.makedirs(save_folder, exist_ok=True)

# 1. Load data
### Load scRNA data
adata = sc.read_h5ad('/thinker/3.tangjiale/shinanxi/celloracle/output_20000.h5ad')
### Load Base GRN data
base_GRN = co.data.load_human_promoter_base_GRN()
base_GRN.head()

# 2. Create Oracle object
# After reading in the two required data, we first construct an Oracle object.
# When importing the single-cell transcriptome matrix, you need to provide the clustering name and the dimension reduction name.
# The clustering name is stored in the obs attribute of the AnnData, which can be viewed via adata.obs.columns.
# The dimension reduction name is stored in the obsm attribute, which can be viewed via adata.obsm.keys().
oracle = co.Oracle()
oracle.import_anndata_as_raw_count(adata=adata,
                                   cluster_column_name="celltype",
                                   embedding_name="X_umap")
oracle.import_TF_data(TF_info_matrix=base_GRN)

# CellOracle uses the same strategy as velocyto to visualize cell migration, which requires KNN imputation in advance.
# Then, cluster-specific GRNs are built for all clusters. The Links object stores the inferred GRN and corresponding metadata. Most network structure analyses are performed through the links object.
### KNN imputation
oracle.perform_PCA()
n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]

n_comps = min(n_comps, 50)
print(n_comps)
n_cell = oracle.adata.shape[0]
print(f"cell number is :{n_cell}")
k = int(0.025*n_cell)
print(f"Auto-selected k is :{k}")
oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                      b_maxl=k*4, n_jobs=4)

### Save object
oracle.to_hdf5("/thinker/3.tangjiale/shinanxi/celloracle/wound.celloracle.oracle")

### GRN calculation
links = oracle.get_links(cluster_name_for_GRN_unit="celltype", alpha=10, verbose_level=10)
links.links_dict.keys()
links.to_hdf5(file_path="/thinker/3.tangjiale/shinanxi/celloracle/wound.celloracle.links")

oracle = co.load_hdf5("wound.celloracle.oracle")
links = co.load_hdf5(file_path="wound.celloracle.links")

# Here we first filter the links and re-fit the ridge regression model.
### Make predictive models for simulation
links.filter_links()
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)

# Next, we simulate the effect of TF perturbation on cell state to study its potential function and regulatory mechanism. We select CD44 as the transcription factor to be studied.
# In perturbation simulation calculations, we can use any expression value, but avoid extremely high values far from the natural gene expression range; the allowed upper limit is twice the maximum gene expression value. The paper used 0 and 2 as thresholds for studying knockout and overexpression. We take knockout as an example for code demonstration.
### Knock-out

##### ZEB1 Knock-out 
oracle.simulate_shift(perturb_condition={'ZEB1': 0.0}, n_propagation=3)
oracle.estimate_transition_prob(n_neighbors=200,
                                knn_random=True,
                                sampled_fraction=1)
# Calculate embedding
oracle.calculate_embedding_shift(sigma_corr=0.05)

# Next, we need to filter parameters to achieve optimal visualization. We visualize the simulation results as a vector field on a digitized grid, where the transition vectors of single cells are grouped by grid points.
# First, we need to determine the n_grid and min_mass parameters. n_grid is the number of grid points, and min_mass is the threshold for cell density. The official recommendation is to start with n_grid=40. The function oracle.suggest_mass_thresholds() displays a series of min_mass parameter values for selection.
n_grid = 40
oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
oracle.suggest_mass_thresholds(n_suggestion=12)
plt.savefig('./figures/kO_ZEB1_vector_fields.png', bbox_inches='tight')
plt.savefig('./figures/kO_ZEB1_vector_fields.pdf', bbox_inches='tight')
min_mass=7.4
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)

# Next, we need to find an optimal scale parameter. The scale parameter controls the length of the vectors. If the vectors are not visible, try reducing the parameter to enlarge the vector length. However, if large vectors appear in the randomized results (right side), it indicates the scale parameter is too small.
import matplotlib
matplotlib.use('Agg')  # Set backend to Agg before importing plt
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=[13, 6])
scale_simulation = 12
# Show quiver plot
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: ZEB1 KO")
fig.savefig("./figures/kO_ZEB1_Simulated_cell_identity_shift_vector_ZEB1_KO.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/kO_ZEB1_Simulated_cell_identity_shift_vector_ZEB1_KO.pdf", dpi=600, bbox_inches='tight') 
# Show quiver plot that was calculated with randomized graph.
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")
fig.savefig("./figures/kO_ZEB1_Randomized_simulation_vector.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/kO_ZEB1_Randomized_simulation_vector.pdf", dpi=600, bbox_inches='tight') 

# Next, we can simulate how TF perturbation affects cell state as in the previous steps, and visualize the results as a vector field plot. Compare the simulated perturbation vectors with the developmental gradient vectors.
# The pseudotime data is first transformed into an n x n digitized grid, then the 2D gradient of pseudotime is calculated to obtain the vector field, and finally the two vector fields are compared by calculating the inner product of these vectors. Through comparison, we can intuitively understand how transcription factors influence cell fate decisions during development.
from celloracle.applications import Gradient_calculator
# Instantiate Gradient calculator object
gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")
gradient.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
gradient.calculate_mass_filter(min_mass=min_mass, plot=True)
plt.savefig('./figures/kO_ZEB1_Rmass_filter.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/kO_ZEB1_Rmass_filter.pdf', dpi=600, bbox_inches='tight')
gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly":3}, plot=True)
plt.savefig('./figures/kO_ZEB1_Rdata_grid.png', dpi=600, transparent=True)
plt.savefig('./figures/kO_ZEB1_Rdata_grid.pdf', dpi=600, transparent=True)
# Calculate gradient
gradient.calculate_gradient()

# Show results
scale_dev = 40
gradient.visualize_results(scale=scale_dev, s=5)
plt.savefig('./figures/kO_ZEB1_Rgradient_results.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/kO_ZEB1_Rgradient_results.pdf', dpi=600, bbox_inches='tight')

# The length of the vectors also reflects the magnitude of the inner product value, defined as the perturbation score (PS). Negative PS means that TF perturbation hinders differentiation, while positive PS means that TF perturbation promotes differentiation.
from celloracle.applications import Oracle_development_module

# Make Oracle_development_module to compare two vector fields
dev = Oracle_development_module()

# Load development flow
dev.load_differentiation_reference_data(gradient_object=gradient)

# Load simulation result
dev.load_perturb_simulation_data(oracle_object=oracle)

# Calculate inner product scores
dev.calculate_inner_product()
dev.calculate_digitized_ip(n_bins=10)

# Here, we still need to adjust the vm parameter to achieve optimal visualization of the PS score. If no color is visible on the left, you can use a smaller vm parameter to enlarge the visual scale. If color appears on the right, it indicates the vm parameter is too small.
vm = 0.1
fig, ax = plt.subplots(1, 2, figsize=[12, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax[0])
ax[0].set_title(f"PS")
fig.savefig('./figures/kO_ZEB1_RPS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/kO_ZEB1_RPS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)
dev.plot_inner_product_random_on_grid(vm=vm, s=50, ax=ax[1])
ax[1].set_title(f"PS calculated with Randomized simulation vector")
fig.savefig('./figures/kO_ZEB1_Randomized_PS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/kO_ZEB1_Randomized_PS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)

# Show perturbation scores with perturbation simulation vector field
fig, ax = plt.subplots(figsize=[6, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax)
dev.plot_simulation_flow_on_grid(scale=scale_simulation, show_background=False, ax=ax)
fig.savefig('./figures/kO_ZEB1_final.png', dpi=300, bbox_inches='tight')
fig.savefig('./figures/kO_ZEB1_final.pdf', dpi=300, bbox_inches='tight')


##### ZEB1 Overexpression 1.5x 
oracle.simulate_shift(perturb_condition={'ZEB1': 1.5}, n_propagation=3)
oracle.estimate_transition_prob(n_neighbors=200,
                                knn_random=True,
                                sampled_fraction=1)
# Calculate embedding
oracle.calculate_embedding_shift(sigma_corr=0.05)

# Next, we need to filter parameters to achieve optimal visualization...
n_grid = 40
oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
oracle.suggest_mass_thresholds(n_suggestion=12)
plt.savefig('./figures/OE15_ZEB1_vector_fields.png', bbox_inches='tight')
plt.savefig('./figures/OE15_ZEB1_vector_fields.pdf', bbox_inches='tight')
min_mass=7.4
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=[13, 6])
scale_simulation = 12
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: ZEB1 OE15")
fig.savefig("./figures/OE15_ZEB1_Simulated_cell_identity_shift_vector_ZEB1_OE15.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/OE15_ZEB1_Simulated_cell_identity_shift_vector_ZEB1_OE15.pdf", dpi=600, bbox_inches='tight') 
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")
fig.savefig("./figures/OE15_ZEB1_Randomized_simulation_vector.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/OE15_ZEB1_Randomized_simulation_vector.pdf", dpi=600, bbox_inches='tight') 

from celloracle.applications import Gradient_calculator
gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")
gradient.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
gradient.calculate_mass_filter(min_mass=min_mass, plot=True)
plt.savefig('./figures/OE15_ZEB1_Rmass_filter.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/OE15_ZEB1_Rmass_filter.pdf', dpi=600, bbox_inches='tight')
gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly":3}, plot=True)
plt.savefig('./figures/OE15_ZEB1_Rdata_grid.png', dpi=600, transparent=True)
plt.savefig('./figures/OE15_ZEB1_Rdata_grid.pdf', dpi=600, transparent=True)
gradient.calculate_gradient()

scale_dev = 40
gradient.visualize_results(scale=scale_dev, s=5)
plt.savefig('./figures/OE15_ZEB1_Rgradient_results.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/OE15_ZEB1_Rgradient_results.pdf', dpi=600, bbox_inches='tight')

from celloracle.applications import Oracle_development_module
dev = Oracle_development_module()
dev.load_differentiation_reference_data(gradient_object=gradient)
dev.load_perturb_simulation_data(oracle_object=oracle)
dev.calculate_inner_product()
dev.calculate_digitized_ip(n_bins=10)

vm = 0.1
fig, ax = plt.subplots(1, 2, figsize=[12, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax[0])
ax[0].set_title(f"PS")
fig.savefig('./figures/OE15_ZEB1_RPS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/OE15_ZEB1_RPS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)
dev.plot_inner_product_random_on_grid(vm=vm, s=50, ax=ax[1])
ax[1].set_title(f"PS calculated with Randomized simulation vector")
fig.savefig('./figures/OE15_ZEB1_Randomized_PS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/OE15_ZEB1_Randomized_PS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)

fig, ax = plt.subplots(figsize=[6, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax)
dev.plot_simulation_flow_on_grid(scale=scale_simulation, show_background=False, ax=ax)
fig.savefig('./figures/OE15_ZEB1_final.png', dpi=300, bbox_inches='tight')
fig.savefig('./figures/OE15_ZEB1_final.pdf', dpi=300, bbox_inches="tight")


##### FOSB Knock-out
oracle.simulate_shift(perturb_condition={'FOSB': 0.0}, n_propagation=3)
oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
oracle.calculate_embedding_shift(sigma_corr=0.05)

n_grid = 40
oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
oracle.suggest_mass_thresholds(n_suggestion=12)
plt.savefig('./figures/kO_FOSB_vector_fields.png', bbox_inches='tight')
plt.savefig('./figures/kO_FOSB_vector_fields.pdf', bbox_inches='tight')
min_mass=7.4
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=[13, 6])
scale_simulation = 12
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: FOSB KO")
fig.savefig("./figures/kO_FOSB_Simulated_cell_identity_shift_vector_FOSB_KO.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/kO_FOSB_Simulated_cell_identity_shift_vector_FOSB_KO.pdf", dpi=600, bbox_inches='tight') 
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")
fig.savefig("./figures/kO_FOSB_Randomized_simulation_vector.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/kO_FOSB_Randomized_simulation_vector.pdf", dpi=600, bbox_inches='tight') 

from celloracle.applications import Gradient_calculator
gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")
gradient.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
gradient.calculate_mass_filter(min_mass=min_mass, plot=True)
plt.savefig('./figures/kO_FOSB_Rmass_filter.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/kO_FOSB_Rmass_filter.pdf', dpi=600, bbox_inches='tight')
gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly":3}, plot=True)
plt.savefig('./figures/kO_FOSB_Rdata_grid.png', dpi=600, transparent=True)
plt.savefig('./figures/kO_FOSB_Rdata_grid.pdf', dpi=600, transparent=True)
gradient.calculate_gradient()

scale_dev = 40
gradient.visualize_results(scale=scale_dev, s=5)
plt.savefig('./figures/kO_FOSB_Rgradient_results.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/kO_FOSB_Rgradient_results.pdf', dpi=600, bbox_inches='tight')

from celloracle.applications import Oracle_development_module
dev = Oracle_development_module()
dev.load_differentiation_reference_data(gradient_object=gradient)
dev.load_perturb_simulation_data(oracle_object=oracle)
dev.calculate_inner_product()
dev.calculate_digitized_ip(n_bins=10)

vm = 0.1
fig, ax = plt.subplots(1, 2, figsize=[12, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax[0])
ax[0].set_title(f"PS")
fig.savefig('./figures/kO_FOSB_RPS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/kO_FOSB_RPS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)
dev.plot_inner_product_random_on_grid(vm=vm, s=50, ax=ax[1])
ax[1].set_title(f"PS calculated with Randomized simulation vector")
fig.savefig('./figures/kO_FOSB_Randomized_PS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/kO_FOSB_Randomized_PS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)

fig, ax = plt.subplots(figsize=[6, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax)
dev.plot_simulation_flow_on_grid(scale=scale_simulation, show_background=False, ax=ax)
fig.savefig('./figures/kO_FOSB_final.png', dpi=300, bbox_inches='tight')
fig.savefig('./figures/kO_FOSB_final.pdf', dpi=300, bbox_inches='tight')


##### FOSB Overexpression 1x 
oracle.simulate_shift(perturb_condition={'FOSB': 1}, n_propagation=3)
oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
oracle.calculate_embedding_shift(sigma_corr=0.05)

n_grid = 40
oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
oracle.suggest_mass_thresholds(n_suggestion=12)
plt.savefig('./figures/OE10_FOSB_vector_fields.png', bbox_inches='tight')
plt.savefig('./figures/OE10_FOSB_vector_fields.pdf', bbox_inches='tight')
min_mass=7.4
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=[13, 6])
scale_simulation = 12
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: FOSB OE10")
fig.savefig("./figures/OE10_FOSB_Simulated_cell_identity_shift_vector_FOSB_OE10.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/OE10_FOSB_Simulated_cell_identity_shift_vector_FOSB_OE10.pdf", dpi=600, bbox_inches='tight') 
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")
fig.savefig("./figures/OE10_FOSB_Randomized_simulation_vector.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/OE10_FOSB_Randomized_simulation_vector.pdf", dpi=600, bbox_inches='tight') 

from celloracle.applications import Gradient_calculator
gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")
gradient.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
gradient.calculate_mass_filter(min_mass=min_mass, plot=True)
plt.savefig('./figures/OE10_FOSB_Rmass_filter.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/OE10_FOSB_Rmass_filter.pdf', dpi=600, bbox_inches='tight')
gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly":3}, plot=True)
plt.savefig('./figures/OE10_FOSB_Rdata_grid.png', dpi=600, transparent=True)
plt.savefig('./figures/OE10_FOSB_Rdata_grid.pdf', dpi=600, transparent=True)
gradient.calculate_gradient()

scale_dev = 40
gradient.visualize_results(scale=scale_dev, s=5)
plt.savefig('./figures/OE10_FOSB_Rgradient_results.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/OE10_FOSB_Rgradient_results.pdf', dpi=600, bbox_inches='tight')

from celloracle.applications import Oracle_development_module
dev = Oracle_development_module()
dev.load_differentiation_reference_data(gradient_object=gradient)
dev.load_perturb_simulation_data(oracle_object=oracle)
dev.calculate_inner_product()
dev.calculate_digitized_ip(n_bins=10)

vm = 0.1
fig, ax = plt.subplots(1, 2, figsize=[12, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax[0])
ax[0].set_title(f"PS")
fig.savefig('./figures/OE10_FOSB_RPS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/OE10_FOSB_RPS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)
dev.plot_inner_product_random_on_grid(vm=vm, s=50, ax=ax[1])
ax[1].set_title(f"PS calculated with Randomized simulation vector")
fig.savefig('./figures/OE10_FOSB_Randomized_PS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/OE10_FOSB_Randomized_PS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)

fig, ax = plt.subplots(figsize=[6, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax)
dev.plot_simulation_flow_on_grid(scale=scale_simulation, show_background=False, ax=ax)
fig.savefig('./figures/OE10_FOSB_final.png', dpi=300, bbox_inches='tight')
fig.savefig('./figures/OE10_FOSB_final.pdf', dpi=300, bbox_inches="tight")


##### JUN Knock-out 
oracle.simulate_shift(perturb_condition={'JUN': 0.0}, n_propagation=3)
oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
oracle.calculate_embedding_shift(sigma_corr=0.05)

n_grid = 40
oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
oracle.suggest_mass_thresholds(n_suggestion=12)
plt.savefig('./figures/kO_JUN_vector_fields.png', bbox_inches='tight')
plt.savefig('./figures/kO_JUN_vector_fields.pdf', bbox_inches='tight')
min_mass=7.4
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=[13, 6])
scale_simulation = 12
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: JUN KO")
fig.savefig("./figures/kO_JUN_Simulated_cell_identity_shift_vector_JUN_KO.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/kO_JUN_Simulated_cell_identity_shift_vector_JUN_KO.pdf", dpi=600, bbox_inches='tight') 
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")
fig.savefig("./figures/kO_JUN_Randomized_simulation_vector.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/kO_JUN_Randomized_simulation_vector.pdf", dpi=600, bbox_inches='tight') 

from celloracle.applications import Gradient_calculator
gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")
gradient.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
gradient.calculate_mass_filter(min_mass=min_mass, plot=True)
plt.savefig('./figures/kO_JUN_Rmass_filter.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/kO_JUN_Rmass_filter.pdf', dpi=600, bbox_inches='tight')
gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly":3}, plot=True)
plt.savefig('./figures/kO_JUN_Rdata_grid.png', dpi=600, transparent=True)
plt.savefig('./figures/kO_JUN_Rdata_grid.pdf', dpi=600, transparent=True)
gradient.calculate_gradient()

scale_dev = 40
gradient.visualize_results(scale=scale_dev, s=5)
plt.savefig('./figures/kO_JUN_Rgradient_results.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/kO_JUN_Rgradient_results.pdf', dpi=600, bbox_inches='tight')

from celloracle.applications import Oracle_development_module
dev = Oracle_development_module()
dev.load_differentiation_reference_data(gradient_object=gradient)
dev.load_perturb_simulation_data(oracle_object=oracle)
dev.calculate_inner_product()
dev.calculate_digitized_ip(n_bins=10)

vm = 0.1
fig, ax = plt.subplots(1, 2, figsize=[12, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax[0])
ax[0].set_title(f"PS")
fig.savefig('./figures/kO_JUN_RPS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/kO_JUN_RPS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)
dev.plot_inner_product_random_on_grid(vm=vm, s=50, ax=ax[1])
ax[1].set_title(f"PS calculated with Randomized simulation vector")
fig.savefig('./figures/kO_JUN_Randomized_PS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/kO_JUN_Randomized_PS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)

fig, ax = plt.subplots(figsize=[6, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax)
dev.plot_simulation_flow_on_grid(scale=scale_simulation, show_background=False, ax=ax)
fig.savefig('./figures/kO_JUN_final.png', dpi=300, bbox_inches='tight')
fig.savefig('./figures/kO_JUN_final.pdf', dpi=300, bbox_inches='tight')


##### JUN Overexpression 1.5x
oracle.simulate_shift(perturb_condition={'JUN': 1.5}, n_propagation=3)
oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
oracle.calculate_embedding_shift(sigma_corr=0.05)

n_grid = 40
oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
oracle.suggest_mass_thresholds(n_suggestion=12)
plt.savefig('./figures/OE15_JUN_vector_fields.png', bbox_inches='tight')
plt.savefig('./figures/OE15_JUN_vector_fields.pdf', bbox_inches='tight')
min_mass=7.4
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=[13, 6])
scale_simulation = 12
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: JUN OE15")
fig.savefig("./figures/OE15_JUN_Simulated_cell_identity_shift_vector_JUN_OE15.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/OE15_JUN_Simulated_cell_identity_shift_vector_JUN_OE15.pdf", dpi=600, bbox_inches='tight') 
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")
fig.savefig("./figures/OE15_JUN_Randomized_simulation_vector.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/OE15_JUN_Randomized_simulation_vector.pdf", dpi=600, bbox_inches='tight') 

from celloracle.applications import Gradient_calculator
gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")
gradient.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
gradient.calculate_mass_filter(min_mass=min_mass, plot=True)
plt.savefig('./figures/OE15_JUN_Rmass_filter.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/OE15_JUN_Rmass_filter.pdf', dpi=600, bbox_inches='tight')
gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly":3}, plot=True)
plt.savefig('./figures/OE15_JUN_Rdata_grid.png', dpi=600, transparent=True)
plt.savefig('./figures/OE15_JUN_Rdata_grid.pdf', dpi=600, transparent=True)
gradient.calculate_gradient()

scale_dev = 40
gradient.visualize_results(scale=scale_dev, s=5)
plt.savefig('./figures/OE15_JUN_Rgradient_results.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/OE15_JUN_Rgradient_results.pdf', dpi=600, bbox_inches='tight')

from celloracle.applications import Oracle_development_module
dev = Oracle_development_module()
dev.load_differentiation_reference_data(gradient_object=gradient)
dev.load_perturb_simulation_data(oracle_object=oracle)
dev.calculate_inner_product()
dev.calculate_digitized_ip(n_bins=10)

vm = 0.1
fig, ax = plt.subplots(1, 2, figsize=[12, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax[0])
ax[0].set_title(f"PS")
fig.savefig('./figures/OE15_JUN_RPS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/OE15_JUN_RPS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)
dev.plot_inner_product_random_on_grid(vm=vm, s=50, ax=ax[1])
ax[1].set_title(f"PS calculated with Randomized simulation vector")
fig.savefig('./figures/OE15_JUN_Randomized_PS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/OE15_JUN_Randomized_PS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)

fig, ax = plt.subplots(figsize=[6, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax)
dev.plot_simulation_flow_on_grid(scale=scale_simulation, show_background=False, ax=ax)
fig.savefig('./figures/OE15_JUN_final.png', dpi=300, bbox_inches='tight')
fig.savefig('./figures/OE15_JUN_final.pdf', dpi=300, bbox_inches="tight")


##### NFKB1 Knock-out 
oracle.simulate_shift(perturb_condition={'NFKB1': 0.0}, n_propagation=3)
oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
oracle.calculate_embedding_shift(sigma_corr=0.05)

n_grid = 40
oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
oracle.suggest_mass_thresholds(n_suggestion=12)
plt.savefig('./figures/kO_NFKB1_vector_fields.png', bbox_inches='tight')
plt.savefig('./figures/kO_NFKB1_vector_fields.pdf', bbox_inches='tight')
min_mass=7.4
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=[13, 6])
scale_simulation = 12
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: NFKB1 KO")
fig.savefig("./figures/kO_NFKB1_Simulated_cell_identity_shift_vector_NFKB1_KO.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/kO_NFKB1_Simulated_cell_identity_shift_vector_NFKB1_KO.pdf", dpi=600, bbox_inches='tight') 
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")
fig.savefig("./figures/kO_NFKB1_Randomized_simulation_vector.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/kO_NFKB1_Randomized_simulation_vector.pdf", dpi=600, bbox_inches='tight') 

from celloracle.applications import Gradient_calculator
gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")
gradient.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
gradient.calculate_mass_filter(min_mass=min_mass, plot=True)
plt.savefig('./figures/kO_NFKB1_Rmass_filter.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/kO_NFKB1_Rmass_filter.pdf', dpi=600, bbox_inches='tight')
gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly":3}, plot=True)
plt.savefig('./figures/kO_NFKB1_Rdata_grid.png', dpi=600, transparent=True)
plt.savefig('./figures/kO_NFKB1_Rdata_grid.pdf', dpi=600, transparent=True)
gradient.calculate_gradient()

scale_dev = 40
gradient.visualize_results(scale=scale_dev, s=5)
plt.savefig('./figures/kO_NFKB1_Rgradient_results.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/kO_NFKB1_Rgradient_results.pdf', dpi=600, bbox_inches='tight')

from celloracle.applications import Oracle_development_module
dev = Oracle_development_module()
dev.load_differentiation_reference_data(gradient_object=gradient)
dev.load_perturb_simulation_data(oracle_object=oracle)
dev.calculate_inner_product()
dev.calculate_digitized_ip(n_bins=10)

vm = 0.1
fig, ax = plt.subplots(1, 2, figsize=[12, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax[0])
ax[0].set_title(f"PS")
fig.savefig('./figures/kO_NFKB1_RPS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/kO_NFKB1_RPS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)
dev.plot_inner_product_random_on_grid(vm=vm, s=50, ax=ax[1])
ax[1].set_title(f"PS calculated with Randomized simulation vector")
fig.savefig('./figures/kO_NFKB1_Randomized_PS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/kO_NFKB1_Randomized_PS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)

fig, ax = plt.subplots(figsize=[6, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax)
dev.plot_simulation_flow_on_grid(scale=scale_simulation, show_background=False, ax=ax)
fig.savefig('./figures/kO_NFKB1_final.png', dpi=300, bbox_inches='tight')
fig.savefig('./figures/kO_NFKB1_final.pdf', dpi=300, bbox_inches='tight')


##### NFKB1 Overexpression 1.5x
oracle.simulate_shift(perturb_condition={'NFKB1': 1.5}, n_propagation=3)
oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
oracle.calculate_embedding_shift(sigma_corr=0.05)

n_grid = 40
oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
oracle.suggest_mass_thresholds(n_suggestion=12)
plt.savefig('./figures/OE15_NFKB1_vector_fields.png', bbox_inches='tight')
plt.savefig('./figures/OE15_NFKB1_vector_fields.pdf', bbox_inches='tight')
min_mass=7.4
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=[13, 6])
scale_simulation = 12
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: NFKB1 OE15")
fig.savefig("./figures/OE15_NFKB1_Simulated_cell_identity_shift_vector_NFKB1_OE15.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/OE15_NFKB1_Simulated_cell_identity_shift_vector_NFKB1_OE15.pdf", dpi=600, bbox_inches='tight') 
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")
fig.savefig("./figures/OE15_NFKB1_Randomized_simulation_vector.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/OE15_NFKB1_Randomized_simulation_vector.pdf", dpi=600, bbox_inches='tight') 

from celloracle.applications import Gradient_calculator
gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")
gradient.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
gradient.calculate_mass_filter(min_mass=min_mass, plot=True)
plt.savefig('./figures/OE15_NFKB1_Rmass_filter.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/OE15_NFKB1_Rmass_filter.pdf', dpi=600, bbox_inches='tight')
gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly":3}, plot=True)
plt.savefig('./figures/OE15_NFKB1_Rdata_grid.png', dpi=600, transparent=True)
plt.savefig('./figures/OE15_NFKB1_Rdata_grid.pdf', dpi=600, transparent=True)
gradient.calculate_gradient()

scale_dev = 40
gradient.visualize_results(scale=scale_dev, s=5)
plt.savefig('./figures/OE15_NFKB1_Rgradient_results.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/OE15_NFKB1_Rgradient_results.pdf', dpi=600, bbox_inches='tight')

from celloracle.applications import Oracle_development_module
dev = Oracle_development_module()
dev.load_differentiation_reference_data(gradient_object=gradient)
dev.load_perturb_simulation_data(oracle_object=oracle)
dev.calculate_inner_product()
dev.calculate_digitized_ip(n_bins=10)

vm = 0.1
fig, ax = plt.subplots(1, 2, figsize=[12, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax[0])
ax[0].set_title(f"PS")
fig.savefig('./figures/OE15_NFKB1_RPS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/OE15_NFKB1_RPS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)
dev.plot_inner_product_random_on_grid(vm=vm, s=50, ax=ax[1])
ax[1].set_title(f"PS calculated with Randomized simulation vector")
fig.savefig('./figures/OE15_NFKB1_Randomized_PS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/OE15_NFKB1_Randomized_PS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)

fig, ax = plt.subplots(figsize=[6, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax)
dev.plot_simulation_flow_on_grid(scale=scale_simulation, show_background=False, ax=ax)
fig.savefig('./figures/OE15_NFKB1_final.png', dpi=300, bbox_inches='tight')
fig.savefig('./figures/OE15_NFKB1_final.pdf', dpi=300, bbox_inches="tight")


##### STAT1 Knock-out #####
oracle.simulate_shift(perturb_condition={'STAT1': 0.0}, n_propagation=3)
oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
oracle.calculate_embedding_shift(sigma_corr=0.05)

n_grid = 40
oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
oracle.suggest_mass_thresholds(n_suggestion=12)
plt.savefig('./figures/kO_STAT1_vector_fields.png', bbox_inches='tight')
plt.savefig('./figures/kO_STAT1_vector_fields.pdf', bbox_inches='tight')
min_mass=7.4
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=[13, 6])
scale_simulation = 12
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: STAT1 KO")
fig.savefig("./figures/kO_STAT1_Simulated_cell_identity_shift_vector_STAT1_KO.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/kO_STAT1_Simulated_cell_identity_shift_vector_STAT1_KO.pdf", dpi=600, bbox_inches='tight') 
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")
fig.savefig("./figures/kO_STAT1_Randomized_simulation_vector.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/kO_STAT1_Randomized_simulation_vector.pdf", dpi=600, bbox_inches='tight') 

from celloracle.applications import Gradient_calculator
gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")
gradient.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
gradient.calculate_mass_filter(min_mass=min_mass, plot=True)
plt.savefig('./figures/kO_STAT1_Rmass_filter.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/kO_STAT1_Rmass_filter.pdf', dpi=600, bbox_inches='tight')
gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly":3}, plot=True)
plt.savefig('./figures/kO_STAT1_Rdata_grid.png', dpi=600, transparent=True)
plt.savefig('./figures/kO_STAT1_Rdata_grid.pdf', dpi=600, transparent=True)
gradient.calculate_gradient()

scale_dev = 40
gradient.visualize_results(scale=scale_dev, s=5)
plt.savefig('./figures/kO_STAT1_Rgradient_results.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/kO_STAT1_Rgradient_results.pdf', dpi=600, bbox_inches='tight')

from celloracle.applications import Oracle_development_module
dev = Oracle_development_module()
dev.load_differentiation_reference_data(gradient_object=gradient)
dev.load_perturb_simulation_data(oracle_object=oracle)
dev.calculate_inner_product()
dev.calculate_digitized_ip(n_bins=10)

vm = 0.1
fig, ax = plt.subplots(1, 2, figsize=[12, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax[0])
ax[0].set_title(f"PS")
fig.savefig('./figures/kO_STAT1_RPS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/kO_STAT1_RPS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)
dev.plot_inner_product_random_on_grid(vm=vm, s=50, ax=ax[1])
ax[1].set_title(f"PS calculated with Randomized simulation vector")
fig.savefig('./figures/kO_STAT1_Randomized_PS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/kO_STAT1_Randomized_PS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)

fig, ax = plt.subplots(figsize=[6, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax)
dev.plot_simulation_flow_on_grid(scale=scale_simulation, show_background=False, ax=ax)
fig.savefig('./figures/kO_STAT1_final.png', dpi=300, bbox_inches='tight')
fig.savefig('./figures/kO_STAT1_final.pdf', dpi=300, bbox_inches='tight')


##### STAT1 Overexpression 1.1x #####
oracle.simulate_shift(perturb_condition={'STAT1': 1.1}, n_propagation=3)
oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
oracle.calculate_embedding_shift(sigma_corr=0.05)

n_grid = 40
oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
oracle.suggest_mass_thresholds(n_suggestion=12)
plt.savefig('./figures/OE11_STAT1_vector_fields.png', bbox_inches='tight')
plt.savefig('./figures/OE11_STAT1_vector_fields.pdf', bbox_inches='tight')
min_mass=7.4
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=[13, 6])
scale_simulation = 12
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: STAT1 OE11")
fig.savefig("./figures/OE11_STAT1_Simulated_cell_identity_shift_vector_STAT1_OE11.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/OE11_STAT1_Simulated_cell_identity_shift_vector_STAT1_OE11.pdf", dpi=600, bbox_inches='tight') 
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")
fig.savefig("./figures/OE11_STAT1_Randomized_simulation_vector.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/OE11_STAT1_Randomized_simulation_vector.pdf", dpi=600, bbox_inches='tight') 

from celloracle.applications import Gradient_calculator
gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")
gradient.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
gradient.calculate_mass_filter(min_mass=min_mass, plot=True)
plt.savefig('./figures/OE11_STAT1_Rmass_filter.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/OE11_STAT1_Rmass_filter.pdf', dpi=600, bbox_inches='tight')
gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly":3}, plot=True)
plt.savefig('./figures/OE11_STAT1_Rdata_grid.png', dpi=600, transparent=True)
plt.savefig('./figures/OE11_STAT1_Rdata_grid.pdf', dpi=600, transparent=True)
gradient.calculate_gradient()

scale_dev = 40
gradient.visualize_results(scale=scale_dev, s=5)
plt.savefig('./figures/OE11_STAT1_Rgradient_results.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/OE11_STAT1_Rgradient_results.pdf', dpi=600, bbox_inches='tight')

from celloracle.applications import Oracle_development_module
dev = Oracle_development_module()
dev.load_differentiation_reference_data(gradient_object=gradient)
dev.load_perturb_simulation_data(oracle_object=oracle)
dev.calculate_inner_product()
dev.calculate_digitized_ip(n_bins=10)

vm = 0.1
fig, ax = plt.subplots(1, 2, figsize=[12, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax[0])
ax[0].set_title(f"PS")
fig.savefig('./figures/OE11_STAT1_RPS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/OE11_STAT1_RPS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)
dev.plot_inner_product_random_on_grid(vm=vm, s=50, ax=ax[1])
ax[1].set_title(f"PS calculated with Randomized simulation vector")
fig.savefig('./figures/OE11_STAT1_Randomized_PS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/OE11_STAT1_Randomized_PS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)

fig, ax = plt.subplots(figsize=[6, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax)
dev.plot_simulation_flow_on_grid(scale=scale_simulation, show_background=False, ax=ax)
fig.savefig('./figures/OE11_STAT1_final.png', dpi=300, bbox_inches='tight')
fig.savefig('./figures/OE11_STAT1_final.pdf', dpi=300, bbox_inches="tight")


##### FOXP1 Knock-out #####
oracle.simulate_shift(perturb_condition={'FOXP1': 0.0}, n_propagation=3)
oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
oracle.calculate_embedding_shift(sigma_corr=0.05)

n_grid = 40
oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
oracle.suggest_mass_thresholds(n_suggestion=12)
plt.savefig('./figures/kO_FOXP1_vector_fields.png', bbox_inches='tight')
plt.savefig('./figures/kO_FOXP1_vector_fields.pdf', bbox_inches='tight')
min_mass=7.4
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=[13, 6])
scale_simulation = 12
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: FOXP1 KO")
fig.savefig("./figures/kO_FOXP1_Simulated_cell_identity_shift_vector_FOXP1_KO.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/kO_FOXP1_Simulated_cell_identity_shift_vector_FOXP1_KO.pdf", dpi=600, bbox_inches='tight') 
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")
fig.savefig("./figures/kO_FOXP1_Randomized_simulation_vector.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/kO_FOXP1_Randomized_simulation_vector.pdf", dpi=600, bbox_inches='tight') 

from celloracle.applications import Gradient_calculator
gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")
gradient.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
gradient.calculate_mass_filter(min_mass=min_mass, plot=True)
plt.savefig('./figures/kO_FOXP1_Rmass_filter.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/kO_FOXP1_Rmass_filter.pdf', dpi=600, bbox_inches='tight')
gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly":3}, plot=True)
plt.savefig('./figures/kO_FOXP1_Rdata_grid.png', dpi=600, transparent=True)
plt.savefig('./figures/kO_FOXP1_Rdata_grid.pdf', dpi=600, transparent=True)
gradient.calculate_gradient()

scale_dev = 40
gradient.visualize_results(scale=scale_dev, s=5)
plt.savefig('./figures/kO_FOXP1_Rgradient_results.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/kO_FOXP1_Rgradient_results.pdf', dpi=600, bbox_inches='tight')

from celloracle.applications import Oracle_development_module
dev = Oracle_development_module()
dev.load_differentiation_reference_data(gradient_object=gradient)
dev.load_perturb_simulation_data(oracle_object=oracle)
dev.calculate_inner_product()
dev.calculate_digitized_ip(n_bins=10)

vm = 0.1
fig, ax = plt.subplots(1, 2, figsize=[12, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax[0])
ax[0].set_title(f"PS")
fig.savefig('./figures/kO_FOXP1_RPS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/kO_FOXP1_RPS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)
dev.plot_inner_product_random_on_grid(vm=vm, s=50, ax=ax[1])
ax[1].set_title(f"PS calculated with Randomized simulation vector")
fig.savefig('./figures/kO_FOXP1_Randomized_PS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/kO_FOXP1_Randomized_PS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)

fig, ax = plt.subplots(figsize=[6, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax)
dev.plot_simulation_flow_on_grid(scale=scale_simulation, show_background=False, ax=ax)
fig.savefig('./figures/kO_FOXP1_final.png', dpi=300, bbox_inches='tight')
fig.savefig('./figures/kO_FOXP1_final.pdf', dpi=300, bbox_inches='tight')


##### FOXP1 Overexpression 1.5x #####
oracle.simulate_shift(perturb_condition={'FOXP1': 1.5}, n_propagation=3)
oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
oracle.calculate_embedding_shift(sigma_corr=0.05)

n_grid = 40
oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
oracle.suggest_mass_thresholds(n_suggestion=12)
plt.savefig('./figures/OE15_FOXP1_vector_fields.png', bbox_inches='tight')
plt.savefig('./figures/OE15_FOXP1_vector_fields.pdf', bbox_inches='tight')
min_mass=7.4
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=[13, 6])
scale_simulation = 12
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: FOXP1 OE15")
fig.savefig("./figures/OE15_FOXP1_Simulated_cell_identity_shift_vector_FOXP1_OE15.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/OE15_FOXP1_Simulated_cell_identity_shift_vector_FOXP1_OE15.pdf", dpi=600, bbox_inches='tight') 
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")
fig.savefig("./figures/OE15_FOXP1_Randomized_simulation_vector.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/OE15_FOXP1_Randomized_simulation_vector.pdf", dpi=600, bbox_inches='tight') 

from celloracle.applications import Gradient_calculator
gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")
gradient.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
gradient.calculate_mass_filter(min_mass=min_mass, plot=True)
plt.savefig('./figures/OE15_FOXP1_Rmass_filter.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/OE15_FOXP1_Rmass_filter.pdf', dpi=600, bbox_inches='tight')
gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly":3}, plot=True)
plt.savefig('./figures/OE15_FOXP1_Rdata_grid.png', dpi=600, transparent=True)
plt.savefig('./figures/OE15_FOXP1_Rdata_grid.pdf', dpi=600, transparent=True)
gradient.calculate_gradient()

scale_dev = 40
gradient.visualize_results(scale=scale_dev, s=5)
plt.savefig('./figures/OE15_FOXP1_Rgradient_results.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/OE15_FOXP1_Rgradient_results.pdf', dpi=600, bbox_inches='tight')

from celloracle.applications import Oracle_development_module
dev = Oracle_development_module()
dev.load_differentiation_reference_data(gradient_object=gradient)
dev.load_perturb_simulation_data(oracle_object=oracle)
dev.calculate_inner_product()
dev.calculate_digitized_ip(n_bins=10)

vm = 0.1
fig, ax = plt.subplots(1, 2, figsize=[12, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax[0])
ax[0].set_title(f"PS")
fig.savefig('./figures/OE15_FOXP1_RPS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/OE15_FOXP1_RPS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)
dev.plot_inner_product_random_on_grid(vm=vm, s=50, ax=ax[1])
ax[1].set_title(f"PS calculated with Randomized simulation vector")
fig.savefig('./figures/OE15_FOXP1_Randomized_PS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/OE15_FOXP1_Randomized_PS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)

fig, ax = plt.subplots(figsize=[6, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax)
dev.plot_simulation_flow_on_grid(scale=scale_simulation, show_background=False, ax=ax)
fig.savefig('./figures/OE15_FOXP1_final.png', dpi=300, bbox_inches='tight')
fig.savefig('./figures/OE15_FOXP1_final.pdf', dpi=300, bbox_inches="tight")


##### SOX4 Knock-out 
oracle.simulate_shift(perturb_condition={'SOX4': 0.0}, n_propagation=3)
oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
oracle.calculate_embedding_shift(sigma_corr=0.05)

n_grid = 40
oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
oracle.suggest_mass_thresholds(n_suggestion=12)
plt.savefig('./figures/kO_SOX4_vector_fields.png', bbox_inches='tight')
plt.savefig('./figures/kO_SOX4_vector_fields.pdf', bbox_inches='tight')
min_mass=7.4
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=[13, 6])
scale_simulation = 12
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: SOX4 KO")
fig.savefig("./figures/kO_SOX4_Simulated_cell_identity_shift_vector_SOX4_KO.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/kO_SOX4_Simulated_cell_identity_shift_vector_SOX4_KO.pdf", dpi=600, bbox_inches='tight') 
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")
fig.savefig("./figures/kO_SOX4_Randomized_simulation_vector.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/kO_SOX4_Randomized_simulation_vector.pdf", dpi=600, bbox_inches='tight') 

from celloracle.applications import Gradient_calculator
gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")
gradient.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
gradient.calculate_mass_filter(min_mass=min_mass, plot=True)
plt.savefig('./figures/kO_SOX4_Rmass_filter.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/kO_SOX4_Rmass_filter.pdf', dpi=600, bbox_inches='tight')
gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly":3}, plot=True)
plt.savefig('./figures/kO_SOX4_Rdata_grid.png', dpi=600, transparent=True)
plt.savefig('./figures/kO_SOX4_Rdata_grid.pdf', dpi=600, transparent=True)
gradient.calculate_gradient()

scale_dev = 40
gradient.visualize_results(scale=scale_dev, s=5)
plt.savefig('./figures/kO_SOX4_Rgradient_results.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/kO_SOX4_Rgradient_results.pdf', dpi=600, bbox_inches='tight')

from celloracle.applications import Oracle_development_module
dev = Oracle_development_module()
dev.load_differentiation_reference_data(gradient_object=gradient)
dev.load_perturb_simulation_data(oracle_object=oracle)
dev.calculate_inner_product()
dev.calculate_digitized_ip(n_bins=10)

vm = 0.1
fig, ax = plt.subplots(1, 2, figsize=[12, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax[0])
ax[0].set_title(f"PS")
fig.savefig('./figures/kO_SOX4_RPS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/kO_SOX4_RPS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)
dev.plot_inner_product_random_on_grid(vm=vm, s=50, ax=ax[1])
ax[1].set_title(f"PS calculated with Randomized simulation vector")
fig.savefig('./figures/kO_SOX4_Randomized_PS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/kO_SOX4_Randomized_PS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)

fig, ax = plt.subplots(figsize=[6, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax)
dev.plot_simulation_flow_on_grid(scale=scale_simulation, show_background=False, ax=ax)
fig.savefig('./figures/kO_SOX4_final.png', dpi=300, bbox_inches='tight')
fig.savefig('./figures/kO_SOX4_final.pdf', dpi=300, bbox_inches='tight')


##### ELF1 Knock-out
oracle.simulate_shift(perturb_condition={'ELF1': 0.0}, n_propagation=3)
oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
oracle.calculate_embedding_shift(sigma_corr=0.05)

n_grid = 40
oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
oracle.suggest_mass_thresholds(n_suggestion=12)
plt.savefig('./figures/kO_ELF1_vector_fields.png', bbox_inches='tight')
plt.savefig('./figures/kO_ELF1_vector_fields.pdf', bbox_inches='tight')
min_mass=7.4
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=[13, 6])
scale_simulation = 12
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: ELF1 KO")
fig.savefig("./figures/kO_ELF1_Simulated_cell_identity_shift_vector_ELF1_KO.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/kO_ELF1_Simulated_cell_identity_shift_vector_ELF1_KO.pdf", dpi=600, bbox_inches='tight') 
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")
fig.savefig("./figures/kO_ELF1_Randomized_simulation_vector.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/kO_ELF1_Randomized_simulation_vector.pdf", dpi=600, bbox_inches='tight') 

from celloracle.applications import Gradient_calculator
gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")
gradient.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
gradient.calculate_mass_filter(min_mass=min_mass, plot=True)
plt.savefig('./figures/kO_ELF1_Rmass_filter.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/kO_ELF1_Rmass_filter.pdf', dpi=600, bbox_inches='tight')
gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly":3}, plot=True)
plt.savefig('./figures/kO_ELF1_Rdata_grid.png', dpi=600, transparent=True)
plt.savefig('./figures/kO_ELF1_Rdata_grid.pdf', dpi=600, transparent=True)
gradient.calculate_gradient()

scale_dev = 40
gradient.visualize_results(scale=scale_dev, s=5)
plt.savefig('./figures/kO_ELF1_Rgradient_results.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/kO_ELF1_Rgradient_results.pdf', dpi=600, bbox_inches='tight')

from celloracle.applications import Oracle_development_module
dev = Oracle_development_module()
dev.load_differentiation_reference_data(gradient_object=gradient)
dev.load_perturb_simulation_data(oracle_object=oracle)
dev.calculate_inner_product()
dev.calculate_digitized_ip(n_bins=10)

vm = 0.1
fig, ax = plt.subplots(1, 2, figsize=[12, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax[0])
ax[0].set_title(f"PS")
fig.savefig('./figures/kO_ELF1_RPS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/kO_ELF1_RPS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)
dev.plot_inner_product_random_on_grid(vm=vm, s=50, ax=ax[1])
ax[1].set_title(f"PS calculated with Randomized simulation vector")
fig.savefig('./figures/kO_ELF1_Randomized_PS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/kO_ELF1_Randomized_PS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)

fig, ax = plt.subplots(figsize=[6, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax)
dev.plot_simulation_flow_on_grid(scale=scale_simulation, show_background=False, ax=ax)
fig.savefig('./figures/kO_ELF1_final.png', dpi=300, bbox_inches='tight')
fig.savefig('./figures/kO_ELF1_final.pdf', dpi=300, bbox_inches='tight')


##### ELF1 Overexpression 1.5x 
oracle.simulate_shift(perturb_condition={'ELF1': 1.5}, n_propagation=3)
oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
oracle.calculate_embedding_shift(sigma_corr=0.05)

n_grid = 40
oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
oracle.suggest_mass_thresholds(n_suggestion=12)
plt.savefig('./figures/OE15_ELF1_vector_fields.png', bbox_inches='tight')
plt.savefig('./figures/OE15_ELF1_vector_fields.pdf', bbox_inches='tight')
min_mass=7.4
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=[13, 6])
scale_simulation = 12
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: ELF1 OE15")
fig.savefig("./figures/OE15_ELF1_Simulated_cell_identity_shift_vector_ELF1_OE15.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/OE15_ELF1_Simulated_cell_identity_shift_vector_ELF1_OE15.pdf", dpi=600, bbox_inches='tight') 
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")
fig.savefig("./figures/OE15_ELF1_Randomized_simulation_vector.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/OE15_ELF1_Randomized_simulation_vector.pdf", dpi=600, bbox_inches='tight') 

from celloracle.applications import Gradient_calculator
gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")
gradient.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
gradient.calculate_mass_filter(min_mass=min_mass, plot=True)
plt.savefig('./figures/OE15_ELF1_Rmass_filter.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/OE15_ELF1_Rmass_filter.pdf', dpi=600, bbox_inches='tight')
gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly":3}, plot=True)
plt.savefig('./figures/OE15_ELF1_Rdata_grid.png', dpi=600, transparent=True)
plt.savefig('./figures/OE15_ELF1_Rdata_grid.pdf', dpi=600, transparent=True)
gradient.calculate_gradient()

scale_dev = 40
gradient.visualize_results(scale=scale_dev, s=5)
plt.savefig('./figures/OE15_ELF1_Rgradient_results.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/OE15_ELF1_Rgradient_results.pdf', dpi=600, bbox_inches='tight')

from celloracle.applications import Oracle_development_module
dev = Oracle_development_module()
dev.load_differentiation_reference_data(gradient_object=gradient)
dev.load_perturb_simulation_data(oracle_object=oracle)
dev.calculate_inner_product()
dev.calculate_digitized_ip(n_bins=10)

vm = 0.1
fig, ax = plt.subplots(1, 2, figsize=[12, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax[0])
ax[0].set_title(f"PS")
fig.savefig('./figures/OE15_ELF1_RPS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/OE15_ELF1_RPS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)
dev.plot_inner_product_random_on_grid(vm=vm, s=50, ax=ax[1])
ax[1].set_title(f"PS calculated with Randomized simulation vector")
fig.savefig('./figures/OE15_ELF1_Randomized_PS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/OE15_ELF1_Randomized_PS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)

fig, ax = plt.subplots(figsize=[6, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax)
dev.plot_simulation_flow_on_grid(scale=scale_simulation, show_background=False, ax=ax)
fig.savefig('./figures/OE15_ELF1_final.png', dpi=300, bbox_inches='tight')
fig.savefig('./figures/OE15_ELF1_final.pdf', dpi=300, bbox_inches="tight")


##### ETS1 Knock-out
oracle.simulate_shift(perturb_condition={'ETS1': 0.0}, n_propagation=3)
oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
oracle.calculate_embedding_shift(sigma_corr=0.05)

n_grid = 40
oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
oracle.suggest_mass_thresholds(n_suggestion=12)
plt.savefig('./figures/kO_ETS1_vector_fields.png', bbox_inches='tight')
plt.savefig('./figures/kO_ETS1_vector_fields.pdf', bbox_inches='tight')
min_mass=7.4
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=[13, 6])
scale_simulation = 12
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: ETS1 KO")
fig.savefig("./figures/kO_ETS1_Simulated_cell_identity_shift_vector_ETS1_KO.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/kO_ETS1_Simulated_cell_identity_shift_vector_ETS1_KO.pdf", dpi=600, bbox_inches='tight') 
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")
fig.savefig("./figures/kO_ETS1_Randomized_simulation_vector.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/kO_ETS1_Randomized_simulation_vector.pdf", dpi=600, bbox_inches='tight') 

from celloracle.applications import Gradient_calculator
gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")
gradient.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
gradient.calculate_mass_filter(min_mass=min_mass, plot=True)
plt.savefig('./figures/kO_ETS1_Rmass_filter.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/kO_ETS1_Rmass_filter.pdf', dpi=600, bbox_inches='tight')
gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly":3}, plot=True)
plt.savefig('./figures/kO_ETS1_Rdata_grid.png', dpi=600, transparent=True)
plt.savefig('./figures/kO_ETS1_Rdata_grid.pdf', dpi=600, transparent=True)
gradient.calculate_gradient()

scale_dev = 40
gradient.visualize_results(scale=scale_dev, s=5)
plt.savefig('./figures/kO_ETS1_Rgradient_results.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/kO_ETS1_Rgradient_results.pdf', dpi=600, bbox_inches='tight')

from celloracle.applications import Oracle_development_module
dev = Oracle_development_module()
dev.load_differentiation_reference_data(gradient_object=gradient)
dev.load_perturb_simulation_data(oracle_object=oracle)
dev.calculate_inner_product()
dev.calculate_digitized_ip(n_bins=10)

vm = 0.1
fig, ax = plt.subplots(1, 2, figsize=[12, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax[0])
ax[0].set_title(f"PS")
fig.savefig('./figures/kO_ETS1_RPS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/kO_ETS1_RPS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)
dev.plot_inner_product_random_on_grid(vm=vm, s=50, ax=ax[1])
ax[1].set_title(f"PS calculated with Randomized simulation vector")
fig.savefig('./figures/kO_ETS1_Randomized_PS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/kO_ETS1_Randomized_PS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)

fig, ax = plt.subplots(figsize=[6, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax)
dev.plot_simulation_flow_on_grid(scale=scale_simulation, show_background=False, ax=ax)
fig.savefig('./figures/kO_ETS1_final.png', dpi=300, bbox_inches='tight')
fig.savefig('./figures/kO_ETS1_final.pdf', dpi=300, bbox_inches='tight')


##### ETS1 Overexpression 1.5x
oracle.simulate_shift(perturb_condition={'ETS1': 1.5}, n_propagation=3)
oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
oracle.calculate_embedding_shift(sigma_corr=0.05)

n_grid = 40
oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
oracle.suggest_mass_thresholds(n_suggestion=12)
plt.savefig('./figures/OE15_ETS1_vector_fields.png', bbox_inches='tight')
plt.savefig('./figures/OE15_ETS1_vector_fields.pdf', bbox_inches='tight')
min_mass=7.4
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=[13, 6])
scale_simulation = 12
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: ETS1 OE15")
fig.savefig("./figures/OE15_ETS1_Simulated_cell_identity_shift_vector_ETS1_OE15.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/OE15_ETS1_Simulated_cell_identity_shift_vector_ETS1_OE15.pdf", dpi=600, bbox_inches='tight') 
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")
fig.savefig("./figures/OE15_ETS1_Randomized_simulation_vector.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/OE15_ETS1_Randomized_simulation_vector.pdf", dpi=600, bbox_inches='tight') 

from celloracle.applications import Gradient_calculator
gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")
gradient.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
gradient.calculate_mass_filter(min_mass=min_mass, plot=True)
plt.savefig('./figures/OE15_ETS1_Rmass_filter.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/OE15_ETS1_Rmass_filter.pdf', dpi=600, bbox_inches='tight')
gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly":3}, plot=True)
plt.savefig('./figures/OE15_ETS1_Rdata_grid.png', dpi=600, transparent=True)
plt.savefig('./figures/OE15_ETS1_Rdata_grid.pdf', dpi=600, transparent=True)
gradient.calculate_gradient()

scale_dev = 40
gradient.visualize_results(scale=scale_dev, s=5)
plt.savefig('./figures/OE15_ETS1_Rgradient_results.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/OE15_ETS1_Rgradient_results.pdf', dpi=600, bbox_inches='tight')

from celloracle.applications import Oracle_development_module
dev = Oracle_development_module()
dev.load_differentiation_reference_data(gradient_object=gradient)
dev.load_perturb_simulation_data(oracle_object=oracle)
dev.calculate_inner_product()
dev.calculate_digitized_ip(n_bins=10)

vm = 0.1
fig, ax = plt.subplots(1, 2, figsize=[12, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax[0])
ax[0].set_title(f"PS")
fig.savefig('./figures/OE15_ETS1_RPS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/OE15_ETS1_RPS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)
dev.plot_inner_product_random_on_grid(vm=vm, s=50, ax=ax[1])
ax[1].set_title(f"PS calculated with Randomized simulation vector")
fig.savefig('./figures/OE15_ETS1_Randomized_PS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/OE15_ETS1_Randomized_PS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)

fig, ax = plt.subplots(figsize=[6, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax)
dev.plot_simulation_flow_on_grid(scale=scale_simulation, show_background=False, ax=ax)
fig.savefig('./figures/OE15_ETS1_final.png', dpi=300, bbox_inches='tight')
fig.savefig('./figures/OE15_ETS1_final.pdf', dpi=300, bbox_inches="tight")


##### NFATC3 Knock-out #####
oracle.simulate_shift(perturb_condition={'NFATC3': 0.0}, n_propagation=3)
oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
oracle.calculate_embedding_shift(sigma_corr=0.05)

n_grid = 40
oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
oracle.suggest_mass_thresholds(n_suggestion=12)
plt.savefig('./figures/kO_NFATC3_vector_fields.png', bbox_inches='tight')
plt.savefig('./figures/kO_NFATC3_vector_fields.pdf', bbox_inches='tight')
min_mass=7.4
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=[13, 6])
scale_simulation = 12
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: NFATC3 KO")
fig.savefig("./figures/kO_NFATC3_Simulated_cell_identity_shift_vector_NFATC3_KO.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/kO_NFATC3_Simulated_cell_identity_shift_vector_NFATC3_KO.pdf", dpi=600, bbox_inches='tight') 
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")
fig.savefig("./figures/kO_NFATC3_Randomized_simulation_vector.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/kO_NFATC3_Randomized_simulation_vector.pdf", dpi=600, bbox_inches='tight') 

from celloracle.applications import Gradient_calculator
gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")
gradient.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
gradient.calculate_mass_filter(min_mass=min_mass, plot=True)
plt.savefig('./figures/kO_NFATC3_Rmass_filter.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/kO_NFATC3_Rmass_filter.pdf', dpi=600, bbox_inches='tight')
gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly":3}, plot=True)
plt.savefig('./figures/kO_NFATC3_Rdata_grid.png', dpi=600, transparent=True)
plt.savefig('./figures/kO_NFATC3_Rdata_grid.pdf', dpi=600, transparent=True)
gradient.calculate_gradient()

scale_dev = 40
gradient.visualize_results(scale=scale_dev, s=5)
plt.savefig('./figures/kO_NFATC3_Rgradient_results.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/kO_NFATC3_Rgradient_results.pdf', dpi=600, bbox_inches='tight')

from celloracle.applications import Oracle_development_module
dev = Oracle_development_module()
dev.load_differentiation_reference_data(gradient_object=gradient)
dev.load_perturb_simulation_data(oracle_object=oracle)
dev.calculate_inner_product()
dev.calculate_digitized_ip(n_bins=10)

vm = 0.1
fig, ax = plt.subplots(1, 2, figsize=[12, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax[0])
ax[0].set_title(f"PS")
fig.savefig('./figures/kO_NFATC3_RPS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/kO_NFATC3_RPS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)
dev.plot_inner_product_random_on_grid(vm=vm, s=50, ax=ax[1])
ax[1].set_title(f"PS calculated with Randomized simulation vector")
fig.savefig('./figures/kO_NFATC3_Randomized_PS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/kO_NFATC3_Randomized_PS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)

fig, ax = plt.subplots(figsize=[6, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax)
dev.plot_simulation_flow_on_grid(scale=scale_simulation, show_background=False, ax=ax)
fig.savefig('./figures/kO_NFATC3_final.png', dpi=300, bbox_inches='tight')
fig.savefig('./figures/kO_NFATC3_final.pdf', dpi=300, bbox_inches='tight')


##### NFATC3 Overexpression 1.5x #####
oracle.simulate_shift(perturb_condition={'NFATC3': 1.5}, n_propagation=3)
oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
oracle.calculate_embedding_shift(sigma_corr=0.05)

n_grid = 40
oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
oracle.suggest_mass_thresholds(n_suggestion=12)
plt.savefig('./figures/OE15_NFATC3_vector_fields.png', bbox_inches='tight')
plt.savefig('./figures/OE15_NFATC3_vector_fields.pdf', bbox_inches='tight')
min_mass=7.4
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=[13, 6])
scale_simulation = 12
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: NFATC3 OE15")
fig.savefig("./figures/OE15_NFATC3_Simulated_cell_identity_shift_vector_NFATC3_OE15.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/OE15_NFATC3_Simulated_cell_identity_shift_vector_NFATC3_OE15.pdf", dpi=600, bbox_inches='tight') 
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")
fig.savefig("./figures/OE15_NFATC3_Randomized_simulation_vector.png", dpi=600, bbox_inches='tight') 
fig.savefig("./figures/OE15_NFATC3_Randomized_simulation_vector.pdf", dpi=600, bbox_inches='tight') 

from celloracle.applications import Gradient_calculator
gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")
gradient.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
gradient.calculate_mass_filter(min_mass=min_mass, plot=True)
plt.savefig('./figures/OE15_NFATC3_Rmass_filter.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/OE15_NFATC3_Rmass_filter.pdf', dpi=600, bbox_inches='tight')
gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly":3}, plot=True)
plt.savefig('./figures/OE15_NFATC3_Rdata_grid.png', dpi=600, transparent=True)
plt.savefig('./figures/OE15_NFATC3_Rdata_grid.pdf', dpi=600, transparent=True)
gradient.calculate_gradient()

scale_dev = 40
gradient.visualize_results(scale=scale_dev, s=5)
plt.savefig('./figures/OE15_NFATC3_Rgradient_results.png', dpi=600, bbox_inches='tight')
plt.savefig('./figures/OE15_NFATC3_Rgradient_results.pdf', dpi=600, bbox_inches='tight')

from celloracle.applications import Oracle_development_module
dev = Oracle_development_module()
dev.load_differentiation_reference_data(gradient_object=gradient)
dev.load_perturb_simulation_data(oracle_object=oracle)
dev.calculate_inner_product()
dev.calculate_digitized_ip(n_bins=10)

vm = 0.1
fig, ax = plt.subplots(1, 2, figsize=[12, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax[0])
ax[0].set_title(f"PS")
fig.savefig('./figures/OE15_NFATC3_RPS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/OE15_NFATC3_RPS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)
dev.plot_inner_product_random_on_grid(vm=vm, s=50, ax=ax[1])
ax[1].set_title(f"PS calculated with Randomized simulation vector")
fig.savefig('./figures/OE15_NFATC3_Randomized_PS.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
fig.savefig('./figures/OE15_NFATC3_Randomized_PS.pdf', dpi=600, bbox_inches='tight', pad_inches=0.5)

fig, ax = plt.subplots(figsize=[6, 6])
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax)
dev.plot_simulation_flow_on_grid(scale=scale_simulation, show_background=False, ax=ax)
fig.savefig('./figures/OE15_NFATC3_final.png', dpi=300, bbox_inches='tight')
fig.savefig('./figures/OE15_NFATC3_final.pdf', dpi=300, bbox_inches="tight")

