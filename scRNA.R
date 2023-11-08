library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
#library(ggunchull)
library(tidydr)
library(ggsci)
library(Cairo)
library(harmony)
library(dplyr)
library(jsonlite)
library(lubridate)
library(tidyverse)
##scRNA sequencing analysis for GSE125970
##extract scRNA-seq for Rectum tissue from all tissues
data<-read.table("E:\\drug\\pQTL\\HO\\GSE125970_raw_UMIcounts.txt",row.names=1,fill=T)
colnames(data)<-data[1,]
meta<-data.frame(colnames(data))
meta$V1<-meta$colnames.data.
meta<-meta %>%
  separate(V1, c("cells", "sample","cell_type"), "_")
data11<-subset(meta,V3==c("Rectum-2"))
data12<-subset(meta,V3=="Rectum-1")
data0<-rbind(data11,data12)
data3<-data[,colnames(data)%in%data0$colnames.data.]  ##subset the Rectum samples from all samples
data31<-data.matrix(data3[-1,])

##split the raw meta into sample, cell types, and individuals cells
meta1<-data.frame(colnames(data3))
meta1$V1<-meta1$colnames.data3.
meta1<-meta1 %>%
  separate(V1, c("cells", "sample","cell_type"), "_")
rownames(meta1)<-meta1$colnames.data3.
scRNA_HEM <- CreateSeuratObject(counts = data31, meta.data = meta1,min.cells = 3, project = "WT")
scRNA_HEM@meta.data$orig.ident<-scRNA_HEM@meta.data$V3
scRNA_HEM <- SplitObject(scRNA_HEM, split.by = "orig.ident")

for (i in 1:length(scRNA_HEM)) {
  scRNA_HEM[[i]] <- NormalizeData(scRNA_HEM[[i]])
  scRNA_HEM[[i]] <- FindVariableFeatures(scRNA_HEM[[i]], selection.method = "vst")
}
##
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNA_HEM)

##
scRNA_HEM <- IntegrateData(anchorset = scRNA.anchors)
dim(scRNA_HEM)
#2000 3898
table(scRNA_HEM@meta.data$orig.ident)
#Rectum-1 Rectum-2 
#2842     1056

##integrated assay
DefaultAssay(scRNA_HEM) <- "integrated"
scRNA_HEM <- ScaleData(scRNA_HEM, verbose = FALSE, vars.to.regress = c("nCount_RNA"))
scRNA_HEM <- RunPCA(scRNA_HEM, verbose = FALSE,npcs = 100)

###cluster
scRNA_HEM <- FindNeighbors(object = scRNA_HEM, dims = 1:20)
scRNA_HEM <- FindClusters(object = scRNA_HEM, resolution = 0.14) 

###umap
scRNA_HEM <- RunUMAP(scRNA_HEM, reduction = "pca", dims = 1:20)
DimPlot(object = scRNA_HEM, reduction = 'umap',label = TRUE,group.by = 'cell_type',pt.size = 1) +
  scale_color_manual(values = pal <- c("#009292","#ff6db6","#db6d00",
                                                "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
                                                "#920000","#924900"))
DefaultAssay(scRNA_HEM) <- "RNA"
DotPlot(scRNA_HEM, features = c("ERLEC1"),group.by = "cell_type") + RotatedAxis() + scale_size(range = c(1, 8))+theme(legend.position="top")
FeaturePlot(scRNA_HEM,features =c("ERLEC1"),pt.size = 1) & scale_color_viridis_c()
library(ggsignif)
library(ggpubr)
my_comparisons <- list( c("Enteriendocrine", "TA"),c("Stem", "TA"), c("Enterocyte", "TA"), c("Goblet", "TA"),c("Paneth-like", "TA"), c("Progenitor", "TA"))
VlnPlot(scRNA_HEM,features = "ERLEC1",group.by = "cell_type",pt.size = 0)+
  scale_fill_manual(values = pal <- c("#009292","#ff6db6","#db6d00",
                                               "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
                                               "#920000","#924900"))+NoLegend()+
                                                 geom_signif(comparisons = my_comparisons)+ylim(0, 4)