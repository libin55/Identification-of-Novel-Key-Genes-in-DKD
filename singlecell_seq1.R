#=======================================================
#set the working files and load the packages
#=======================================================

rm(list=ls())
library(dplyr)
library(tidyr)
library(Biobase)
library(limma)
library(gdata)
library(limma)
library(edgeR)
library(dplyr)
library(tidyr)
library(Biobase)
library(limma)
library(Seurat)
library(gridExtra)
library(ggplot2)
library(reshape2)
library("R.utils")

setwd('D:\\SCIwork\\F49DKD\\s0single\\GSE131882_RAW')



gunzip("GSM3823939_control.s1.dgecounts.rds.gz")
gunzip("GSM3823940_control.s2.dgecounts.rds.gz")
gunzip("GSM3823941_control.s3.dgecounts.rds.gz")
gunzip("GSM3823942_diabetes.s1.dgecounts.rds.gz")
gunzip("GSM3823943_diabetes.s2.dgecounts.rds.gz")
gunzip("GSM3823944_diabetes.s3.dgecounts.rds.gz")


#=======================================================
#
#=======================================================


con1 <- readRDS('GSM3823939_control.s1.dgecounts.rds')
con1 <- con1[["readcount"]][["exon"]][["all"]]
sobj1 <- CreateSeuratObject(counts = con1)
dim(sobj1)
sobj1 <- subset(sobj1, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 )
dim(sobj1)




con2 <- readRDS('GSM3823940_control.s2.dgecounts.rds')
con2 <- con2[["readcount"]][["exon"]][["all"]]
sobj2 <- CreateSeuratObject(counts = con2)
dim(sobj2)
sobj2 <- subset(sobj2, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 )
dim(sobj2)



con3 <- readRDS('GSM3823941_control.s3.dgecounts.rds')
con3 <- con3[["readcount"]][["exon"]][["all"]]
sobj3 <- CreateSeuratObject(counts = con3)
dim(sobj3)
sobj3 <- subset(sobj3, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 )
dim(sobj3)


con4 <- readRDS('GSM3823942_diabetes.s1.dgecounts.rds')
con4 <- con4[["readcount"]][["exon"]][["all"]]
sobj4 <- CreateSeuratObject(counts = con4)
dim(sobj4)
sobj4 <- subset(sobj4, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 )
dim(sobj4)




con5 <- readRDS('GSM3823943_diabetes.s2.dgecounts.rds')
con5 <- con5[["readcount"]][["exon"]][["all"]]
sobj5 <- CreateSeuratObject(counts = con5)
dim(sobj5)
sobj5 <- subset(sobj5, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 )
dim(sobj2)



con6 <- readRDS('GSM3823944_diabetes.s3.dgecounts.rds')
con6 <- con6[["readcount"]][["exon"]][["all"]]
sobj6 <- CreateSeuratObject(counts = con6)
dim(sobj6)
sobj6 <- subset(sobj6, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 )
dim(sobj6)



#=======================================================
#
#=======================================================


sobj <- merge(sobj1, y = c(sobj2,  sobj3,  sobj4,  sobj5,  sobj6),
              add.cell.ids = c("Co1", "Co2", "Co3", "DN1", "DN2", "DN3"),
              project = "DKD")
sobj



dim(sobj)

expr_data <- as.data.frame(sobj$RNA@data)

expr_data[1:5,1:5]

#=======================================================
#
#=======================================================




setwd("D:\\SCIwork\\F49DKD\\s0single\\00analysis")

load("gtf_df.Rda")

expr_data$gene_id <- rownames(expr_data)

expr_data <- merge(gtf_df, expr_data,  by='gene_id')

expr_data$gene_id <- NULL
expr_data$gene_biotype <- NULL
expr_data <- expr_data[which(!duplicated(expr_data$gene_name)),]
rownames(expr_data)  <- expr_data$gene_name
expr_data$gene_name<- NULL
expr_data[1:5,1:5]


sobj <- CreateSeuratObject(counts = expr_data)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")

#################################################################################
#step4 QC stats
#################################################################################

# Visualize QC metrics as a violin plot
VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pdf(file = 'pre1.pdf', height = 5, width = 10)
plot1+plot2
dev.off()



#################################################################################
#step5 remove low quality cells
#################################################################################


dim(sobj)
sobj <- subset(sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 10)
dim(sobj)


# Visualize QC metrics as a violin plot
plot12 <- VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pdf(file = 'pre2.pdf', height = 5, width = 10)
plot12
dev.off()



#################################################################################
#step5 normalize the data matirx
#################################################################################

hist(colSums(sobj$RNA@data),
     breaks = 100,
     main = "Total expression before normalisation",
     xlab = "Sum of expression")



pdf(file = 'pre3.pdf', height = 5, width = 10)
hist(colSums(sobj$RNA@data),
     breaks = 100,
     main = "Total expression before normalisation",
     xlab = "Sum of expression")
dev.off()



#################################################################################
#step4 QC stats
#################################################################################

# Visualize QC metrics as a violin plot
VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pdf(file = 'pre1.pdf', height = 5, width = 10)
plot1+plot2
dev.off()



#################################################################################
#step5 remove low quality cells
#################################################################################


dim(sobj)
sobj <- subset(sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)
dim(sobj)


# Visualize QC metrics as a violin plot
plot12 <- VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pdf(file = 'pre2.pdf', height = 5, width = 10)
plot12
dev.off()



#################################################################################
#step5 normalize the data matirx
#################################################################################

hist(colSums(sobj$RNA@data),
     breaks = 100,
     main = "Total expression before normalisation",
     xlab = "Sum of expression")



pdf(file = 'pre3.pdf', height = 5, width = 10)
hist(colSums(sobj$RNA@data),
     breaks = 100,
     main = "Total expression before normalisation",
     xlab = "Sum of expression")
dev.off()



sobj <- NormalizeData(sobj,
                      normalization.method = "LogNormalize",
                      scale.factor = 10000)



hist(colSums(sobj$RNA@data),
     breaks = 100,
     main = "Total expression after normalisation",
     xlab = "Sum of expression")  


pdf(file = 'pre4.pdf', height = 5, width = 10)
hist(colSums(sobj$RNA@data),
     breaks = 100,
     main = "Total expression before normalisation",
     xlab = "Sum of expression")
dev.off()


#################################################################################
##step6 identity the cluster-related genes
#################################################################################

sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(sobj), 10) #head(pbmc$RNA@var.features,10)

plot1 <- VariableFeaturePlot(sobj)

plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

pdf(file = 'pre5.pdf', height = 4, width = 8)
plot1 + plot2
dev.off()


#################################################################################
##step7 ????????????  ????????????????????????,ScaleData()??????????????????????????????
#???????????????????????????0,?????????1???
#???????????????pbmc[["RNA"]]@scale.data
#################################################################################

all.genes <- rownames(sobj)
sobj <- ScaleData(sobj, features = all.genes)
#pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")


#################################################################################
##step8 PCA ?????????
#################################################################################

sobj <- RunPCA(sobj,   npcs = 50,   verbose = TRUE,  seed.use = 1234,
               features = VariableFeatures(object = sobj))

print(sobj[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(sobj, dims = 1:2, reduction = "pca")
DimPlot(sobj, reduction = "pca",split.by = 'ident')
DimHeatmap(sobj, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(sobj, dims = 1:15, cells = 500, balanced = TRUE)


#################################################################################
##step9 PCA ?????????????????????????????????????????? Determine the 'dimensionality' of the dataset
#################################################################################
#To overcome the extensive technical noise in any single feature for scRNA-seq data,
#Seurat clusters cells based on their PCA scores, with each PC essentially representing a
#'metafeature' that combines information across a correlated feature set.
#The top principal components therefore represent a robust compression of the dataset.
#However, how many components should we choose to include? 10? 20? 100?
#In Macosko et al, we implemented a resampling test inspired by the JackStraw procedure.
#We randomly permute a subset of the data (1% by default) and rerun PCA, constructing a 'null distribution' of feature scores,
#and repeat this procedure. We identify 'significant' PCs as those who have a strong enrichment of low p-value
# ??????????????????????????????????????????
sobj <- JackStraw(sobj, num.replicate = 100)
sobj <- ScoreJackStraw(sobj, dims = 1:20)

pdf(file = 'jackstraw.pdf', height = 5, width = 5)
JackStrawPlot(sobj, dims = 1:20)
dev.off()

pdf(file = 'elbow.pdf', height = 5, width = 5)
ElbowPlot(sobj)
dev.off()




#################################################################################
#step10 ??????
#################################################################################
#
# remove.packages("Matrix")
# packageurl <- "https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.3-2.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")
#

#?????????17???PCA??????????????????
sobj <- FindNeighbors(sobj, dims = 1:14)
#resolution?????????,??????????????????cluster????????????,????????????0.4???1.2?????????
sobj <- FindClusters(sobj, resolution = 0.4)
# Look at cluster IDs of the first 5 cells
head(Idents(sobj), 5)
cluster_assig <- as.data.frame(Idents(sobj))
head(cluster_assig)
names(cluster_assig ) <- 'cluster'
table(cluster_assig$cluster)


#################################################################################
##step11 ??????????????????????????????
#################################################################################
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
sobj <- RunUMAP(sobj, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
sobj <- RunTSNE(sobj, dims = 1:10)


library(ggsci)
plot1<-DimPlot(sobj, reduction = "umap",label = TRUE)
plot2<-DimPlot(sobj, reduction = "tsne",label = TRUE)
plot1 + plot2


pdf(file = 'cluster.pdf', height = 5, width = 10)
plot1 + plot2
dev.off()



save(sobj, file = "sobj.Rdata")
head(sobj@reductions$tsne@cell.embeddings)
load("sobj.Rdata")



#################################################################################
#step12 ???????????????????????? (cluster biomarkers)
#################################################################################

# find all markers of cluster 1
cluster1.markers <- FindMarkers(sobj, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers)

# find markers for every cluster compared to all remaining cells, report only the positive ones
sobj.markers <- FindAllMarkers(sobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sobj.markers %>% group_by(cluster) %>% top_n(n = 100)

cluster1.markers <- FindMarkers(sobj, ident.1 = 0,
                                logfc.threshold = 0.25,
                                test.use = "roc",
                                only.pos = TRUE)
head(cluster1.markers)

top10 <- sobj.markers %>% group_by(cluster) %>% top_n(n = 100)

write.csv(top10, file = 'top10.csv')

write.table(top10, file = 'top10.txt')

