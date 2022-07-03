#=======================================================
#set the working files and load the packages
#=======================================================

rm(list=ls())
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
library(SingleR)
library(celldex)


setwd("D:\\SCIwork\\F49DKD\\s0single\\00analysis")


load("sobj.Rdata")

dim(sobj)


#################################################################################


#################################################################################



table(sobj$seurat_clusters)

new.cluster.ids <- c("unknown",
                     "unknown", 
                     "unknown", 
                     "unknown",
                     "unknown",
                     "unknown",
                     "unknown", 
                     "unknown", 
                     "unknown",
                     "immune",
                     "unknown",
                     "unknown")

names(new.cluster.ids) <- levels(sobj)

sobj <- RenameIdents(sobj, new.cluster.ids)

table(Idents(sobj))


cd_genes <- c('CD3G', 'CD8A', 'FCGR3A', 'NCAM1')
DotPlot(object =sobj, features = cd_genes)

sobj@meta.data$subtype <-  Idents(sobj)

sobj@meta.data$subtype

pdf(file = 'immune.pdf', height = 6, width = 6)
DotPlot(object =sobj, features = cd_genes)
dev.off()



#################################################################################


#################################################################################




sobj_immune = sobj[,Idents(sobj)=="immune"]

sobj_immune@meta.data$group <- substr(x=sobj_immune@meta.data$orig.ident, start = 1, stop = 2)

table(sobj_immune$group)

data_group <- as.data.frame(sobj_immune@meta.data)

write.csv(data_group, file = "data_group.csv")


pdf(file = "s1.pdf", heigh=5, width=5)
DimPlot(sobj_immune, reduction = "tsne", 
        group.by = 'group', label = TRUE, 
        pt.size = 0.5) 
dev.off()


