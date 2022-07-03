#=======================================================
#set the working files and load the packages
#=======================================================

library(GEOquery)
rm(list=ls())
library(dplyr)
library(tidyr)
library(Biobase)
library(limma)
library(gdata)
library(limma)
library(edgeR)



setwd("D:\\SCIwork\\F49DKD\\s1download\\GSE142025\\GSE142025")

#=======================================================
#
#=======================================================

surv <- read.csv('group.csv', header = T, row.names = 1)

table(surv$group)

surv$group[surv$group=='Advanced_DN'] <- 'AdDKD'
surv$group[surv$group=='Early_DN'] <- 'EaDKD'
surv$group[surv$group=='Control'] <- 'Normal'  

table(surv$group)

surv <- subset(surv, surv$group != 'Normal')

table(surv$group)


##########################################################################################
## 
###########################################################################################


exprset <- read.csv(file = 'exprset.csv', header = T, row.names = 1)

exprset <-  as.data.frame(t(exprset))

exprset[1:5,1:5] 

normalize <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}

exprset <- as.data.frame(apply(exprset, 2, normalize))

exprset[1:5,1:5] 

exprset <-  as.data.frame(t(exprset))

exprset[1:5,1:5] 

fpkmToTpm <- function(fpkm)
{exp(log(fpkm) - log(sum(fpkm)) + log(1e6))}

exprset <- as.data.frame (apply(exprset, 2, fpkmToTpm))

exprset[1:5,1:5] 

exprset <- as.data.frame(t(exprset))

exprset[1:5,1:5] 

exprset$sample <- rownames(exprset)

str(exprset)



##########################################################################################
## 
###########################################################################################


data <- merge(surv,exprset, by='sample')

data[1:5,1:5]

rownames(data) <- data$sample

data$sample <- NULL

group <- subset(data, select=c('group'))

head(group)

group <- group$group

data$group <- NULL

data <- as.data.frame(t(data))

data[1:5,1:5]

exprSet <- data

table(group)


setwd("D:\\SCIwork\\F49DKD\\s3EaAd")



# ===================================================


# ===================================================


design <- model.matrix(~group)
y <- DGEList(counts=exprSet,group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y,pair = c("EaDKD","AdDKD"))
topTags(et)
ordered_tags <- topTags(et, n=100000)

allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$PValue)==FALSE,]
diff=allDiff
newData=y$pseudo.counts

foldChange=1.5
padj=0.05



# ===================================================


# ===================================================



write.csv(diff, file="edger_diff.csv")

diffSig = diff[(diff$PValue < padj & (diff$logFC>foldChange | diff$logFC<(-foldChange))),]
write.csv(diffSig, file="edger_diffSig.csv")

diffUp = diff[(diff$PValue < padj & (diff$logFC>foldChange)),]
write.csv(diffUp, file="edger_up.csv")

diffDown = diff[(diff$PValue < padj & (diff$logFC<(-foldChange))),]
write.csv(diffDown, file="edger_down.csv")

dim(diffSig)
dim(diffUp)
dim(diffDown)




