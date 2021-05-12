library(dplyr)
library(tidyverse)
library(DESeq2)

target = "targetFile.txt"

phenoData <- read.table(file = target, header = T)

################################
################################
# read in all count files 

for (i in 1:length(phenoData$label)) {
  
  assign(x = paste0(phenoData$label[i], "_raw"), 
         value = read.table(as.character(phenoData$files)[i]),
         envir = .GlobalEnv )
  
}

raw_count_pattern <- ls(pattern = "_raw")

for (i in 1:length(raw_count_pattern)) {
  
  assign(x = raw_count_pattern[i], 
         value = get(raw_count_pattern[i]) %>% 
           rename(., V1 = "geneID")  %>% 
           rename(.,  V2 = "count"))  
}

################################
################################
# populate countMatrix

# rownames(countMatrix) <- get(raw_count_pattern[1])$geneID

countMatrix <- as.data.frame(mget(raw_count_pattern))

rownames(countMatrix) <- countMatrix[,1]

countMatrix <- countMatrix %>% dplyr::select(-matches("geneID"))

rm(list = ls(pattern = "_raw"))
################################
################################
# DESEQ2

library("DESeq2")
## -------------------------------------------------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = phenoData,
                              design = ~ group)
## -------------------------------------------------------------------------------------------------------------------
dds <- DESeq(dds)
resultsNames(dds)

## -------------------------------------------------------------------------------------------------------------------
vsd <- varianceStabilizingTransformation(dds, blind=T)
colnames(vsd) <- gsub("_raw.count","",colnames(vsd))

## -------------------------------------------------------------------------------------------------------------------
rv <- rowVars(assay(dds))
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(dds)[select,]))

# the contribution to the total variance for each component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
cond <- phenoData$group
intgroup <- "group"
intgroup.df <- as.data.frame(colData(dds)[, intgroup, drop=FALSE])
rownames(intgroup.df) <- phenoData$label

# assembly the data for the plot
d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=cond, intgroup.df, name=colnames(dds))
d2 <- data.frame(PC1=pca$x[,1], PC3=pca$x[,3], group=cond, intgroup.df, name=colnames(dds))

# ggplot(data=d2, aes_string(x="PC1", y="PC3", color="group")) + geom_point(size=3) + 
#   xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
#   ylab(paste0("PC3: ",round(percentVar[2] * 100),"% variance")) +
#   coord_fixed()
## -------------------------------------------------------------------------------------------------------------------

PC1= as.data.frame(pca$x[,1])
PC2= as.data.frame(pca$x[,2])
PC3= as.data.frame(pca$x[,3])

pComp.df <- data.frame(PC1, PC2 , PC3)
rownames(pComp.df) <- phenoData$label
pComp.df <- rownames_to_column(pComp.df, var = "samples")
pComp.df <- left_join(pComp.df, phenoData, by = c("samples" = "label"))
pComp.df <- pComp.df %>% mutate_at(vars(matches("pca")), ~ ./100000)



hc <- hclust(dist(t(assay(vsd))), method="ward.D")
# plot(hc, hang=-1, ylab="Height", las=2, 
#     xlab="Method: Euclidean distance - Ward criterion", main="Cluster dendrogram")





