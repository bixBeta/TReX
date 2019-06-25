install.packages("Seurat")
library(Seurat)
library(dplyr)

myPheno <- read.table(file = "file_paths_scRNA.txt", sep = "\t", header = T)

####################################################################################
####################################################################################
# read in filtered matrices for all samples 
for (i in 1:nrow(myPheno)) {
      assign(x = paste0("RD10x_",myPheno$sample_id[i]), 
      Read10X(data.dir = as.character(myPheno$filtered_matrix[i])), 
      envir = .GlobalEnv)
}

####################################################################################
####################################################################################
# create seurat objects from imported matrices
ten_x_pattern <- objects(pattern = "RD10x_")

for (i in 1:length(ten_x_pattern)) {
      assign(x = paste0("SObj_", strsplit(ten_x_pattern[i], "_")[[1]][2]), 
      value = CreateSeuratObject(counts = get(ten_x_pattern[i]), 
      project = "MECFS_PBMC", min.cells = 3, min.features = 200),
      envir = .GlobalEnv)
}

####################################################################################
####################################################################################
# merge all seurat objects into one main object
# n <- length(sobj_pattern)
# n2 <- seq(0,n,2)
# 
# z <- capture.output( for (i in 2:length(sobj_pattern)) {
#   cat(sobj_pattern[i], ",")
# })

sobj_pattern <- noquote(objects(pattern = "SObj_"))

# test <- merge(SObj_1744, 
#                   y = mget(objects(pattern = "SObj_")) , project = "MECFS_PBMC")

all_pbmcs <- merge(get(sobj_pattern[1]), 
                   y = mget(sobj_pattern[2:length(sobj_pattern)]) , project = "MECFS_PBMC")

