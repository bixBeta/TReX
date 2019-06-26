library(Seurat)
load("MECFS_merged_seurat_object.RData")

####################################################################################
####################################################################################
# store mitochondrial percentage in object meta data
all_pbmcs <- PercentageFeatureSet(all_pbmcs, pattern = "^MT-", col.name = "percent.mt")

head(all_pbmcs@meta.data)

# Visualize QC metrics as a violin plot
VlnPlot(all_pbmcs, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(all_pbmcs, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all_pbmcs, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

####################################################################################
####################################################################################
# filter test 1
pbmc <- subset(all_pbmcs, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 20)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)


####################################################################################
####################################################################################
# dplyr mutate to change typeof library_id
batchInfo <- 
  myPheno %>%
  select(1,4,5)

batchInfo <- 
  batchInfo %>%
  mutate_each(funs = (as.character), library_id)

head(pbmc@meta.data)
tail(pbmc@meta.data)
localMeta <- left_join(x = pbmc@meta.data, y = batchInfo, by = c("orig.ident" = "library_id"))

# add pool and batch info
pbmc <- AddMetaData(object = pbmc, metadata = localMeta$pool, col.name = "pool")
pbmc <- AddMetaData(object = pbmc, metadata = localMeta$batch, col.name = "batch")

