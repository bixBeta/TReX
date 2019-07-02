# ####################################################################################
# ####################################################################################
# # Normalization
# system("date")
# pbmc <- SCTransform(pbmc, vars.to.regress = c("percent.mt", "batch", "nCount_RNA"), verbose = T)
# system("date")
# 
# saveRDS(pbmc , "filter-2-normalized_scTransformed_pbmc.rds")
# 
# pbmc <- RunPCA(pbmc, verbose = F)
# pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = T)
# pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
# pbmc <- FindClusters(pbmc, verbose = FALSE)
# 
# DimPlot(pbmc, label = TRUE, group.by = "batch") 
# pbmc <- RunTSNE(pbmc, dims = 1:30, nthreads = 48, max_iter = 2000)
# DimPlot(pbmc, label = TRUE, group.by = "batch", reduction = "tsne") 
# 

####################################################################################
####################################################################################
# use CCA w/ Seurat::FindIntegrationAnchors()

# update meta on all_pbmcs

meta.all.pbmcs <- all_pbmcs@meta.data

newPheno <- 
  newPheno %>%
  mutate_each(funs = (as.character), library_id)

newPheno <- left_join(x = meta.all.pbmcs, y = newPheno, by = c("orig.ident" = "library_id"))
nrow(meta.all.pbmcs) == nrow(newPheno)

newPheno$pool[newPheno$pool=="no"]<- "N/A"

all_pbmcs <- AddMetaData(object = all_pbmcs, metadata = newPheno$pool, col.name = "pool")
all_pbmcs <- AddMetaData(object = all_pbmcs, metadata = newPheno$batch, col.name = "batch")
all_pbmcs <- AddMetaData(object = all_pbmcs, metadata = newPheno$ENID.., col.name = "ENID")
all_pbmcs <- AddMetaData(object = all_pbmcs, metadata = newPheno$CPET.Day, col.name = "CPET.Day")
all_pbmcs <- AddMetaData(object = all_pbmcs, metadata = newPheno$Specimen.ID, col.name = "Specimen.ID")
all_pbmcs <- AddMetaData(object = all_pbmcs, metadata = newPheno$Analysis.ID, col.name = "Analysis.ID")
all_pbmcs <- AddMetaData(object = all_pbmcs, metadata = newPheno$Phenotype, col.name = "Phenotype")

# split objects by batch 
all.pbmcs.list <- SplitObject(all_pbmcs, split.by = "batch")

for (i in 1:length(all.pbmcs.list)) {
  all.pbmcs.list[[i]] <- NormalizeData(all.pbmcs.list[[i]], verbose = FALSE)
  all.pbmcs.list[[i]] <- FindVariableFeatures(all.pbmcs.list[[i]], selection.method = "vst", nfeatures = 2000, 
                                             verbose = FALSE)
}


reference.list <- all.pbmcs.list[c("batch1", "batch2", "batch3", "batch4", "batch5", "batch6", "batch7")]
all.pbmcs.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)

integrated.pbmcs <- IntegrateData(anchorset = all.pbmcs.anchors, dims = 1:30)

saveRDS(all.pbmcs.anchors, "all.pbmc.anchors.rds")
saveRDS(integrated.pbmcs, "integrated.pbmcs.rds")


####################################################################################
####################################################################################
# switch to integrated assay. The variable features of this assay 
# are automatically set during
# IntegrateData
DefaultAssay(integrated.pbmcs) <- "integrated"

# Run the standard workflow for visualization and clustering
integrated.pbmcs <- ScaleData(integrated.pbmcs, verbose = FALSE)
integrated.pbmcs <- RunPCA(integrated.pbmcs, npcs = 30, verbose = FALSE)
integrated.pbmcs <- RunUMAP(integrated.pbmcs, reduction = "pca", dims = 1:30)

DimPlot(integrated.pbmcs, label = F, group.by = "batch") 
DimPlot(integrated.pbmcs, label = F, group.by = "ENID") 

integrated.pbmcs <- RunTSNE(integrated.pbmcs, dims = 1:30, nthreads = 48, max_iter = 2000)
DimPlot(pbmc, label = TRUE, group.by = "batch", reduction = "tsne") 


table(integrated.pbmcs@meta.data)

umIdent <- StashIdent(integrated.pbmcs, save.name = "UMAP.Clusters")
table(umIdent@meta.data$UMAP.Clusters, umIdent@meta.data$orig.ident)

# or
## integrated.pbmcs[["UMAP.Clusters.Id"]] <- Idents(object = integrated.pbmcs)

####################################################################################
####################################################################################
# clustering
integrated.pbmcs <- FindNeighbors(integrated.pbmcs, dims = 1:30)
integrated.pbmcs <- FindClusters(integrated.pbmcs, resolution = 0.5)

DimPlot(integrated.pbmcs, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + 
  NoLegend()
DimPlot(integrated.pbmcs, reduction = "umap", group.by = "Specimen.ID", label = F, repel = TRUE)

saveRDS(integrated.pbmcs, "integrated.pbmcs.cluster.0.5.res.rds")

# test <- FindClusters(integrated.pbmcs, resolution = 0.9)
# 
# DimPlot(test, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + 
#   NoLegend()

BuildClusterTree(test)

test_DGE <- FindMarkers(test, ident.1 = "1", ident.2 = "2")
test_DGE <- FindMarkers(test, ident.1 = "batch3", group.by = 'batch', subset.ident = "1")


ms <- FindAllMarkers(test, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ms %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)



DimPlot(integrated.pbmcs, reduction = "umap", split.by = "batch")
