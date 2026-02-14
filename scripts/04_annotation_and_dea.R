library(Seurat)
library(SingleR)
library(celldex)
library(SingleCellExperiment)

scc <- readRDS(snakemake@input[[1]])
ref <- celldex::ImmGenData()

# Annotation
pred_clust <- SingleR(test = as.SingleCellExperiment(scc), ref = ref,
                      clusters = scc$seurat_clusters, labels = ref$label.fine)

scc$immgen_clust <- scc$seurat_clusters
levels(scc$immgen_clust) <- pred_clust$labels

# Save Plots
pdf(snakemake@output[["plots"]])
print(DimPlot(scc, reduction = "umap", group.by = "immgen_clust"))
print(FeaturePlot(scc, features = "BCR", split.by = "treatment"))
print(FeaturePlot(scc, features = "IGKC", split.by = "treatment"))
dev.off()

# DEA
clusterN.markers <- FindMarkers(scc, ident.1 = 3, ident.2 = c(0, 1), min.pct = 0.25)
write.csv(clusterN.markers, file = snakemake@output[["markers"]])

markers.tto <- FindMarkers(scc, ident.1 = "post", ident.2 = "pre", group.by = 'treatment')
write.csv(markers.tto, file = snakemake@output[["dea"]])

png(snakemake@output[["immgen_png"]], width = 800, height = 600)
DimPlot(scc, reduction = "umap", group.by = "immgen_clust")
dev.off()

png(snakemake@output[["bcr_png"]], width = 1000, height = 600)
FeaturePlot(scc, features = "BCR", split.by = "treatment")
dev.off()

png(snakemake@output[["igkc_png"]], width = 1000, height = 600)
FeaturePlot(scc, features = "IGKC", split.by = "treatment")
dev.off()

saveRDS(scc, file = snakemake@output[["rds"]])