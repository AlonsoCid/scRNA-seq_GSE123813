library(Seurat)

scc <- readRDS(snakemake@input[[1]])

# Clustering
scc <- FindNeighbors(scc, dims = 1:16)
scc <- FindClusters(scc, resolution = 0.3) 
scc <- RunUMAP(scc, n.neighbors = 20L, min.dist = 0.3, dims = 1:16)

png(snakemake@output[["elbow_png"]], width = 600, height = 600)
ElbowPlot(scc)
dev.off()

png(snakemake@output[["umap_png"]], width = 1000, height = 600)
DimPlot(scc, reduction = "umap", split.by = "treatment", group.by = "cluster")
dev.off()

# Visualization
pdf(snakemake@output[["plots"]])
print(DimPlot(scc, reduction = "umap", group.by = "cluster"))
print(DimPlot(scc, reduction = "umap", split.by = "patient", group.by = "cluster"))
print(DimPlot(scc, reduction = "umap", split.by = "treatment", group.by = "cluster"))
dev.off()

saveRDS(scc, file = snakemake@output[["rds"]])