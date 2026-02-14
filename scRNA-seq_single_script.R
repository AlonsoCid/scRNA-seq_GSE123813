wd <- getwd()

# Libraries
library(readxl)
library(tidyverse)
library(GEOquery)
library(knitr)
library(kableExtra)
library(Seurat)
library(SingleR)
library(SingleCellExperiment)
library(scrapper)
library(celldex)

# 1. Sample data
## Get the data from GEO
gseid="GSE123813"
GEO <- getGEOSuppFiles(gseid, fetch_files = F)
getGEOSuppFiles(gseid, makeDirectory = T, filter_regex="_scc_metadata")
getGEOSuppFiles(gseid, makeDirectory = T, filter_regex="_scc_scRNA_counts")

## Create a Seurat object
metadata <- data.table::fread(file.path(wd,gseid,"GSE123813_scc_metadata.txt.gz"))
expression_data <- read.table(file.path(wd,gseid,"GSE123813_scc_scRNA_counts.txt.gz"),
                              header=TRUE, row.names=1, sep="\t", check.names=FALSE)
## Reorder expression matrix
expression_data <- expression_data[,metadata$cell.id]

## length(intersect(metadata$cell.id,colnames(expression_data)))
all.equal(metadata$cell.id,colnames(expression_data))

## The same but not in the same order
expression_matrix <- as.matrix(expression_data)
scc <- CreateSeuratObject(counts = expression_matrix, project = "GSE123813")
scc <- AddMetaData(object = scc, metadata = metadata)

scc@meta.data$cluster <- metadata$cluster
scc@meta.data$patient <- metadata$patient
scc@meta.data$treatment <- metadata$treatment

save(scc,file="GSE123813.raw.seurat.RData")

# 2. The Seurat object
## Basic data exploration
load(file="GSE123813.raw.seurat.RData")
scc
scc@assays
head(scc@meta.data)
table(scc$patient, scc$treatment)

# 4. QC
## Identify mitochondrial genes
grep("MT-",rownames(scc), value=T)

scc[["percent.mt"]] <- PercentageFeatureSet(scc, pattern = "^MT-")
VlnPlot(scc, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
VlnPlot(scc, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "treatment" )

plot2 <- FeatureScatter(scc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2

# 5. Preprocessing

scc <- NormalizeData(scc, normalization.method = "LogNormalize", scale.factor = 10000)
scc <- FindVariableFeatures(scc, selection.method = "vst", nfeatures = 2000)
scc <- ScaleData(scc)
scc <- RunPCA(scc, features = VariableFeatures(object = scc))

ElbowPlot(scc)

# 6. Clustering
scc <- FindNeighbors(scc, dims = 1:16)
scc <- FindClusters(scc, resolution = 0.3, ident="cluster") #resolution=0.3, reported by authors
scc <- RunUMAP(scc, n.neighbors=20L, min.dist=0.3, dims = 1:16)
DimPlot(scc, reduction = "umap", group.by="cluster")
## By patient
DimPlot(scc, reduction = "umap", split.by="patient", group.by="cluster")
## By treatment
DimPlot(scc, reduction = "umap", split.by="treatment", group.by="cluster")

# 7. Cluster annotation
ref = celldex::ImmGenData()

pred_clust <- SingleR(test=as.SingleCellExperiment(scc), ref=ref,
                      clusters = scc$seurat_clusters, labels=ref$label.fine)
table(pred_clust$labels)

scc$immgen_clust <- scc$seurat_clusters
levels(scc$immgen_clust) <- pred_clust$labels
DimPlot(scc, reduction = "umap", group.by="immgen_clust")

# 8. Gene expression
FeaturePlot(scc, features = "BCR", split.by="treatment")

# 9. Differential Expression Analysis
clusterN.markers <- FindMarkers(scc, ident.1 = 3, ident.2 = c(0, 1), min.pct = 0.25)
head(clusterN.markers, n = 10)

markers.tto <- FindMarkers(scc, ident.1 = "post", ident.2 ="pre", group.by = 'treatment')
head(markers.tto, n=10)

FeaturePlot(scc, features = "IGKC", split.by="treatment")
