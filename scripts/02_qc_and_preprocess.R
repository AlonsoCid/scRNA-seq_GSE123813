library(Seurat)

scc <- readRDS(snakemake@input[[1]])

# QC
scc[["percent.mt"]] <- PercentageFeatureSet(scc, pattern = "^MT-")

pdf(snakemake@output[["plot"]])
print(VlnPlot(scc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt")))
print(FeatureScatter(scc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
dev.off()

png(snakemake@output[["vln_png"]], width = 800, height = 600)
VlnPlot(scc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "treatment")
dev.off()

png(snakemake@output[["scatter_png"]], width = 600, height = 600)
FeatureScatter(scc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

# Preprocessing
scc <- NormalizeData(scc, normalization.method = "LogNormalize", scale.factor = 10000)
scc <- FindVariableFeatures(scc, selection.method = "vst", nfeatures = 2000)
scc <- ScaleData(scc)
scc <- RunPCA(scc, features = VariableFeatures(object = scc))

saveRDS(scc, file = snakemake@output[["rds"]])