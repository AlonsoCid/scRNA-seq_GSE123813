library(Seurat)
library(GEOquery)
library(data.table)

gseid <- "GSE123813"
# Create directories if they don't exist
dir.create("data/external", recursive = TRUE, showWarnings = FALSE)

# Fetch files
getGEOSuppFiles(gseid, makeDirectory = TRUE, filter_regex = "_scc_metadata")
getGEOSuppFiles(gseid, makeDirectory = TRUE, filter_regex = "_scc_scRNA_counts")

# Load data
metadata <- fread(file.path(gseid, "GSE123813_scc_metadata.txt.gz"))
expression_data <- read.table(file.path(gseid, "GSE123813_scc_scRNA_counts.txt.gz"),
                              header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# Reorder and Create Seurat Object
expression_data <- expression_data[, metadata$cell.id]
scc <- CreateSeuratObject(counts = as.matrix(expression_data), project = "GSE123813")
scc <- AddMetaData(object = scc, metadata = metadata)

# Assign specific metadata
scc$cluster <- metadata$cluster
scc$patient <- metadata$patient
scc$treatment <- metadata$treatment

saveRDS(scc, file = snakemake@output[[1]])