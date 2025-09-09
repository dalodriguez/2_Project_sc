library(Seurat)

data_dir <- "data/SCO-42-2_S21_analysis_ch10/outs/filtered_feature_bc_matrix"

expression_matrix <- Read10X(data.dir = data_dir)
seurat <- CreateSeuratObject(counts = expression_matrix)

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

seurat <- subset(seurat, subset =  nFeature_RNA < 8000 & percent.mt < 15)

seurat <- NormalizeData(seurat, normalization.method = "LogNormalize")

seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)

seurat <- ScaleData(seurat)

seurat <- RunPCA(seurat, features = VariableFeatures(seurat))

seurat <- FindNeighbors(seurat, dims = 1:15)

seurat <- FindClusters(seurat)

seurat <- RunUMAP(seurat, dims = 1:15)

saveRDS(seurat, file = "data/SCO_analysis_seurat.rds")