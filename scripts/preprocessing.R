library(Seurat)
library(dplyr)
library(ggplot2)

# set the working directory
wd <- "/home/alperen/PycharmProjects/single_cell_assesment/"
setwd(wd)

# define the output directory
output.dir <- paste0(wd, "output/")

# read the data
adata <- ReadMtx(mtx = "resource/GSE266577_counts_raw.mtx.gz",
                 cells = "resource/GSE266577_barcodes.txt.gz",
                 features = "resource/GSE266577_seurat_features.txt.gz",
                 feature.column = 1)

# create the seurat object
seurat.obj <- CreateSeuratObject(counts = adata,
                                 project = "alperen", 
                                 min.cells = 500,
                                 min.features = 1000)

# Add percentage of mitochondrial genes and total counts metadata
seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "^MT-")

# Visualize QC metrics as violin plots
violin.qc <- VlnPlot(seurat.obj,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                  ncol = 3,
                  raster = F)

print(violin.qc)

# Scatter plots for QC metrics
scatter.qc <- FeatureScatter(seurat.obj,
                             feature1 = "nCount_RNA",
                             feature2 = "percent.mt",
                             raster = F) +
  FeatureScatter(seurat.obj,
                 feature1 = "nCount_RNA",
                 feature2 = "nFeature_RNA",
                 raster = F)

print(scatter.qc)

# filter the data with the values obtained with plots
seurat.obj <- subset(seurat.obj,
                     subset = nFeature_RNA > 200 &
                       nFeature_RNA < 7000 & 
                       nCount_RNA > 1000 &
                       nCount_RNA < 50000 & 
                       percent.mt < 5)


# normalize the data
seurat.obj <- NormalizeData(seurat.obj,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000)
seurat.obj <- FindVariableFeatures(seurat.obj,
                                   selection.method = "vst",
                                   nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat.obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat.obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# PCA analysis
# scale the data for PCA and UMAP purposes
all.genes <- rownames(seurat.obj)
seurat.obj <- ScaleData(seurat.obj, features = all.genes)

seurat.obj <- RunPCA(seurat.obj, 
                     features = VariableFeatures(
                       object = seurat.obj)
                     )

VizDimLoadings(object = seurat.obj)
DimPlot(seurat.obj, reduction = "pca")
DimHeatmap(seurat.obj)
ElbowPlot(seurat.obj)

seurat.obj <- FindNeighbors(seurat.obj, dims = 1:17)
seurat.obj <- FindClusters(seurat.obj, resolution = 0.5)
seurat.obj <- RunUMAP(seurat.obj, dims = 1:19, )
DimPlot(seurat.obj, reduction = "umap")

# save the data to load it later 
saveRDS(seurat.obj, file = "output/assesment_singlecell.rds")


# find markers for every cluster compared to all remaining cells, report only the positive
# ones
all.markers <- FindAllMarkers(seurat.obj, only.pos = TRUE, logfc.threshold = 3)
all.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 3)

marker_list <- list()
idents <- unique(seurat.obj@meta.data$RNA_snn_res.0.5)
for (i in idents) {
  name <- paste0("cluster", i)
  marker_list[[name]] <- FindMarkers(seurat.obj, 
                                     ident.1 = i, logfc.threshold = 3, test.use = "roc", only.pos = TRUE)
}
