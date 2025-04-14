setwd( )

source("ks_scAverExp.R")
source("ks_scDotplot.R")
source("dimplot_new.R")

## Original data loading
library(Seurat)

# Read expression matrix
scRNA_data <- fread(file="GSE149614_HCC.scRNAseq.S71915.count.txt",
                    quote="", header = "auto", sep='\t')
scRNA_data <- column_to_rownames(scRNA_data, var = "V1")
scRNA <- CreateSeuratObject(scRNA_data)

# Load metadata
metadata <- read.delim("GSE149614_HCC.metadata.updated.txt", row.names=1)

# Add metadata to Seurat object
scRNA <- AddMetaData(scRNA, metadata)
head(scRNA@meta.data)

## Seurat Pipeline
nfeatures = 2000
ndim = 15
neigh = 20
dist = 0.5
res = 0.6

# Subset tumor samples
scRNA_tumor <- subset(scRNA, subset = site %in% c("Tumor","Lymph","PVTT"))

# Normalize data
scRNA_tumor <- NormalizeData(scRNA_tumor, scale.factor = 10000,
                             normalization.method = "LogNormalize")

# Identify highly variable features
scRNA_tumor <- FindVariableFeatures(scRNA_tumor, nfeatures = nfeatures, 
                                    selection.method = "vst")

# Scale data
scRNA_tumor <- ScaleData(scRNA_tumor, features = VariableFeatures(scRNA_tumor))

# Run PCA
scRNA_tumor <- RunPCA(scRNA_tumor, assay = 'RNA', slot = 'scale.data')

# Run Harmony for batch effect correction
scRNA_tumor <- RunHarmony(scRNA_tumor, group.by.vars = "sample", dims.use = 1:50,
                          assay.use = "RNA")

# Find neighbors
scRNA_tumor <- FindNeighbors(scRNA_tumor, k.param = neigh,
                             dims = 1:ndim, reduction = "harmony")

# Find clusters
scRNA_tumor <- FindClusters(scRNA_tumor, resolution = 0.4, n.iter = 50)

# Run UMAP
scRNA_tumor <- RunUMAP(scRNA_tumor, dims = 1:ndim,
                       n.neighbors = neigh, min.dist = dist, 
                       reduction = "harmony", reduction.name = "umap_harmony")

# Run tSNE
scRNA_tumor <- RunTSNE(scRNA_tumor, dims = 1:ndim,
                       n.neighbors = neigh, min.dist = dist, 
                       reduction = "harmony", reduction.name = "tsne_harmony")

## Quick QC with FeaturePlot
my_colors <- c("#F0F0F0",'#EDD1D8', '#f4a3a8', '#e38191', '#cc607d', '#ad466c', '#8b3058', '#672044')
FeaturePlot(scRNA_tumor, features = "ALB", reduction = "tsne_harmony") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(colour = "black", size = 15),
        plot.title = element_text(size = 17, hjust = 0.5),
        panel.border = element_blank()) + 
  scale_color_gradientn(colors = my_colors) + 
  labs(x = ' ', y = ' ', title = "ALB")

## Visualize by sample or cluster

plot <- dimplot_new(scRNA_tumor, reduction = "tsne_harmony", pt.size = 0.05, label = F, group.by = c("site"))
plot

plot <- dimplot_new(scRNA_tumor, reduction = "tsne_harmony", pt.size = 0.05, label = F, group.by = c("seurat_clusters"))
plot

## Marker gene identification
scRNA.markers <- FindAllMarkers(scRNA_tumor, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10_markers <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10_table = unstack(top10_markers, gene ~ cluster)
names(top10_table) = gsub("X","cluster",names(top10_table))

## Cell type annotation
scRNA_tumor@meta.data <- scRNA_tumor@meta.data %>%
  mutate(celltype_major_steven = recode(seurat_clusters,
                                        `0` = "c01_Malignant", `1` = "c02_T/NK_cells", `2` = "c03_Macrophages",
                                        `3` = "c01_Malignant", `4` = "c04_Monocytes", `5` = "c01_Malignant",
                                        `6` = "c05_Cycling_Malignant", `7` = "c06_Endothelial", `8` = "c07_Fibroblasts",
                                        `9` = "c08_Plasma_cells", `10` = "c01_Malignant", `11` = "c09_Cycling_T/NK_cells",
                                        `12` = "c10_B_cells", `13` = "c11_Cycling_Myeloid_cells", `14` = "c01_Malignant",
                                        `15` = "c12_Pericytes", `16` = "c01_Malignant", `17` = "c01_Malignant"
  ))

# Re-plot cluster or cell type
plot <- dimplot_new(scRNA_tumor, reduction = "tsne_harmony", pt.size = 0.01, label = F, group.by = c("celltype_major_steven"))
plot
plot <- dimplot_new(scRNA_tumor, reduction = "tsne_harmony", pt.size = 0.05, label = F, group.by = c("sample"))
plot
plot <- dimplot_new(scRNA_tumor, reduction = "tsne_harmony", pt.size = 0.05, label = F, group.by = c("site"))
plot

## Plot selected marker genes
select <- c("EPCAM","ALB","CD3D","KLRB1","CD8A","CD4","SPP1","C1QA","LYZ","S100A9","MKI67","PECAM1","COL1A1","IGHG1","MS4A1","RGS5")
for (i in select) {
  setwd("/scRNA-seq")
  fig2 = FeaturePlot(scRNA_tumor, features = i, reduction = "tsne_harmony") + 
    theme_bw() +
    theme(panel.grid = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
          axis.title = element_text(colour = "black", size = 15), plot.title = element_text(size = 17, hjust = 0.5),
          panel.border = element_blank(), legend.position = "right") + 
    scale_color_gradientn(colors = my_colors) + 
    labs(x = ' ', y = ' ', title = i)
  ggsave(paste0(i, "tsne.png"), fig2, height = 3.5, width = 4.7)
}

# Average expression and dotplot
featureSets <- list(
  "c01_Malignant" = c("EPCAM","KRT7"),
  "c02_T/NK_cells" = c("NKG7","KLRB1","CD3D"),
  "c03_Macrophages" = c("C1QA","C1QB","CD68"),
  "c04_Monocytes" = c("S100A8", "S100A9","IL1B"),
  'c05_Cycling_Malignant' = c("MKI67"),
  "c06_Endothelial" = c("PECAM1", "VWF1"),
  'c07_Fibroblasts' = c("COL1A1","COL1A2","ACTA2"),
  "c08_Plasma_cells" = c("IGHG1","IGHA1","JCHAIN","SDC1"),
  "c09_Cycling_T/NK_cells" = c("GZMA","TOP2A"),
  "c10_B_cells" = c("MS4A1","CD79A","CD27"),
  "c11_Cycling_Myeloid_cells" = c("FCER1G"),
  "c12_Pericytes" = c("RGS5","CLDN5")
)

scRNA_tumor$celltype <- scRNA_tumor$celltype_major_steven
Exp_scRNA <- ks_scAverExp(obj = scRNA_tumor, features = featureSets, CellTypes = unique(scRNA_tumor$celltype_major_steven))

library(dittoSeq)
cols <- c("c01_Malignant" = "#A6D719", "c02_T/NK_cells" = "#176EBF", "c03_Macrophages" = "#00A8DE",
          "c04_Monocytes" = "#AEE0E8", 'c05_Cycling_Malignant' = "#00A9A3", "c06_Endothelial" = "#FBD324",
          'c07_Fibroblasts' = "#F28A24", "c08_Plasma_cells" = "#A52828", "c09_Cycling_T/NK_cells" = "#A37CB7",
          "c10_B_cells" = "#F2D7EE", "c11_Cycling_Myeloid_cells" = "#CD6981", "c12_Pericytes" = "#FBD0C0")

ss = ks_scDotplot(obj = scRNA_tumor, features = featureSets, CellTypes = unique(scRNA_tumor$celltype_major_steven),
                  avgPctMat = Exp_scRNA, pal = cols)
ss
