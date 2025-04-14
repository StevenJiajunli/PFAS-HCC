
## Prepare matrix for raincloud plot

HCC_cohorts_diffexp <- readRDS("HCC_cohorts_diffexp.rds")

data <- HCC_cohorts_diffexp[["TCGA"]][,c("ESR1","APOA1","IGF1","PPARGC1A","SERPINE1","PON1")]

# Add tissue labels
data$tissue[substr(rownames(data),14,15) == "01"] <- "tumor"
data$tissue[substr(rownames(data),14,15) != "01"] <- "normal"

library(ggpubr)

data$tissue <- as.factor(data$tissue)

# Draw raincloud plot

source("Rainclouds_plot.R")

P5 <- ggplot(data, aes(x = tissue, y = IGF1, fill = tissue)) +
  geom_flat_violin(aes(fill = tissue), position = position_nudge(x = 0.1, y = 0), trim = TRUE, alpha = .5, colour = NA) +
  geom_point(aes(x = .55, y = IGF1, colour = tissue), position = position_jitter(width = .05), size = 1, shape = 20) +
  geom_boxplot(aes(x = tissue, y = IGF1, fill = tissue), outlier.shape = NA, alpha = .5, width = .1, colour = "black") +
  scale_colour_manual(values = c("#85BBDA", "#1C60AC")) +
  scale_fill_manual(values = c("#85BBDA", "#1C60AC")) +
  theme_classic() +
  geom_signif(
    comparisons = list(c("tumor", "normal")),
    map_signif_level = TRUE,
    test = "wilcox.test",
    y_position = c(5)
  ) +
  theme(
    axis.title = element_text(size = 15),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 12)
  ) +
  ylab('IGF1 Expression')

P5

## scRNA-seq 数据 FeaturePlot / FeaturePlot for scRNA-seq
my_colors <- c("#F0F0F0",'#AED4E5', '#81B5D5', '#5795C7', '#3371B3', '#345D82', '#1E4C9C')

FeaturePlot(scRNA_tumor, features = "PON1", reduction = "tsne_harmony") + 
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_text(colour = "black", size = 15),
    plot.title = element_text(size = 17, hjust = 0.5),
    panel.border = element_blank()
  ) + 
  scale_color_gradientn(colors = my_colors) + 
  labs(x = ' ', y = ' ', title = "PON1")

## Spatial transcriptomics

## Alternative deconvolution with SpaCET
library(Seurat)
library(hdf5r)
library(SpaCET)

visiumPath <- "/HCC-3L"
SpaCET_obj <- create.SpaCET.object.10X(visiumPath = visiumPath)
SpaCET_obj <- SpaCET.quality.control(SpaCET_obj, min.genes = 1)

SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "QualityControl", spatialFeatures = c("UMI", "Gene"), imageBg = TRUE)
SpaCET_obj <- SpaCET.deconvolution(SpaCET_obj, cancerType = "LIHC", coreNo = 1)
SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "CellFraction", spatialFeatures = c("Malignant"))
